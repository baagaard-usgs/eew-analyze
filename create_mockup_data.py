#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import sys
import logging
from importlib import import_module
import numpy

import matplotlib.pyplot as pyplot
import matplotlib_extras

import h5py # avoid gdal loading incompatible HDF5 library
from eewperformance import analysisdb
from eewperformance import shakemap
from eewperformance import mockup_maps
from eewperformance import analysis_utils
from eewperformance import gdalraster

DEFAULTS = """
[event]
id = ci15481673
alerts = [0,4,21]

[shakemap]
projection = EPSG:3311

[shaking_time]
function = eewperformance.userdisplay.shaking_time_vs
#vs_kmps = 3.55 ; User display
#vs_kmps = 3.4 ; From NC record section
vs_kmps = 3.5 ; Avg NC/SC

[mmi_predicted]
function = eewperformance.shakemap.mmi_via_gmpe_gmice
gmpe = ASK2014
gmice = default

[alerts]
mmi_threshold = 3.5
magnitude_threshold = 3.95001

[maps]
tiler = cartopy_extra_tiles.esri_tiles.ESRI
tiler_style = streetmap
tiler_cache_dir = ~/data_scratch/images/tiles
zoom_level = 8
width_in = 5.0
height_in = 5.3

warning_time_contour_interval = 2.0

[files]
event_dir = ./data/[EVENTID]/
plots_dir = ./data/mockup/

analysis_db = ./data/analysisdb.sqlite
"""

# ----------------------------------------------------------------------
class Event(object):
    """Earthquake information for early warning system analysis.
    """

    def __init__(self, steps, config, eq_id):
        """Constructor.

        :type steps: ArgumentParser
        :param steps: Processing steps to perform.

        :type config: ConfigParser
        :param config: Analysis configuration.

        :type eq_id: str
        :param eq_id: ComCat earthquake id.
        """
        self.steps = steps
        self.config = config
        self.eqId = eq_id

        if steps.show_progress:
            self.showProgress = True

        self.db = None
        self.shakemap = None
        self.alerts = None
        self.event = None
        self.shakingTime = None
        return

    def process(self):
        """Perform processing steps for earthquake.

        :type multiprocessing: bool
        :param multiprocessing: True if using multiple processes.
        """
        self._load_data()

        if self.steps.process_event or self.steps.all:
            self._process_event()

        if self.steps.plot_event_maps or self.steps.all:
            self._plot_maps()

        return

    def _load_data(self):
        """Load data needed for processing event.
        """
        if self.showProgress:
            print("Loading data for {}...".format(self.eqId))

        # Database
        self.db = analysisdb.AnalysisData(self.config.get("files", "analysis_db"))
        
        # ShakeMap
        dataDir = analysis_utils.get_dir(self.config, "event_dir").replace("[EVENTID]", self.eqId)
        filename = os.path.join(dataDir, "custom_grid.xml.gz")
        if not os.path.isfile(filename):
            filename = os.path.join(dataDir, "grid.xml.gz")
        gmice = self.config.get("mmi_predicted", "gmice")
        if gmice == "WaldEtal1999":
            gmiceFn = shakemap.mmi_WaldEtal1999
            gmiceGrid = None
        elif gmice == "WordenEtal2012":
            gmiceFn = shakemap.mmi_WordenEtal2012
            gmiceGrid = None
        elif gmice == "default":
            gmiceFn = None
            shakemapInfo = self.db.comcat_shakemap(self.eqId)
            if "Wald99" in shakemapInfo["pgm2mi"]:
                gmiceGrid = "WaldEtal1999"
            elif "WGRW11" in shakemapInfo["pgm2mi"]:
                gmiceGrid = "WordenEtal2012"
            else:
                gmiceGrid = None
        else:
            raise ValueError("Unknown GMICE '{}'.".format(gmice))
        self.shakemap = shakemap.ShakeMap(gmiceFn)
        self.shakemap.load(filename)
        self.shakemap.gmiceGrid = gmiceGrid

        # Analsis DB data (event and alerts)
        self.event = self.db.comcat_event(self.eqId)
        server = self.config.get("shakealert.production", "server")
        alerts = self.db.alerts(self.eqId, server)
        self.alerts = [alerts[index] for index in map(int, analysis_utils.config_get_list(self.config.get("event", "alerts")))]
        
        # Shaking time
        functionPath = self.config.get("shaking_time", "function").split(".")
        fn = getattr(import_module(".".join(functionPath[:-1])), functionPath[-1])
        self.shakingTime = fn(self.event, self.shakemap.data, dict(self.config.items("shaking_time")))

        return

    def _process_event(self):
        """For given event, fetch data, process data, generate plots, and generate report.
        
        :type plotAlertMaps: bool
        :param plotAlertMaps: If true, plot map with predicted MMI and warning time contours for each alert.
        """
        magAlertThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        mmiAlertThreshold = self.config.getfloat("alerts", "mmi_threshold")
        if self.showProgress:
            print("Processing event {event[event_id]} with alert thresholds M{mag} and MMI {mmi} ...".format(event=self.event, mag=magAlertThreshold, mmi=mmiAlertThreshold))        

        functionPath = self.config.get("mmi_predicted", "function").split(".")
        fn = getattr(import_module(".".join(functionPath[:-1])), functionPath[-1])
            
        shape = self.shakemap.data["mmi"].shape
        warningTimeZero = numpy.zeros((1,), dtype="timedelta64[us]")
        warningTime = gdalraster.NO_DATA_VALUE * 1.0e+6 * numpy.ones(shape, dtype="timedelta64[us]")
        mmiPred = gdalraster.NO_DATA_VALUE * numpy.ones(shape, numpy.float32)
        alertVersion = gdalraster.NO_DATA_VALUE * numpy.ones(shape, numpy.int32)
        alertTime = gdalraster.NO_DATA_VALUE * 1.0e+6 * numpy.ones(shape, numpy.float32)

        gmpe = self.config.get("mmi_predicted", "gmpe")
        gmice = self.config.get("mmi_predicted", "gmice")
        if gmice == "default":
            gmice = self.shakemap.gmiceGrid
        
        thresholdReached = False
        for alert in self.alerts:
            if numpy.datetime64(alert["timestamp"]) > numpy.max(self.shakingTime):
                # Skip alerts with no positive warning times in
                # domain. Changes in estimated earthquake location
                # could result in later alerts having positive warning
                # times.
                logging.getLogger(__name__).debug("Skipping alert version {ver} with no positive warning times.".format(ver=alert["version"]))
                continue
            if not thresholdReached and alert["magnitude"] < magAlertThreshold:
                continue
            else:
                if not thresholdReached:
                    thresholdReached = True
                
            mmiPredCur = fn(alert, self.shakemap.data, gmpe, gmice)
            warningTimeCur = self.shakingTime - numpy.datetime64(alert["timestamp"])
            
            # Update alert time if greater than previous
            maskAlert = numpy.bitwise_and(warningTimeCur > warningTime, mmiPredCur >= mmiAlertThreshold)
            warningTime[maskAlert] = warningTimeCur[maskAlert]
            alertVersion[maskAlert] = alert["version"]
            alertTime[maskAlert] = numpy.datetime64(alert["timestamp"]) - numpy.datetime64(self.event["origin_time"])

            # Update predicted MMI if greater than previous AND
            # positive warning time. Assumes action will be taken if
            # alert threshold is reached (cannot be undone if later
            # updates reduce predicted MMI).
            maskMMI = numpy.bitwise_and(mmiPredCur > mmiPred, warningTimeCur >= warningTimeZero)
            mmiPred[maskMMI] = mmiPredCur[maskMMI]


            
        mmiObs = self.shakemap.data["mmi"]
        values = [
            ("mmi_obs", mmiObs,),
            ("mmi_pred", mmiPred,),
            ("alert_version", alertVersion,),
            ("alert_time", analysis_utils.timedelta_to_seconds(alertTime),),
            ("warning_time", analysis_utils.timedelta_to_seconds(warningTime),),
            ]
        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)

        server = self.config.get("shakealert.production", "server")
        gmpe = self.config.get("mmi_predicted", "gmpe")
        if magAlertThreshold is None:
            magAlertThreshold = params.getfloat("alerts", "magnitude_threshold")
        if mmiAlertThreshold is None:
            mmiAlertThreshold = params.getfloat("alerts", "mmi_threshold")
        label = "{eqId}-{server}-{gmpe}-M{magAlertThreshold:.1f}-MMI{mmiAlertThreshold:.1f}".format(
            eqId=self.event["event_id"], server=server, gmpe=gmpe, magAlertThreshold=magAlertThreshold, mmiAlertThreshold=mmiAlertThreshold)
        filename = "mockup_" + label + ".tiff"
        gdalraster.write(os.path.join(plotsDir, filename), values, self.shakemap.num_lon(), self.shakemap.num_lat(), self.shakemap.spatial_ref(), self.shakemap.geo_transform())
        return
    
    def _plot_maps(self):
        """Plot maps with analysis results for event.
        """
        if self.showProgress:
            print("Plotting maps for event {event[event_id]}...".format(event=self.event))

        selection = self.steps.plot_event_maps or "all"

        mapPanels = mockup_maps.EventMaps(self.config, self.eqId, self.event, self.alerts)
        mapPanels.load_data()
        if "mmi" in selection or "all" == selection:
            mapPanels.mmi_observed()
            mapPanels.mmi_predicted()
            mapPanels.mmi_residual()
        if "alert" in selection or "all" == selection:
            mapPanels.alert_version()
            mapPanels.alert_time()
            mapPanels.alert_warning_time()
        return

# ----------------------------------------------------------------------
class MockUpApp(object):
    """
    Analyze ShakeAlert performance using ShakeMap.
    """

    def __init__(self):
        """Constructor.
        """
        self.config = None
        self.showProgress = False
        return

    def main(self):
        """Main entry point
        """
        # Initialization
        args = self._parse_command_line()
        logLevel = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=logLevel, filename="analyzer.log")
        if args.show_progress:
            self.showProgress = True
        self.initialize(args.config)

        # Show parameters
        if args.show_parameters or args.all:
            self.show_parameters()

        # Show parameters
        if args.show_alerts or args.all:
            self.show_alerts()

        pyplot.style.use("size-presentation")
        pyplot.style.use("color-lightbg")
        matplotlib_extras.colors.add_general()
            
        # Event processing
        if args.process_event or args.plot_event_maps or args.all:
            eqId = self.config.get("event", "id")
            event = Event(args, self.config, eqId)
            event.process()
        return

    def initialize(self, config_filenames):
        """Set parameters from config file and DEFAULTS.

        :type config_filename: str
        :param config_filename: Name of configuration (INI) file with parameters.
        """
        import ConfigParser
        import io
        config = ConfigParser.SafeConfigParser()
        config.readfp(io.BytesIO(DEFAULTS))
        for filename in config_filenames.split(","):
            if self.showProgress:
                print("Fetching parameters from {}...".format(filename))
            config.read(filename)

        self.config = config

        return
    
    def show_parameters(self):
        """Write parameters to stdout.
        """
        self.config.write(sys.stdout)
        return

    def show_alerts(self):
        """Write parameters to stdout.
        """
        db = analysisdb.AnalysisData(self.config.get("files", "analysis_db"))
        server = self.config.get("shakealert.production", "server")
        eqId = self.config.get("event", "id")
        alerts = db.alerts(eqId, server)

        for alert in alerts:
            print("{alert[version]:3d} {alert[magnitude]:.1f} {alert[timestamp]}".format(alert=alert))
        return

    def _parse_command_line(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config", required=True)
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--show-alerts", action="store_true", dest="show_alerts")
        parser.add_argument("--process-event", action="store_true", dest="process_event", default=True)
        parser.add_argument("--plot-event-maps", action="store", dest="plot_event_maps", default=None, choices=[None, "all", "mmi", "alert"])
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug", default=True)
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    MockUpApp().main()


# End of file
