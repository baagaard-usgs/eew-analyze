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
import multiprocessing
import numpy

import h5py # avoid gdal loading incompatible HDF5 library
from eewperformance import analysisdb
from eewperformance import perfmetrics
from eewperformance import shakemap
from eewperformance import maps
from eewperformance import plotsxy
from eewperformance import reports
from eewperformance import analysis_utils
from eewperformance import gdalraster

DEFAULTS = """
[events]
# Example:
# nc72923380 = Mw 4.6 Paicines, 2017-11-13

[shakemap]
projection = EPSG:3311

[shaking_time]
function = eewperformance.userdisplay.shaking_time_vs
#vs_kmps = 3.55 ; User display
#vs_kmps = 3.4 ; From NC record section
vs_kmps = 3.5 ; From NC record section

[mmi_predicted]
function = eewperformance.shakemap.mmi_via_gmpe_gmice
gmpe = ASK2014
gmice = default

[alerts]
#mmi_threshold = 0.0
#magnitude_threshold = 2.95

mmi_threshold = 2.0
magnitude_threshold = 4.45

[fragility_curves]
object = eewperformance.fragility_curves.PublicFearAvoidance
cost_action = 0.1
damage_low_mmi = 2.5
damage_high_mmi = 5.5

#object = eewperformance.fragility_curves.PublicInjury
#cost_action = 0.1
#damage_low_mmi = 4.5
#damage_high_mmi = 7.5

#object = eewperformance.fragility_curves.StepDamage
#cost_action = 0.1
#damage_mmi = 3.5

[optimize]
mmi_threshold_min = 2.0
mmi_threshold_max = 4.0
mmi_threshold_step = 0.5

magnitude_threshold_min = 2.95
magnitude_threshold_max = 4.45
magnitude_threshold_step = 0.50

[qgis]
prefix_path = None
#prefix_path = /Applications/QGIS.app/Contents/MacOS
# PYTHONPATH=/Applications/QGIS.app/Contents/Resources/python:$PYTHONPATH

[maps]
tiler = cartopy_extra_tiles.esri_tiles.ESRI
tiler_style = streetmap
tiler_cache_dir = ~/data_scratch/images/tiles
#projection = EPSG:3857
zoom_level = 8
width_in = 8.0
height_in = 8.3

[files]
event_dir = ./data/[EVENTID]/
analysis_cache_dir = ./data/cache/
plots_dir = ./data/plots/
report = report.pdf

analysis_db = ./data/analysisdb.sqlite
population_density = ~/data/gis/populationdensity.tiff
"""

# ----------------------------------------------------------------------
def event_worker(event):
    """Perform EEW analysis of earthquake.

    :type event: Event
    :param event: Earthquake information
    """
    event.process()
    return


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
        self.populationDensity = None
        return

    def process(self):
        """Perform processing steps for earthquake.

        :type multiprocessing: bool
        :param multiprocessing: True if using multiple processes.
        """
        self._load_data()

        if self.steps.process_events or self.steps.all:
            self._process_event(plot_alert_maps=self.steps.plot_alert_maps)

        if self.steps.optimize_events or self.steps.all:
            self._optimize_thresholds()

        if self.steps.plot_maps or self.steps.all:
            self._plot_maps()

        if self.steps.plot_figures or self.steps.all:
            self._plot_figures()

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
            from shakemap import mmi_WaldEtal1999
            gmiceFn = mmi_WaldEtal1999
        elif gmice == "WordenEtal2012":
            from shakemap import mmi_WordenEtal2012
            gmiceFn = mmi_WordenEtal2012
        elif gmice == "default":
            gmiceFn = None
        else:
            raise ValueError("Unknown GMICE '{}'.".format(gmice))
        self.shakemap = shakemap.ShakeMap(gmiceFn)
        self.shakemap.load(filename)

        # Analsis DB data (event and alerts)
        self.event = self.db.comcat_event(self.eqId)
        server = self.config.get("shakealert.production", "server")
        if not server.startswith("catalog-magnitude"):
            self.alerts = self.db.alerts(self.eqId, server)
        else:
            self.alerts = [{
                "event_id": -111,
                "longitude": self.event["longitude"],
                "latitude": self.event["latitude"],
                "depth_km": self.event["depth_km"],
                "origin_time": self.event["origin_time"],
                "magnitude": self.event["magnitude"],
                "timestamp": self.event["origin_time"],
                }]
            if server == "catalog-magnitude-bias":
                bias = self.db.comcat_shakemap(self.eqId)["mmi_bias"]
                self.alerts[0]["magnitude"] += bias
                self.alerts[0]["event_id"] = -222

        # Shaking time
        functionPath = self.config.get("shaking_time", "function").split(".")
        fn = getattr(import_module(".".join(functionPath[:-1])), functionPath[-1])
        self.shakingTime = fn(self.event, self.shakemap.data, dict(self.config.items("shaking_time")))

        # Population density
        filename = analysis_utils.get_dir(self.config, "population_density")
        self.populationDensity = gdalraster.resample(filename, self.shakemap.num_lon(), self.shakemap.num_lat(), self.shakemap.spatial_ref(), self.shakemap.geo_transform())

        return
    
    def _process_event(self, plot_alert_maps=False):
        """For given event, fetch data, process data, generate plots, and generate report.
        
        :type plotAlertMaps: bool
        :param plotAlertMaps: If true, plot map with predicted MMI and warning time contours for each alert.
        """
        magAlertThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        mmiAlertThreshold = self.config.getfloat("alerts", "mmi_threshold")
        if self.showProgress:
            print("Processing event {event[event_id]} with alert thresholds M{mag} and MMI {mmi} ...".format(event=self.event, mag=magAlertThreshold, mmi=mmiAlertThreshold))
            
        costSavings = perfmetrics.CostSavings(self.config)
        stats = costSavings.compute(self.event, self.shakemap, self.alerts, self.shakingTime, self.populationDensity, magAlertThreshold, mmiAlertThreshold, plot_alert_maps)
        stats.update({
            "comcat_id": self.event["event_id"],
            "eew_server": self.config.get("shakealert.production", "server"),
            "dm_id": self.alerts[0]["event_id"] if len(self.alerts) > 0 else -1,
            "dm_timestamp": self.alerts[0]["timestamp"] if len(self.alerts) > 0 else "",
            "gmpe": self.config.get("mmi_predicted", "gmpe"),
            "fragility": self.config.get("fragility_curves", "object").split(".")[-1],
            "magnitude_threshold": magAlertThreshold,
            "mmi_threshold": mmiAlertThreshold,
            })
        self.db.add_performance(stats, replace=True)

        return

    def _optimize_thresholds(self):
        """Determine optimum threshold by looping over range of alert
        thresholds for earthquake magnitude and MMI.

        Note: Results are added to analysis database for extraction
        and determination of the optimum value later.
        """
        thresholdStart = self.config.getfloat("optimize", "mmi_threshold_min")
        thresholdStop = self.config.getfloat("optimize", "mmi_threshold_max")
        thresholdStep = self.config.getfloat("optimize", "mmi_threshold_step")
        mmiThresholds = numpy.arange(thresholdStart, thresholdStop+0.1*thresholdStep, thresholdStep)

        thresholdStart = self.config.getfloat("optimize", "magnitude_threshold_min")
        thresholdStop = self.config.getfloat("optimize", "magnitude_threshold_max")
        thresholdStep = self.config.getfloat("optimize", "magnitude_threshold_step")
        magThresholds = numpy.arange(thresholdStart, thresholdStop+0.1*thresholdStep, thresholdStep)

        statsExtra = {
            "comcat_id": self.event["event_id"],
            "eew_server": self.config.get("shakealert.production", "server"),
            "dm_id": self.alerts[0]["event_id"],
            "dm_timestamp": self.alerts[0]["timestamp"],
            "gmpe": self.config.get("mmi_predicted", "gmpe"),
            "fragility": self.config.get("fragility_curves", "object").split(".")[-1],
            }
        
        costSavings = perfmetrics.CostSavings(self.config)
        for magnitude in magThresholds:
            for mmi in mmiThresholds:
                if self.showProgress:
                    print("Processing event {event[event_id]} with alert thresholds M{mag} and MMI {mmi} ...".format(event=self.event, mag=magnitude, mmi=mmi))

                stats = costSavings.compute(self.event, self.shakemap, self.alerts, self.shakingTime, self.populationDensity, magnitude, mmi)
                stats.update(statsExtra)
                stats["magnitude_threshold"] = magnitude
                stats["mmi_threshold"] = mmi
                self.db.add_performance(stats, replace=True)
        return

    def _plot_maps(self):
        """Plot maps with analysis results for event.
        """
        if self.showProgress:
            print("Plotting maps for event {event[event_id]}...".format(event=self.event))

        selection = self.steps.plot_maps if self.steps.plot_maps else "all"

        mapPanels = maps.MapPanels(self.config, self.eqId, self.event, self.alerts)
        mapPanels.load_data()
        if selection == "mmi" or selection == "all":
            mapPanels.mmi_observed()
            mapPanels.mmi_predicted()
            mapPanels.mmi_residual()
        if selection == "alert" or selection == "all":
            mapPanels.alert_category()
        return

    def _plot_figures(self):
        """Plot xy figures with analysis results for event.
        """
        if self.showProgress:
            print("Plotting figures for event {event[event_id]}...".format(event=self.event))

        selection = self.steps.plot_figures if self.steps.plot_figures else "all"
        figures = plotsxy.Figures(self.config, self.event)
        if selection == "alert_error" or selection == "all":
            figures.alert_error(self.alerts)
        return


# ----------------------------------------------------------------------
class EEWAnalyzeApp(object):
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

        # Event processing
        if args.process_events or args.optimize_events or args.plot_maps or args.plot_figures or args.all:
            if args.nthreads <= 0:
                for eqId in self.config.options("events"):
                    event = Event(args, self.config, eqId)
                    event.process()
            else:
                pool = multiprocessing.Pool(args.nthreads)
                for eqId in self.config.options("events"):
                    event = Event(args, self.config, eqId)
                    pool.apply_async(event_worker, args=(event,))
                pool.close()
                pool.join()

        # Generate report
        if args.generate_report or args.all:
            self.generate_report()
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

    def generate_report(self):
        """Assemble plots, etc into PDF file.
        """
        if self.showProgress:
            print("Generating report...")

        summary = reports.AnalysisSummary(self.config)
        summary.generate(self.config.options("events"))
        return
        
    def _parse_command_line(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config", required=True)
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--process-events", action="store_true", dest="process_events")
        parser.add_argument("--optimize-events", action="store_true", dest="optimize_events")
        parser.add_argument("--plot-alert-maps", action="store_true", dest="plot_alert_maps")
        parser.add_argument("--plot-maps", action="store", dest="plot_maps", default=None, choices=[None, "all", "mmi", "alert"])
        parser.add_argument("--plot-figures", action="store", dest="plot_figures", default=None, choices=[None, "all", "alert_error"])
        parser.add_argument("--generate-report", action="store_true", dest="generate_report")
        parser.add_argument("--num-threads", action="store", type=int, dest="nthreads", default=1)
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug", default=True)
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    EEWAnalyzeApp().main()


# End of file
