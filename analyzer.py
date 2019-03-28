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
import argparse
from importlib import import_module
import multiprocessing
import numpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pyplot
import matplotlib_extras.colors

import h5py # avoid gdal loading incompatible HDF5 library
from eewperformance import analysisdb
from eewperformance import perfmetrics
from eewperformance import shakemap
from eewperformance import maps
from eewperformance import plotsxy
from eewperformance import reports
from eewperformance import analysis_utils
from eewperformance import gdalraster
from eewperformance import local_color

DEFAULTS = u"""
[events]
# Example:
# nc72923380 = Mw 4.6 Paicines, 2017-11-13

[shakealert.production]
login_url = None
log_url = None
server = eew-bk-prod1
username = None
password = None

[shakealert.demonstration]
login_url = None
log_url = None
server = eew2demo
username = None
password = None

[shakemap]
projection = EPSG:3311

[shaking_time]
function = eewperformance.userdisplay.shaking_time_vs
vs_kmps = 3.5
vp_kmps = 6.1

[mmi_predicted]
function = eewperformance.shakemap.mmi_via_gmpe_gmice
gmpe = ASK2014
gmice = default

[alerts]
# Current User Display?
#mmi_threshold = 0.0
#magnitude_threshold = 2.95001

# Proposed public
#mmi_threshold = 2.0
#magnitude_threshold = 4.45001

# Optimum thresholds
mmi_threshold = 3.5
magnitude_threshold = 3.95001

[fragility_curves]
object = eewperformance.fragility_curves.LinearDamage
label = FearAvoidanceLinear

[optimize]
mmi_threshold_min = 2.0
mmi_threshold_max = 4.0
mmi_threshold_step = 0.5

magnitude_threshold_min = 3.45001
magnitude_threshold_max = 4.45001
magnitude_threshold_step = 0.50

[maps]
tiler = cartopy_extra_tiles.esri_tiles.ESRI
tiler_style = streetmap
tiler_cache_dir = ~/data_scratch/images/tiles
zoom_level = 8
width_in = 5.0
height_in = 5.3

warning_time_contour_interval = 2.0

[plots]
raster = False

[files]
event_dir = ./data/[EVENTID]/
analysis_cache_dir = ./data/cache/
plots_dir = ./data/plots/
report = report.pdf

analysis_db = ./data/analysisdb.sqlite
population_density = ~/data/gis/census/populationdensity.tiff
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

        if self.steps.plot_event_maps or self.steps.all:
            self._plot_maps()

        if self.steps.plot_event_figures or self.steps.all:
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

        if not self.config.has_section("theoretical"):
            self.alerts = self.db.alerts(self.eqId, server)
        else:
            event_id = self.config.get("theoretical", "event_id")
            use_first_alert_time = self.config.getboolean("theoretical", "use_first_alert_time")
            latency = self.config.getfloat("theoretical", "alert_offset")
            add_bias = self.config.getboolean("theoretical", "add_magnitude_bias")
            bias = self.db.comcat_shakemap(self.eqId)["mmi_bias"] if add_bias else 0.0
            if use_first_alert_time:
                alerts = self.db.alerts(self.eqId, server)
                if len(alerts) == 0:
                    self.alerts = []
                else:
                    alert_time = self.db.alerts(self.eqId, server)[0]["timestamp"]
                    self.alerts = [{
                        "event_id": event_id,
                        "longitude": self.event["longitude"],
                        "latitude": self.event["latitude"],
                        "depth_km": self.event["depth_km"],
                        "origin_time": self.event["origin_time"],
                        "magnitude": self.event["magnitude"] + bias,
                        "timestamp": alert_time,
                    }]
            else:
                self.alerts = [{
                    "event_id": event_id,
                    "longitude": self.event["longitude"],
                    "latitude": self.event["latitude"],
                    "depth_km": self.event["depth_km"],
                    "origin_time": self.event["origin_time"],
                    "magnitude": self.event["magnitude"] + bias,
                    "timestamp": self.event["origin_time"],
                }]

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
            "fragility": self.config.get("fragility_curves", "label"),
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
            "dm_id": self.alerts[0]["event_id"] if len(self.alerts) > 0 else -1,
            "dm_timestamp": self.alerts[0]["timestamp"] if len(self.alerts) > 0 else "",
            "gmpe": self.config.get("mmi_predicted", "gmpe"),
            "fragility": self.config.get("fragility_curves", "label"),
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

        selection = self.steps.plot_event_maps or "all"

        mapPanels = maps.EventMaps(self.config, self.eqId, self.event, self.alerts)
        mapPanels.load_data()
        if "mmi" in selection or "all" == selection:
            mapPanels.mmi_observed()
            mapPanels.mmi_predicted()
            mapPanels.mmi_residual()
        if "alert" in selection or "all" == selection:
            mapPanels.alert_category()
            mapPanels.cost_savings()
        return

    def _plot_figures(self):
        """Plot xy figures with analysis results for event.
        """
        if self.showProgress:
            print("Plotting figures for event {event[event_id]}...".format(event=self.event))

        selection = self.steps.plot_event_figures or "all"
        figures = plotsxy.EventFigures(self.config, self.event)
        if "alert_error" in selection or "all" == selection:
            mmi_bias = self.db.comcat_shakemap(self.eqId)["mmi_bias"]
            figures.alert_error(self.alerts, mmi_bias)
        if "mmi_correlation" in selection or "all" == selection:
            figures.mmi_correlation()
        if "warning_time_mmi" in selection or "all" == selection:
            figures.warning_time_mmi()
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

    def main(self, **kwargs):
        """Main entry point
        """
        # Initialization
        args = argparse.Namespace(**kwargs) if kwargs else self._parse_command_line()
        logLevel = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=logLevel, filename="analyzer.log")
        if args.show_progress:
            self.showProgress = True
        self.initialize(args.config)

        # Show parameters
        if args.show_parameters or args.all:
            self.show_parameters()

        pyplot.style.use("size-presentation")
        pyplot.style.use("color-"+args.color_style)
        matplotlib_extras.colors.add_general()
        if args.color_style == "lightbg":
            local_color.lightbg()
        else:
            local_color.darkbg()
            
        # Event processing
        if args.process_events or args.optimize_events or args.plot_event_maps or args.plot_event_figures or args.all:
            if args.nthreads <= 0:
                for eqId in self.config.options("events"):
                    event = Event(args, self.config, eqId)
                    event.process()
            else:
                pool = multiprocessing.Pool(args.nthreads)
                result = []
                for eqId in self.config.options("events"):
                    event = Event(args, self.config, eqId)
                    r = pool.apply_async(event_worker, args=(event,))
                    result.append(r)
                for r in result:
                    r.get()
                pool.close()
                pool.join()

        # Summary maps and figures
        if args.plot_summary_maps or args.all:
            self.plot_summary_maps(args.plot_summary_maps or args.all)
        if args.plot_summary_figures or args.all:
            self.plot_summary_figures(args.plot_summary_figures or args.all)

        # Generate report
        if args.generate_report or args.all:
            self.generate_report("True" if args.generate_report == "summary" else False)
        return

    def initialize(self, config_filenames):
        """Set parameters from config file and DEFAULTS.

        :type config_filename: str
        :param config_filename: Name of configuration (INI) file with parameters.
        """
        import io
        import six
        if six.PY2:
            import ConfigParser
            config = ConfigParser.SafeConfigParser()
        else:    
            import configparser
            config = configparser.SafeConfigParser()
        config.readfp(io.StringIO(DEFAULTS))
        for filename in config_filenames.split(","):
            if not os.path.isfile(filename):
                raise IOError("Could not find configuration file '{}'.".format(filename))
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

    def plot_summary_maps(self, selection):
        """Plot summary maps.
        """
        if self.showProgress:
            print("Plotting summary maps...")

        db = analysisdb.AnalysisData(self.config.get("files", "analysis_db"))
        events = self.config.options("events")
        figures = maps.SummaryMaps(self.config, events, db)

        if "events" in selection or "all" == selection:
            figures.earthquakes()
        if "performance" in selection or "all" == selection:
            figures.cost_savings("area_costsavings_eew")
            figures.cost_savings("population_costsavings_eew")
        return
    
    def plot_summary_figures(self, selection):
        """Plot summary figures.
        """
        if self.showProgress:
            print("Plotting summary figures...")

        db = analysisdb.AnalysisData(self.config.get("files", "analysis_db"))
        events = self.config.options("events")
        figures = plotsxy.SummaryFigures(self.config, events, db)

        if "magnitude_time" in selection or "all" == selection:
            figures.magnitude_versus_time()
        if "optimum_thresholds" in selection or "all" == selection:
            figures.optimal_mmithresholds()
        if "metric_time" in selection or "all" == selection:
            figures.costsavings_versus_time()
        if "metric_magnitude" in selection or "all" == selection:
            figures.costsavings_versus_magnitude()
        if "cost_functions" in selection or "all" == selection:
            figures.cost_functions()
        if "metric_cost_functions" in selection or "all" == selection:
            figures.metric_cost_functions()
        if "metric_theoretical" in selection or "all" == selection:
            figures.metric_theoretical()
        return
    
    def generate_report(self, summary_only=True):
        """Assemble plots, etc into PDF file.
        """
        if self.showProgress:
            print("Generating report...")

        summary = reports.AnalysisSummary(self.config, summary_only)
        summary.generate(self.config.options("events"))
        return
        
    def _parse_command_line(self):
        """Parse command line arguments.
        """
        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config", required=True)
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--process-events", action="store_true", dest="process_events")
        parser.add_argument("--optimize-events", action="store_true", dest="optimize_events")
        parser.add_argument("--plot-alert-maps", action="store_true", dest="plot_alert_maps")
        parser.add_argument("--plot-event-maps", action="store", dest="plot_event_maps", default=None, choices=[None, "all", "mmi", "alert"])
        parser.add_argument("--plot-event-figures", action="store", dest="plot_event_figures", default=None, choices=[None, "all", "alert_error", "mmi_correlation", "warning_time_mmi"])
        parser.add_argument("--plot-summary-maps", action="store", dest="plot_summary_maps", default=None, choices=[None, "all", "events", "performance"])
        parser.add_argument("--plot-summary-figures", action="store", dest="plot_summary_figures", default=None, choices=[None, "all", "magnitude_time", "optimum_thresholds", "metric_time", "metric_magnitude", "cost_functions", "metric_cost_functions", "metric_theoretical"])
        parser.add_argument("--generate-report", action="store", dest="generate_report", default=None, choices=[None,"summary", "full"])
        parser.add_argument("--num-threads", action="store", type=int, dest="nthreads", default=0)
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug", default=True)
        parser.add_argument("--color-style", action="store", dest="color_style", default="lightbg", choices=["darkbg", "lightbg"])
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    EEWAnalyzeApp().main()


# End of file
