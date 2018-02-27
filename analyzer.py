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
from multiprocessing import Pool
import numpy

from shakemap import ShakeMap # openquake before osgeo
from analysisdb import AnalysisData
from perfmetrics import CostSavings
from maps import MapPanels
from plotsxy import Figures
from reports import AnalysisSummary
import analysis_utils
import gdalraster

DEFAULTS = """
[events]
# Example:
# nc72923380 = Mw 4.6 Paicines, 2017-11-13

[shakemap]
projection = EPSG:3311

[shaking_time]
function = userdisplay.shaking_time_vs
#vs_kmps = 3.55 ; User display
#vs_kmps = 3.4 ; From NC record section
vs_kmps = 3.5 ; From NC record section

[mmi_predicted]
function = shakemap.mmi_via_gmpe_gmice
gmpe = ASK2014
#gmice = WordenEtal2012
gmice = WaldEtal1999

[alerts]
#mmi_threshold = 0.0
mmi_threshold = 2.0
magnitude_threshold = 2.95

magnitude_threshold = 3.95
#magnitude_threshold = 4.45

[fragility_curves]
object = fragility_curves.PublicFearAvoidance
cost_action = 0.1
damage_low_mmi = 2.5
damage_high_mmi = 5.5

[optimize]
mmi_threshold_min = 2.0
mmi_threshold_max = 4.0
mmi_threshold_step = 0.5

[qgis]
prefix_path = None
#prefix_path = /Applications/QGIS.app/Contents/MacOS
# PYTHONPATH=/Applications/QGIS.app/Contents/Resources/python:$PYTHONPATH

[maps]
projection = EPSG:3857
width_pixels = 1280
height_pixels = 1024
bg_color = white
basemap = esri_streetmap.xml

[files]
event_dir = ./data/[EVENTID]/
analysis_cache_dir = ./data/cache/
plots_dir = ./data/plots/
report = report.pdf

analysis_db = ./data/analysisdb_NEW.sqlite
population_density = ~/data/gis/populationdensity.tiff
"""

# ----------------------------------------------------------------------
class EEWAnalyzeApp(object):
    """
    Analyze ShakeAlert performance using ShakeMap.
    """

    def __init__(self):
        """Constructor.
        """
        self.config = None

        self.event = None
        self.shakemap = None
        self.alerts = None
        self.shakingTime = None
        self.populationDensity = None
        self.maps = None
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

        if args.process_events or args.all:
            self.maps = MapPanels(self.config)
            self.maps.load_basemap()
            for eqId in self.config.options("events"):
                self.load_data(eqId)
                self.process_event(plotAlertMaps=args.plot_alert_maps)

        if args.optimize_events or args.all:
            self.maps = None
            for eqId in self.config.options("events"):
                self.load_data(eqId)
                self.optimize_threshold()

        if args.plot_maps or args.all:
            if self.maps is None:
                self.maps = MapPanels(self.config)
                self.maps.load_basemap()
            for eqId in self.config.options("events"):
                maps = args.plot_maps if args.plot_maps else "all"
                self.plot_maps(eqId, maps)

        if args.plot_figures or args.all:
            for eqId in self.config.options("events"):
                self.load_data(eqId)
                selection = args.plot_figures if args.plot_figures else "all"
                self.plot_figures(eqId, selection)

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

    def load_data(self, eqId):
        """Load ShakeMap and ShakeAlert data for event.

        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).
        """
        if self.showProgress:
            print("Loading data for {}...".format(eqId))

        # ShakeMap
        dataDir = analysis_utils.get_dir(self.config, "event_dir").replace("[EVENTID]", eqId)
        filename = os.path.join(dataDir, "custom_grid.xml.gz")
        if not os.path.isfile(filename):
            filename = os.path.join(dataDir, "grid.xml.gz")
        self.shakemap = ShakeMap()
        self.shakemap.load(filename)

        # Analsis DB data (event and alerts)
        db = AnalysisData(self.config.get("files", "analysis_db"))
        self.event = db.comcat_event(eqId)
        server = self.config.get("shakealert.production", "server")
        self.alerts = db.alerts(eqId, server)

        # Shaking time
        functionPath = self.config.get("shaking_time", "function").split(".")
        fn = getattr(import_module(".".join(functionPath[:-1])), functionPath[-1])
        self.shakingTime = fn(self.event, self.shakemap.data, dict(self.config.items("shaking_time")))

        # Population density
        filename = analysis_utils.get_dir(self.config, "population_density")
        self.populationDensity = gdalraster.resample(filename, self.shakemap.num_lon(), self.shakemap.num_lat(), self.shakemap.spatial_ref(), self.shakemap.geo_transform())
        return
    
    def process_event(self, plotAlertMaps=False):
        """For given event, fetch data, process data, generate plots, and generate report.
        
        :type plotAlertMaps: bool
        :param plotAlertMaps: If true, plot map with predicted MMI and warning time contours for each alert.
        """
        magAlertThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        mmiAlertThreshold = self.config.getfloat("alerts", "mmi_threshold")
        if self.showProgress:
            print("Processing event {event[event_id]} with alert thresholds M{mag} and MMI {mmi} ...".format(event=self.event, mag=magAlertThreshold, mmi=mmiAlertThreshold))
            
        costSavings = CostSavings(self.config, self.maps)
        stats = costSavings.compute(self.event, self.shakemap, self.alerts, self.shakingTime, self.populationDensity, mmiAlertThreshold, plotAlertMaps)
        stats.update({
            "comcat_id": self.event["event_id"],
            "eew_server": self.config.get("shakealert.production", "server"),
            "dm_id": self.alerts[0]["event_id"],
            "dm_timestamp": self.alerts[0]["timestamp"],
            "gmpe": self.config.get("mmi_predicted", "gmpe"),
            "fragility": self.config.get("fragility_curves", "object").split(".")[-1],
            "magnitude_threshold": magAlertThreshold,
            "mmi_threshold": mmiAlertThreshold,
            })
        db = AnalysisData(self.config.get("files", "analysis_db"))            
        db.add_performance(stats, replace=True)
        return

    def optimize_threshold(self):
        thresholdStart = self.config.getfloat("optimize", "mmi_threshold_min")
        thresholdStop = self.config.getfloat("optimize", "mmi_threshold_max")
        thresholdStep = self.config.getfloat("optimize", "mmi_threshold_step")
        thresholds = numpy.arange(thresholdStart, thresholdStop+0.1*thresholdStep, thresholdStep)

        statsExtra = {
            "comcat_id": self.event["event_id"],
            "eew_server": self.config.get("shakealert.production", "server"),
            "dm_id": self.alerts[0]["event_id"],
            "dm_timestamp": self.alerts[0]["timestamp"],
            "gmpe": self.config.get("mmi_predicted", "gmpe"),
            "fragility": self.config.get("fragility_curves", "object").split(".")[-1],
            "magnitude_threshold": self.config.getfloat("alerts", "magnitude_threshold"),
            }
        
        db = AnalysisData(self.config.get("files", "analysis_db"))            
        
        costSavings = CostSavings(self.config, self.maps)
        #threadPool = Pool(self.maxthreads)
        for iperf, threshold in enumerate(thresholds):
            #[threadPool.apply_async(self._optimize_worker, (threshold, statsExtra, db) for threshold in thresholds)]
            stats = costSavings.compute(self.event, self.shakemap, self.alerts, self.shakingTime, self.populationDensity, threshold, plotAlertMaps=False)
            stats.update(statsExtra)
            stats["mmi_threshold"] = threshold
            db.add_performance(stats, replace=True)

        return

    def _optimize_worker(self, threshold, statsExtra, db):
        stats = costSavings.compute(self.event, self.shakemap, self.alerts, self.shakingTime, self.populationDensity, threshold, plotAlertMaps=False)
        stats.update(statsExtra)
        stats["mmi_threshold"] = threshold
        db.add_performance(stats, replace=True)
        return
            
    
    def plot_maps(self, eqId, maps):
        """Plot maps with ShakeAlert performance information.

        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).
        """
        if self.maps is None:
            self.maps = MapPanels(self.config)
        self.maps.load_basemap()
        self.maps.load_data(eqId)
        if maps == "mmi" or maps == "all":
            self.maps.mmi_observed()
            self.maps.mmi_predicted()
            self.maps.mmi_residual()
            self.maps.alert_category()
        return
    
    def plot_figures(self, eqId, selection):
        """Plot figures with ShakeAlert performance information.

        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).
        """
        figures = Figures(self.config, self.event)
        if selection == "alert_error" or selection == "all":
            figures.alert_error(self.alerts)
        return
    
    def generate_report(self):
        """Assemble plots, etc into PDF file.
        """
        if self.showProgress:
            print("Generating report...")

        summary = AnalysisSummary(self.config)
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
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug", default=True)
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    EEWAnalyzeApp().main()


# End of file
