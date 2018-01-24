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
import numpy
import datetime
import dateutil.parser
from importlib import import_module

from comcat import DetailEvent
from analysisdb import AnalysisData
from shakemap import ShakeMap
from perfmetrics import CostSavings
from maps import MapPanels
from plotsxy import Figures
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
vs_kmps = 3.4 ; From NC record section

[mmi_predicted]
function = shakemap.gmpe
gmpe = ASK2014

[alerts]
#mmi_threshold = 0
#magnitude_threshold = 2.95
mmi_threshold = 1.5
magnitude_threshold = 4.45

[fragility_curves]
object = fragility_curves.PublicAnxiety
cost_action = 0.1
damage_low_mmi = 2.5
damata_high_mmi = 5.5

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
plots_dir = ./data/plots/

analysis_db = ./data/analysisdb.sqlite
analysis_event = ./data/[EVENTID]/analysis_data.tiff
population_density = ~/data/gis/populationdensity.tiff
"""

# ----------------------------------------------------------------------


def config_get_list(list_string):
    """Convert list as string to list.

    :type list_string: list
    :param list_string: List as string.
    :returns: List of strings.
    """
    l = [f.strip() for f in list_string[1:-1].split(",")]
    return l


def get_dir(params, name):
    """Get expanded directory name in [files] section.

    :type params: ConfigParser
    :param params: Configuration options

    :type name: str
    :param name: Option in [files] section.
    """
    return os.path.expanduser(params.get("files", name))

# ----------------------------------------------------------------------
class EEWAnalyzeApp(object):
    """
    Analyze ShakeAlert performance using ShakeMap.
    """
    
    def __init__(self):
        """Constructor.
        """
        self.params = None
        
        self.event = None
        self.shakemap = None
        self.alerts = None
        self.shakingTime = None
        self.populationDensity = None
        self.maps = None
        return

    def main(self):
        """Main entry point
        """
        # Initialization
        args = self._parse_command_line()
        logLevel = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=logLevel, filename="analyze_events.log")
        if args.show_progress:
            self.showProgress = True
        self.initialize(args.config)

        # Show parameters
        if args.show_parameters or args.all:
            self.show_parameters()

        if args.process_events or args.all:
            self.maps = MapPanels(self.params)
            self.maps.load_basemap()
            for eqId in self.params.options("events"):
                self.load_data(eqId)
                self.process_event(plotAlertMaps=args.plot_alert_maps)

        if args.plot_maps or args.all:
            if self.maps is None:
                self.maps = MapPanels(self.params)
                self.maps.load_basemap()
            for eqId in self.params.options("events"):
                maps = args.plot_maps if args.plot_maps else "all"
                self.plot_maps(eqId, maps)

        if args.plot_figures or args.all:
            for eqId in self.params.options("events"):
                self.load_data(eqId)
                selection = args.plot_figures if args.plot_figures else "all"
                self.plot_figures(eqId, selection)
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

        self.params = config

        return
    
    def show_parameters(self):
        """Write parameters to stdout.
        """
        self.params.write(sys.stdout)
        return

    def load_data(self, eqId):
        """Load ShakeMap and ShakeAlert data for event.

        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).
        """
        if self.showProgress:
            print("Loading data for {}...".format(eqId))

        # ShakeMap
        dataDir = get_dir(self.params, "event_dir").replace("[EVENTID]", eqId)
        filename = os.path.join(dataDir, "grid.xml.gz")
        self.shakemap = ShakeMap()
        self.shakemap.load(filename)

        # Analsis DB data (event and alerts)
        db = AnalysisData(self.params.get("files", "analysis_db"))
        self.event = db.comcat_event(eqId)
        self.alerts = db.alerts(eqId)

        # Shaking time
        functionPath = self.params.get("shaking_time", "function").split(".")
        fn = getattr(import_module(".".join(functionPath[:-1])), functionPath[-1])
        self.shakingTime = fn(self.event, self.shakemap.data, dict(self.params.items("shaking_time")))

        # Population density
        filename = get_dir(self.params, "population_density")
        self.populationDensity = gdalraster.resample(filename, self.shakemap.num_lon(), self.shakemap.num_lat(), self.shakemap.spatial_ref(), self.shakemap.geo_transform())
        return
    
    def process_event(self, mmiAlertThreshold=None, plotAlertMaps=False):
        """For given event, fetch data, process data, generate plots, and generate report.
        
        :type mmiAlertThreshold: float
        :param mmiAlertThreshold:
            MMI threshold for sending alert. Regions with predicted
            MMI above this threshold would receive an alert.
        """
        if mmiAlertThreshold is None:
            mmiAlertThreshold = self.params.getfloat("alerts","mmi_threshold")
        if self.showProgress:
            print("Processing event {event[event_id]} with MMI alert={alert} ...".format(event=self.event, alert=mmiAlertThreshold))
            
        costSavings = CostSavings(self.params, self.maps)
        stats = costSavings.compute(self.event, self.shakemap, self.alerts, self.shakingTime, self.populationDensity, mmiAlertThreshold, plotAlertMaps)
        print stats
        return

    def plot_maps(self, eqId, maps):
        """Plot maps with ShakeAlert performance information.

        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).
        """
        if self.maps is None:
            self.maps = MapPanels(self.params)
        self.maps.load_basemap()
        self.maps.load_data(eqId)
        if maps == "mmi" or maps =="all":
            self.maps.mmi_observed()
            self.maps.mmi_predicted()
            self.maps.mmi_warning_time()
            self.maps.mmi_residual()
        return
    
    def plot_figures(self, eqId, selection):
        """Plot figures with ShakeAlert performance information.

        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).
        """
        figures = Figures(self.params, self.event)
        if selection == "alert_error" or selection =="all":
            figures.alert_error(self.alerts)
        return
    
    def generate_report(self):
        """Assemble plots, etc into PDF file.
        """
        if self.showProgress:
            print("Generating report...")
        raise NotImplementedError(":TODO: @brad")
        return
        
    def _parse_command_line(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config", required=True)
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--process-events", action="store_true", dest="process_events")
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
