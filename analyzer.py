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
import datetime
import dateutil.parser
from importlib import import_module

from comcat import DetailEvent
from analysisdb import AnalysisData
from shakemap import ShakeMap

DEFAULTS = """
[events]
# Example:
# nc72923380 = Mw 4.6 Paicines, 2017-11-13

[shaking_time]
function = userdisplay.shaking_time_vs
vs = 3.0e+3

[mmi_predicted]
function = userdisplay.gmpe

[alerts]
mmi = 2.5

[map]

[files]
event_dir = ./data/[EVENTID]/
analysis_db = ./data/analysisdb.sqlite
"""

# ----------------------------------------------------------------------


def _config_get_list(list_string):
    """Convert list as string to list.

    :type list_string: list
    :param list_string: List as string.
    :returns: List of strings.
    """
    l = [f.strip() for f in list_string[1:-1].split(",")]
    return l


# ----------------------------------------------------------------------
class EEWAnalyzeApp(object):
    """
    Analyze ShakeAlert performance using ShakeMap.
    """
    
    def __init__(self):
        """Constructor.
        """
        self.params = None
        self.shakemap = None
        self.alerts = None
        self.warningTime = None
        self.alertRegion = None
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
            for eqId in self.params.options("events"):
                self.load_data(eqId)
                self.process_event(mmiAlert=self.params.get("alerts","mmi"))
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
                print("Fetching parameters from %s..." % filename)
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
        dataDir = self.params.get("files", "event_dir").replace("[EVENTID]", eqId)
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
        return
    
    def process_event(self, mmiAlert):
        """For given event, fetch data, process data, generate plots, and generate report.
        
        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).
        """
        if self.showProgress:
            print("Processing event {event[event_id]} with MMI alert={alert} ...".format(event=self.event, alert=mmiAlert))

        raise NotImplementedError(":TODO: @brad")
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
        parser.add_argument("--plot-map", action="store", dest="plot", default=None, choices=[None, "all", "map-mmi", "map-alert"])
        parser.add_argument("--generate-report", action="store_true", dest="generate_report")
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug")
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    EEWAnalyzeApp().main()


# End of file
