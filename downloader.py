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

from comcat import DetailEvent
from shakealert import EEWServer, DMLogXML
from analysisdb import AnalysisData

DEFAULTS = """
[events]
# Example: nc72923380 = Mw 4.6 Paicines, 2017-11-13

[shakealert.production]
login_url = None
log_url = None
server = None
username = None
password = None

[shakealert.demonstration]
login_url = None
log_url = None
server = None
username = None
password = None

[files]
dmlogs_dir = ./data/dmlogs/
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


def _data_filename(params, option, eqId, args=None, makeDir=False):
    """Construct relative path for file.

    :type params: ConfigParser
    :param params: Application parameters.
    :type option: str
    :param option: Option in application parameters for filename.
    :type eqId: str
    :param ComCat event id (e.g., nc72923380).
    :param params: Application parameters.
    :type args: tuple
    :param args: Tuple for arguments for substitution in name of parameter file.
    :type makeDir: bool
    :param makeDir: Create directory for file if True.
    """
    eventDir = os.path.join("data", eqId)
    if makeDir or not os.path.isdir(eventDir):
        os.makedirs(eventDir)

    filename = params.get("files", option) if args is None else params.get("files", option) % args
    return os.path.join(eventDir, filename)
    

# ----------------------------------------------------------------------
class DownloaderApp(object):
    """
    Download data (DM logs, ShakeMaps, ComCat events) needed for analysis of ShakeAlert performance.
    """
    
    def __init__(self):
        """Constructor.
        """
        self.params = None
        return

    def main(self):
        """Main entry point
        """
        # Initialization
        args = self._parseCommandLine()
        logLevel = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=logLevel, filename="downloader.log")
        if args.show_progress:
            self.showProgress = True
        self.initialize(args.config)

        # Show parameters
        if args.show_parameters or args.all:
            self.show_parameters()

        if args.fetch_eewalerts:
            beginStr,endStr = args.fetch_eewalerts.split(",")
            dateBegin = dateutil.parser.parse(beginStr).date()
            dateEnd = dateutil.parser.parse(endStr).date()
            self._fetch_eewalerts(dateBegin, dateEnd)

        if args.fetch_events or args.all:
            self._fetch_comcat_events()

        if args.fetch_shakemaps or args.all:
            self._fetch_shakemaps()

        if args.initdb or args.all:
            self._initdb()

        if args.dbstatus or args.all:
            self._dbstatus()

        if args.updatedb or args.all:
            if args.updatedb == "eew_alerts" or args.all:
                self._updatedb_eewalerts()
            if args.updatedb == "comcat_events" or args.all:
                self._updatedb_events()
            if args.updatedb == "matches" or args.all:
                self._updatedb_matches()

        if args.show_matches or args.all:
            self._show_matches()
        return

    def _fetch_eewalerts(self, dateBegin, dateEnd):
        if self.showProgress:
            print("Fetching EEW alerts...")
            
        self.eewserver = EEWServer(self.params)
        self.eewserver.login()

        dmlog = DMLogXML(config=self.params)
        logsDir = self.params.get("files", "dmlogs_dir")
        if not os.path.isdir(logsDir):
            os.makedirs(logsDir)
        date = dateBegin
        dt = datetime.timedelta(days=1)
        while date <= dateEnd:
            dmlog.fetch(self.eewserver, date, logsDir)
            date += dt
        self.eewserver.logout()
        return

    def _fetch_comcat_events(self):
        """Fetch geojson event file from USGS ComCat using web services.
        
        :type eqId: string
        :param eqId: ComCat event id (e.g., nc72923380).
        """
        if self.showProgress:
            print("Fetching earthquakes from ComCat database...")

        event = DetailEvent()
        dirTemplate = self.params.get("files", "event_dir")
        for eqId in self.params.options("events"):
            dataDir = dirTemplate.replace("[EVENTID]", eqId)
            event.fetch(eqId, dataDir)
        return

    def _fetch_shakemaps(self):
        """Fetch geojson event file from USGS ComCat using web services.
        """
        if self.showProgress:
            print("Fetching ShakeMaps...")

        event = DetailEvent()
        dirTemplate = self.params.get("files", "event_dir")
        for eqId in self.params.options("events"):
            dataDir = dirTemplate.replace("[EVENTID]", eqId)
            event.load(os.path.join(dataDir, eqId+".geojson"))
            shakemap = event.get_product("shakemap")
            if len(shakemap) > 1:
                raise ValueError("Expected to get preferred ShakeMap.")
            shakemap[0].fetch("grid.xml", dataDir)
        return

    def _initdb(self):
        """Create analysis database with ShakeAlert DM alerts and ComCat events.
        """
        if self.showProgress:
            print("Setting up analysis database...")
            
        db = AnalysisData(self.params.get("files", "analysis_db"))
        db.init()
        logging.getLogger(__name__).info(db.summary())
        return

    def _dbstatus(self):
        """Show summary of analysis database contents.
        """
        db = AnalysisData(self.params.get("files", "analysis_db"))
        print(db.summary())
        return
    
    def _updatedb_events(self):
        if self.showProgress:
            print("Updating ComCat events in analysis database...")

        db = AnalysisData(self.params.get("files", "analysis_db"))

        event = DetailEvent()
        dirTemplate = self.params.get("files", "event_dir")
        events = self.params.options("events")
        numEvents = len(events)
        for iEvent,eqId in enumerate(events):
            if self.showProgress:
                sys.stdout.write("\rProcessing ComCat events...{:d}%%".format(((iEvent+1)*100)/numEvents))
                sys.stdout.flush()

            dataDir = dirTemplate.replace("[EVENTID]", eqId)
            event.load(os.path.join(dataDir, eqId+".geojson"))
            db.add_event(event)
        sys.stdout.write("\n")
        return
    
    def _updatedb_eewalerts(self, all=False):
        """
        """
        import glob
        
        if self.showProgress:
            print("Updating ShakeAlert DM alerts in analysis database...")

        db = AnalysisData(self.params.get("files", "analysis_db"))

        logsDir = self.params.get("files", "dmlogs_dir")
        files = sorted(glob.glob(os.path.join(logsDir, "dmevent*.log.gz")))
        if not all:
            # Get most recent entry and use files starting on that day.
            pass

        # Read DM logs
        dmlog = DMLogXML(config=self.params)
        numFiles = len(files)
        for iFile,filename in enumerate(files):
            if self.showProgress:
                sys.stdout.write("\rProcessing DM logs...{:d}%%".format(((iFile+1)*100)/numFiles))
                sys.stdout.flush()
            alerts = dmlog.load(filename)
            db.add_alerts(alerts)
        sys.stdout.write("\n")
        return

    def _updatedb_matches(self):
        if self.showProgress:
            print("Finding EEW alerts matchings ComCat events...")

        db = AnalysisData(self.params.get("files", "analysis_db"))
        db.find_matches()
        return 
    
    def _show_matches(self):
        if self.showProgress:
            print("Showing matches between ComCat and ShakeAlert...")

        db = AnalysisData(self.params.get("files", "analysis_db"))
        db.show_matches()
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
        if config_filenames:
            for filename in config_filenames.split(","):
                if self.showProgress:
                    print("Fetching parameters from {}...".format(filename))
                config.read(filename)

        self.params = config

        return
    
    def show_parameters(self):
        """Write parameters to stdout.
        """
        import sys
        self.params.write(sys.stdout)
        return

    def _parseCommandLine(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config")
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--fetch-eewalerts", action="store", dest="fetch_eewalerts", default=None)
        parser.add_argument("--fetch-events", action="store_true", dest="fetch_events")
        parser.add_argument("--fetch-shakemaps", action="store_true", dest="fetch_shakemaps")
        parser.add_argument("--initdb", action="store_true", dest="initdb")
        parser.add_argument("--dbstatus", action="store_true", dest="dbstatus")
        parser.add_argument("--updatedb", action="store", dest="updatedb", choices=["eew_alerts", "comcat_events", "matches"])
        parser.add_argument("--show-matches", action="store_true", dest="show_matches")
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug", default=True)
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    DownloaderApp().main()


# End of file
