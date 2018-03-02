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
import gzip

from eewperformance import comcat
from eewperformance import shakealert
from eewperformance import analysisdb

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
dmlogs_dir = ./data/dmlogs/[SERVER]/
event_dir = ./data/[EVENTID]/
analysis_db = ./data/analysisdb.sqlite
"""

# ----------------------------------------------------------------------
class DownloaderApp(object):
    """
    Download data (DM logs, ShakeMaps, ComCat events) needed for analysis of ShakeAlert performance.
    """
    
    def __init__(self):
        """Constructor.
        """
        self.config = None
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

        self.db = analysisdb.AnalysisData(self.config.get("files", "analysis_db"))

        if args.db_init or args.all:
            self._db_init(args.db_init)

        if args.db_summary or args.all:
            self._db_summary(args.db_summary)

        if args.db_populate or args.all:
            if "eew_alerts" in args.db_populate or args.db_populate == "all" or args.all:
                self._db_populate_eewalerts(all=args.db_populate != "new_eew_alerts", replace=args.db_replace_rows)
            if args.db_populate == "comcat_events" or args.db_populate == "all" or args.all:
                self._db_populate_events(replace=args.db_replace_rows)
            if args.db_populate == "comcat_shakemaps" or args.db_populate == "all" or args.all:
                self._db_populate_shakemaps(replace=args.db_replace_rows)

        if args.show_matches or args.all:
            self._show_matches()
        return

    def _fetch_eewalerts(self, dateBegin, dateEnd):
        if self.showProgress:
            print("Fetching EEW alerts...")

        self.eewserver = shakealert.EEWServer(self.config)
        self.eewserver.login()
        
        dmlog = shakealert.DMLogXML(config=self.config)
        server = self.config.get("shakealert.production", "server")    
        logsDir = self.config.get("files", "dmlogs_dir").replace("[SERVER]", server)
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

        event = comcat.DetailEvent()
        dirTemplate = self.config.get("files", "event_dir")
        for eqId in self.config.options("events"):
            dataDir = dirTemplate.replace("[EVENTID]", eqId)
            event.fetch(eqId, dataDir)
        return

    def _fetch_shakemaps(self):
        """Fetch geojson event file from USGS ComCat using web services.
        """
        if self.showProgress:
            print("Fetching ShakeMaps...")

        event = comcat.DetailEvent()
        dirTemplate = self.config.get("files", "event_dir")
        for eqId in self.config.options("events"):
            dataDir = dirTemplate.replace("[EVENTID]", eqId)
            event.load(os.path.join(dataDir, eqId+".geojson"))
            shakemap = event.get_product("shakemap")
            if len(shakemap) > 1:
                raise ValueError("Expected to get preferred ShakeMap.")
            shakemap[0].fetch("grid.xml", dataDir)
            shakemap[0].fetch("info.json", dataDir)
            if not os.path.isfile(os.path.join(dataDir, "info.json.gz")):
                shakemap[0].fetch("info.xml", dataDir)
                if not os.path.isfile(os.path.join(dataDir, "info.xml.gz")):
                    logging.getLogger(__name__).error("Could not retrieve ShakeMap info JSON or XML file for {}.".format(eqId))
                else:
                    logging.getLogger(__name__).info("Performing minimal conversion of ShakeMap info XML file to JSON for {}.".format(eqId))
                    #self._extract_shakemap_info(dataDir)
                
        return


    def _extract_shakemap_info(self, dataDir):
        """Perform minimal conversion of ShakeMap info XML file to JSON.

        Get MMI bias.
        """
        with gzip.open(os.path.join(dataDir, "info.xml.gz"), "r") as fxml:

            bytes = fxml.read()
            elRoot = etree.fromstring(bytes) # info

            # Values
            mmiBias = float(elRoot.xpath("tag[@name='mi_bias']")[0].get("value"))
            pgmBias = tuple(map(float, elRoot.xpath("tag[@name='bias']")[0].get("value").split()))
            (pgaBias, pgvBias, psa03Bias, psa10Bias, psa30Bias) = pgmBias

            mmiMax = float(elRoot.xpath("tag[@name='mi_max']")[0].get("value"))
            pgvMax = float(elRoot.xpath("tag[@name='pgv_max']")[0].get("value"))
            pgaMax = float(elRoot.xpath("tag[@name='pga_max']")[0].get("value"))
            psa03Max = float(elRoot.xpath("tag[@name='psa03_max']")[0].get("value"))
            psa10Max = float(elRoot.xpath("tag[@name='psa10_max']")[0].get("value"))
            psa30Max = float(elRoot.xpath("tag[@name='psa30_max']")[0].get("value"))
    
            gmpeStr = elRoot.xpath("tag[@name='GMPE']")[0].get("value")
            pgm2miStr = elRoot.xpath("tag[@name='pgm2mi']")[0].get("value")
            shakemapVer = elRoot.xpath("tag[@name='ShakeMap revision']")[0].get("value")

            data = {
                "output": {
                    "ground_motions": {
                        "intensity": {
                            "bias": mmiBias,
                            "max": mmiMax,
                            "units": "intensity",
                        },
                        "pga": {
                            "bias": pgaBias,
                            "max": pgaMax,
                            "units": "%g",
                        },
                        "pgv": {
                            "bias": pgvBias,
                            "max": pgvMax,
                            "units": "cm/s",
                        },
                        "psa03": {
                            "bias": psa03Bias,
                            "max": psa03Max,
                            "units": "%g",
                        },
                        "psa10": {
                            "bias": psa10Bias,
                            "max": psa10Max,
                            "units": "%g",
                        },
                        "psa30": {
                            "bias": psa30Bias,
                            "max": psa30Max,
                            "units": "%g",
                        },
                    },
                },
                "processing": {
                    "ground_motion_modules": {
                        "gmpe": {
                            "module": gmpeStr,
                            "reference": gmpeStr,
                        },
                        "pgm2mi": {
                            "module": pgm2miStr,
                            "reference": pgm2miStr,
                        },
                    },
                    "shakemap_versions": {
                        "shakemap_revision": shakemapVer,
                    },
                },
            }
            with gzip.open(os.path.join(dataDir, "info.json.gz"), "w") as fjson:
                import json
                json.dump(data, fjson)
        return
    
    def _db_init(self, tables):
        """Create analysis database with ShakeAlert DM alerts and ComCat events.
        """
        if self.showProgress:
            print("Setting up analysis database...")
            
        self.db.init(tables)
        logging.getLogger(__name__).info(self.db.tables_info())
        return

    def _db_summary(self, style="tables_info"):
        """Show summary of analysis database contents.
        """
        if style == "tables_info":
            print(self.db.tables_info())
        else:
            print(self.db.summary())
        return

    def _db_populate_events(self, replace=False):
        if self.showProgress:
            print("Updating ComCat events in analysis database...")

        event = comcat.DetailEvent()
        dirTemplate = self.config.get("files", "event_dir")
        events = self.config.options("events")
        numEvents = len(events)
        for iEvent,eqId in enumerate(events):
            if self.showProgress:
                sys.stdout.write("\rProcessing ComCat events...{:d}%%".format(((iEvent+1)*100)/numEvents))
                sys.stdout.flush()

            dataDir = dirTemplate.replace("[EVENTID]", eqId)
            event.load(os.path.join(dataDir, eqId+".geojson"))
            self.db.add_event(event, replace)
        if self.showProgress:
            sys.stdout.write("\n")
        return
    
    def _db_populate_shakemaps(self, replace=False):
        import json
        if self.showProgress:
            print("Updating ComCat ShakeMap info in analysis database...")

        dirTemplate = self.config.get("files", "event_dir")
        events = self.config.options("events")
        numEvents = len(events)
        for iEvent,eqId in enumerate(events):
            if self.showProgress:
                sys.stdout.write("\rProcessing ComCat events...{:d}%%".format(((iEvent+1)*100)/numEvents))
                sys.stdout.flush()

            dataDir = dirTemplate.replace("[EVENTID]", eqId)
            with gzip.open(os.path.join(dataDir, "info.json.gz"), "r") as fh:
                info = json.load(fh)
                info["event_id"] = eqId
                self.db.add_shakemap_info(info, replace)
        if self.showProgress:
            sys.stdout.write("\n")
        return
    
    def _db_populate_eewalerts(self, all=False, replace=False):
        """
        """
        import glob
        import re
        
        if self.showProgress:
            print("Updating ShakeAlert DM alerts in analysis database...")

        server = self.config.get("shakealert.production", "server")    
        logsDir = self.config.get("files", "dmlogs_dir").replace("[SERVER]", server)
        files = sorted(glob.glob(os.path.join(logsDir, "dmevent_*.log.gz")))
        if not all:
            # Get most recent entry and use files starting on that day.
            alert = self.db.most_recent_alert(server)
            dateMostRecent = dateutil.parser.parse(alert["timestamp"]).date()
            pattern = [
                "dmevent_",
                "([0-9]{4})",
                "([0-9]{2})",
                "([0-9]{2})",
                ".log*",
                ]
            remove = []
            for filename in files:
                year, month, day = map(int, re.search("".join(pattern), filename).groups())
                fileDate = datetime.date(year=year, month=month, day=day)
                if fileDate < dateMostRecent:
                    remove.append(filename)
            for filename in remove:
                files.remove(filename)

        # Read DM logs
        dmlog = shakealert.DMLogXML(config=self.config)
        numFiles = len(files)
        logging.getLogger(__name__).info("Processing {:d} DM logs starting with {:s}.".format(numFiles, files[0]))
        for iFile,filename in enumerate(files):
            if self.showProgress:
                sys.stdout.write("\rProcessing DM logs...{:d}%%".format(((iFile+1)*100)/numFiles))
                sys.stdout.flush()
            alerts = dmlog.load(filename)
            self.db.add_alerts(alerts, replace)
        if self.showProgress:
            sys.stdout.write("\n")
        return

    def _show_matches(self):
        if self.showProgress:
            print("Showing matches between ComCat and ShakeAlert...")

        self.db = analysisdb.AnalysisData(self.config.get("files", "analysis_db"))
        self.db.show_matches(self.config.get("shakealert.production", "server"))
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

        self.config = config
        return
    
    def show_parameters(self):
        """Write parameters to stdout.
        """
        import sys
        self.config.write(sys.stdout)
        return

    def _parseCommandLine(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config")
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--fetch-eewalerts", action="store", dest="fetch_eewalerts", default=None, metavar="DATE_BEGIN,DATE_END")
        parser.add_argument("--fetch-events", action="store_true", dest="fetch_events")
        parser.add_argument("--fetch-shakemaps", action="store_true", dest="fetch_shakemaps")
        parser.add_argument("--db-init", action="store", dest="db_init", choices=["eew_alerts", "comcat_events", "comcat_shakemaps", "performance", "all"])
        parser.add_argument("--db-summary", action="store", dest="db_summary", default=None, choices=[None, "tables_only", "summary"])
        parser.add_argument("--db-populate", action="store", dest="db_populate", choices=["all_eew_alerts", "new_eew_alerts", "comcat_events", "comcat_shakemaps", "all"])
        parser.add_argument("--db-replace-rows", action="store_true", dest="db_replace_rows")
        parser.add_argument("--show-matches", action="store_true", dest="show_matches")
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug", default=True)
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    DownloaderApp().main()


# End of file
