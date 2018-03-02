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
analysis_db = ./data/analysisdb_NEW.sqlite
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

        if args.db_status or args.all:
            self._db_status()

        if args.db_populate or args.all:
            if "eew_alerts" in args.db_populate or args.db_populate == "all" or args.all:
                self._db_populate_eewalerts(all=args.db_populate != "new_eew_alerts", replace=args.db_replace_rows)
            if args.db_populate == "comcat_events" or args.db_populate == "all" or args.all:
                self._db_populate_events(replace=args.db_replace_rows)

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
        # XML
        # <tag name="mi_bias" value="-0.31" desc="magnitude bias for Intensity" />
        # <tag name="bias" value="-0.04 -0.14 -0.30 -0.27 -0.45 " desc="magnitude bias (pga pgv psa03 psa10 psa30 )" />
        #
        # <tag name="pga_max" value="83.83" desc="Max value of grid" />
        # <tag name="pgv_max" value="52.53" desc="Max value of grid" />
        # <tag name="mi_max" value="7.89" desc="Max value of grid" />
        # <tag name="psa03_max" value="112.13" desc="Max value of grid" />
        # <tag name="psa10_max" value="54.42" desc="Max value of grid" />
        # <tag name="psa30_max" value="5.72" desc="Max value of grid" />
        #
        # <tag name="GMPE" value="GMPE::BA08" desc="GMPE type" />
        # <tag name="pgm2mi" value="GMICE::Wald99 - Wald, et al.; 1999" desc="Intensity Function" />
        #
        # <tag name="ShakeMap revision" value="3.5.687" desc="ShakeMap source code revision number" />
        #
        #
        # JSON
        # output
        #    ground_motions
        #      intensity
        #          bias
        #          max
        #          units: "intensity"
        #      pga
        #          bias
        #          max
        #          units: "%g"
        #      pgv
        #          bias
        #          max
        #          units: "cm/s"
        #      psa03
        #          bias
        #          max
        #          units: "%g"
        #      psa10
        #          bias
        #          max
        #          units: "%g"
        #      pds30
        #          bias
        #          max
        #          units: "%g"
        # processing
        #     ground_motion_modules
        #         gmpe
        #             module
        #             reference
        #         pgm2mi
        #             module
        #             reference
        #     shakemap_versions
        #         shakemap_revision
        with gzip.open(os.path.join(dataDir, "info.xml.gz"), "r") as fxml:
            import pdb
            pdb.set_trace()
            
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
        logging.getLogger(__name__).info(db.summary())
        return

    def _db_status(self):
        """Show summary of analysis database contents.
        """
        print(db.summary())
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
        parser.add_argument("--db-init", action="store", dest="db_init", choices=["eew_alerts", "comcat_events", "performance", "all"])
        parser.add_argument("--db-status", action="store_true", dest="db_status")
        parser.add_argument("--db-populate", action="store", dest="db_populate", choices=["all_eew_alerts", "new_eew_alerts", "comcat_events", "all"])
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
