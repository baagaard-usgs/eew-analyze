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
import logging


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
            import dateutil.parser
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

        if args.updatedb or args.all:
            self._updatedb()
        return

    def _fetch_eewalerts(self, dateBegin, dateEnd):
        if self.showProgress:
            print("Fetching EEW alerts...")
            
        from shakealert import EEWServer, DMLogXML
        import datetime

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
        from comcat import DetailEvent

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

        shakemap = self.event.get_product("shakemap")
        if len(shakemap) > 1:
            raise ValueError("Expected to get preferred ShakeMap.")
        filename = _data_filename(self.params, "shakemap", self.event.id)
        shakemap[0].fetch("shape.zip", filename)
        return

    def _initdb(self):
        from analysisdb import AnalysisData

        if self.showProgress:
            print("Setting up analysis database...")
            
        db = AnalysisData()
        db.init()
        return

    def _updatedb(self):
        from analysisdb import AnalysisData

        if self.showProgress:
            print("Updating analysis database...")

        db = AnalysisData()
        db.open()
        db.update()
        db.close()
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
        import sys
        self.params.write(sys.stdout)
        return

    def _parseCommandLine(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config", required=True)
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--fetch-eewalerts", action="store", dest="fetch_eewalerts", default=None)
        parser.add_argument("--fetch-events", action="store_true", dest="fetch_events")
        parser.add_argument("--fetch-shakemaps", action="store_true", dest="fetch_shakemaps")
        parser.add_argument("--initdb", action="store_true", dest="initdb")
        parser.add_argument("--updatedb", action="store_true", dest="updatedb")
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug")
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    DownloaderApp().main()


# End of file
