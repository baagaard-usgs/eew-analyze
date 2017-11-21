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
import datetime
import requests
import gzip

DEFAULTS = """
[shakealert]
#servers = [eew-bk-prod1,eew2demo]
servers = [eew-bk-prod1]
# username = USERNAME
# password = PASSWORD

[eew2demo]
date_start = 2012-01-01
date_end = 2016-01-31

[eew-bk-prod1]
#date_start = 2016-02-01
date_start = 2017-11-01
date_end = 2017-11-17
"""

TIMEOUT_SECS = 30

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
class DownloadApp(object):
    """
    Download DM logs from ShakeAlert servers.
    """
    
    def __init__(self):
        """Constructor.
        """
        self.params = None
        self.connection = None
        return

    def main(self):
        """Main entry point
        """
        # Initialization
        args = self._parseCommandLine()
        logLevel = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=logLevel, filename="downlaod_dmlogs.log")
        if args.show_progress:
            self.showProgress = True
        self.initialize(args.config)

        # Show parameters
        if args.show_parameters or args.all:
            self.show_parameters()

        # Download logs
        if args.fetch or args.all:
            self.fetch()

        return

    def fetch(self):
        """Download DM logs from ShakeAlert server.
        """

        servers = _config_get_list(self.params.get("shakealert", "servers"))
        for server in servers:
            self.logDir = os.path.join("data", "dmlogs", server)
            if not os.path.isdir(self.logDir):
                os.makedirs(self.logDir)

            connection = self.login(server)
            startYear, startMonth, startDay = map(int, self.params.get(server, "date_start").split("-"))
            endYear, endMonth, endDay = map(int, self.params.get(server, "date_end").split("-"))
            dateStart = datetime.date(year=startYear, month=startMonth, day=startDay)
            dateEnd = datetime.date(year=endYear, month=endMonth, day=endDay)
            date = dateStart
            dt = datetime.timedelta(days=1)
            while date <= dateEnd:
                self.fetch_log(connection, server, date)
                date += dt
            connection.close()
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

    def login(self, server):
        """Log in to ShakeAlert server.
        """
        URL_TEMPLATE = "https://eew.geo.berkeley.edu:8443/[SERVER]/dmreview/"
        urlCookie = URL_TEMPLATE.replace("[SERVER]", server)
        urlLogin = urlCookie + "j_security_check"

        connection = requests.session()
        connection.headers["User-Agent"] = "Mozilla/5.0"
        response = connection.post(urlCookie, timeout=TIMEOUT_SECS)
        response.raise_for_status()
        
        payload = {
            "j_username": self.params.get("shakealert", "username"),
            "j_password": self.params.get("shakealert", "password"),
        }
        response = connection.post(urlLogin, cookies=response.cookies, data=payload, timeout=TIMEOUT_SECS)
        response.raise_for_status()
        return connection
    
    def fetch_log(self, connection, server, date):
        """Fetch DM log from ShakeAlert server.
        
        :type date: datetime.date
        :param date: Date of log to download in the form YYYYMM.
        """
        URL_TEMPLATE = "https://eew.geo.berkeley.edu:8443/[SERVER]/dmreview/dmlogs/dm_[YYYYMMDD].txt"
        tstamp = "%d%02d%02d" % (date.year, date.month, date.day,)

        if self.showProgress:
            print("Fetching DM log 'dm_%s.txt'..." % tstamp)

        url = URL_TEMPLATE.replace("[SERVER]", server).replace("[YYYYMMDD]", tstamp)
        filename = os.path.join(self.logDir, "dm_%s.txt.gz" % tstamp)
        
        try:
            response = connection.get(url, timeout=TIMEOUT_SECS)
            response.raise_for_status()
            data = response.text.decode("utf-8")
        except requests.exceptions.RequestException as htpe:
            try:
                response = connection.get(url, timeout=TIMEOUT_SECS)
                response.raise_for_status()
                data = response.text.decode("utf-8")
            except requests.exceptions.RequestException as msg:
                raise Exception("Could not connect to server - %s." % url)

        with gzip.open(filename, "w") as fh:
            fh.write(data)
        return
    
    def _parseCommandLine(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config", required=True)
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--fetch", action="store_true", dest="fetch", default=False)
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug")
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    DownloadApp().main()


# End of file
