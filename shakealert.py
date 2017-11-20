# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import requests
import re
import numpy

TIMEOUT_SECS = 30 # How many seconds to wait for download

class EEWServer(object):
    """EEW ShakeAlert server holding status and logs.
    """

    def __init__(self, name, baseURL="https://eew.geo.berkeley.edu:8443"):
        """Constructor with name and URL.

        :type name: str
        :param name: Name of EEW ShakeAlert server.

        :type baseURL: str
        :param baseURL: Base URL where server data is stored.
        """
        self.name = name
        self.baseURL = baseURL

        self.connection = None
        return

    def login(self, username, password):
        """Log in to server and open connection.
        """
        URL_TEMPLATE = self.baseURL + "/[SERVER]/dmreview/"
        urlCookie = URL_TEMPLATE.replace("[SERVER]", self.name)
        urlLogin = urlCookie + "j_security_check"

        payload = {
            "j_username": username,
            "j_password": password,
            }
        
        try:
            connection = requests.session()
            connection.headers["User-Agent"] = "Mozilla/5.0"
            response = connection.post(urlCookie, timeout=TIMEOUT_SECS)
            response.raise_for_status()
        
            response = connection.post(urlLogin, cookies=response.cookies, data=payload, timeout=TIMEOUT_SECS)
            response.raise_for_status()
        except requests.exceptions.RequestException as htpe:
            try:
                response = connection.post(urlCookie, timeout=TIMEOUT_SECS)
                response.raise_for_status()
        
                response = connection.post(urlLogin, cookies=response.cookies, data=payload, timeout=TIMEOUT_SECS)
                response.raise_for_status()
            except requests.exceptions.RequestException as msg:
                raise Exception("Could not connect to server - %s." % self.baseURL)

        self.connection = connection
        return

    def logout(self):
        """Close connection to server.
        """
        if self.connection:
            self.connection.close()
        return
        
    
class DMLog(object):
    """DM log fetcher and parser.
    """

    def __init__(self, filename=None):
        """Constructor.

        :type filename: str
        :param filename: Name of file for local storage of DM log.
        """
        if filename:
            self.load(filename)
        return

    def fetch(self, server, event, filename):
        """Fetch DM log from ShakeAlert system.
        
        :type server: EEWServer
        :param server: EEW ShakeAlert server.
        :type event: DetailEvent
        :param event: ComCat detailed event.
        :type filename: str
        :param filename: Name of file for locally storing DM log.
        """
        URL_TEMPLATE = server.baseURL + "/[SERVER]/dmreview/dmlogs/dm_[YEARMMDD].txt"
        DEMONSTRATION_SERVER = "eew2demo"
        
        t = event.time
        tstamp = "%d%02d%02d" % (t.year, t.month, t.day,)
        serverName = server.name if tstamp >= 20160201 else DEMONSTRATION_SERVER
        url = URL_TEMPLATE.replace("[SERVER]", serverName).replace("[YEARMMDD]", tstamp)

        try:
            response = server.connection.get(url, timeout=TIMEOUT_SECS)
            response.raise_for_status()
        except requests.exceptions.RequestException as htpe:
            try:
                response = connection.get(url, timeout=TIMEOUT_SECS)
                response.raise_for_status()
            except requests.exceptions.RequestException as htpe:
                raise Exception("Could not connect to server - %s." % self.baseURL)

        lines = response.text.decode("utf-8").split("\n")
        self._fix_timestamp(lines)
        buffer = "\n".join(lines)
        with open(filename, "w") as fh:
            fh.write(buffer)
        self.data = self._parse(buffer)
        return

    def load(self, filename):
        """Load DM log file.

        :type filename: str
        :param filename: Name of local file with DM log.
        """
        with open(filename, "r") as fh:
            self.data = self._parse(fh)
        return

    def alerts(self, event):
        """Extract alerts for event from DM log.

        :type event: DetailEvent
        :param event: ComCat detailed event.
        """
        originDiff = numpy.abs(self.data["origin_time"] - numpy.datetime64(event.time))
        indicesSorted = numpy.argsort(originDiff)
        latitudeDiff = numpy.abs(self.data["latitude"] - event.latitude)
        longitudeDiff = numpy.abs(self.data["longitude"] - event.longitude)

        import pdb
        pdb.set_trace()
        
        return
    
    def _fix_timestamp(self, lines):
        """Fix time stamp in DM log files. Add date to time stamp and
        change ':' between seconds and milliseconds to '.'.

        :type lines: List of str
        :param lines: List of lines in file.
        """
        pattern = "^(?P<timestamp>[0-9]{2}:[0-9]{2}:[0-9]{2}:[0-9]{3})\|"
        tstampRe = re.compile(pattern)
        linesNew = []
        for index,line in enumerate(lines):
            m = tstampRe.match(line)
            if m:
                told = m.group("timestamp")
                rindex = told.rfind(":")
                tnew = told[:rindex] + "." + told[rindex+1:]
                lines[index] = line.replace(told, "%d-%02d-%02d %s" % (t.year, t.month, t.day, tnew))
        return
                
    def _parse(self, fh):
        """Parse DM log file.
        
        :type lines: str
        :param lines: Lines from DM log file.
        :returns: Numpy structured array with DM log contents.
        """
        pattern = [
            "(?P<timestamp>[0-9]{2}-[0-9]+-[0-9]+ [0-9]{2}:[0-9]{2}:[0-9]{2}.[0-9]{3})\|",
            "[ ]*(?P<action>[a-zA-Z ]+):",
            "[ ]*(?P<system>[a-zA-Z]+)",
            "[ ]+(?P<id>[0-9]+)",
            "[ ]+(?P<version>[\-0-9]+)",
            "[ ]+(?P<category>[a-zA-Z]+)",
            "[ ]+(?P<type>[a-zA-Z]+)",
            "[ ]+(?P<magnitude>[0-9]+\.[0-9]+)",
            "[ ]+(?P<magnitude_uncertainty>[-0-9]+\.[0-9]+)",
            "[ ]+(?P<latitude>[\-0-9]+\.[0-9]+)",
            "[ ]+(?P<latitude_uncertainty>[-0-9]+\.[0-9]+)",
            "[ ]+(?P<longitude>[\-0-9]+\.[0-9]+)",
            "[ ]+(?P<longitude_uncertainty>[0-9]+\.[0-9]+)",
            "[ ]+(?P<depth>[0-9]+\.[0-9]+)",
            "[ ]+(?P<depth_uncertainty>[0-9]+\.[0-9]+)",
            "[ ]+(?P<origin_time>[0-9]{2}-[0-9]+-[0-9]+ [0-9]{2}:[0-9]{2}:[0-9]{2}.[0-9]{3})",
            ]
        dtype = [
            ("timestamp", "datetime64[ms]",),
            ("action", "S10",),
            ("system", "S7",),
            ("id", "int64",),
            ("version", "int32",),
            ("category", "S5",),
            ("type", "S7",),
            ("magnitude", "float32",),
            ("magnitude_uncertainty", "float32",),
            ("latitude", "float32",),
            ("latitude_uncertainty", "float32",),
            ("longitude", "float32",),
            ("longitude_uncertainty", "float32",),
            ("depth", "float32",),
            ("depth_uncertainty", "float32",),
            ("origin_time", "datetime64[ms]",),
            ]

        return numpy.fromregex(fh, "".join(pattern), dtype=dtype)
    
# End of file
