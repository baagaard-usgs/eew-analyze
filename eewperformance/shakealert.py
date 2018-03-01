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
import datetime
import os
import gzip
import logging
import dateutil
from lxml import etree

TIMEOUT_SECS = 30 # How many seconds to wait for download
DEMONSTRATION_BEGIN = datetime.date(year=2012, month=1, day=27)
DEMONSTRATION_END = datetime.date(year=2016, month=6, day=4)
PRODUCTION_BEGIN = datetime.date(year=2016, month=6, day=5)
PRODUCTION_END = datetime.date(year=3000, month=1, day=1)

class EEWServer(object):
    """EEW ShakeAlert server holding status and logs.
    """

    def __init__(self, config):
        """Constructor with configuraton.

        :type config: ConfigParser
        :param config: Configuration for application.

        """
        self.config = config
        self.connectionProd = None
        self.connectionDemo = None
        return

    def login(self):
        """Log in to server and open connection.
        """
        self._loginProd()
        self._loginDemo()
        return

    def logout(self):
        """Close connection to server.
        """
        if self.connectionProd:
            self.connectionProd.close()
        if self.connectionDemo:
            self.connectionDemo.close()
        return
        

    def _loginDemo(self):
        """Log in to server with DM log for demonstration system.
        """
        urlTemplate = self.config.get("shakealert.demonstration", "login_url")
        name = self.config.get("shakealert.demonstration", "server")
        username = self.config.get("shakealert.demonstration", "username")
        password = self.config.get("shakealert.demonstration", "password")
        urlCookie = urlTemplate.replace("[SERVER]", name)
        urlLogin = urlCookie + "/j_security_check"

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
                raise Exception("Could not connect to server - %s." % urlCookie)

        self.connectionDemo = connection
        return

    def _loginProd(self):
        """Log in to server with DM log for production system.
        """
        urlTemplate = self.config.get("shakealert.production", "login_url")
        name = self.config.get("shakealert.production", "server")
        url = urlTemplate.replace("[SERVER]", name)

        try:
            connection = requests.session()
            connection.headers["User-Agent"] = "Mozilla/5.0"
            response = connection.post(url, timeout=TIMEOUT_SECS)
            response.raise_for_status()
        except requests.exceptions.RequestException as htpe:
            try:
                response = connection.post(url, timeout=TIMEOUT_SECS)
                response.raise_for_status()
            except requests.exceptions.RequestException as msg:
                raise Exception("Could not connect to server - %s." % url)

        self.connectionProd = connection
        return


class DMLogXML(object):
    """XML DM log fetcher and parser.
    """

    def __init__(self, config=None, filename=None):
        """Constructor.

        :type filename: str
        :param filename: Name of file for local storage of DM log.
        """
        self.config = config
        if filename:
            self.load(filename)
        return

    def fetch(self, server, date, logsDir):
        """Fetch DM XML log from ShakeAlert system.
        
        :type server: EEWServer
        :param server: EEW ShakeAlert server.

        :type date: datetime.date
        :param event: Date of log to download.

        :type logsDir: str
        :param logsDir: Directory in which to store file locally.
        """
        tstamp = "%d%02d%02d" % (date.year, date.month, date.day,)
        if date >= PRODUCTION_BEGIN and date <= PRODUCTION_END:
            urlTemplate = self.config.get("shakealert.production", "log_url")
            logServer = self.config.get("shakealert.production", "server")
            connection = server.connectionProd
        elif date >= DEMONSTRATION_BEGIN and date <= DEMONSTRATION_END:
            urlTemplate = self.config.get("shakealert.demonstration", "log_url")
            logServer = self.config.get("shakealert.demonstration", "server")
            connection = server.connectionDemo
        else:
            raise ValueError("DM XML logs not available for event on %d-%02d-%02d." % (date.year, date.month, date.day,))

        url = urlTemplate.replace("[SERVER]", logServer).replace("[YYYYMMDD]", tstamp)
        filename = os.path.join(logsDir, os.path.split(url)[1])
        try:
            response = connection.get(url, timeout=TIMEOUT_SECS)
            response.raise_for_status()
        except requests.exceptions.RequestException as htpe:
            try:
                response = connection.get(url, timeout=TIMEOUT_SECS)
                response.raise_for_status()
            except requests.exceptions.RequestException as htpe:
                logging.getLogger(__name__).info("Could not download log %s. Log may not exist." % url)
                return
                
        suffix = ""
        if not filename.endswith(".gz"):
            suffix = ".gz"
        with gzip.open(filename+suffix, "w") as fh:
            fh.write(response.text.decode("utf-8"))
        return

    def load(self, filename):
        """Load DM log file.

        :type filename: str
        :param filename: Name of local file with DM log.
        """
        suffix = ""
        if not filename.endswith(".gz"):
            suffix = ".gz"
        with gzip.open(filename+suffix, "r") as fh:
            alerts = self._parse(fh)
        return alerts

    def _parse(self, fh):
        """Parse DM XML log file.
        
        :type lines: str
        :param lines: Lines from DM log file.
        :returns: Numpy structured array with DM log contents.
        """
        bytes = fh.read()
        fh.close()
        pattern = [
            "(?P<timestamp>[0-9]{4}-[0-9]+-[0-9]+T[0-9]{2}:[0-9]{2}:[0-9]{2}.[0-9]+Z)[\s]*",
            "(?P<xml>\<\?xml[\s\S]+?</event_message>)"
        ]
        xmlBlocks = re.findall("".join(pattern), bytes)
        alerts = []
        for timestamp,xmlBlock in xmlBlocks:
            elMsg = etree.fromstring(xmlBlock.replace("UTF-16", "utf-8"))
            elCoreInfo = self._getChild(elMsg, "core_info")
            elMag = self._getChild(elCoreInfo, "mag")

            # Use time step in message, fallback to time stamp in log.
            if elMsg.get("timestamp"):
                timestamp = elMsg.get("timestamp")
            if elCoreInfo.get("id").lower().startswith("test"):
                continue
            try:
                alerts.append({
                    "category": elMsg.get("category"),
                    "server": elMsg.get("instance").replace("dm@","") if elMsg.get("instance") else "unknown",
                    "message_type": elMsg.get("message_type"),
                    "timestamp": dateutil.parser.parse(timestamp),
                    "version": int(elMsg.get("version")),
                    
                    "event_id": int(elCoreInfo.get("id")),
                    
                    "magnitude": float(elMag.text),
                    "magnitude_type": elMag.get("units"),
                    
                    "latitude": self._getChildValue(elCoreInfo, "lat", vtype=float),
                    "longitude": self._getChildValue(elCoreInfo, "lon", vtype=float),
                    "depth_km": self._getChildValue(elCoreInfo, "depth", vtype=float),
                    "origin_time": dateutil.parser.parse(self._getChildValue(elCoreInfo, "orig_time")),
                    "likelihood": self._getChildValue(elCoreInfo, "likelihood", vtype=float, default=-1),
                    "num_stations": self._getChildValue(elCoreInfo, "num_stations", vtype=int, default=-1),
                })
            except Exception as ex:
                import pdb
                pdb.set_trace()
        return alerts
        

    def _getChild(self, el, name):
        """Get child element 'name'. Raises IOError if more than one child element with name is found.

        :type el: XML element
        :param el: XML element

        :type name: str
        :param name: Name of child element.
        """
        elList = el.xpath(name)
        if len(elList) != 1:
            raise IOError("Expect one '{0}' child element in XML element. Found {1} elements in {3}.".format(name, len(elList), el.tag))
        return elList[0]

    def _getChildValue(self, el, name, vtype=str, default=None):
        """Get text value of child element 'name'. Returns default value if exactly one child is not found.

        :type el: XML element
        :param el: XML element

        :type name: str
        :param name: Name of child element.
        """
        elList = el.xpath(name)
        if len(elList) != 1:
            return default
        value =  elList[0].text
        if vtype == str:
            return value
        value = float(value)
        if vtype == int:
            value = int(value)
        return value


    
class DMLogASCII(object):
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
        URL_TEMPLATE = server.productionURL + "/[SERVER]/dmreview/dmlogs/dm_[YEARMMDD].txt"
        
        t = event.time
        tstamp = "%d%02d%02d" % (t.year, t.month, t.day,)
        if int(tstamp) >= DEMONSTRAION_END:
            serverName = server.name
            connection = server.connection
        else:
            serverName = DEMONSTRATION_SERVER
            connection = server.connectionDemo
        url = URL_TEMPLATE.replace("[SERVER]", serverName).replace("[YEARMMDD]", tstamp)
        try:
            response = connection.get(url, timeout=TIMEOUT_SECS)
            response.raise_for_status()
        except requests.exceptions.RequestException as htpe:
            try:
                response = connection.get(url, timeout=TIMEOUT_SECS)
                response.raise_for_status()
            except requests.exceptions.RequestException as htpe:
                print("Could not get log %s." % url)

        lines = response.text.decode("utf-8").split("\n")
        self._fix_timestamp(lines, event)
        buffer = "\n".join(lines)
        with open(filename, "w") as fh:
            fh.write(buffer)
        import StringIO
        self.data = self._parse(StringIO.StringIO(buffer))
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
        maskDM = numpy.bitwise_and(self.data["system"] == "dm", self.data["action"] == "Created")
        maskDM = numpy.bitwise_and(maskDM, self.data["type"] == "new")
        dataDM = self.data[maskDM]
        
        originDiff = numpy.abs(dataDM["origin_time"] - numpy.datetime64(event.time))
        latitudeDiff = numpy.abs(dataDM["latitude"] - event.latitude)
        longitudeDiff = numpy.abs(dataDM["longitude"] - event.longitude)
        indicesSorted = numpy.argsort(originDiff)

        alertIndex = None
        for index in indicesSorted:
            if originDiff[index] < numpy.timedelta64(1, "s") and latitudeDiff[index] < 0.3 and longitudeDiff[index] < 0.3:
                alertIndex = index
                break
        if not alertIndex:
            # Add closest match to log
            return None

        alertId = dataDM["id"][alertIndex]
        maskDM = numpy.bitwise_and(self.data["system"] == "dm", self.data["id"] == alertId)
        maskDM = numpy.bitwise_and(maskDM, self.data["version"] >= 0)
        maskDM = numpy.bitwise_and(maskDM, self.data["action"] == "Published")
        alerts = self.data[maskDM]
        alerts.sort(order="version")
        return alerts
    
    def _fix_timestamp(self, lines, event):
        """Fix time stamp in DM log files. Add date to time stamp and
        change ':' between seconds and milliseconds to '.'.

        :type lines: List of str
        :param lines: List of lines in file.
        """
        t = event.time
        patternTimeStamp = "^(?P<timestamp>[0-9]{2}:[0-9]{2}:[0-9]{2}:[0-9]{3})\|"
        patternOriginTime = "(?P<origin_time>[0-9]{2}-[0-9]+-[0-9]+ [0-9]{2}:[0-9]{2}:[0-9]{2}.[0-9]{3})"
        tstampRe = re.compile(patternTimeStamp)
        originTRe = re.compile(patternOriginTime)
        for index,line in enumerate(lines):
            matchTimeStamp = tstampRe.match(line)
            matchOriginTime = originTRe.search(line)
            lineNew = line
            if matchTimeStamp:
                told = matchTimeStamp.group("timestamp")
                rindex = told.rfind(":")
                tnew = told[:rindex] + "." + told[rindex+1:]
                lineNew = lineNew.replace(told, "%d-%02d-%02d %s" % (t.year, t.month, t.day, tnew))
                lines[index] = lineNew
            if matchOriginTime:
                told = matchOriginTime.group("origin_time")
                tnew = "20"+told
                lineNew = lineNew.replace(told, tnew)
                lines[index] = lineNew
        return
                
    def _parse(self, fh):
        """Parse DM log file.
        
        :type lines: str
        :param lines: Lines from DM log file.
        :returns: Numpy structured array with DM log contents.
        """
        pattern = [
            "(?P<timestamp>[0-9]{4}-[0-9]+-[0-9]+ [0-9]{2}:[0-9]{2}:[0-9]{2}.[0-9]{3})\|",
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
            "[ ]+(?P<origin_time>[0-9]{4}-[0-9]+-[0-9]+ [0-9]{2}:[0-9]{2}:[0-9]{2}.[0-9]{3})",
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
