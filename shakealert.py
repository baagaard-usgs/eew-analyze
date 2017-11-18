# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import urllib2

TIMEOUT_SECS = 30 # How many seconds to wait for download
WAIT_SECS = 3

class DMLog(object):
    """DM log parser.
    """

    def __init__(self, filename=None):
        """Constructor.

        :type filename: str
        :param filename: Name of file for local storage of DM log.
        """
        self.alerts = None
        
        if filename:
            self.load(filename)
        return

    def fetch(self, event, server, filename):
        """Fetch DM log from ShakeAlert production system.
        
        :type event: DetailEvent
        :param event: ComCat detailed event.
        :type server: str
        :param server: Name of EEW production server with DM log.
        :type filename: str
        :param filename: Name of file for locally storing DM log.
        """
        
        URL_TEMPLATE = "https://eew.geo.berkeley.edu:8443/[SERVER]/dmreview/dmlogs/dm_[YEARMMDD].txt"

        t = event.time
        tstamp = "%d%02d%02d" % (t.year, t.month, t.day,)
        url = URL_TEMPLATE.replace("[SERVER]", server).replace("[YEARMMDD]", tstamp)

        import pdb
        pdb.set_trace()
        
        try:
            fh = urllib2.urlopen(url, timeout=TIMEOUT_SECS)
            data = fh.read().decode("utf-8")
            fh.close()
        except urllib2.HTTPError as htpe:
            time.wait(WAIT_SECS)
            try:
                fh = urllib.request.urlopen(url,timeout=TIMEOUT_SECS)
                data = fh.read().decode("utf-8")
                fh.close()
            except Exception as msg:
                raise Exception("Could not connect to server - %s." % url).with_traceback(msg2.__traceback__)

        with open(filename, "w") as fh:
            fh.write(data)
        self._extract_alerts(event)
        return

    def load(self, filename):
        """Load DM log file.

        :type filename: str
        :param filename: Name of local file with DM log.
        """
        with open(filename, "r") as fh:
            self.data = fh.read().decode("utf-8")
        return

    def get_num_dmalerts(self):
        """Get number of DM alerts.
        :returns: Number of alerts
        """
        raise NotImplementedError(":TODO: @brad")
        return
    
    def get_dmalert(self, index):
        """Get DM alert information.
        """
        raise NotImplementedError(":TODO: @brad")
        return

    def _extract_alerts(self, event):
        """Extract alert information.
        """
        raise NotImplementedError(":TODO: @brad")
        return
    

# End of file
