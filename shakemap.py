# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import gzip
import cStringIO
import numpy
import logging
from lxml import etree

class ShakeMap(object):
    """ShakeMap reader.
    """

    def __init__(self):
        """Constructor.
        """
        self.data = None
        self.grid = None
        return

    def load(self, filename):
        """Load ShakeMap grid file.

        :type filename: str
        :param filename: Name of local file with ShakeMap grid data.
        """
        suffix = ""
        if not filename.endswith(".gz"):
            suffix = ".gz"
        with gzip.open(filename+suffix, "r") as fh:
            self._parse(fh)
        return

    def _parse(self, fh):
        """Parse ShakeMap grid XML file.
        
        :type fh: File handle
        :param fh: File handle to ShakeMap grid XML file.
        :returns: Numpy structured array with ShakeMap grid data.
        """
        NS = "http://earthquake.usgs.gov/eqcenter/shakemap"
        namespaces = {"ns": NS}
        
        bytes = fh.read()
        elRoot = etree.fromstring(bytes)

        # Grid specification
        elGrid = elRoot.xpath("ns:grid_specification", namespaces=namespaces)[0]
        self.grid = {
            "longitude_min": float(elGrid.get("lon_min")),
            "longitude_max": float(elGrid.get("lon_max")),
            "latitude_min": float(elGrid.get("lat_min")),
            "latitude_max": float(elGrid.get("lat_max")),
            "longitude_spacing": float(elGrid.get("nominal_lon_spacing")),
            "latitude_spacing": float(elGrid.get("nominal_lat_spacing")),
            "num_longitude": int(elGrid.get("nlon")),
            "num_latitude": int(elGrid.get("nlat")),
        }
        
        # Grid values
        lonIndex = int(elRoot.xpath("ns:grid_field[@name='LON']", namespaces=namespaces)[0].get("index"))-1
        latIndex = int(elRoot.xpath("ns:grid_field[@name='LAT']", namespaces=namespaces)[0].get("index"))-1
        mmiIndex = int(elRoot.xpath("ns:grid_field[@name='MMI']", namespaces=namespaces)[0].get("index"))-1
        #pgaIndex = int(elRoot.xpath("ns:grid_field[@name='PGA']", namespaces=namespaces)[0].get("index"))-1
        #pgvIndex = int(elRoot.xpath("ns:grid_field[@name='PGV']", namespaces=namespaces)[0].get("index"))-1
        colIndices = (lonIndex, latIndex, mmiIndex,)
        colNames = "longitude, latitude, mmi"
        colFormats = "float32, float32, float32"

        data = numpy.loadtxt(cStringIO.StringIO(elRoot.xpath("ns:grid_data", namespaces=namespaces)[0].text), usecols=colIndices)
        self.data = numpy.core.records.fromarrays(data.transpose(), names=colNames, formats=colFormats)
        return
        

# End of file
