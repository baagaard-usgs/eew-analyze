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

from openquake_gmpe import OpenQuakeGMPE

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
        pgaIndex = int(elRoot.xpath("ns:grid_field[@name='PGA']", namespaces=namespaces)[0].get("index"))-1
        pgvIndex = int(elRoot.xpath("ns:grid_field[@name='PGV']", namespaces=namespaces)[0].get("index"))-1
        vs30Index = int(elRoot.xpath("ns:grid_field[@name='SVEL']", namespaces=namespaces)[0].get("index"))-1
        colIndices = (lonIndex, latIndex, mmiIndex, pgaIndex, pgvIndex, vs30Index)
        colNames = "longitude, latitude, mmi, pga, pgv, vs30"
        colFormats = "float32, float32, float32, float32, float32, float32"

        data = numpy.loadtxt(cStringIO.StringIO(elRoot.xpath("ns:grid_data", namespaces=namespaces)[0].text), usecols=colIndices)
        self.data = numpy.core.records.fromarrays(data.transpose(), names=colNames, formats=colFormats)
        self.data["mmi"] = mmi_WordenEtal2012(self.data["pga"], self.data["pgv"])
        return
        

def mmi_WordenEtal2012(pga, pgv):
    """Use Worden et al. (2012) to convert PGA and PGV to MMI.
    
    :type pga: Numpy array
    :param pga: PGA in percent g.
    
    :type pgv: Numpy array
    :param pgv: PGV in cm/s.
    """
    MIN_FLOAT = 1.0e-20
    G_ACC = 9.80665
    
    mmiPGA = numpy.zeros(pga.shape)
    mmiPGV = numpy.zeros(pgv.shape)
    
    logY = numpy.log10(MIN_FLOAT + pga*G_ACC)
    maskLower = logY <= 1.57
    mmiPGA = maskLower*(1.78 + 1.55*logY) + ~maskLower*(-1.60 + 3.70*logY)
    mmiPGA = numpy.maximum(1.0, mmiPGA)
    mmiPGA = numpy.minimum(10.0, mmiPGA)

    
    logY = numpy.log10(MIN_FLOAT + pgv)
    maskLower = logY <= 0.53
    mmiPGV = maskLower*(3.78 + 1.47*logY) + ~maskLower*(2.89 + 3.16*logY)
    mmiPGV = numpy.maximum(1.0, mmiPGV)
    mmiPGV = numpy.minimum(10.0, mmiPGV)
    
    maskPGA = mmiPGA > 0.0
    maskPGV = mmiPGV > 0.0
    mmi = 0.5*(mmiPGA+mmiPGV)*(maskPGA*maskPGV) + mmiPGA*(maskPGA*~maskPGV) + mmiPGV*(~maskPGA*maskPGV)

    return mmi


def gmpe(alert, points, options):
    """Get predicted MMI for given alert.

    RotD50 to "larger" (max) PGA/PGV from Beyer and Bommer, BSSA (2006) doi: 10.1785/0120050210.

    :type alert: dict
    :param alert: ShakeAlert alert dictionary (from AnalysisData).

    :type points: Numpy structured array
    :param points: Array with 'longitude' and 'latitude' point locations.

    :type options: dict
    :param options: Config options for GMPE.
    """
    ROTD50_TO_PGA_LARGER = 1.1
    ROTD50_TO_PGV_LARGER = 1.0
    
    oqGMPE = OpenQuakeGMPE(options["gmpe"])
    values = oqGMPE.computeMean(alert, points)

    values["pgaG"] *= ROTD50_TO_PGA_LARGER
    values["pgvCmps"] *= ROTD50_TO_PGV_LARGER
    
    return mmi_WordenEtal2012(values["pgaG"]*100.0, values["pgvCmps"])



# End of file
