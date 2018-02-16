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

    def num_lon(self):
        """Get number of points along longitude direction.
        """
        return self.grid["num_longitude"]
    
    def num_lat(self):
        """Get number of points along latitude direction.
        """
        return self.grid["num_latitude"]
    
    def spatial_ref(self):
        """Get spatial reference system for ShakeMap.
        
        :returns: OSR SpatialReference
        """
        from osgeo import osr
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326) # WGS84
        return srs

    def geo_transform(self):
        """Get GDAL geometric transformation.

        The geometric transformation describes the coordinates of the
        points associated with the lon/lat grid.
        
        We use a grid origin that is offset 1/2 cell in lon/lat from
        the cell centers given by the longitude/latitude coordinates
        of the ShakeMap.

        :returns: GDAL GeoTransform for ShakeMap.

        """
        numLon = self.grid["num_longitude"]
        numLat = self.grid["num_latitude"]
        dLon = (self.grid["longitude_max"]-self.grid["longitude_min"]) / (self.grid["num_longitude"]-1)
        dLat = -(self.grid["latitude_max"]-self.grid["latitude_min"]) / (self.grid["num_latitude"]-1)
        originLon = self.grid["longitude_min"] - 0.5*dLon
        originLat = self.grid["latitude_max"] + 0.5*dLat
        return (originLon, dLon, 0, originLat, 0, dLat,)

    def pixel_area(self, projection):
        """Get pixel area of points in lon/lat grid.

        We project the points and compute the area assume the grid is
        a rhomus in the projected coordinate system.
        
        :type projection: str
        :param projection: Name of projection in the form EPSF:XXXX.
        """
        from osgeo import osr
        import greatcircle
        
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(projection.replace("EPSG:","")))
        return greatcircle.area_km2(self.data["longitude"], self.data["latitude"], self.grid["longitude_spacing"], self.grid["latitude_spacing"], srs)        
    
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
    
    logY = numpy.log10(MIN_FLOAT + pga*G_ACC) # acceleration in cm/s**2
    maskLower = logY <= 1.57
    mmiPGA = maskLower*(1.78 + 1.55*logY) + ~maskLower*(-1.60 + 3.70*logY)
    
    logY = numpy.log10(MIN_FLOAT + pgv) # velocity in cm/s
    maskLower = logY <= 0.53
    mmiPGV = maskLower*(3.78 + 1.47*logY) + ~maskLower*(2.89 + 3.16*logY)
    
    maskPGA = pga > 0.0
    maskPGV = pgv > 0.0
    mmi = 0.5*(mmiPGA+mmiPGV)*(maskPGA*maskPGV) + mmiPGA*(maskPGA*~maskPGV) + mmiPGV*(~maskPGA*maskPGV)
    mmi = numpy.maximum(1.0, mmi)
    mmi = numpy.minimum(10.0, mmi)
    return mmi


def mmi_WaldEtal1999(pga, pgv):
    """Use Wald et al. (1999) to convert PGA and PGV to MMI.
    
    :type pga: Numpy array
    :param pga: PGA in percent g.
    
    :type pgv: Numpy array
    :param pgv: PGV in cm/s.
    """
    MIN_FLOAT = 1.0e-20
    G_ACC = 9.80665
    
    mmiPGA = numpy.zeros(pga.shape)
    mmiPGV = numpy.zeros(pgv.shape)
    
    logY = numpy.log10(MIN_FLOAT + pga*G_ACC) # acceleration in cm/s**2
    maskLower = logY <= 1.82
    mmiPGA = maskLower*(1.00 + 2.20*logY) + ~maskLower*(-1.66 + 3.66*logY)
    
    logY = numpy.log10(MIN_FLOAT + pgv) # velocity in cm/s
    maskLower = logY <= 0.76
    mmiPGV = maskLower*(3.40 + 2.10*logY) + ~maskLower*(2.35 + 3.47*logY)
    
    maskPGA = pga > 0.0
    maskPGV = pgv > 0.0
    wtPGA = 1.0*(mmiPGA <= 5.0) + (7.0 - mmiPGA)/(7.0 - 5.0)*numpy.bitwise_and(mmiPGA > 5.0, mmiPGA < 7.0)
    wtPGV = 1.0*(mmiPGV >= 7.0) + (mmiPGA - 5.0)/(7.0 - 5.0)*numpy.bitwise_and(mmiPGA > 5.0, mmiPGA < 7.0)

    mmi = (wtPGA*mmiPGA + wtPGV*mmiPGV)*(maskPGA*maskPGV) + mmiPGA*(maskPGA*~maskPGV) + mmiPGV*(~maskPGA*maskPGV)
    mmi = numpy.maximum(1.0, mmi)
    mmi = numpy.minimum(10.0, mmi)
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

    if options["gmice"] == "WordenEtal2012":
        gmice = mmi_WordenEtal2012
    elif options["gmice"] == "WaldEtal1999":
        gmice = mmi_WaldEtal1999
    
    return gmice(values["pgaG"]*100.0, values["pgvCmps"])


if __name__ == "__main__":
    event = {
        "magnitude": 6.0,
        "longitude": -120.0,
        "latitude": 37.00,
        "depth_km": 8.0,
    }
    
    cols = [
        ("longitude", "float32",),
        ("latitude", "float32",),
        ("vs30", "float32",),
    ]
    longitude = event["longitude"] + numpy.arange(0.2, 25.0, 0.005)
    points = numpy.zeros(longitude.shape[-1], dtype=cols)
    points["longitude"] = longitude
    points["latitude"] = event["latitude"]
    points["vs30"] = 400.0 # m/s

    options = {
        "gmpe": "CB2014",
        "gmice": "WordenEtal2012",
    }

    
    mmiWorden = gmpe(event, points, options)
    options["gmice"] = "WaldEtal1999"
    mmiWald = gmpe(event, points, options)

    import greatcircle
    dist = greatcircle.distance(event["longitude"], event["latitude"], points["longitude"], points["latitude"])
    import matplotlib.pyplot as pyplot
    pyplot.semilogx(dist/1.0e+3, mmiWorden, 'r-', dist/1.0e+3, mmiWald, 'b--')

    mmiThreshold = 2.0
    pyplot.axhline(mmiThreshold, color="black", linestyle=":")
    pyplot.xlabel("Distance (km)")
    pyplot.ylabel("MMI")
    pyplot.legend(["Worden et al. (2012)", "Wald et al. (1999)"])
    pyplot.show()

        
# End of file
