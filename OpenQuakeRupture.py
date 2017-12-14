# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#


import logging
import numpy

import openquake.hazardlib.gsim.base
from openquake.hazardlib.geo.mesh import Mesh
from openquake.hazardlib.geo.surface.planar import PlanarSurface
from openquake.hazardlib.geo.point import Point

# ----------------------------------------------------------------------
class OpenQuakeRupture(object):
    """
    Earthquake rupture interface to OpenQuake hazardlib.
    """
    def __init__(self, magnitude, strikeDeg, dipDeg, rakeDeg, lengthKm, widthKm):
        """
        Constructor.
        """
        self.eqMagnitude = magnitude
        self.strikeDeg = strikeDeg
        self.dipDeg = dipDeg
        self.rakeDeg = rakeDeg
        self.lengthKm = lengthKm
        self.widthKm = widthKm

        self.strikeRad = strikeDeg*numpy.pi/180.0
        self.dipRad = dipDeg*numpy.pi/180.0
        return

    def initialize(self):
        """
        Initialize rupture geometry.
        """
        ORIGIN_LONLAT = (-122.4149, 37.7749) # San Francisco
        SPACING_KM = 5.0
        
        import pyproj
        self.proj = pyproj.Proj(proj="tmerc", lon_0=ORIGIN_LONLAT[0], lat_0=ORIGIN_LONLAT[1], datum="WGS84")
        x0,y0 = self.proj(ORIGIN_LONLAT[0], ORIGIN_LONLAT[1])
        self.x0 = x0
        self.y0 = y0

        # Top-left corner
        topLx = x0 - 1.0e+3*0.5*self.lengthKm*numpy.sin(self.strikeRad)
        topLy = y0 - 1.0e+3*0.5*self.lengthKm*numpy.cos(self.strikeRad)
        lon, lat = self.proj(topLx, topLy, inverse=True)
        topL = Point(lon, lat)

        # Top-right corner
        topRx = x0 + 1.0e+3*0.5*self.lengthKm*numpy.sin(self.strikeRad)
        topRy = y0 + 1.0e+3*0.5*self.lengthKm*numpy.cos(self.strikeRad)
        lon, lat = self.proj(topRx, topRy, inverse=True)
        topR = Point(lon, lat)

        # Bottom-left corner
        botLx = topLx + 1.0e+3*self.widthKm*numpy.cos(self.strikeRad)*numpy.cos(self.dipRad)
        botLy = topLy - 1.0e+3*self.widthKm*numpy.sin(self.strikeRad)*numpy.cos(self.dipRad)
        botLz = 1.0e+3*self.widthKm*numpy.sin(self.dipRad)
        lon, lat = self.proj(botLx, botLy, inverse=True)
        botL = Point(lon, lat, 1.0e-3*botLz)
        
        # Bottom-right corner
        botRx = topRx + 1.0e+3*self.widthKm*numpy.cos(self.strikeRad)*numpy.cos(self.dipRad)
        botRy = topRy - 1.0e+3*self.widthKm*numpy.sin(self.strikeRad)*numpy.cos(self.dipRad)
        botRz = 1.0e+3*self.widthKm*numpy.sin(self.dipRad)
        lon, lat = self.proj(botRx, botRy, inverse=True)
        botR = Point(lon, lat, 1.0e-3*botRz)

        self.rupture = PlanarSurface(SPACING_KM, self.strikeDeg, self.dipDeg, topL, topR, botR, botL)
        return

    def profileStrike(self, lengthKm=400.0, dxKm=2.5):
        """
        Create Mesh object for profile along fault strike.
        """
        d = 1.0e+3*numpy.arange(-0.5*lengthKm, 0.5*lengthKm+0.1*dxKm, dxKm)
        x = self.x0 + d*numpy.sin(self.strikeRad)
        y = self.y0 + d*numpy.cos(self.strikeRad)

        import matplotlib.pyplot as pyplot
        lon,lat = self.proj(x, y, inverse=True)
        return d, Mesh(lon, lat)

    def profilePerpStrike(self, lengthKm=400.0, dxKm=2.5):
        """
        Create Mesh object for profile perpendicular to fault strike.
        """
        d = 1.0e+3*numpy.arange(-0.5*lengthKm, 0.5*lengthKm+0.1*dxKm, dxKm)
        x = self.x0 + d*numpy.cos(self.strikeRad)
        y = self.y0 - d*numpy.sin(self.strikeRad)
        lon,lat = self.proj(x, y, inverse=True)
        return d, Mesh(lon, lat)

# End of file
