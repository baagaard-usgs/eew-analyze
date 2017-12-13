# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import numpy

EARTH_MEAN_RADIUS_M = 6371.0e+3
DEG_TO_RAD = numpy.pi / 180.0

def distance(refLon, refLat, ptsLon, ptsLat):
    """Get great circle distance in meters from reference point to points.

    Source: https://en.wikipedia.org/wiki/Great-circle_distance

    :type refLon: float
    :param refLon: Longitude of reference point in degrees.

    :type refLat: float
    :param refLat: Latitude of reference point in degrees.

    :type ptsLon: Numpy array
    :param ptsLon: Longitude of points in degrees.

    :type ptsLat: Numpy array
    :param ptsLat: Latitude of points in degrees.
    """
    refLonR = refLon * DEG_TO_RAD
    refLatR = refLat * DEG_TO_RAD
    ptsLonR = ptsLon * DEG_TO_RAD
    ptsLatR = ptsLat * DEG_TO_RAD

    p = numpy.sin(0.5*(ptsLatR-refLatR))**2 \
        + numpy.cos(refLatR)*numpy.cos(ptsLatR)*numpy.sin(0.5*(ptsLonR-refLonR))**2
    return EARTH_MEAN_RADIUS_M * 2.0*numpy.arcsin(p**0.5)


# ======================================================================
# Test
if __name__ == "__main__":
    import pyproj
    import math

    refLon = -122.56
    refLat = 37.75
    
    utmZone = int(math.floor(refLon + 180.0)/6.0) + 1
    proj = pyproj.Proj(proj="utm", zone=utmZone, ellps="WGS84")
    originX, originY = proj(refLon, refLat)

    points = numpy.array([
        [-122.60, 37.85],
        [-120.00, 37.45],
        [-121.45, 36.56],
        [-115.00, 35.00],
        ], dtype=numpy.float64)

    x, y = proj(points[:,0], points[:,1])
    dist = ((x-originX)**2 + (y-originY)**2)**0.5

    distGR = distance(refLon, refLat, points[:,0], points[:,1])

    print(1.0 - distGR / dist)
    

# End of file
