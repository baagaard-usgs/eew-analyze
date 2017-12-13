# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import numpy
import pyproj
import math

def shaking_time_vs(event, points, options):
    """Get shaking time based on origin time and shear wave speed.

    :type event: dict
    :param event: ComCat event dictionary (from AnalysisData).

    :type points: Numpy structured array
    :param points: Array with 'longitude' and 'latitude' point locations.

    :type options: dict
    :param options: Config options for shaking time.
    """
    vs = float(options["vs"])
    originTime = numpy.datetime64(event["origin_time"])
    originLon = event["longitude"]
    originLat = event["latitude"]
    originDepth = event["depth_km"]*1.0e+3
    
    utmZone = int(math.floor(originLon + 180.0)/6.0) + 1
    proj = pyproj.Proj(proj="utm", zone=utmZone, ellps="WGS84")
    originX, originY = proj(originLon, originLat)
    x, y = proj(points["longitude"], points["latitude"])
    dist = ((x-originX)**2 + (y-originY)**2 + originDepth**2)**0.5
    shakingTime = originTime + numpy.array(1000.0*dist/vs, dtype="timedelta64[ms]")
    return shakingTime


# End of file
