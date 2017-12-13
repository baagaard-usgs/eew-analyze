# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import numpy

import greatcircle

def shaking_time_vs(event, points, options):
    """Get shaking time based on origin time and shear wave speed.

    :type event: dict
    :param event: ComCat event dictionary (from AnalysisData).

    :type points: Numpy structured array
    :param points: Array with 'longitude' and 'latitude' point locations.

    :type options: dict
    :param options: Config options for shaking time.
    """
    SEC_TO_MSEC = 1000.0
    
    vs = float(options["vs_mps"])
    originTime = numpy.datetime64(event["origin_time"])
    originLon = event["longitude"]
    originLat = event["latitude"]
    originDepth = event["depth_km"]*1.0e+3

    distHoriz = greatcircle.distance(originLon, originLat, points["longitude"], points["latitude"])
    dist = (distHoriz**2 + originDepth**2)**0.5
    shakingTime = originTime + numpy.array(SEC_TO_MSEC*dist/vs, dtype="timedelta64[ms]")
    return shakingTime


# End of file
