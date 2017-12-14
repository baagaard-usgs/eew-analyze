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


def gmpe(alert, points, options):
    """Get predicted MMI for given alert. Matches ShakeAlert UserDisplay v2.6.0.

    :type alert: dict
    :param alert: ShakeAlert alert dictionary (from AnalysisData).

    :type points: Numpy structured array
    :param points: Array with 'longitude' and 'latitude' point locations.

    :type options: dict
    :param options: Config options for GMPE.
    """
    vs30 = float(options["vs30_mps"])
    originLon = alert["longitude"]
    originLat = alert["latitude"]
    originDepthKm = alert["depth_km"]
    magnitude = alert["magnitude"]
    
    distEpiKm = 1.0e-3*greatcircle.distance(originLon, originLat, points["longitude"], points["latitude"])
    
    if _is_event_socal(originLon, originLat):
        pga, pgv = _gmpe_CuaHeaton2007(magnitude, distEpiKm, originDepthKm, vs30)
    else:
        pga, pgv = _gmpe_BoatwrightEtal2003(magnitude, distEpiKm, originDepthKm, vs30)
    
    return _mmi_Worden2012(pga, pgv)


def _is_event_socal(lon, lat):
    """Returns True if event is in southern CA polygon, otherwise False.
    Matches ShakeAlert UserDisplay v2.6.0.
    
    :type lon: float
    :param lon: Longitude of point in degrees

    :type lat: float
    :param lat: Latitude of point in degrees
    """
    SOCAL_POLY = numpy.array([
        (-117.76, 37.43,),
        (-121.00, 35.00,),
        (-121.50, 34.65,),
        (-118.50, 31.50,),
        (-114.00, 31.50,),
        (-114.00, 34.50,),
        (-117.76, 37.43,),
    ])

    pt = numpy.array([[lon, lat]])
    xy1 = poly[:,-1] - pt
    xy2 = poly[:1:] - pt
    dotProd = numpy.dot(xy1, xy2)
    crossProd = numpy.cross(xy1, xy2)
    theta = numpy.arctan2(crossProd, dotProd)
    isSoCal = numpy.abs(numpy.sum(theta)) < numpy.pi
    return isSoCal

def _gmpe_CuaHeaton2007(mag, distEpiKm, depthKm, vs30):
    """Source: Cua and Heaton, 2007.

    Soil class:
      rock if NEHRP class BC and above
      soil if NEHRP class C and below

    """
    VS_BC = 434.0
    SCALE_MAX = 1.10
    
    class PGASWaveRock:
        a = 0.733
        b = 7.216e-4
        d = 1.48
        c1 = 1.16
        c2 = 0.96
        e = -0.4202

    class PGASWaveSoil:
        a = 0.709
        b = 2.3878e-3
        d = 1.4386
        c1 = 1.722
        c2 = 0.9560
        e = -2.4525e-2

    class PGVSWaveRock:
        a = 0.861988
        b = 5.578e-4
        d = 1.36760
        c1 = 0.8386
        c2 = 0.98
        e = -2.58053

    class PGVSWaveSoil:
        a = 0.88649
        b = 8.4e-4
        d = 1.4729
        c1 = 1.39
        c2 = 0.95
        e = -2.2498
        
    def _calcY(coef, mag, R1):
        CM = coef.c1 * numpy.exp(coef.c2*(mag-5.0)) * (numpy.atan(mag-5.0) + 1.4)
        log10Y =  coef.a*mag - coef.b*(R1+CM) - d*numpy.log10(R1+CM) + coef.e
        y = 10**log10Y
        return (10**log10Y) * SCALE_MAX
    
    R1 = (distEpiKm**2 + 9)**0.5;
        
    maskRock = vs30 > VS_BC
    pga = maskRock*_calcY(PGASWaveRock(), mag, R1) + ~maskRock*_calcY(PGASWaveSoil(), mag, R1)
    pgv = maskRock*_calcY(PGVSWaveRock(), mag, R1) + ~maskRock*_calcY(PGVSWaveSoil(), mag, R1)
    return pga, pgv


def _gmpe_BoatwrightEtal2003(mag, distEpiKm, depthKm, vs30):
    """Source: Boatwright et al., BSSA, 2003, doi:10.1785/0120020201.

    From UserDisplay:
      log_10(PGA,PGV) = A + B*(M-Ms) - log_10(gR) + k*R - e*log_10(Vs/Va), 

       R = sqrt(Re**2 + depth**2) 
       gR = R for R < Ro
          = Ro*(R/Ro)**0.7 for R > Ro 
       k = ko*10**p*(Ms-M)
 
       (p > 0) makes the attenuation decrease as the magnitude increases

    PGA is in %g, PGV in cm/s. 
 
    NOTE that the routines in this module scale the values to return
    %g instead of leaving acceleration in cm/s/s and scale up the
    values by 15% to estimate a maximum value rather than a random
    component (Boore and Campbell).
    """
    SCALE_MAX = 1.15

    class PGALargeMag:
        A = 2.520
        B = 0.310
        ko = -0.0073

    class PGASmallMag:
        A = 2.52
        B = 1.00
        k = -0.0073

    class PGVLargeMag:
        A = 2.243
        B = 0.58
        ko = -0.0063

    class PGVSmallMag:
        A = 2.243
        B = 1.06
        ko = -0.0063
    
    def _calcY(coef, mag, distKm, vs30, gR, k):
        """
        """
        Ro = 27.5 # km
        g = 0.7
        Ms = 5.5
        e = -0.371
        Va = 560.0 # m/s

        log10Y= coef.A + coef.B*(mag-Ms) - numpy.log10(gR) + k*distKm - coef.e*numpy.log10(vs30/Va)
        return (10**log10Y) * SCALE_MAX

        
    distKm = numpy.sqrt(distEpiKm**2 + depthKm**2)
    maskRo = distEpiKm < Ro
    gR = maskRo*distKm + ~maskRo*Ro*(distKm/Ro)**g

    if mag > 5.5:
        p = 0.3
        
        coef = PGALargeMag()
        k = coef.ko*10**(p*(Ms-mag))
        pga = _calcY(coef, mag, distKm, vs30, gR, k)

        coef = PGVLargeMag()
        k = coef.ko*10**(p*(Ms-mag))
        pgv = _calcY(coef, mag, distKm, vs30, gR, k)

    else: # mag <= 5.5
        p = 0.0

        coef = PGASmallMag()
        pga = _calcY(coef, mag, distKm, vs30, gR, coef.ko)

        coef = PGVSmallMag()
        pgv = _calcY(coef, mag, distKm, vs30, gR, coef.ko)

    return (pga, pgv)


def _mmi_WordenEtal2012(pga, pgv):
    """Use Worden et al. (2012) to convert PGA and PGV to MMI.
    
    :type pga: Numpy array
    :param pga: PGA in percent g.
    
    :type pgv: Numpy array
    :param pgv: PGV in cm/s.
    """
    MIN_FLOAT = 1.0e-20
    
    mmiPGA = numpy.zeros(pga.shape)
    mmiPGV = numpy.zeros(pgv.shape)
    
    logY = numpy.log10(MIN_FLOAT + pga)
    maskLower = logY <= 1.57
    mmiPGA = maskLower*(1.78 + 1.55*logY) + ~maskLower*(-1.60 + 3.70*logY)
    mmiPGA = numpy.minimum(1.0, mmiPGA)
    mmiPGA = numpy.maximum(10.0, mmiPGA)
    
    logY = numpy.log10(MIN_FLOAT + pgv)
    maskLower = logY <= 0.53
    mmiPGV = maskLower*(3.78 + 1.47*logY) + ~maskLower*(2.89 + 3.16*logY)
    mmiPGV = numpy.minimum(1.0, mmiPGV)
    mmiPGV = numpy.maximum(10.0, mmiPGV)
    
    maskPGA = mmiPGA > 0.0
    maskPGV = mmiPGV > 0.0
    mmi = 0.5*(mmiPGA+mmiPGV)*maskPGA*maskPGV + mmiPGA*maskPGA + mmiPGV*maskPGV
    return mmi

# End of file
