# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#


import logging
import numpy

from . import greatcircle

import openquake.hazardlib.gsim.base
import openquake.hazardlib.imt
import openquake.hazardlib.const


# ----------------------------------------------------------------------
class OpenQuakeGMPE(object):
    """
    Ground-motion prediction equation (GMPE) interface to OpenQuake hazardlib.
    """
    FIELDS = {
        "pgaG": openquake.hazardlib.imt.PGA(),
        "pgvCmps": openquake.hazardlib.imt.PGV(),
    }
    
    def __init__(self, name):
        """
        Constructor.

        :type name: str
        :param name:
            Name for GMPE [BSSA2014, ASK2014, CB2014, CY2014]
        """
        if name == "BSSA2014":
            import openquake.hazardlib.gsim.boore_2014
            self.gmpe = openquake.hazardlib.gsim.boore_2014.BooreEtAl2014()
        elif name == "ASK2014":
            import openquake.hazardlib.gsim.abrahamson_2014
            self.gmpe = openquake.hazardlib.gsim.abrahamson_2014.AbrahamsonEtAl2014()
        elif name == "CB2014":
            import openquake.hazardlib.gsim.campbell_bozorgnia_2014
            self.gmpe = openquake.hazardlib.gsim.campbell_bozorgnia_2014.CampbellBozorgnia2014()
        elif name == "CY2014":
            import openquake.hazardlib.gsim.chiou_youngs_2014
            self.gmpe = openquake.hazardlib.gsim.chiou_youngs_2014.ChiouYoungs2014()
        else:
            raise ValueError("Unknown OpenQuake GMPE '%s'." % name)
        return

    def computeMean(self, event, points):
        """
        Compute mean PGA (g) and PGV (cm/s).

        :type event: dict
        :param event:
            Dictionary with event parameters ["magnitude", "longitude", "latitude", "depth_km"]

        :type points: dict of Numpy arrays or Numpy structured array
        :param points:
            Point locations and metadata ["longitude", "latitude", "vs30"].

        :returns: tuple
           Tuple of PGA and PGV.
        
        """
        ruptureContext = self._ruptureContext(event)
        sitesContext = self._sitesContext(points)
        distContext = self._distanceContext(event, points, ruptureContext)

        fields = OpenQuakeGMPE.FIELDS
        
        numFields = len(fields.keys())
        numPoints = points["longitude"].shape[0]
        dtype = {
            "names": [s for s in sorted(fields.keys())],
            "formats": ["float64"]*numFields,
        }
        data = numpy.zeros(numPoints, dtype=dtype)

        for field in fields:
            imt = self.FIELDS[field]
            (values,stddev) = self.gmpe.get_mean_and_stddevs(sitesContext, ruptureContext, distContext, imt, [openquake.hazardlib.const.StdDev.TOTAL])
            data[field] = self.gmpe.to_imt_unit_values(values)
        return data

    def _ruptureContext(self, event):
        """Set rupture parameters. 

        Assumptions:
            1. Vertical, strike-slip fault.
            2. Width from a circular rupture with a stress drop of 3.0 MPa.

        :type event: dict
        :param event:
            Dictionary with event parameters ["magnitude", "depth_km"]
        """
        STRESS_DROP = 3.0e+6 # Pa
        seismicMoment = 10**(1.5*(event["magnitude"]+10.7)-7.0)
        radiusKm = 1.0e-3 * (7.0/16.0*seismicMoment/STRESS_DROP)**(1/3.0)
        widthKm = 2.0 * numpy.minimum(event["depth_km"], radiusKm)
        
        context = openquake.hazardlib.gsim.base.RuptureContext()
        context.mag = event["magnitude"]
        context.dip = 90.0
        context.rake = 0.0
        context.ztor = numpy.maximum(0.0, event["depth_km"] - 0.5 * widthKm)
        context.hypo_depth = event["depth_km"]
        context.width = widthKm # only used in CB2014 and ASK2014 hanging wall terms
        return context

    def _sitesContext(self, points):
        """Set site parameters.

        Assumptions:
            1. Z1.0 from vs30 using relation given in ASK2014.
            2. Z2.5 from vs30 using relation given in CB2014.

        :type points: dict of Numpy arrays or Numpy structured array
        :param points:
            Point locations and metadata ["vs30"].
        """
        # Compute Z1.0 using ASK2014 reference relation (similar to CY2014)
        z1pt0Km = 1.0e-3 * numpy.exp((-7.67 / 4.) * numpy.log((points["vs30"]**4 + 610.**4) / (1360.**4 + 610.**4)))

        # Compute Z2.5 using CB2014 California equation.
        z2pt5Km = numpy.exp(7.089 - 1.144 * numpy.log(points["vs30"]))
        
        context = openquake.hazardlib.gsim.base.SitesContext()
        context.vs30 = points["vs30"]
        context.z1pt0 = z1pt0Km
        context.z2pt5 = z2pt5Km
        context.vs30measured = False*numpy.ones(points["vs30"].shape, dtype=numpy.bool)
        return context

    def _distanceContext(self, event, points, ruptureContext):
        """Get rupture distance information using great circle path.

        Assumptions:
            1. Vertical, strike-slip fault.
            2. Rx and Ry0 = 0.0 (neglect hanging wall effects in CY2014, CB2014, ASK2014)

        :type event: dict
        :param event:
            Dictionary with event parameters ["longitude", "latitude", "depth_km"]

        :type points: dict of Numpy arrays or Numpy structured array
        :param points:
            Point locations and metadata ["longitude", "latitude"].
        """
        context = openquake.hazardlib.gsim.base.DistancesContext()
        distEpiKm = 1.0e-3*greatcircle.distance(event["longitude"], event["latitude"], points["longitude"], points["latitude"])
        distRupKm = (distEpiKm**2 + ruptureContext.ztor**2)**0.5 # Assumes vertical fault and point source (ignore strike)
        context.rjb = distEpiKm
        context.rrup = distRupKm
        context.rx = numpy.zeros(distEpiKm.shape) # Neglect hanging wall effects
        context.ry0 = numpy.zeros(distEpiKm.shape) # Neglect hanging wall effects
        return context

if __name__ == "__main__":
    # Test
    event = {
        "magnitude": 4.38, # overwritten
        "depth_km": 12.3,
        "longitude": -122.0,
        "latitude": 38.0,
    }
    colNames = "longitude, latitude, vs30"
    colFormats = "float32, float32, float32"
    longitude = event["longitude"] + numpy.arange(0.0, 4.01, 0.05)
    latitude = event["latitude"]*numpy.ones(longitude.shape)
    npts = longitude.shape[-1]
    points = numpy.zeros(npts, dtype=[("longitude", numpy.float32,), ("latitude", numpy.float32,), ("vs30", numpy.float32,)])
    points["longitude"] = longitude
    points["latitude"] = latitude
    points["vs30"] = 540.0

    import matplotlib.pyplot as pyplot
    import matplotlib_extras
    figure = pyplot.figure(figsize=(12.0, 5.5))
    nrows = 1
    ncols = 2
    rectFactory = matplotlib_extras.axes.RectFactory(figure, nrows, ncols, margins=((0.7, 1.0, 0.1), (0.50, 0, 0.25)))
    
    axPGA = figure.add_axes(rectFactory.rect(col=1))
    axPGA.autoscale(tight=True)
    axPGA.set_xlabel("Joyner-Boore Distance (km)")
    axPGA.set_ylabel("PGA (g)")

    axPGV = figure.add_axes(rectFactory.rect(col=2))
    axPGV.autoscale(tight=True)
    axPGV.set_xlabel("Joyner-Boore Distance (km)")
    axPGV.set_ylabel("PGV (cm/s)")

    gmpes = {
        "BSSA2014": "red",
        "ASK2014": "blue",
        "CB2014": "orange",
        "CY2014": "purple",
    }

    magnitudes = numpy.arange(4.0, 7.01, 1.0)
    for label, color in gmpes.items():
        gmpe = OpenQuakeGMPE(label)
        for magnitude in magnitudes:
            event["magnitude"] = magnitude
            data = gmpe.computeMean(event, points)
            rupContext = gmpe._ruptureContext(event)
            distContext = gmpe._distanceContext(event, points, rupContext)
            axPGA.loglog(distContext.rjb, data["pgaG"], color=color)
            axPGV.loglog(distContext.rjb, data["pgvCmps"], color=color)

    import matplotlib.pyplot as pyplot
    pyplot.show()
        

    # Multiple events, fixed distance

    
# End of file
