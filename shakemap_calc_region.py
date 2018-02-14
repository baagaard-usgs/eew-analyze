#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import sys
import logging
from importlib import import_module
import numpy

from shakemap import ShakeMap
from osgeo import osr

import greatcircle
from comcat import DetailEvent

DEFAULTS = """
[events]
# Example:
# nc72923380 = Mw 4.6 Paicines, 2017-11-13

[mmi_predicted]
function = shakemap.gmpe
gmpe = ASK2014

[region]
mmi_threshold = 2.0
resolution_longitude = 0.008333333
resolution_latitude = 0.008333333
buffer_km = 25.0

[map]
scale = 1.0e+6
zoom_level = 8

[files]
event_dir = ./data/[EVENTID]/
"""

# ----------------------------------------------------------------------
class ShakeMapRegionApp(object):
    """Calculate target rectangular region for customized ShakeMap that encapsulates
    region with MMI above some threshold.
    """
    SPEC_COLS = [
        ("eqid", "|S12",),
        ("longitude_min", "float32",),
        ("longitude_max", "float32",),
        ("latitude_min", "float32",),
        ("latitude_max", "float32",),
        ("longitude_npoints", "int32",),
        ("latitude_npoints", "int32",),
        ("new", "int32",),
    ]
    
    def __init__(self):
        """Constructor.
        """
        self.config = None
        return

    def main(self):
        """Main entry point
        """
        # Initialization
        args = self._parseCommandLine()
        logLevel = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=logLevel, filename="shakemap_calc_region.log")
        if args.show_progress:
            self.showProgress = True
        self.initialize(args.config)

        # Show parameters
        if args.show_parameters or args.all:
            self.show_parameters()

        if args.calc_regions or args.all:
            self.calc_regions(args.filename, args.show_map)

        return

    def calc_regions(self, filename, showMap=False):
        """Compute target rectangular region for customized ShakeMap.
        """
        event = DetailEvent()
        shakemap = ShakeMap()
        eventDir = self.config.get("files", "event_dir")

        functionPath = self.config.get("mmi_predicted", "function").split(".")
        gmpe = getattr(import_module(".".join(functionPath[:-1])), functionPath[-1])

        buffer_km = self.config.getfloat("region", "buffer_km")
        
        events = self.config.options("events")
        numEvents = len(events)

        gridSpecs = numpy.zeros(numEvents, dtype=ShakeMapRegionApp.SPEC_COLS)
        
        for iEvent,eqId in enumerate(events):
            if self.showProgress:
                sys.stdout.write("\rProcessing events...{:d}%%".format(((iEvent+1)*100)/numEvents))
                sys.stdout.flush()

            dataDir = eventDir.replace("[EVENTID]", eqId)
            event.load(os.path.join(dataDir, eqId+".geojson"))

            distance = self._calc_distance(event, gmpe)
            gridSpec = self._calc_grid(event, distance + buffer_km*1000.0)

            shakemap.load(os.path.join(dataDir, "grid.xml.gz"))
            if shakemap.num_lon() >= gridSpec["longitude_npoints"] and shakemap.num_lat() >= gridSpec["latitude_npoints"]:
                gridSpec["new"] = 0
            if showMap:
                self._show_map(event, gridSpec)
                
            gridSpecs[iEvent] = tuple(gridSpec[col[0]] for col in ShakeMapRegionApp.SPEC_COLS)
            
        sys.stdout.write("\n")
            
        header = ", ".join(gridSpecs.dtype.names)
        numpy.savetxt(filename, gridSpecs, header=header, fmt="%s  %11.6f %11.6f %10.6f %10.6f  %4d %4d  %d")
        return

    def _calc_distance(self, event, gmpe):
        """
        :returns: Distance in m
        """
        eventDict = {
            "magnitude": event.magnitude,
            "longitude": event.longitude,
            "latitude": event.latitude,
            "depth_km": event.depth,
            }

        cols = [
            ("longitude", "float32",),
            ("latitude", "float32",),
            ("vs30", "float32",),
        ]
        longitude = event.longitude + numpy.arange(0.2, 25.0, 0.005)
        points = numpy.zeros(longitude.shape[-1], dtype=cols)
        points["longitude"] = longitude
        points["latitude"] = event.latitude
        points["vs30"] = 400.0 # m/s

        mmi = gmpe(eventDict, points, dict(self.config.items("mmi_predicted")))
        mmiThreshold = self.config.getfloat("region", "mmi_threshold")
        index = numpy.where(mmi < mmiThreshold)[0][0]

        if False: # :DEBUG:
            dist = greatcircle.distance(event.longitude, event.latitude, points["longitude"], points["latitude"])
            import matplotlib.pyplot as pyplot
            pyplot.plot(dist/1.0e+3, mmi, 'r')
            pyplot.axhline(mmiThreshold, color="black", linestyle="--")
            pyplot.show()
        
        return greatcircle.distance(event.longitude, event.latitude, points["longitude"][index], points["latitude"][index])


    def _calc_grid(self, event, distance):
        """Compute grid specifications based on distance from event.

        Compute extent of grid in x and y directions in projected
        coordinate system, and then compute lon/lat of grid extent.
        """
        x,y = self._project(event.longitude, event.latitude)
        
        # Compute extent of grid from epicenterXY
        lonMin = self._project(x - distance, y, inverse=True)[0]
        latMin = self._project(x, y - distance, inverse=True)[1]

        # Compute number of points along lon/lat directions.
        dLon = self.config.getfloat("region", "resolution_longitude")
        dLat = self.config.getfloat("region", "resolution_latitude")

        numLonHalf = 1 + int((event.longitude - lonMin) / dLon + 0.5)
        numLatHalf = 1 + int((event.latitude - latMin) / dLat + 0.5)

        lonMin = event.longitude - dLon * (numLonHalf - 1)
        latMin = event.latitude - dLat * (numLatHalf - 1)
        
        lonMax = lonMin + dLon * (2*numLonHalf - 2)
        latMax = latMin + dLat * (2*numLatHalf - 2)

        gridSpec = {
            "eqid": event.id,
            "longitude_min": lonMin,
            "longitude_max": lonMax,
            "latitude_min": latMin,
            "latitude_max": latMax,
            "longitude_npoints": 2*numLonHalf-1,
            "latitude_npoints": 2*numLatHalf-1,
            "new": 1,
            }
        return gridSpec

    def _project(self, xIn, yIn, inverse=False):
        """
        """
        srsWGS84 = osr.SpatialReference()
        srsWGS84.ImportFromEPSG(4326)

        srsProj = osr.SpatialReference()
        srsProj.ImportFromEPSG(3311)

        if not inverse:
            transf = osr.CoordinateTransformation(srsWGS84, srsProj)
        else:
            transf = osr.CoordinateTransformation(srsProj, srsWGS84)
        xyOut = numpy.array(transf.TransformPoints(numpy.vstack((xIn, yIn)).transpose()))
        return (xyOut[:,0], xyOut[:,1])

    def _show_map(self, event, gridSpec):
        """
        """
        import matplotlib.pyplot as pyplot
        import cartopy.crs
        from cartopy_extra_tiles.cached_tiler import CachedTiler
        from cartopy_extra_tiles.esri_tiles import ESRI
        
        INCH_PER_METER = 0.0254

        tiler = CachedTiler(ESRI(style="streetmap"), cache_dir="~/data_scratch/images/tiles")
        
        lonMin = gridSpec["longitude_min"]
        lonMax = gridSpec["longitude_max"]
        latMin = gridSpec["latitude_min"]
        latMax = gridSpec["latitude_max"]

        refDist = greatcircle.distance(event.longitude, event.latitude, lonMin, event.latitude)
        mapScale = self.config.getfloat("map", "scale")
        figWidthIn = 2 * refDist / mapScale / INCH_PER_METER
        figHeightIn = figWidthIn
        
        figure = pyplot.figure(figsize=(figWidthIn, figHeightIn))
        ax = pyplot.axes(projection=tiler.crs)
        ax.set_extent([lonMin, lonMax, latMin, latMax])
        ax.add_image(tiler, self.config.getint("map", "zoom_level"))
        ax.plot(event.longitude, event.latitude, "r*", ms=12, transform=cartopy.crs.Geodetic())
        pyplot.show()
        return
    
    def initialize(self, config_filenames):
        """Set parameters from config file and DEFAULTS.

        :type config_filename: str
        :param config_filename: Name of configuration (INI) file with parameters.
        """
        import ConfigParser
        import io
        config = ConfigParser.SafeConfigParser()
        config.readfp(io.BytesIO(DEFAULTS))
        if config_filenames:
            for filename in config_filenames.split(","):
                if self.showProgress:
                    print("Fetching parameters from {}...".format(filename))
                config.read(filename)

        self.config = config

        return
    
    def show_parameters(self):
        """Write parameters to stdout.
        """
        import sys
        self.config.write(sys.stdout)
        return

    def _parseCommandLine(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config")
        parser.add_argument("--filename", action="store", dest="filename", required=True)
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--calc-regions", action="store_true", dest="calc_regions")
        parser.add_argument("--show-map", action="store_true", dest="show_map")
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug")
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    ShakeMapRegionApp().main()


# End of file
