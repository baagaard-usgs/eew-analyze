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
import logging
import numpy
from osgeo import ogr, gdal

gdal.UseExceptions()

DEFAULTS = u"""
[census.geometry]
shapefile = tl_2010_06_tabblock10.shp
field_block = GEOID10
field_area = ALAND10

[census.population]
shapefile = tabblock2010_06_pophu.shp
field_block = BLOCKID10
field_population = POP10

[raster]
resolution_deg = 0.005
filename = population_density.tiff
"""

# ----------------------------------------------------------------------

class PopulationDensityApp(object):
    """
    Create GeoTIFF raster with population density in people per square
    kilometer from census block shapefiles.
    """

    def __init__(self):
        """Constructor.

        self.params = None
        """
        return

    def main(self):
        """Main entry point.
        """
        # Initialization
        args = self._parseCommandLine()
        logLevel = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=logLevel, filename="generate_popdensity.log")
        if args.show_progress:
            self.showProgress = True
        self.initialize(args.config)

        if args.show_parameters or args.all:
            self.show_parameters()

        if args.generate or args.all:
            popDensityLayer = self._joinShapefiles()
            self._rasterize(popDensityLayer)
        return

    def initialize(self, config_filenames):
        """Set parameters from config file and DEFAULTS.

        :type config_filename: str
        :param config_filename: Name of configuration (INI) file with parameters.
        """
        import configparser
        import io
        config = configparser.SafeConfigParser()
        config.readfp(io.StringIO(DEFAULTS))
        if config_filenames:
            for filename in config_filenames.split(","):
                if self.showProgress:
                    print("Fetching parameters from %s..." % filename)
                config.read(filename)

        self.params = config

        return
    
    def show_parameters(self):
        """Write parameters to stdout.
        """
        import sys
        self.params.write(sys.stdout)
        return

    def _joinShapefiles(self):
        driverShapefile = ogr.GetDriverByName("ESRI Shapefile")

        # Get id and area from Census block file with geometry
        filename = self.params.get("census.geometry", "shapefile")
        fieldBlock = self.params.get("census.geometry", "field_block")
        fieldArea = self.params.get("census.geometry", "field_area")
        geometryDataSrc = driverShapefile.Open(filename, 0)
        geometryLayer = geometryDataSrc.GetLayer()
        geometryDict = dict([(feature[fieldBlock], feature[fieldArea],) for feature in geometryLayer])

        filename = self.params.get("census.population", "shapefile")
        fieldBlock = self.params.get("census.population", "field_block")
        fieldPopulation = self.params.get("census.population", "field_population")
        populationDataSrc = driverShapefile.Open(filename, 0)
        populationLayer = populationDataSrc.GetLayer()
        populationLayerDefn = populationLayer.GetLayerDefn()

        driverMemory = ogr.GetDriverByName("MEMORY")
        popDensityDataSrc = driverMemory.CreateDataSource("temp")
        driverMemory.Open("temp", gdal.GA_Update)
        popDensityLayer = popDensityDataSrc.CreateLayer("temp", populationLayer.GetSpatialRef(), populationLayerDefn.GetGeomType())
        popDensityLayerDefn = popDensityLayer.GetLayerDefn()

        fieldId = ogr.FieldDefn("blockid", ogr.OFTString)
        fieldId.SetWidth(15)
        popDensityLayer.CreateField(fieldId)
        popDensityLayer.CreateField(ogr.FieldDefn("area", ogr.OFTReal))
        popDensityLayer.CreateField(ogr.FieldDefn("population", ogr.OFTReal))
        popDensityLayer.CreateField(ogr.FieldDefn("density", ogr.OFTReal))

        # Extract, merge, and compute density
        for populationFeature in populationLayer:
            tempFeature = ogr.Feature(popDensityLayerDefn)
            tempFeature.SetGeometry(populationFeature.GetGeometryRef().Clone())

            blockId = populationFeature[fieldBlock]
            population = float(populationFeature[fieldPopulation])
            if not blockId in geometryDict:
                raise ValueError("Could not find census block '%d' in geometry dictionary." % blockId)
            area = max(float(geometryDict[blockId]), 1.0)
            density = population / area * 1.0e+6

            tempFeature.SetField("blockid", blockId)
            tempFeature.SetField("area", area)
            tempFeature.SetField("population", population)
            tempFeature.SetField("density", density)
            popDensityLayer.CreateFeature(tempFeature)
            del tempFeature
        return popDensityDataSrc

    def _rasterize(self, popDensityDataSrc):
        NODATA_VALUE = -999

        cellSizeDeg = self.params.getfloat("raster", "resolution_deg")

        popDensityLayer = popDensityDataSrc.GetLayer()
        xMin, xMax, yMin, yMax = popDensityLayer.GetExtent()

        # Create the destination data source
        numX = int(1 + (xMax - xMin) / cellSizeDeg)
        numY = int(1 + (yMax - yMin) / cellSizeDeg)
        rasterDataSrc = gdal.GetDriverByName('GTiff').Create("populationdensity.tiff", numX, numY, gdal.GA_Update, gdal.GDT_Float32)
        rasterDataSrc.SetGeoTransform((xMin, cellSizeDeg, 0, yMax, 0, -cellSizeDeg))
        rasterDataSrc.SetProjection(popDensityLayer.GetSpatialRef().ExportToWkt())
        band = rasterDataSrc.GetRasterBand(1)
        band.SetNoDataValue(NODATA_VALUE)

        # Rasterize
        gdal.RasterizeLayer(rasterDataSrc, [1], popDensityLayer, burn_values=[0], options=["ATTRIBUTE=density"])
        return
    
    def _parseCommandLine(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config")
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--generate", action="store_true", dest="generate")
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug")
        return parser.parse_args()


# ======================================================================
if __name__ == "__main__":
    PopulationDensityApp().main()


# End of file
