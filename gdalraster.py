# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import numpy
from lxml import etree

from osgeo import gdal, ogr, osr

gdal.UseExceptions()

NO_DATA_VALUE = -999.0

def resample(filename, destNumX, destNumY, destSRS, destGeoTransform):
    """Resample and crop raster to match specified grid.

    :type filename: str
    :param filename: 
        Filename with raster

    :type destNumX: int
    :param destNumX: Number of points in destination grid in x direction.

    :type destNumY: int
    :param destNumY: Number of points in destination grid in y direction.

    :type destSRS: OSR SpatialReference
    :param destSRS: Spatial reference in destination grid.

    :type destGeoTransform: GDAL GeoTransform
    :param destGeoTransform: Geometric transformation of destination grid.
    """
    src = gdal.Open(filename, gdal.GA_ReadOnly)
    nbands = src.RasterCount

    dest = gdal.GetDriverByName("MEM").Create("temp", destNumX, destNumY, nbands, gdal.GDT_Float32)
    dest.SetGeoTransform(destGeoTransform)
    dest.SetProjection(destSRS.ExportToWkt())
    gdal.ReprojectImage(src, dest, src.GetProjection(), dest.GetProjection(), gdal.GRA_Bilinear)
    dest.FlushCache()

    del src
    
    return numpy.array(dest.GetRasterBand(1).ReadAsArray()).ravel()


def write(filename, values, numX, numY, spatialRef, geoTransform):
    """Write values to GeoTiff raster file with given spatial reference and geometric transformation.

    :type filename: str
    :param filename:
        Filename for GeoTiff raster.

    :type values: List
    :param values:
        List of tuples with value name and Numpy array.

    :type spatialRef: OSR SpatialReference
    :param spatialRef: Spatial reference associated with values.

    :type geoTranform: GDAL GeoTransform
    :param geoTransform: Geometric transformation associated with values.
    """
    nbands = len(values)
    
    dest = gdal.GetDriverByName("GTiff").Create(filename, numX, numY, nbands, gdal.GDT_Float32)
    dest.SetGeoTransform(geoTransform)
    dest.SetProjection(spatialRef.ExportToWkt())

    for ivalue,(name,value,) in enumerate(values):
        band = dest.GetRasterBand(1+ivalue)
        band.SetDescription(name)
        band.SetNoDataValue(NO_DATA_VALUE)
        band.WriteArray(value.reshape((numY,numX,)))
        band.FlushCache()
    dest.FlushCache()
    return


def clone_new_data(name, values, src, noDataValue=None):
    """Create GDAL raster in memory with values as data with grid specs from src raster.

    :type name: str
    :param name:
        Name for raster layer (description).

    :type values: Numpy array
    :param values: Array of values.

    :type src: GDAL raster
    :param src: GDAL raster with from which values are derived used to specific grid information.

    """
    NBANDS = 1
    
    dest = gdal.GetDriverByName("MEM").Create("", src.RasterXSize, src.RasterYSize, NBANDS, gdal.GDT_Float32)
    dest.SetGeoTransform(src.GetGeoTransform())
    dest.SetProjection(src.GetProjection())

    band = dest.GetRasterBand(1)
    band.SetDescription(name)
    if not noDataValue is None:
        band.SetNoDataValue(noDataValue)
    band.WriteArray(values)
    band.FlushCache()
    dest.FlushCache()
    return dest


def contours_from_raster(raster, bandName, cstart=0.0, cinterval=10.0, clevels=[]):
    """Create GDAL raster in memory with values as data with grid specs from src raster.

    :type raster: GDAL raster
    :param raster: 
        GDAL raster.

    :type band: int
    :param band:
        Band number to use for contours. =1 for first band.

    :type cstart: float
    :param cstart:
        Starting value for contours

    :type cinterval: float
    :param cinterval:
        Contour interval.

    :type levels: list
    :param levels:
        List of contour values. cstart and cstep are ignored if levels is given.
    """
    cband = None
    for iband in range(raster.RasterCount):
        band = raster.GetRasterBand(1+iband)
        if band.GetDescription() == bandName:
            cband = band
    if not cband:
        filename = raster.GetFileList()[0]
        raise ValueError("Could not find band '{band}' in raster '{filename}'.".format(band=bandName, filename=filename))
    
    driverMemory = ogr.GetDriverByName("MEMORY")
    contours = driverMemory.CreateDataSource("temp")
    driverMemory.Open("temp", gdal.GA_Update)

    spatialRef = osr.SpatialReference()
    spatialRef.ImportFromWkt(raster.GetProjectionRef())
    
    ogrLayer = contours.CreateLayer(bandName+"_contour", spatialRef, ogr.wkbLineString25D)
    ogrLayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    ogrLayer.CreateField(ogr.FieldDefn("value", ogr.OFTReal))

    gdal.ContourGenerate(cband, cinterval, cstart, clevels, 1, NO_DATA_VALUE, ogrLayer, 0, 1)

    return contours


# End of file
