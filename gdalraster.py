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

from osgeo import gdal, osr

gdal.UseExceptions()

def resample(filename, gridSpecs):
    """Resample and crop raster to match specified grid.

    The grid specs come from ShakeMap which uses a vertex-based
    grid. We are writing a cell-base grid, so the ShakeMap vertices
    become the cell centers. The grid origin is 1/2 cell in lon/lat
    from the cell center.

    :type filename: str
    :param filename: 
        Filename with raster

    :type gridSpecs: dict
    :param gridSpecs: 
        Grid specifications to match.
    """
    src = gdal.Open(filename, gdal.GA_ReadOnly)

    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    numLon = gridSpecs["num_longitude"]
    numLat = gridSpecs["num_latitude"]
    dLon = (gridSpecs["longitude_max"]-gridSpecs["longitude_min"]) / gridSpecs["num_longitude"]
    dLat = -(gridSpecs["latitude_max"]-gridSpecs["latitude_min"]) / gridSpecs["num_latitude"]
    originLon = gridSpecs["longitude_min"] - 0.5*dLon
    originLat = gridSpecs["latitude_max"] + 0.5*dLat
    nbands = 1

    dest = gdal.GetDriverByName("MEM").Create("temp", numLon, numLat, nbands, gdal.GDT_Float32)
    dest.SetGeoTransform((originLon, dLon, 0, originLat, 0, dLat,))
    dest.SetProjection(wgs84.ExportToWkt())
    gdal.ReprojectImage(src, dest, src.GetProjection(), dest.GetProjection(), gdal.GRA_Bilinear)
    dest.FlushCache()
    
    return numpy.array(dest.GetRasterBand(1).ReadAsArray()).ravel()


def write(filename, values, gridSpecs):
    """Write values to GeoTiff raster file with grid specs from dictionary.

    The grid specs come from ShakeMap which uses a vertex-based
    grid. We are writing a cell-base grid, so the ShakeMap vertices
    become the cell centers. The grid origin is 1/2 cell in lon/lat
    from the cell center.

    :type filename: str
    :param filename:
        Filename for GeoTiff raster.

    :type values: List
    :param values:
        List of tuples with value name and Numpy array.

    :type gridSpects: dict
    :param gridSpecs: 
        Grid specification.

    """
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    numLon = gridSpecs["num_longitude"]
    numLat = gridSpecs["num_latitude"]
    dLon = (gridSpecs["longitude_max"]-gridSpecs["longitude_min"]) / (gridSpecs["num_longitude"]-1)
    dLat = -(gridSpecs["latitude_max"]-gridSpecs["latitude_min"]) / (gridSpecs["num_latitude"]-1)
    originLon = gridSpecs["longitude_min"] - 0.5*dLon
    originLat = gridSpecs["latitude_max"] + 0.5*dLat
    nbands = len(values)
    
    dest = gdal.GetDriverByName("GTiff").Create(filename, numLon, numLat, nbands, gdal.GDT_Float32)
    dest.SetGeoTransform((originLon, dLon, 0, originLat, 0, dLat,))
    dest.SetProjection(wgs84.ExportToWkt())

    for ivalue,(name,value,) in enumerate(values):
        band = dest.GetRasterBand(1+ivalue)
        band.SetDescription(name)
        band.WriteArray(value.reshape((numLat,numLon,)))
        band.FlushCache()
    dest.FlushCache()
    return


def clone(filename, name, values, src):
    """Write values to GeoTiff raster file with grid specs from src raster.

    :type filename: str
    :param filename:
        Filename for GeoTiff raster.

    :type values: Numpy array
    :param values:
        Array of values.

    :type src: GDAL rastter
    :param src: 
        GDAL raster with from which values are derived used to specific grid information.

    """
    numLon = ":TODO:"
    numLat = ":TODO:"
    nbands = 1
    
    dest = gdal.GetDriverByName("GTiff").Create(filename, numLon, numLat, nbands, gdal.GDT_Float32)
    dest.SetGeoTransform(src.GetGeoTransform())
    dest.SetProjection(src.GetProjection())

    band = dest.GetRasterBand(1)
    band.SetDescription(name)
    band.WriteArray(values)
    band.FlushCache()
    dest.FlushCache()
    return

# End of file
