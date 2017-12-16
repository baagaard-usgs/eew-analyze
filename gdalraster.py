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
    """Resample/crop raster to match specified grid.

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
    """Write Get shaking time based on origin time and shear wave speed.

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

    write_bands_virtual(filename, dest)
    return

def write_bands_virtual(filename, raster):

    dataset = etree.Element("VRTDataset")
    dataset.set("rasterXSize", "{:d}".format(raster.RasterXSize))
    dataset.set("rasterYSize", "{:d}".format(raster.RasterYSize))

    spatialref = etree.SubElement(dataset, "SRS")
    spatialref.text = raster.GetProjection()
    
    geotrans = etree.SubElement(dataset, "GeoTransform")
    geotrans.text = ", ".join(map(str, raster.GetGeoTransform()))

    band = etree.SubElement(dataset, "VRTRasterBand")
    band.set("dataType", "Float32")
    band.set("band", "1")

    colorinterp = etree.SubElement(band, "ColorInterp")
    colorinterp.text = "Gray"

    src = etree.SubElement(band, "SimpleSource")

    srcFilename = etree.SubElement(src, "SourceFilename")
    srcFilename.set("relativeToVRT", "1")
    srcFilename.text = os.path.split(filename)[1]

    srcBand = etree.SubElement(src, "SourceBand")

    srcRect = etree.SubElement(src, "SrcRect")
    srcRect.set("xOff", "0")
    srcRect.set("yOff", "0")
    srcRect.set("xSize", "{:d}".format(raster.RasterXSize))
    srcRect.set("ySize", "{:d}".format(raster.RasterYSize))
    
    dstRect = etree.SubElement(src, "DstRect")
    dstRect.set("xOff", "0")
    dstRect.set("yOff", "0")
    dstRect.set("xSize", "{:d}".format(raster.RasterXSize))
    dstRect.set("ySize", "{:d}".format(raster.RasterYSize))

    
    for iband in range(raster.RasterCount):
        srcBand.text = "{:d}".format(1+iband)
        description = raster.GetRasterBand(1+iband).GetDescription()

        tree = etree.ElementTree(dataset)
        fileroot = filename[:filename.rfind(".")]
        vfilename = "{0}-{1}.vrt".format(fileroot, description)
        tree.write(vfilename)

# End of file
