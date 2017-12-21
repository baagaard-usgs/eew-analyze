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

def gdal_to_qgisraster(gdalRaster):
    """Convert GDAL raster layer to QGIS raster layer.

    QGIS 2.x does not support setting raster values, so we write to a
    GeoTiff file using GDAL and then read it using QGIS.
    """
    numLon = ":TODO:"
    numLat = ":TODO:"
    dLon = ":TODO:"
    dLat = ":TODO:"
    originLon = ":TODO:"
    originLat = ":TODO:"
    nbands = 1
    

    # Write to temporary file
    filename = ":TODO:"
    temp = gdal.GetDriverByName("Tiff").Create(filename, numLon, numLat, nbands, gdal.GDT_Float32)
    temp.SetGeoTransform((originLon, dLon, 0, originLat, 0, dLat,))
    temp.SetProjection(gdalRaster.GetProjection())
    band = temp.GetRasterBand(1)
    band.SetDescription(gdalRaster)
    band.WriteArray(gdalRaster.GetRasterBand(1).ReadFromArray())
    band.FlushCache()
    temp.FlushCache()

    qgisRaster = qgis.core.QgsRasterLayer(filename)
    if not qgisRaster.isValid():
        raise IOError("Could not read raster layer from '{}'.".format(filename))
    return qgisRaster

def ogr_to_qgisvector(ogrVector):
    """Convert OGR vector layer to QGIS vector layer.

    QGIS 2.x does not support constructing a QGIS layer directory from
    an OGR layer, so we create the QGIS a vector layer in memory and
    transfer the features.
    """
    return
                    

def extract_raster_bands(gdalRaster):
    """Extract a GDAL raster layer into a QGIS raster layer for each band.

    QGIS 2.x does not support setting raster values, so we write a
    virtual raster layer XML file and then read it using QGIS.
    """
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

    return bands
        
# End of file
