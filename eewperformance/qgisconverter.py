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

import qgis.core
import PyQt4.QtCore
from osgeo import gdal, osr, ogr


gdal.UseExceptions()

def gdal_to_qgisraster(gdalRaster, filename="raster.tiff"):
    """Convert GDAL raster layer to QGIS raster layer.

    QGIS 2.x does not support setting raster values, so we write to a
    GeoTiff file using GDAL and then read it using QGIS.

    :type gdalRaster: GDAL raster layer
    :param gdalRaster:
        GDAL raster layer

    :type filename: str
    :param filename:
        Filename for temporary file.

    :returns: QGIS raster layer.
    """
    nbands = gdalRaster.RasterCount
    temp = gdal.GetDriverByName("GTiff").Create(filename, gdalRaster.RasterXSize, gdalRaster.RasterYSize, nbands, gdal.GDT_Float32)
    temp.SetGeoTransform(gdalRaster.GetGeoTransform())
    temp.SetProjection(gdalRaster.GetProjection())
    for iband in range(gdalRaster.RasterCount):
        gdalBand = gdalRaster.GetRasterBand(1+iband)
        tempBand = temp.GetRasterBand(1+iband)
        tempBand.SetDescription(gdalBand.GetDescription())
        if not gdalBand.GetNoDataValue() is None:
            tempBand.SetNoDataValue(gdalBand.GetNoDataValue())
        tempBand.WriteArray(gdalBand.ReadAsArray())
        tempBand.FlushCache()
    temp.FlushCache()

    qgisRaster = qgis.core.QgsRasterLayer(filename)
    if not qgisRaster.isValid():
        raise IOError("Could not read raster layer from '{}'.".format(filename))
    return qgisRaster


def ogr_to_qgisvector(ogrVector, filename="vector.shp"):
    """Convert OGR vector layer to QGIS vector layer.

    QGIS 2.x does not support constructing a QGIS layer directory from
    an OGR layer, so we write to an ESRI Shapefile and then read it using QGIS.

    :type ogrVector: OGR vector data source.
    :param ogrVector:
        OGR vector data source.

    :type filename: str
    :param filename:
        Filename for temporary file.

    :returns: QGIS vector layer.
    """
    ogrLayer = ogrVector.GetLayer()
    name = ogrLayer.GetDescription()

    dest = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(filename)
    destLayer = dest.CreateLayer("temp", ogrLayer.GetSpatialRef(), ogrLayer.GetGeomType())

    ogrLayerDefn = ogrLayer.GetLayerDefn()
    for ifield in range(ogrLayerDefn.GetFieldCount()):
        destLayer.CreateField(ogrLayerDefn.GetFieldDefn(ifield))

    # Add features to the dest layer
    destLayerDefn = destLayer.GetLayerDefn()
    for ifeature in range(ogrLayer.GetFeatureCount()):
        ogrFeature = ogrLayer.GetFeature(ifeature)
        destFeature = ogr.Feature(destLayerDefn)
        for ifield in range(destLayerDefn.GetFieldCount()):
            destFeature.SetField(destLayerDefn.GetFieldDefn(ifield).GetNameRef(), ogrFeature.GetField(ifield))
        destFeature.SetGeometry(ogrFeature.GetGeometryRef())
        destLayer.CreateFeature(destFeature)
        del destFeature
    dest.FlushCache()
    del destLayer
    del dest

    # Create QgsVectorLayer.
    qgisVector = qgis.core.QgsVectorLayer(filename, name, "ogr")
    if not qgisVector.isValid():
        raise IOError("Could not read vector layer from '{}'.".format(filename))
    return qgisVector
                    

def extract_raster_bands(gdalRaster):
    """Extract a GDAL raster layer into a QGIS raster layer for each band.

    QGIS 2.x does not support setting raster values, so we write a
    virtual raster layer XML file and then read it using QGIS.

    :type gdalRaster: GDAL raster layer
    :param gdalRaster:
        GDAL raster layer

    :returns: Dictionary (layer name -> QGIS raster layer).
    """
    filenameGDAL = gdalRaster.GetFileList()[0]

    dataset = etree.Element("VRTDataset")
    dataset.set("rasterXSize", "{:d}".format(gdalRaster.RasterXSize))
    dataset.set("rasterYSize", "{:d}".format(gdalRaster.RasterYSize))

    spatialref = etree.SubElement(dataset, "SRS")
    spatialref.text = gdalRaster.GetProjection()
    
    geotrans = etree.SubElement(dataset, "GeoTransform")
    geotrans.text = ", ".join(map(str, gdalRaster.GetGeoTransform()))

    band = etree.SubElement(dataset, "VRTRasterBand")
    band.set("dataType", "Float32")
    band.set("band", "1")

    colorinterp = etree.SubElement(band, "ColorInterp")
    colorinterp.text = "Gray"

    src = etree.SubElement(band, "SimpleSource")

    srcFilename = etree.SubElement(src, "SourceFilename")
    srcFilename.set("relativeToVRT", "1")
    srcFilename.text = os.path.split(filenameGDAL)[1]

    srcBand = etree.SubElement(src, "SourceBand")

    srcRect = etree.SubElement(src, "SrcRect")
    srcRect.set("xOff", "0")
    srcRect.set("yOff", "0")
    srcRect.set("xSize", "{:d}".format(gdalRaster.RasterXSize))
    srcRect.set("ySize", "{:d}".format(gdalRaster.RasterYSize))
    
    dstRect = etree.SubElement(src, "DstRect")
    dstRect.set("xOff", "0")
    dstRect.set("yOff", "0")
    dstRect.set("xSize", "{:d}".format(gdalRaster.RasterXSize))
    dstRect.set("ySize", "{:d}".format(gdalRaster.RasterYSize))

    bands = {}
    for iband in range(gdalRaster.RasterCount):
        srcBand.text = "{:d}".format(1+iband)
        description = gdalRaster.GetRasterBand(1+iband).GetDescription()

        tree = etree.ElementTree(dataset)
        fileroot = filenameGDAL[:filenameGDAL.rfind(".")]
        vfilename = "{0}-{1}.vrt".format(fileroot, description)
        tree.write(vfilename)

        layer = qgis.core.QgsRasterLayer(vfilename)
        bands[description] = layer

    return bands
        
# End of file
