import qgis.core
import qgis.gui
import PyQt4.QtCore
import PyQt4.QtGui

from osgeo import gdal, osr
gdal.UseExceptions()


standalone = True

if standalone:
    app = qgis.core.QgsApplication([], True)
    #qgis.core.QgsApplication.setPrefixPath("/Applications/QGIS.app/Contents/MacOS", True)
    app.initQgis()

values = ["warning_time", "mmi_obs", "mmi_pred", "population_density",]
renderValues = ["mmi_pred", "population_density", "basemap"]

layers = {}
for value in values:
    filename = "./data/nc72923380/analysis_data-{}.vrt".format(value)
    layer = qgis.core.QgsRasterLayer(filename)
    if not layer.isValid():
        print("Could not load raster from '{}'".format(filename))
    else:
        layer.setName(value)
        layers[value] = layer
        
layers["mmi_obs"].setDrawingStyle("SingleBandPsuedoColor")
layers["mmi_obs"].loadNamedStyle("mmi_style.qml")

layers["mmi_pred"].setDrawingStyle("SingleBandPsuedoColor")
layers["mmi_pred"].loadNamedStyle("mmi_style.qml")

layers["population_density"].setDrawingStyle("SingleBandGray")
layers["population_density"].loadNamedStyle("pop_style.qml")

layers["basemap"] = qgis.core.QgsRasterLayer("esri_streetmap.xml")
if not layers["basemap"].isValid():
    print("Could not load ESRI streetmap raster.")
layers["basemap"].setDrawingStyle('MultiBandSingleBandGray')
layers["basemap"].renderer().setGrayBand(2)

layerRegistry = qgis.core.QgsMapLayerRegistry.instance()
layerRegistry.addMapLayers([layers[value] for value in values+["basemap"]])

def contours_from_raster(filename, contourStart=1.5, contourInterval=0.5):
    src = gdal.Open(filename, gdal.GA_ReadOnly)

    dest = gdal.GetDriverByName("MEM").Create("temp", numLon, numLat, nbands, gdal.GDT_Float32)

    driverMemory = ogr.GetDriverByName("MEMORY")
    contourDataSrc = driverMemory.CreateDataSource("temp")
    driverMemory.Open("temp", gdal.GA_Update)

    contourLayer = contourDataSrc.CreateLayer("temp", src.GetSpatialRef(), populationLayerDefn.GetGeomType())
    contourLayerDefn = popDensityLayer.GetLayerDefn()

    contourLayer.CreateField(ogr.FieldDefn("id", ogr.OFInteger))
    contourLayer.CreateField(ogr.FieldDefn("mmi_obs", ogr.OFTReal))
    gdal.ContourGenerate(src.GetRasterBand(1), contourInterval, contourStart, 0, [], 0, 0, contourLayer, 0, 1)
    del contourDataSrc

    # Create QgsVectorLayer (memory) and transfer geometry
    
    
    return contourLayer

# contour mmi_obs
#mmiObsContours = contours_from_raster("./data//nc72923380/analysis_data-mmi_obs.vrt", 1.5, 0.5)



canvas = qgis.gui.QgsMapCanvas()
canvas.setCanvasColor(PyQt4.QtCore.Qt.white)
canvas.enableAntiAliasing(True)
canvas.setExtent(layers[renderValues[0]].extent())
canvas.setLayerSet([qgis.gui.QgsMapCanvasLayer(layers[value]) for value in renderValues])

if standalone:
    canvas.window().resize(1920, 1080)
    
renderer = canvas.mapRenderer()
renderer.setDestinationCrs(qgis.core.QgsCoordinateReferenceSystem("EPSG:3311"))
renderer.setProjectionsEnabled(True)
#canvas.show()

image = PyQt4.QtGui.QImage(PyQt4.QtCore.QSize(1920, 1080), PyQt4.QtGui.QImage.Format_ARGB32_Premultiplied)
image.fill(PyQt4.QtCore.Qt.white)

painter = PyQt4.QtGui.QPainter()
painter.begin(image)
painter.setRenderHint(PyQt4.QtGui.QPainter.Antialiasing)


# set output size
renderer.setOutputSize(image.size(), image.logicalDpiX())

renderer.render(painter)
painter.end()
image.save("test.png","png")
