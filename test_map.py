import qgis.core
import qgis.gui
import PyQt4.QtCore
import PyQt4.QtGui

standalone = True

if standalone:
    #qgis.core.QgsApplication.setPrefixPath("/usr", True)
    app = qgis.core.QgsApplication([], True)
    app.initQgis()

values = ["warning_time", "mmi_obs", "mmi_pred", "population_density",]
renderValues = ["mmi_obs", "population_density"]

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

layerRegistry = qgis.core.QgsMapLayerRegistry.instance()
layerRegistry.addMapLayers([layers[value] for value in values])

#    gdal.ContourGenerate(ds.GetRasterBand(1), 10, 0, [], 0, 0, ogr_lyr, 0, 1)





#canvas = qgis.gui.QgsMapCanvas()
#canvas.setCanvasColor(PyQt4.QtCore.Qt.white)
#canvas.enableAntiAliasing(True)
#canvas.setExtent(layers[renderLayers[0].extent())
#canvas.setLayerSet([qgis.gui.QgsMapCanvasLayer(layers[value]) for value in renderValues])
#
#if standalone:
#    canvas.window().resize(1920, 1080)
#    
#renderer = canvas.mapRenderer()
#renderer.setProjectionsEnabled(True)
#renderer.setDestinationCrs(qgis.core.QgsCoordinateReferenceSystem("EPSG:3311"))
#canvas.show()
    
image = PyQt4.QtGui.QImage(PyQt4.QtCore.QSize(1920, 1080), PyQt4.QtGui.QImage.Format_ARGB32_Premultiplied)
image.fill(PyQt4.QtCore.Qt.white)

painter = PyQt4.QtGui.QPainter()
painter.begin(image)
painter.setRenderHint(PyQt4.QtGui.QPainter.Antialiasing)

renderer = qgis.core.QgsMapRenderer()

renderer.setLayerSet([layers[value].id() for value in renderValues])
renderer.setProjectionsEnabled(True)
renderer.setDestinationCrs(qgis.core.QgsCoordinateReferenceSystem("EPSG:3311"))

# set extent
rect = qgis.core.QgsRectangle(renderer.fullExtent())
rect.scale(1.1)
renderer.setExtent(rect)

# set output size
renderer.setOutputSize(image.size(), image.logicalDpiX())

renderer.render(painter)
painter.end()
image.save("test.png","png")
