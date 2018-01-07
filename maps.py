# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import numpy

import qgis.core
import qgis.gui
import PyQt4.QtCore
import PyQt4.QtGui
from osgeo import gdal, osr

import gdalraster
import qgisconverter

gdal.UseExceptions()

class MapPanels(object):
    """Maps of MMI, warning time, etc.
    """
    def __init__(self, config):
        """
        :type config: ConfigParser
        :param config: Configuration for application.
        """
        self.config = config        

        self.baseLayers = {}
        self.dataLayers = {}

        self._initQgis()
        return
    
    def __del__(self):
        """
        """
        self._exitQgis()
        return

    def load_basemap(self):
        if self.baseLayers.has_key("basemap"):
            return
        
        filename = self.config.get("maps", "basemap")
        basemap = qgis.core.QgsRasterLayer(filename)
        if not basemap.isValid():
            raise IOError("Could not load basemap from '{}'.".format(filename))
        basemap.loadNamedStyle("basemap.qml")
        self.baseLayers["basemap"] = basemap

        layerRegistry = qgis.core.QgsMapLayerRegistry.instance()
        layerRegistry.addMapLayers([layer for layer in self.baseLayers.values()])
        return

    def load_data(self, eqId, alertVersion=None):
        """
        :type eqId: str
        :param eqId: ComCat earthquake id.

        :type config: ConfigParser
        :param config: Configuration for application.

        :type filename: str
        :param config: Filename for analysis data in GDAL raster file.
        """
        self.eqId = eqId
        self.alertVersion = alertVersion

        layerRegistry = qgis.core.QgsMapLayerRegistry.instance()
        layerRegistry.removeMapLayers([layer for layer in self.dataLayers.values()])
        self.dataLayers = {}

        # Generate temporary virtual raster bands
        filename = self.config.get("files", "analysis_event").replace("[EVENTID]", self.eqId)
        src = gdal.Open(filename, gdal.GA_ReadOnly)
        for name,layer in qgisconverter.extract_raster_bands(src).items():
            self.dataLayers[name] = layer

        # Contours from temporary virtual raster bands.
        values = ("mmi_obs", "mmi_pred", "warning_time",) if self.alertVersion is None else ("mmi_pred", "warning_time",)
        for value in values:
            if value in ("mmi_obs", "mmi_pred",):
                cstart = 1.5
                cinterval = 0.5
                clevels = []
            else:
                cstart = 0
                cinterval = 2.0
                clevels = []
            cvalue = value+"_contour"
            ogrLayer = gdalraster.contours_from_raster(src, value, cstart, cinterval, clevels)
            fileSuffix = "_{}.shp".format(cvalue)
            filenameLayer = filename.replace(".tiff", fileSuffix)
            self.dataLayers[cvalue] = qgisconverter.ogr_to_qgisvector(ogrLayer, filenameLayer)

        if self.dataLayers.has_key("mmi_obs") and self.dataLayers.has_key("mmi_pred"):
            # MMI residual (observed - predicted)
            mmiObs = numpy.array(src.GetRasterBand(1).ReadAsArray())
            mmiPred = numpy.array(src.GetRasterBand(2).ReadAsArray())
            mmiResidual = mmiObs - mmiPred
            filenameResidual = filename.replace(".tiff", "-mmi_residual.tiff")
            gdalLayer = gdalraster.clone_new_data("mmi_residual", mmiResidual, src)
            self.dataLayers["mmi_residual"] = qgisconverter.gdal_to_qgisraster(gdalLayer, filenameResidual)
            del gdalLayer

        del src

        layerRegistry = qgis.core.QgsMapLayerRegistry.instance()
        layerRegistry.addMapLayers([layer for layer in self.dataLayers.values()])
        return # TEMPORARY
    
        # Epicenter
        options = "delimiter={delimiter}&crs=EPSG:4326&xField={x}&yField={}".format(delimiter="|", x="longitude", y="latitude")
        uri = "file://{filename}&{options}".format(filename, options)
        self.dataLayers["epicenter_obs"] = qgis.core.QgsVectorLayer(uri, "epicenter_obs", "delimitedtext")

        uri = "file://{filename}&{options}".format(filename, options)
        self.dataLayers["epicenter_pred"] = qgis.core.QgsVectorLayer(uri, "epicenter_pred", "delimitedtext")
        
        return    

    def mmi_observed(self):
        """Create map with observed MMI with contours.
        """
        basemap = self.baseLayers["basemap"]
        
        mmi = self.dataLayers["mmi_obs"]
        mmi.loadNamedStyle("mmi.qml")

        mmiContours = self.dataLayers["mmi_obs_contour"]
        self._set_contour_style(mmiContours)
        self._add_labels(mmiContours)
            
        #epicenter = self.dataLayers["epicenter_obs"]
        # :TODO: set style
        
        self._render([mmiContours, mmi, basemap], mmi.extent(), "mmi_obs", "png")
        return

    def mmi_predicted(self, showZeroWarningTime=False):
        """Create map with observed MMI with contours.
        """
        basemap = self.baseLayers["basemap"]
        
        mmi = self.dataLayers["mmi_pred"]
        mmi.loadNamedStyle("mmi.qml")

        mmiContours = self.dataLayers["mmi_pred_contour"]
        self._set_contour_style(mmiContours)
        self._add_labels(mmiContours)

        if showZeroWarningTime:
            # :TODO: Add red contour showing zero warning time
            pass
        
        #epicenter = self.dataLayers["epicenter_obs"]
        # :TODO: set style
        
        self._render([mmiContours, mmi, basemap], mmi.extent(), "mmi_pred", "png")
        return

    def mmi_residual(self):
        basemap = self.baseLayers["basemap"]
        
        mmiResidual = self.dataLayers["mmi_residual"]
        mmiResidual.loadNamedStyle("mmi_residual.qml")

        #epicenter = self.dataLayers["epicenter_obs"]
        # :TODO: set style
        
        self._render([mmiResidual, basemap], mmiResidual.extent(), "mmi_residual", "png")
        return

    def mmi_warning_time(self):
        """Create map with observed MMI with contours.
        """
        basemap = self.baseLayers["basemap"]

        if self.alertVersion is None:
            mmi = self.dataLayers["mmi_obs"]
        else:
            mmi = self.dataLayers["mmi_pred"]
        mmi.loadNamedStyle("mmi.qml")

        warningContours = self.dataLayers["warning_time_contour"]
        self._set_contour_style(warningContours)
        self._add_labels(warningContours)
            
        #epicenter = self.dataLayers["epicenter_obs"]
        # :TODO: set style
        
        self._render([warningContours, mmi, basemap], mmi.extent(), "mmi_warning", "png")
        return

    def alert_categories(self):
        return

    def _add_labels(self, layer):
        labels = qgis.core.QgsPalLayerSettings()
        labels.readFromLayer(layer)
        labels.enabled = True
        labels.fieldName = "mmi_obs"
        labels.fontSizeInMapUnits = False
        labels.fontLimitPixelSize = True
        labels.fontMinPixelSize = 6
        labels.textFont.setFamily("Arial Narrow")
        labels.textFont.setPointSize(8)
        labels.textFont.setWeight(labels.textFont.Normal)
        labels.textColor = PyQt4.QtGui.QColor("#000000")
        labels.placement = qgis.core.QgsPalLayerSettings.AboveLine
        labels.writeToLayer(layer)
        return

    def _set_contour_style(self, layer):
        symbolLayer = layer.rendererV2().symbol().symbolLayer(0)
        symbolLayer.setColor(PyQt4.QtGui.QColor("#000000")) # works
        symbolLayer.setWidth(2.0)
        symbolLayer.setWidthUnit(qgis.core.QgsSymbolV2.Pixel)
        return
    
    def _render(self, layers, extent, name, format):
        """Render image to file in given format.
        """
        imageWidth = self.config.getint("maps", "width_pixels")
        imageHeight = self.config.getint("maps", "height_pixels")
        canvas = qgis.gui.QgsMapCanvas()
        canvas.setCanvasColor(PyQt4.QtCore.Qt.white)
        canvas.enableAntiAliasing(True)
        canvas.setExtent(extent)
        canvas.setLayerSet([qgis.gui.QgsMapCanvasLayer(layer) for layer in layers])
        canvas.window().resize(imageWidth, imageHeight)
        renderer = canvas.mapRenderer()
        
        projection = self.config.get("maps", "projection")
        renderer.setDestinationCrs(qgis.core.QgsCoordinateReferenceSystem(projection))
        renderer.setProjectionsEnabled(True)

        bgColor = self.config.get("maps", "bg_color")

        composer = qgis.core.QgsComposition(renderer)
        composer.setPlotStyle(qgis.core.QgsComposition.Print)
        dpmm = 150 / 25.4
        composer.setPaperSize(imageWidth/dpmm, imageHeight/dpmm)
        composer.setPrintResolution(dpmm*25.4)
        
        x, y = 0,0
        w, h = composer.paperWidth(), composer.paperHeight()
        composerMap = qgis.core.QgsComposerMap(composer, x, y, w, h)
        composer.addItem(composerMap)

        item = qgis.core.QgsComposerScaleBar(composer)
        item.setStyle('Single box') # optionally modify the style
        item.setComposerMap(composerMap)
        item.applyDefaultSize()
        composer.addItem(item)

        # create output image and initialize it
        image = PyQt4.QtGui.QImage(PyQt4.QtCore.QSize(imageWidth, imageHeight), PyQt4.QtGui.QImage.Format_ARGB32)
        image.setDotsPerMeterX(dpmm * 1000)
        image.setDotsPerMeterY(dpmm * 1000)
        image.fill(PyQt4.QtGui.QColor(bgColor))

        # render the composition
        painter = PyQt4.QtGui.QPainter(image)
        composer.renderPage(painter, 0)
        painter.end()

        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        if self.alertVersion is None:
            filename = "{eq}-map_{name}.{format}".format(eq=self.eqId, name=name, format=format)
        else:
            filename = "{eq}-alert{ver:03d}_map_{name}.{format}".format(eq=self.eqId, ver=self.alertVersion, name=name, format=format)
        image.save(os.path.join(plotsDir, filename), format)
        return

    def _initQgis(self):
        self.qgis = qgis.core.QgsApplication([], True)
        qgisPrefixPath = self.config.get("qgis", "prefix_path")
        if qgisPrefixPath != "None":
            # "/Applications/QGIS.app/Contents/MacOS"
            qgis.core.QgsApplication.setPrefixPath(qgisPrefixPath, True)
        self.qgis.initQgis()
        return

    def _exitQgis(self):
        TEMPORARY_LAYERS = [
            "mmi_obs",
            "mmi_obs_contours",
            "mmi_pred",
            "mmi_pred_contours",
            "warning_time",
            "population_density",
        ]

        import subprocess
        for name in TEMPORARY_LAYERS:
            if name in self.dataLayers:
                filename = self.dataLayers[name].source()
                subprocess.call("rm {}".format(filename), shell=True)
        self.dataLayers = {}
                
        if self.qgis:
            self.qgis.exitQgis()

        return
