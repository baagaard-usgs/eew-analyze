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

        self.fontFamily = "Arial Narrow"
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
        layerRegistry.addMapLayers([layer for layer in self.baseLayers.values()], False)
        return

    def load_data(self, eqId, alert=None):
        """
        :type eqId: str
        :param eqId: ComCat earthquake id.

        :type config: ConfigParser
        :param config: Configuration for application.

        :type filename: str
        :param config: Filename for analysis data in GDAL raster file.

        :type alert: dict
        :param alert: ShakeAlert alert information from analysis database.
        """
        self.eqId = eqId
        self.alert = alert

        layerRegistry = qgis.core.QgsMapLayerRegistry.instance()
        layerRegistry.removeMapLayers([layer for layer in self.dataLayers.values()])
        self.dataLayers = {}

        # Generate temporary virtual raster bands
        filename = self.config.get("files", "analysis_event").replace("[EVENTID]", self.eqId)
        src = gdal.Open(filename, gdal.GA_ReadOnly)
        for name,layer in qgisconverter.extract_raster_bands(src).items():
            self.dataLayers[name] = layer

        # Contours from temporary virtual raster bands.
        values = ("mmi_obs", "mmi_pred", "warning_time",) if self.alert is None else ("mmi_pred", "warning_time",)
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
            noDataValue = src.GetRasterBand(2).GetNoDataValue()
            mask = mmiPred != noDataValue
            mmiResidual = mask*(mmiObs - mmiPred) + ~mask*noDataValue
            filenameResidual = filename.replace(".tiff", "-mmi_residual.tiff")
            gdalLayer = gdalraster.clone_new_data("mmi_residual", mmiResidual, src, noDataValue=None)
            self.dataLayers["mmi_residual"] = qgisconverter.gdal_to_qgisraster(gdalLayer, filenameResidual)
            del gdalLayer

        del src

        layerRegistry.addMapLayers([layer for layer in self.dataLayers.values()], False)
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
        
        self._render([mmiContours, mmi, basemap], mmi.extent(), "mmi_obs", "jpg", legendLayers=[mmi], legendTitle="MMI", title="Observed Shaking")
        return

    def mmi_predicted(self):
        """Create map with observed MMI with contours.
        """
        basemap = self.baseLayers["basemap"]
        
        mmi = self.dataLayers["mmi_pred"]
        mmi.loadNamedStyle("mmi.qml")

        mmiContours = self.dataLayers["mmi_pred_contour"]
        self._set_contour_style(mmiContours)
        self._add_labels(mmiContours)

        #epicenter = self.dataLayers["epicenter_obs"]
        # :TODO: set style
        
        self._render([mmiContours, mmi, basemap], mmi.extent(), "mmi_pred", "jpg", legendLayers=[mmi], legendTitle="MMI", title="ShakeAlert Predicted Shaking")
        return

    def mmi_residual(self):
        basemap = self.baseLayers["basemap"]
        
        mmiResidual = self.dataLayers["mmi_residual"]
        mmiResidual.loadNamedStyle("mmi_residual.qml")

        #epicenter = self.dataLayers["epicenter_obs"]
        # :TODO: set style
        
        self._render([mmiResidual, basemap], mmiResidual.extent(), "mmi_residual", "jpg", legendLayers=[mmiResidual], legendTitle="Residual (Obs - Pred)", title="MMI Residual")
        return

    def mmi_warning_time(self, t=None):
        """Create map with observed MMI with contours.
        """
        basemap = self.baseLayers["basemap"]

        if self.alert is None:
            mmi = self.dataLayers["mmi_obs"]
        else:
            mmi = self.dataLayers["mmi_pred"]
        mmi.loadNamedStyle("mmi.qml")

        warningContours = self.dataLayers["warning_time_contour"]
        self._set_warning_time_style(warningContours)
        self._add_labels(warningContours)
            
        #epicenter = self.dataLayers["epicenter_obs"]
        # :TODO: set style

        if self.alert:
            import dateutil.parser
            tstamp = dateutil.parser.parse(self.alert["timestamp"])
            title = "{tstamp:%Y-%m-%d %H:%M:%S.%f} ({t:6.3f}s after origin time), Alert {alert[version]:3d}, M{alert[magnitude]:4.2f}".format(alert=self.alert, tstamp=tstamp, t=t)
        else:
            title = "Warning Time (s) and Observed Shaking (MMI)"
        
        self._render([warningContours, mmi, basemap], mmi.extent(), "mmi_warning", "jpg", legendLayers=[mmi], legendTitle="MMI", title=title)
        return

    def alert_categories(self):
        return

    def _add_labels(self, layer):
        labels = qgis.core.QgsPalLayerSettings()
        labels.readFromLayer(layer)
        labels.enabled = True
        labels.fieldName = "value"
        labels.fontSizeInMapUnits = False
        labels.fontLimitPixelSize = True
        labels.fontMinPixelSize = 6
        labels.textFont.setFamily(self.fontFamily)
        labels.textFont.setPointSize(8)
        labels.textFont.setWeight(labels.textFont.Normal)
        labels.textColor = PyQt4.QtGui.QColor("#000000")
        labels.placement = qgis.core.QgsPalLayerSettings.AboveLine
        labels.writeToLayer(layer)
        return

    def _set_contour_style(self, layer):
        symbolLayer = layer.rendererV2().symbol().symbolLayer(0)
        symbolLayer.setColor(PyQt4.QtGui.QColor("#000000")) # works
        symbolLayer.setWidth(1.5)
        symbolLayer.setWidthUnit(qgis.core.QgsSymbolV2.Pixel)
        return

    def _set_warning_time_style(self, layer):
        dashed = [2, 2]
        contour_rules = (
            ("negative", '"value" < 0.0', "#666666", dashed,),
            ("positive", '"value" >= 0.0', "#000000", None,),
        )
        symbol = qgis.core.QgsSymbolV2.defaultSymbol(layer.geometryType())
        renderer = qgis.core.QgsRuleBasedRendererV2(symbol)

        root_rule = renderer.rootRule()

        for label, expression, lineColor, lineStyle in contour_rules:
            rule = root_rule.children()[0].clone()
            rule.setLabel(label)
            rule.setFilterExpression(expression)
            rule.symbol().setColor(PyQt4.QtGui.QColor(lineColor))
            symbolLayer = rule.symbol().symbolLayer(0)
            if lineStyle is None:
                symbolLayer.setPenStyle(PyQt4.QtCore.Qt.SolidLine)
            else:
                symbolLayer.setCustomDashPatternUnit(symbol.MM)
                symbolLayer.setCustomDashVector(lineStyle)
                symbolLayer.setUseCustomDashPattern(True)
            root_rule.appendChild(rule)

        root_rule.removeChildAt(0)
        layer.setRendererV2(renderer)
        return

    
    def _render(self, layers, extent, name, format, legendLayers=None, legendTitle=None, title=None):
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
        canvas.setScaleLocked(True)
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

        # Map
        x, y = 0,0
        w, h = composer.paperWidth(), composer.paperHeight()
        map = qgis.core.QgsComposerMap(composer, x, y, w, h)
        #grid = map.grid()
        #grid.setUnits(grid.MapUnit)
        #grid.setIntervalX(0.5)
        #grid.setIntervalY(0.5)
        #grid.setFrameStyle(grid.Zebra)
        #grid.setStyle(grid.FrameAnnotationsOnly)
        #grid.setAnnotationFormat(grid.Decimal)
        #grid.setAnnotationPrecision(2)
        #grid.setFrameSideFlag(grid.FrameLeft)
        composer.addItem(map)

        # Scale bar
        scale = qgis.core.QgsComposerScaleBar(composer)
        scale.setStyle('Single Box') # optionally modify the style
        scale.setComposerMap(map)
        scale.applyDefaultSize()
        scale.setNumSegmentsLeft(0)
        scale.setNumSegments(2)
        scale.setFont(PyQt4.QtGui.QFont(self.fontFamily, 8))
        scale.setHeight(1.5)
        scale.setLabelBarSpace(1.0)
        scale.setItemPosition(0.0, composer.paperHeight(), scale.LowerLeft)
        composer.addItem(scale)

        # Title/label
        if title:
            label = qgis.core.QgsComposerLabel(composer)
            label.setText(title)
            label.setHAlign(PyQt4.QtCore.Qt.AlignCenter)
            label.setFont(PyQt4.QtGui.QFont(self.fontFamily, 12, PyQt4.QtGui.QFont.Bold))
            label.adjustSizeToText()
            label.setItemPosition(0.5*composer.paperWidth(), 0, label.UpperMiddle)
            composer.addItem(label)

        # Legend
        if legendLayers:
            legend = qgis.core.QgsComposerLegend(composer)
            legend.modelV2().rootGroup().removeAllChildren()
            for layer in legendLayers:
                legend.modelV2().rootGroup().addLayer(layer)
            legend.setBoxSpace(0.5)
            legendSize = legend.paintAndDetermineSize(None)
            legend.setItemPosition(composer.paperWidth()-legendSize.width(), composer.paperHeight()-legendSize.height())
            legend.setStyleFont(qgis.core.QgsComposerLegendStyle.Title, PyQt4.QtGui.QFont(self.fontFamily, 10))
            legend.setStyleFont(qgis.core.QgsComposerLegendStyle.SymbolLabel, PyQt4.QtGui.QFont(self.fontFamily, 8))
            if legendTitle:
                legend.setTitle(legendTitle)
            legend.setBackgroundColor(PyQt4.QtGui.QColor(255, 255, 255, 128))
            composer.addItem(legend)

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
        if self.alert is None:
            filename = "{eq}-map_{name}.{format}".format(eq=self.eqId, name=name, format=format)
        else:
            filename = "{eq}-alert{ver:03d}_map_{name}.{format}".format(eq=self.eqId, ver=self.alert["version"], name=name, format=format)
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
