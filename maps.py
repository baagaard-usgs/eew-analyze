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
import PyQt4.QtCore
import PyQt4.QtGui
from osgeo import gdal, osr

import gdalraster
import qgisconverter

gdal.UseExceptions()

class MapPanels(object):
    """Maps of MMI, warning time, etc.
    """

    def __init__(self, eqId, config):
        """
        :type eqId: str
        :param eqId: ComCat earthquake id.

        :type config: ConfigParser
        :param config: Configuration for application.
        """
        self.eqId = eqId
        self.config = config
        
        self.qgis = None
        self.layers = {}

        self._initialize()
        return

    def __del__(self):
        """
        """
        self._exitQgis()
        return
        
    
    def mmi_observed(self):
        """Create map with observed MMI with contours.
        """
        basemap = self.layers["basemap"]
        basemap.pipe().hueSaturationFilter().setGrayscaleMode(int(qgis.core.QgsImageOperation.GrayscaleLightness))
        basemap.pipe().brightnessFilter().setBrightness(50)
        basemap.pipe().brightnessFilter().setContrast(-20)
        
        mmi = self.layers["mmi_obs"]
        mmi.loadNamedStyle("mmi.qml")

        mmiContours = self.layers["mmi_obs_contour"]
        symbolLayer = mmiContours.rendererV2().symbol().symbolLayer(0)
        symbolLayer.setColor(PyQt4.QtGui.QColor("#000000"))
        symbolLayer.setWidth(2.0)

        #epicenter = self.layers["epicenter_obs"]
        # :TODO: set style
        
        layerRegistry = qgis.core.QgsMapLayerRegistry.instance()
        layerRegistry.addMapLayers([mmiContours, mmi, basemap])
        self._render(layerRegistry, mmi.extent(), "mmi_obs", "png")
        return

    def mmi_predicted(self):
        return

    def mmi_residual(self):
        return

    def mmi_warning_time(self):
        return

    def alert_categories(self):
        return

    def _initialize(self):
        """
        Setup and load layers to be used in plotting.
        """
        self.layers = {}
        self._initQgis()

        # Basemap
        filename = self.config.get("maps", "basemap")
        self.layers["basemap"] = qgis.core.QgsRasterLayer(filename)
        if not self.layers["basemap"].isValid():
            raise IOError("Could not load basemap from '{}'.".format(filename))
        self.layers["basemap"].setDrawingStyle('MultiBandSingleBandGray')
        self.layers["basemap"].renderer().setGrayBand(2)

        # Temporary virtual raster bands
        filename = self.config.get("files", "analysis_event").replace("[EVENTID]", self.eqId)
        src = gdal.Open(filename, gdal.GA_ReadOnly)
        for name,layer in qgisconverter.extract_raster_bands(src).items():
            self.layers[name] = layer

        # MMI residual (predicted - observed)
        mmiObs = numpy.array(src.GetRasterBand(1).ReadAsArray())
        mmiPred = numpy.array(src.GetRasterBand(2).ReadAsArray())
        mmiResidual = mmiPred - mmiObs
        filename = self.config.get("files", "analysis_event").replace("[EVENTID]", self.eqId).replace(".tiff", "-mmi_residual.tiff")
        gdalLayer = gdalraster.clone_new_data("mmi_residual", mmiResidual, src)
        self.layers["mmi_residual"] = qgisconverter.gdal_to_qgisraster(gdalLayer, filename)
        del gdalLayer

        # Contours from temporary virtual raster bands.
        for value in ("mmi_obs", "mmi_pred", "warning_time"):
            if value in ("mmi_obs", "mmi_pred",):
                cstart = 1.5
                cinterval = 0.5
                clevels = []
            else:
                cstart = 0
                cinterval = 0
                clevels = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
            cvalue = value+"_contour"
            ogrLayer = gdalraster.contours_from_raster(src, value, cstart, cinterval, clevels)
            fileSuffix = "_{}.shp".format(cvalue)
            filename = self.config.get("files", "analysis_event").replace("[EVENTID]", self.eqId).replace(".tiff", fileSuffix)
            self.layers[cvalue] = qgisconverter.ogr_to_qgisvector(ogrLayer, filename)
        del src

        return
        # Epicenter
        options = "delimiter={delimiter}&crs=EPSG:4326&xField={x}&yField={}".format(delimiter="|", x="longitude", y="latitude")
        uri = "file://{filename}&{options}".format(filename, options)
        self.layers["epicenter_obs"] = qgis.core.QgsVectorLayer(uri, "epicenter_obs", "delimitedtext")

        uri = "file://{filename}&{options}".format(filename, options)
        self.layers["epicenter_pred"] = qgis.core.QgsVectorLayer(uri, "epicenter_pred", "delimitedtext")

        return
    
    def _initQgis(self):
        self.qgis = qgis.core.QgsApplication([], False)
        qgisPrefixPath = self.config.get("qgis", "prefix_path")
        if qgisPrefixPath != "None":
            # "/Applications/QGIS.app/Contents/MacOS"
            qgis.core.QgsApplication.setPrefixPath(qgisPrefixPath, True)
        self.qgis.initQgis()
        return

    def _render(self, layerRegistry, extent, name, format):
        """Render image to file in given format.
        """
        renderer = qgis.core.QgsMapRenderer()
        renderer.setExtent(extent)
        renderer.setLayerSet([layer.id() for layer in layerRegistry.children()])
        
        projection = self.config.get("maps", "projection")
        renderer.setDestinationCrs(qgis.core.QgsCoordinateReferenceSystem(projection))
        renderer.setProjectionsEnabled(True)

        
        imageWidth = self.config.getint("maps", "width_pixels")
        imageHeight = self.config.getint("maps", "height_pixels")
        bgColor = self.config.get("maps", "bg_color")
        image = PyQt4.QtGui.QImage(PyQt4.QtCore.QSize(imageWidth, imageHeight), PyQt4.QtGui.QImage.Format_ARGB32_Premultiplied)
        image.fill(PyQt4.QtGui.QColor(bgColor))

        painter = PyQt4.QtGui.QPainter()
        painter.begin(image)
        painter.setRenderHint(PyQt4.QtGui.QPainter.Antialiasing)
        
        renderer.setOutputSize(image.size(), image.logicalDpiX())
        renderer.render(painter)

        painter.end()

        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = "{eq}-map_{name}.{format}".format(eq=self.eqId, name=name, format=format)
        image.save(os.path.join(plotsDir, filename), format)
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
            if name in self.layers:
                filename = self.layers[name].source()
                subprocess.call("rm {}".format(filename), shell=True)
        del self.layers
                
        if self.qgis:
            self.qgis.exitQgis()

        return
