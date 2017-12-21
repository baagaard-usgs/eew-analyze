# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import numpy

from osgeo import gdal, osr

import gdalraster
import qgisconverter

gdal.UseExceptions()

class MapPanels(object):

    def __init__(self, eqId):
        """
        :type eqId: str
        :param eqId: ComCat earthquake id.
        """
        self.eqId = eqId
        self.qgis = None
        self.layers = {}
        return

    def initialize(self):
        """
        Setup and load layers to be used in plotting.
        """
        self.layers = {}

        # Basemap
        filename = self.params.get("maps", "basemap")
        self.layers["basemap"] = qgis.core.QgsRasterLayer(filename)
        if not self.layers["basemap"].isValid():
            raise IOError("Could not load basemap from '{}'.".format(filename))

        # Temporary virtual raster bands
        src = gdal.Open(filename, gdal.GA_ReadOnly)
        for name,layer in qgisconverter.extract_bands(src).items():
            self.layers[name] = layer

        # MMI residual (predicted - observed)
        mmiObs = numpy.array(src.getRasterBand(1).ReadAsArray())
        mmiPred = numpy.array(src.getRasterBand(2).ReadAsArray())
        mmiResidual = mmiPred - mmiObs
        gdalLayer = gdalraster.clone(filename, values, src)
        self.layers["mmi_residual"] = qgisconverter.gdal_to_qgisraster(gdalLayer)
        del gdalLayer
        del src
            
        # Contours from temporary virtual raster bands.
        for value in ("mmi_obs", "mmi_pred", "warning_time",):
            ogrLayer = gdalraster.contours_from_raster(self.layers[value])
            self.layers[value+"_contours"] = qgisconverter.ogr_to_qgisvector(self.layers[value], ogrLayer)

        # Epicenter
        options = "delimiter={delimiter}&crs=EPSG:4326&xField={x}&yField={}".format(delimiter="|", x="longitude", y="latitude")
        uri = "file://{filename}&{options}".format(filename, options)
        self.layers["epicenter_obs"] = qgis.core.QgsVectorLayer(uri, "epicenter_obs", "delimitedtext")

        uri = "file://{filename}&{options}".format(filename, options)
        self.layers["epicenter_pred"] = qgis.core.QgsVectorLayer(uri, "epicenter_pred", "delimitedtext")

                        
                          
        self._initQgis()
        return
    
    def mmi_observed(self):
        """Create map with observed MMI with contours.
        """
        basemap = self.layers["basemap"]

        mmi = self.layers["mmi_obs"]
        mmi.loadNamedStyle("mmi.qml")

        mmiContours = self.layers["mmi_obs_contours"]
        # :TODO: set style
        
        epicenter = self.layers["epicenter_obs"]
        # :TODO: set style
        
        layerRegistry = qgis.core.QgsMapLayerRegistry.instance()
        layerRegistry.addMapLayers([mmiContours, mmi, basemap])
        self._render(mmi.extent(), layerRegistry, "mmi_obs", "png")
        return

    def mmi_predicted(self):
        return

    def mmi_residual(self):
        return

    def mmi_warning_time(self):
        return

    def alert_categories(self):
        return

    def _initQgis(self):
        self.qgis = qgis.core.QgsApplication([], True)
        qgisPrefixPath = self.params.get("qgis", "prefix_path")
        if qgisPrefixPath:
            # "/Applications/QGIS.app/Contents/MacOS"
            qgis.core.QgsApplication.setPrefixPath(qgisPrefixPath, True)
        self.qgis.initQgis()
        return

    def _render(self, layerRegistry, extent, name, format):
        """Render image to file in given format.
        """
        renderer = qgis.core.QgsMapRenderer()
        renderer.setExtent(extent)
        renderer.setLayerSet([layer.id() for layer in layerRegistry.mapLayers()])

        imageWidth = self.params.get("maps", "image_width")
        imageHeight = self.params.get("maps", "image_height")
        bgColor = self.params.get("maps", "bg_color")
        image = PyQt4.QtGui.QImage(PyQt4.QtCore.QSize(imageWidth, imageHeight), PyQt4.QtGui.QImage.Format_ARGB32_Premultiplied)
        image.fill(PyQt4.QtCore.QColor(bgColor))

        painter = PyQt4.QtGui.QPainter()
        painter.begin(image)
        painter.setRenderHint(PyQt4.QtGui.QPainter.Antialiasing)
        
        renderer.setOutputSize(image.size(), image.logicalDpiX())
        renderer.render(painter)

        painter.end()

        plotsDir = self.params.get("files", "plots_dir")
        fullFilename = "{eq}-map_{value}.{format}".format(eq=self.eqId, value=value, format=format)
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

        self.qgis.exitQgis()

        import subprocess
        for name in TEMPORARY_LAYERS:
            filename = self.layers[name].source()
            subprocess.call("rm {}".filename)
        return
