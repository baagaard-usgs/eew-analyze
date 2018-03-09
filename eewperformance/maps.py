# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
from importlib import import_module
import numpy

import matplotlib.pyplot as pyplot
import matplotlib.colors as colors
from osgeo import gdal, osr
from cartopy import crs

from cartopy_extra_tiles import cached_tiler

import gdalraster
import analysis_utils

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

        self.eqId = None
        self.alert = None
        self.data = None
        return
    
    def load_data(self, eqId, alert=None):
        """
        :type eqId: str
        :param eqId: ComCat earthquake id.

        :type alert: dict
        :param alert: ShakeAlert alert information from analysis database.
        """
        self.eqId = eqId
        self.alert = alert

        cacheDir = self.config.get("files", "analysis_cache_dir")
        filename = os.path.join(cacheDir, "analysis_" + analysis_utils.analysis_label(self.config, eqId) + ".tiff")
        rasterData = gdal.Open(filename, gdal.GA_ReadOnly)

        srs = osr.SpatialReference()
        srs.ImportFromWkt(rasterData.GetProjection())
        if srs.GetAuthorityCode("GEOGCS") == "4326":
            rasterCRS = crs.PlateCarree()
        else:
            rasterCRS = crs.epsg(srs.GeoAuthorityCode("PROJCS"))

        geot = rasterData.GetGeoTransform()
        extent = (
            geot[0],
            geot[0] + rasterData.RasterXSize * geot[1],
            geot[3] + rasterData.RasterYSize * geot[5],
            geot[3]
        )

        self.data = {
            "crs": rasterCRS,
            "extent": extent,
            "layers": {},
        }
        for iband in range(rasterData.RasterCount):
            band = rasterData.GetRasterBand(1+iband)
            description = band.GetDescription()
            print iband, description
            self.data["layers"][description] = numpy.array(band.ReadAsArray())

        self._mmi_colormap()
        return    

    def mmi_observed(self):
        """Create map with observed MMI with contours.
        """
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        
        mmi = self.data["layers"]["mmi_obs"]
        ax.imshow(mmi, vmin=0.0, vmax=10.0, extent=dataExtent, transform=dataCRS, origin="upper", cmap="MMI", alpha=0.67, zorder=2)

        contourLevels = numpy.arange(1.0, 10.01, 0.5)
        chandle = ax.contour(mmi, levels=contourLevels, zorder=3, colors="black", origin="upper", extent=dataExtent, transform=dataCRS)
        ax.clabel(chandle, inline=True, fmt="%3.1f", fontsize=8)

        ax.set_title("Observed Shaking")

        #colorbar = pyplot.colorbar(cmap="MMI")
        #colorbar.set_label("MMI")

        # Epicenter

        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = analysis_utils.analysis_label(self.config, self.eqId)
        filename += "-map_mmi_obs.jpg"
        figure.savefig(os.path.join(plotsDir, filename), pad_inches=0.02)
        return

    def mmi_predicted(self):
        """Create map with observed MMI with contours.
        """
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        
        mmi = self.data["layers"]["mmi_pred"]
        ax.imshow(mmi, vmin=0.0, vmax=10.0, extent=dataExtent, transform=dataCRS, origin="upper", cmap="MMI", alpha=0.67, zorder=2)

        contourLevels = numpy.arange(1.0, 10.01, 0.5)
        chandle = ax.contour(mmi, levels=contourLevels, zorder=3, colors="black", origin="upper", extent=dataExtent, transform=dataCRS)
        ax.clabel(chandle, inline=True, fmt="%3.1f", fontsize=8)

        ax.set_title("ShakeAlert Predicted Shaking")

        #colorbar = pyplot.colorbar(cmap="MMI")
        #colorbar.set_label("MMI")

        # Epicenter

        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = analysis_utils.analysis_label(self.config, self.eqId)
        filename += "-map_mmi_pred.jpg"
        # if self.alert
        # :TODO:
        figure.savefig(os.path.join(plotsDir, filename), pad_inches=0.02)
        return

    def mmi_residual(self):
        # residual (observed - predicted)
        mmiObs = self.data["mmi_obs"]
        mmiPred = self.data["mmi_pred"]
        noDataValue = gdalraster.NO_DATA_VALUE
        mask = mmiPred != noDataValue
        mmiResidual = mask*(mmiObs - mmiPred) + ~mask*noDataValue

        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        
        ax.imshow(mmiResidual, vmin=-2.0, vmax=2.0, extent=dataExtent, transform=dataCRS, origin="upper", cmap="BlueRed", alpha=0.67, zorder=2)

        contourLevels = numpy.arange(-2.0, 2.01, 0.5)
        chandle = ax.contour(mmi, levels=contourLevels, zorder=3, colors="black", origin="upper", extent=dataExtent, transform=dataCRS)
        ax.clabel(chandle, inline=True, fmt="%3.1f", fontsize=8)

        ax.set_title("MMI Residual")

        #colorbar = pyplot.colorbar(cmap="MMI")
        #colorbar.set_label("MMI")

        # Epicenter

        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = analysis_utils.analysis_label(self.config, self.eqId)
        filename += "-map_mmi_residual.jpg"
        figure.savefig(os.path.join(plotsDir, filename), pad_inches=0.02)
        return

    def mmi_warning_time(self, t):
        """Create map with predicted MMI with warning time contours.
        """
        basemap = self.baseLayers["basemap"]

        mmi = self.dataLayers["mmi_pred"]
        mmi.loadNamedStyle("mmi.qml")

        warningContours = self.dataLayers["warning_time_contour"]
        self._set_warning_time_style(warningContours)
        self._add_labels(warningContours)
            
        #epicenter = self.dataLayers["epicenter_obs"]
        # :TODO: set style

        import dateutil.parser
        tstamp = dateutil.parser.parse(self.alert["timestamp"])
        title = "{tstamp:%Y-%m-%d %H:%M:%S.%f} ({t:6.3f}s after origin time), Alert {alert[version]:3d}, M{alert[magnitude]:4.2f}".format(alert=self.alert, tstamp=tstamp, t=t)
        self._render([warningContours, mmi, basemap], mmi.extent(), "mmi_warning", "jpg", legendLayers=[mmi], legendTitle="MMI", title=title)
        return

    def alert_category(self):
        basemap = self.baseLayers["basemap"]
        
        popDensity = self.dataLayers["population_density"]
        popDensity.loadNamedStyle("pop_density.qml")

        category = self.dataLayers["alert_category"]
        category.loadNamedStyle("alert_category.qml")

        warningContours = self.dataLayers["warning_time_contour"]
        self._set_warning_time_style(warningContours)
        self._add_labels(warningContours)
            
        #epicenter = self.dataLayers["epicenter_obs"]
        # :TODO: set style
        
        self._render([warningContours, category, popDensity, basemap], category.extent(), "alert_category", "jpg", legendLayers=[category], legendTitle="Alert Category", title="Alert Classification and Warning Time (s)") # , Alert Threshold: MMI X.X
        return

    def _mmi_colormap(self):
        cdict = {
            'red': [
                (0.0/10.0, 255.0/255.0, 255.0/255.0),
                (1.0/10.0, 255.0/255.0, 255.0/255.0),
                (2.0/10.0, 191.0/255.0, 191.0/255.0),
                (3.0/10.0, 160.0/255.0, 160.0/255.0),
                (4.0/10.0, 128.0/255.0, 128.0/255.0),
                (5.0/10.0, 122.0/255.0, 122.0/255.0),
                (6.0/10.0, 255.0/255.0, 255.0/255.0),
                (10.0/10.0, 200.0/255.0, 200.0/255.0),
                ],
            'green': [
                (0.0/10.0, 255.0/255.0, 255.0/255.0),
                (1.0/10.0, 255.0/255.0, 255.0/255.0),
                (2.0/10.0, 204.0/255.0, 204.0/255.0),
                (3.0/10.0, 230.0/255.0, 230.0/255.0),
                (4.0/10.0, 255.0/255.0, 255.0/255.0),
                (5.0/10.0, 255.0/255.0, 255.0/255.0),
                (6.0/10.0, 255.0/255.0, 255.0/255.0),
                (7.0/10.0, 200.0/255.0, 200.0/255.0),
                (8.0/10.0, 145.0/255.0, 145.0/255.0),
                (9.0/10.0,   0.0/255.0,   0.0/255.0),
                (10.0/10.0,  0.0/255.0,   0.0/255.0),
                ],
            'blue': [
                (0.0/10.0, 255.0/255.0, 255.0/255.0),
                (4.0/10.0, 255.0/255.0, 255.0/255.0),
                (5.0/10.0, 147.0/255.0, 147.0/255.0),
                (6.0/10.0,   0.0/255.0,   0.0/255.0),
                (10.0/10.0,   0.0/255.0,   0.0/255.0),
                ],
            }
        mmicmap = colors.LinearSegmentedColormap('MMI', cdict)
        pyplot.register_cmap(name='MMI', data=cdict)
        return
    

    def _create_figure(self):
        tilerPath = self.config.get("maps", "tiler").split(".")
        tilerObj = getattr(import_module(".".join(tilerPath[:-1])), tilerPath[-1])
        tilerStyle = self.config.get("maps", "tiler_style")
        tilerZoom = self.config.getint("maps", "zoom_level")
        tilesDir = self.config.get("maps", "tiler_cache_dir")
        tiler = cached_tiler.CachedTiler(tilerObj(desired_tile_form="L", style=tilerStyle), cache_dir=tilesDir)

        figWidthIn = self.config.getfloat("maps", "width_in")
        figHeightIn = self.config.getfloat("maps", "height_in")
        figure = pyplot.figure(figsize=(figWidthIn, figHeightIn))
        ax = pyplot.axes(projection=tiler.crs)
        ax.set_extent(self.data["extent"])
        ax.add_image(tiler, tilerZoom, zorder=0, cmap="gray")        
        return figure
    
    def _set_warning_time_style(self, layer):
        dashed = [2, 2]
        contour_rules = (
            ("negative", '"value" < 0.0', "#666666", dashed,),
            ("positive", '"value" >= 0.0', "#000000", None,),
        )
        return

    
# End of file
