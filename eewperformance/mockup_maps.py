# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
from importlib import import_module
import dateutil.parser
from datetime import datetime
import numpy

import matplotlib.pyplot as pyplot
import matplotlib.colors as colors
import matplotlib.colorbar
import matplotlib.patches as patches
from osgeo import gdal, osr
from cartopy import crs

from cartopy_extra_tiles import cached_tiler
import matplotlib_extras

import gdalraster
import analysis_utils

gdal.UseExceptions()

class EventMaps(object):
    """Maps of MMI, warning time, etc.
    """
    def __init__(self, config, eqId, event, alerts):
        """
        :type config: ConfigParser
        :param config: Configuration for application.
        """
        self.config = config
        self.eqId = eqId
        self.event = event
        self.alerts = alerts

        self.data = None
        return

    def load_data(self):
        """Load data from GDAL raster file.
        """
        plotsDir = self.config.get("files", "plots_dir")

        server = self.config.get("shakealert.production", "server")
        gmpe = self.config.get("mmi_predicted", "gmpe")
        magAlertThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        mmiAlertThreshold = self.config.getfloat("alerts", "mmi_threshold")
        if magAlertThreshold is None:
            magAlertThreshold = params.getfloat("alerts", "magnitude_threshold")
        if mmiAlertThreshold is None:
            mmiAlertThreshold = params.getfloat("alerts", "mmi_threshold")
        label = "{eqId}-{server}-{gmpe}-M{magAlertThreshold:.1f}-MMI{mmiAlertThreshold:.1f}".format(
            eqId=self.event["event_id"], server=server, gmpe=gmpe, magAlertThreshold=magAlertThreshold, mmiAlertThreshold=mmiAlertThreshold)
        filename = "mockup_" + label + ".tiff"
        rasterData = gdal.Open(os.path.join(plotsDir, filename), gdal.GA_ReadOnly)

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
            data = numpy.array(band.ReadAsArray())
            data = numpy.ma.masked_values(data, band.GetNoDataValue())
            self.data["layers"][description] = data

        self._mmi_colormap()
        return    

    def mmi_observed(self):
        """Create map with observed MMI with contours.
        """
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        wgs84CRS = crs.Geodetic()
        
        mmi = self.data["layers"]["mmi_obs"]
        im = ax.imshow(mmi, vmin=0.0, vmax=10.0, extent=dataExtent, transform=dataCRS, origin="upper", cmap="MMI", alpha=0.67, zorder=2)

        contourLevels = numpy.arange(1.0, 10.01, 0.5)
        chandle = ax.contour(mmi, levels=contourLevels, zorder=3, colors="black", origin="upper", extent=dataExtent, transform=dataCRS)
        ax.clabel(chandle, inline=True, fmt="%3.1f", zorder=3)

        ax.plot(self.event["longitude"], self.event["latitude"], transform=wgs84CRS, marker="*", mfc="red", mec="black", c="white", ms=18, zorder=4)

        ax.set_title("ShakeMap - Observed Shaking")
        legend = pyplot.legend(handles=self.mmiPatches, title="MMI", handlelength=0.8, borderpad=0.3, labelspacing=0.2, loc="lower left")

        self._save(figure, "mmi_obs")
        return

    def mmi_predicted(self):
        """Create map with observed MMI with contours.
        """
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        wgs84CRS = crs.Geodetic()

        if len(self.alerts) > 0:
            mmi = self.data["layers"]["mmi_pred"]
            im = ax.imshow(mmi, vmin=0.0, vmax=10.0, extent=dataExtent, transform=dataCRS, origin="upper", cmap="MMI", alpha=0.67, zorder=2)

            contourLevels = numpy.arange(1.0, 10.01, 0.5)
            chandle = ax.contour(mmi, levels=contourLevels, zorder=3, colors="black", origin="upper", extent=dataExtent, transform=dataCRS)
            ax.clabel(chandle, inline=True, fmt="%3.1f", zorder=3)

            cols = [
                ("longitude", "float32",),
                ("latitude", "float32",),
                ]
            epicenters = numpy.zeros(len(self.alerts), dtype=cols)
            for ialert,alert in enumerate(self.alerts):
                epicenters[ialert] = (alert["longitude"], alert["latitude"],)
            ax.plot(epicenters["longitude"], epicenters["latitude"], transform=wgs84CRS, marker="*", mfc="red", mec="black", c="white", ms=18, zorder=4)

        ax.set_title("ShakeAlert - Predicted Shaking")
        pyplot.legend(handles=self.mmiPatches, title="MMI", handlelength=0.8, borderpad=0.3, labelspacing=0.2, loc="lower left")

        self._save(figure, "mmi_pred")
        return

    def mmi_residual(self):
        RESIDUAL_MAX = 2.0
        
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        wgs84CRS = crs.Geodetic()
        
        mmiObs = self.data["layers"]["mmi_obs"]
        mmiPred = self.data["layers"]["mmi_pred"]
        mmiResidual = mmiObs - mmiPred

        cmap = pyplot.cm.get_cmap("RdBu_r")
        cmap._init()
        alpha = numpy.minimum(1.0, numpy.abs(numpy.linspace(-2*RESIDUAL_MAX, 2*RESIDUAL_MAX, cmap.N)))
        alpha[alpha < 0.5] = 0.0
        cmap._lut[:-3,-1] = alpha
        
        im = ax.imshow(mmiResidual, vmin=-RESIDUAL_MAX, vmax=RESIDUAL_MAX, extent=dataExtent, transform=dataCRS, origin="upper", cmap=cmap, alpha=0.67, zorder=2)

        noData = numpy.ma.masked_array(numpy.ones(mmiResidual.shape), ~mmiResidual.mask)
        ax.imshow(noData, extent=dataExtent, transform=dataCRS, origin="upper", cmap="gray_r", vmin=0, vmax=1, alpha=0.3, zorder=3)

        contourLevels = numpy.arange(-RESIDUAL_MAX, RESIDUAL_MAX+0.01, 0.5)
        chandle = ax.contour(mmiResidual, levels=contourLevels, zorder=4, colors="black", origin="upper", extent=dataExtent, transform=dataCRS)
        ax.clabel(chandle, inline=True, fmt="%3.1f", zoerder=4)

        ax.plot(self.event["longitude"], self.event["latitude"], transform=wgs84CRS, marker="*", mfc="red", mec="black", c="white", ms=18, zorder=5)

        ax.set_title("MMI Residual")

        matplotlib_extras.axes.add_background_axes(figure, [0.01, 0.31, 0.16, 0.37])
        cbax = figure.add_axes([0.02, 0.33, 0.02, 0.33])
        colorbar = pyplot.colorbar(im, cax=cbax, format="%4.1f")
        colorbar.set_label("Residual (Obs-Pred)")
        
        self._save(figure, "mmi_residual")
        return

    def alert_warning_time(self):
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        wgs84CRS = crs.Geodetic()
        
        warningTime = self.data["layers"]["warning_time"]
        tmin = numpy.min(warningTime.ravel())
        tmax = numpy.max(warningTime.ravel())
        if tmax > tmin:
            contourLevels = numpy.arange(2.0*numpy.floor(0.5*tmin), numpy.ceil(tmax)+0.01, 2.0)
            chandle = ax.contour(warningTime, levels=contourLevels, zorder=4, colors="black", origin="upper", extent=dataExtent, transform=dataCRS)
            ax.clabel(chandle, inline=True, fmt="%.0f s", color="black", zorder=5)
        
        ax.plot(self.event["longitude"], self.event["latitude"], transform=wgs84CRS, marker="*", mfc="red", mec="black", c="white", ms=18, zorder=6)

        ax.set_title("Warning Time (s)")

        self._save(figure, "alert_warning_time")
        return

    def alert_version(self):
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        wgs84CRS = crs.Geodetic()
        
        version = self.data["layers"]["alert_version"]
        vmin = numpy.min(version.ravel())
        vmax = numpy.max(version.ravel())
        norm = colors.Normalize(vmin, vmax)
        im = ax.imshow(version, extent=dataExtent, transform=dataCRS, origin="upper", cmap="viridis", norm=norm, alpha=0.5, zorder=1)

        ax.plot(self.event["longitude"], self.event["latitude"], transform=wgs84CRS, marker="*", mfc="red", mec="black", c="white", ms=18, zorder=6)

        ax.set_title("Alert Version")

        matplotlib_extras.axes.add_background_axes(figure, [0.01, 0.31, 0.16, 0.37])
        cbax = figure.add_axes([0.02, 0.33, 0.02, 0.33])
        colorbar = pyplot.colorbar(im, cax=cbax, format="%2.0f")
        colorbar.set_label("Version")

        self._save(figure, "alert_version")
        return

    def alert_time(self):
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        wgs84CRS = crs.Geodetic()
        
        alertTime = self.data["layers"]["alert_time"]
        tmin = numpy.min(alertTime.ravel())
        tmax = numpy.max(alertTime.ravel())
        print tmin,tmax
        norm = colors.Normalize(tmin, tmax)
        im = ax.imshow(alertTime, extent=dataExtent, transform=dataCRS, origin="upper", cmap="viridis", norm=norm, alpha=0.5, zorder=1)

        ax.plot(self.event["longitude"], self.event["latitude"], transform=wgs84CRS, marker="*", mfc="red", mec="black", c="white", ms=18, zorder=6)

        ax.set_title("Alert Time after Origin Time (s)")
        matplotlib_extras.axes.add_background_axes(figure, [0.01, 0.31, 0.16, 0.37])
        cbax = figure.add_axes([0.02, 0.33, 0.02, 0.33])
        colorbar = pyplot.colorbar(im, cax=cbax, format="%4.1f")
        colorbar.set_label("Alert Time (s)")

        self._save(figure, "alert_time")
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
                (1.0/10.0, 255.0/255.0, 255.0/255.0),
                (4.0/10.0, 255.0/255.0, 255.0/255.0),
                (5.0/10.0, 147.0/255.0, 147.0/255.0),
                (6.0/10.0,   0.0/255.0,   0.0/255.0),
                (10.0/10.0,   0.0/255.0,   0.0/255.0),
                ],
            'alpha': [
                (0.0/10.0, 0.0/255.0, 0.0/255.0),
                (1.0/10.0, 0.0/255.0, 255.0/255.0),
                (2.0/10.0, 255.0/255.0, 255.0/255.0),
                (3.0/10.0, 255.0/255.0, 255.0/255.0),
                (10.0/10.0, 255.0/255.0, 255.0/255.0),
                ]
            }
        mmicmap = colors.LinearSegmentedColormap("MMI", cdict)
        pyplot.register_cmap(name="MMI", cmap=mmicmap)

        self.mmiPatches = [
            patches.Patch(ec="black", fc=(200.0/255.0,   0.0/255.0,   0.0/255.0), label="X"),
            patches.Patch(ec="black", fc=(255.0/255.0,   0.0/255.0,   0.0/255.0), label="IX"),
            patches.Patch(ec="black", fc=(255.0/255.0, 145.0/255.0,   0.0/255.0), label="VIII"),
            patches.Patch(ec="black", fc=(255.0/255.0, 200.0/255.0,   0.0/255.0), label="VII"),
            patches.Patch(ec="black", fc=(255.0/255.0, 255.0/255.0,   0.0/255.0), label="VI"),
            patches.Patch(ec="black", fc=(122.0/255.0, 255.0/255.0, 147.0/255.0), label="V"),
            patches.Patch(ec="black", fc=(128.0/255.0, 255.0/255.0, 255.0/255.0), label="IV"),
            patches.Patch(ec="black", fc=(160.0/255.0, 230.0/255.0, 255.0/255.0), label="III"),
            patches.Patch(ec="black", fc=(191.0/255.0, 204.0/255.0, 255.0/255.0), label="II"),
            patches.Patch(ec="black", fc=(255.0/255.0, 255.0/255.0, 255.0/255.0), label="I"),
            ]
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

        rectFactory = matplotlib_extras.axes.RectFactory(figure, margins=((0, 0, 0), (0, 0, 0.25)))
        ax = pyplot.axes(rectFactory.rect(), projection=tiler.crs)
        ax.set_extent(self.data["extent"])

        ax.add_image(tiler, tilerZoom, zorder=0, cmap="gray")
        return figure

    def _save(self, figure, label):
        """
        """
        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)

        server = self.config.get("shakealert.production", "server")
        gmpe = self.config.get("mmi_predicted", "gmpe")
        magAlertThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        mmiAlertThreshold = self.config.getfloat("alerts", "mmi_threshold")
        if magAlertThreshold is None:
            magAlertThreshold = params.getfloat("alerts", "magnitude_threshold")
        if mmiAlertThreshold is None:
            mmiAlertThreshold = params.getfloat("alerts", "mmi_threshold")
        labelEv = "{eqId}-{server}-{gmpe}-M{magAlertThreshold:.1f}-MMI{mmiAlertThreshold:.1f}".format(
            eqId=self.event["event_id"], server=server, gmpe=gmpe, magAlertThreshold=magAlertThreshold, mmiAlertThreshold=mmiAlertThreshold)
        filename = "{}-map_{}.jpg".format(labelEv, label)
        figure.savefig(os.path.join(plotsDir, filename), pad_inches=0.02)
        pyplot.close(figure)
        return
    

# End of file
