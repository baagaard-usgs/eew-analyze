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

from . import gdalraster
from . import analysis_utils

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
        cacheDir = self.config.get("files", "analysis_cache_dir")
        filename = "analysis_" + analysis_utils.analysis_event_label(self.config, self.eqId) + ".tiff"
        rasterData = gdal.Open(os.path.join(cacheDir, filename), gdal.GA_ReadOnly)

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
        self._alert_colormap()
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
        ax.clabel(chandle, inline=True, fmt="%3.1f")

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
            ax.clabel(chandle, inline=True, fmt="%3.1f")

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
        ax.clabel(chandle, inline=True, fmt="%3.1f")

        ax.plot(self.event["longitude"], self.event["latitude"], transform=wgs84CRS, marker="*", mfc="red", mec="black", c="white", ms=18, zorder=5)

        ax.set_title("MMI Residual")

        matplotlib_extras.axes.add_background_axes(figure, [0.01, 0.31, 0.16, 0.37])
        cbax = figure.add_axes([0.02, 0.33, 0.02, 0.33])
        colorbar = pyplot.colorbar(im, cax=cbax, format="%4.1f")
        colorbar.set_label("Residual (Obs-Pred)")
        
        self._save(figure, "mmi_residual")
        return

    def mmi_warning_time(self, t):
        """Create map with predicted MMI with warning time contours.
        """
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        wgs84CRS = crs.Geodetic()
        
        mmi = self.data["layers"]["mmi_pred"]
        im = ax.imshow(mmi, vmin=0.0, vmax=10.0, extent=dataExtent, transform=dataCRS, origin="upper", cmap="MMI", alpha=0.67, zorder=2)

        warningTime = self.data["layers"]["warning_time"]
        tmin = numpy.min(warningTime.ravel())
        tmax = numpy.max(warningTime.ravel())
        contourLevels = numpy.arange(2.0*numpy.floor(0.5*tmin), numpy.ceil(tmax)+0.01, 2.0)
        chandle = ax.contour(warningTime, levels=contourLevels, zorder=4, colors="black", origin="upper", extent=dataExtent, transform=dataCRS)
        ax.clabel(chandle, inline=True, fmt="%3.1fs")

        cols = [
            ("longitude", "float32",),
            ("latitude", "float32",),
            ]
        epicenters = numpy.zeros(len(self.alerts), dtype=cols)
        for ialert,alert in enumerate(self.alerts):
            epicenters[ialert] = (alert["longitude"], alert["latitude"],)
        ax.plot(epicenters["longitude"], epicenters["latitude"], transform=wgs84CRS, marker="*", mfc="red", mec="black", c="white", ms=18, zorder=5)

        import dateutil.parser
        tstamp = dateutil.parser.parse(self.alerts[0]["timestamp"])
        title = "{tstamp:%Y-%m-%d %H:%M:%S.%f} ({t:6.3f}s after origin time), Alert {alert[version]:3d}, M{alert[magnitude]:4.2f}".format(alert=self.alerts[0], tstamp=tstamp, t=t)
        ax.set_title(title)
        pyplot.legend(handles=self.mmiPatches, title="MMI", handlelength=0.8, borderpad=0.3, labelspacing=0.2, loc="lower left")

        self._save(figure, "mmi_pred")
        return

    def alert_category(self):
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        wgs84CRS = crs.Geodetic()
        
        popDensity = self.data["layers"]["population_density"]
        norm = colors.LogNorm(vmin=0.01, vmax=1000)
        ax.imshow(popDensity, norm=norm, extent=dataExtent, transform=dataCRS, origin="upper", cmap="gray_r", alpha=0.2, zorder=2)

        category = self.data["layers"]["alert_category"]
        norm = colors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], 4)
        im = ax.imshow(category, extent=dataExtent, transform=dataCRS, origin="upper", cmap="AlertCategory", norm=norm, alpha=0.5, zorder=1)

        warningTime = self.data["layers"]["warning_time"]
        tmin = numpy.min(warningTime.ravel())
        tmax = numpy.max(warningTime.ravel())
        if tmax > tmin:
            contourLevels = numpy.arange(2.0*numpy.floor(0.5*tmin), numpy.ceil(tmax)+0.01, 2.0)
            chandle = ax.contour(warningTime, levels=contourLevels, zorder=4, colors="black", origin="upper", extent=dataExtent, transform=dataCRS)
            ax.clabel(chandle, inline=True, fmt="%.0f s")
        
        ax.plot(self.event["longitude"], self.event["latitude"], transform=wgs84CRS, marker="*", mfc="red", mec="black", c="white", ms=18, zorder=6)

        ax.set_title("Alert Classification and Warning Time")

        pyplot.legend(handles=self.alertPatches, title="Alert Classification", handlelength=0.8, borderpad=0.3, labelspacing=0.2, loc="lower left")
        
        self._save(figure, "alert_category")
        return

    def cost_savings(self):
        """Map of normalized cost savings.
        """
        figure = self._create_figure()
        ax = figure.gca()
        
        dataExtent = self.data["extent"]
        dataCRS = self.data["crs"]
        wgs84CRS = crs.Geodetic()
        
        layers = self.data["layers"]

        popDensity = layers["population_density"]
        norm = colors.LogNorm(vmin=0.01, vmax=1000)
        ax.imshow(popDensity, norm=norm, extent=dataExtent, transform=dataCRS, origin="upper", cmap="gray_r", alpha=0.2, zorder=2)

        categoryTN = 0.0
        categoryFN = 1.0
        maskTN = numpy.ma.masked_values(layers["alert_category"], categoryTN).mask
        maskFN = numpy.ma.masked_values(layers["alert_category"], categoryFN).mask
        maskAlert = layers["alert_category"] > 1.5
        
        costSavings = maskAlert*(layers["cost_no_eew"]-layers["cost_eew"]) - maskFN*(layers["cost_no_eew"]-layers["cost_perfect_eew"])
        cmap = pyplot.cm.get_cmap("PiYG")
        costSavingsColors = colors.Normalize(-1.0, 1.0)(costSavings)
        costSavingsColors = cmap(costSavingsColors)
        alpha = 0.67 * numpy.ones(costSavings.shape)
        alpha[maskTN] = 0.0
        costSavingsColors[..., -1] = alpha
        
        im = ax.imshow(costSavingsColors, extent=dataExtent, transform=dataCRS, origin="upper", zorder=2)
        
        warningTime = self.data["layers"]["warning_time"]
        tmin = numpy.min(warningTime.ravel())
        tmax = numpy.max(warningTime.ravel())
        if tmax > tmin:
            cstep = self.config.getfloat("maps", "warning_time_contour_interval")
            contourLevels = numpy.arange(cstep*numpy.floor(0.5*tmin), numpy.ceil(tmax)+0.01, cstep)
            chandle = ax.contour(warningTime, levels=contourLevels, colors="black", origin="upper", extent=dataExtent, transform=dataCRS, linewidth=0.5, zorder=4)
            ax.clabel(chandle, inline=True, fmt="%.0f s")
        
        ax.plot(self.event["longitude"], self.event["latitude"], transform=wgs84CRS, marker="*", mfc="c_yellow", mec="black", c="white", ms=15, zorder=7)

        ax.set_title("Cost Savings")

        matplotlib_extras.axes.add_background_axes(figure, [0.01, 0.31, 0.16, 0.37])
        cbax = figure.add_axes([0.02, 0.33, 0.02, 0.33])
        colorbar = matplotlib.colorbar.ColorbarBase(cbax, cmap=cmap, norm=colors.Normalize(-1.0, 1.0), format="%4.1f")
        colorbar.set_label("Cost Savings")


        
        self._save(figure, "cost_savings")
        return

    def _alert_colormap(self):
        clist = (
            ("#abd9e9", "True Negative",),
            ("#d7191c", "False Negative",),
            ("#fdae61", "False Positive",),
            ("#2c7bb6", "True Positive",),
            )
        alertcmap = colors.ListedColormap([c[0] for c in clist])
        pyplot.register_cmap(name="AlertCategory", cmap=alertcmap)

        self.alertPatches = [patches.Patch(ec="black", fc=c[0], label=c[1]) for c in clist]
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
        filename = analysis_utils.analysis_event_label(self.config, self.eqId)
        filename += "-map_{}.jpg".format(label)
        figure.savefig(os.path.join(plotsDir, filename), pad_inches=0.02)
        pyplot.close(figure)
        return
    

class SummaryMaps(object):
    """Maps of earthquakes, Q, etc.
    """
    def __init__(self, config, events, db):
        """
        :type config: ConfigParser
        :param config: Configuration for application.
        """
        self.config = config
        self.events = events
        self.db = db
        return

    def _create_figure(self, extent):
        tilerPath = self.config.get("maps", "tiler").split(".")
        tilerObj = getattr(import_module(".".join(tilerPath[:-1])), tilerPath[-1])
        tilerStyle = self.config.get("maps", "tiler_style")
        tilerZoom = 6
        tilesDir = self.config.get("maps", "tiler_cache_dir")
        tiler = cached_tiler.CachedTiler(tilerObj(desired_tile_form="L", style=tilerStyle), cache_dir=tilesDir)

        figWidthIn = 1.5*5.0
        figHeightIn = 1.5*4.2
        figure = pyplot.figure(figsize=(figWidthIn, figHeightIn))

        figure.subplots_adjust(bottom=0.01, top=0.96, left=0.01, right=0.99)
        ax = pyplot.axes(projection=tiler.crs)
        ax.set_extent(extent)

        ax.add_image(tiler, tilerZoom, zorder=0, cmap="gray")
        return figure

    def _save(self, figure, label):
        """
        """
        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = "eqset_" + analysis_utils.analysis_label(self.config)
        filename += "-map_{}.jpg".format(label)
        figure.savefig(os.path.join(plotsDir, filename), pad_inches=0.02)
        pyplot.close(figure)
        return
    

    def earthquakes(self):
        """Create map with earthquakes colored by origin time.
        """
        COLS = [
            ("longitude", "float32",),
            ("latitude", "float32",),
            ("magnitude", "float32",),
            ("origin_time", "datetime64[us]",),
            ]
            
        eqs = numpy.zeros(len(self.events), dtype=COLS)
        for i,eqId in enumerate(self.events):
            event = self.db.comcat_event(eqId)
            ot = numpy.datetime64(dateutil.parser.parse(event["origin_time"]))
            eqs[i] = (event["longitude"], event["latitude"], event["magnitude"], ot)
            
        extent = [
            numpy.min(eqs["longitude"])-2.0, numpy.max(eqs["longitude"])+2.0,
            numpy.min(eqs["latitude"])-0.5, numpy.max(eqs["latitude"])+0.5,
        ]

        from matplotlib.dates import YearLocator,date2num,DateFormatter
            
        figure = self._create_figure(extent)
        ax = figure.gca()
        ms = 0.05 * 10**(0.75*eqs["magnitude"])
        ot = date2num(eqs["origin_time"].astype(datetime))
        sc = ax.scatter(eqs["longitude"], eqs["latitude"], s=ms, c=ot, cmap="viridis", transform=crs.Geodetic(), edgecolors="black", alpha=0.67, zorder=4)

        matplotlib_extras.axes.add_background_axes(figure, [0.025, 0.015, 0.13, 0.36])
        cbax = figure.add_axes([0.03, 0.02, 0.02, 0.33])
        colorbar = pyplot.colorbar(mappable=sc, cax=cbax, ticks=YearLocator(), format=DateFormatter('%Y'))
        colorbar.set_label("Origin Time")
        
        # domains (manual)
        import cartopy.geodesic as geodesic
        import cartopy.feature as feature
        import shapely.geometry as geometry
        
        geod = geodesic.Geodesic()
        circleSF = geod.circle(-122.419, 37.775, 2.0*111e+3)
        circleLA = geod.circle(-117.396, 33.953, 2.5*111e+3)
        geometryDomains = [geometry.Polygon(circleSF), geometry.Polygon(circleLA)]
        ax.add_feature(feature.ShapelyFeature(geometryDomains, crs.PlateCarree(), facecolor="none", edgecolor="blue"), zorder=3)
                            
        self._save(figure, "events")
        return

    def cost_savings(self, metric="area_costsavings_eew"):
        """Create map with earthquakes colored by cost savings with marker size by magnitude.
        """
        COLS = [
            ("longitude", "float32",),
            ("latitude", "float32",),
            ("magnitude", "float32",),
            ("cost_savings", "float32",),
            ]

        server = self.config.get("shakealert.production", "server")
        gmpe = self.config.get("mmi_predicted", "gmpe")
        fragility = self.config.get("fragility_curves", "object").split(".")[-1]
        
        perfs = None
        for eqId in self.events:
            p = self.db.performance_stats(eqId, server, gmpe, fragility)
            if perfs is None:
                perfs = p
            else:
                perfs = numpy.append(perfs, p)

        mmiThreshold = self.config.getfloat("alerts", "mmi_threshold")
        magThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        mask = numpy.logical_and(numpy.ma.masked_values(perfs["mmi_threshold"], mmiThreshold).mask,
                                 numpy.ma.masked_values(perfs["magnitude_threshold"], magThreshold).mask)
        perfs = perfs[mask]
        
        eqs = numpy.zeros(len(self.events), dtype=COLS)
        for i,p in enumerate(perfs):
            event = self.db.comcat_event(p["comcat_id"])
            eqs[i] = (event["longitude"], event["latitude"], event["magnitude"], p[metric])
            
        extent = [
            numpy.min(eqs["longitude"])-2.0, numpy.max(eqs["longitude"])+2.0,
            numpy.min(eqs["latitude"])-0.5, numpy.max(eqs["latitude"])+0.5,
        ]

        figure = self._create_figure(extent)
        ax = figure.gca()
        ms = 0.05 * 10**(0.75*eqs["magnitude"])
        sc = ax.scatter(eqs["longitude"], eqs["latitude"], s=ms, c=eqs["metric"], cmap="plasma", vmin=0.0, vmax=1.0, transform=crs.Geodetic(), edgecolors="black", alpha=0.67, zorder=4)

        matplotlib_extras.axes.add_background_axes(figure, [0.025, 0.015, 0.11, 0.36])
        cbax = figure.add_axes([0.03, 0.02, 0.02, 0.33])
        colorbar = pyplot.colorbar(mappable=sc, cax=cbax)
        label = "Cost Savings: Area" if metric == "area_costsavings_eew" else "Cost Savings: Population"
        colorbar.set_label(label)
        ax.set_title(label)
        
        self._save(figure, metric)
        return

# End of file
