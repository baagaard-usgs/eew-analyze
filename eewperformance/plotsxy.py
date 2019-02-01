# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import dateutil.parser
from datetime import datetime
import numpy

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
import matplotlib.patches as patches

from osgeo import gdal, osr
import matplotlib_extras

from . import analysis_utils
from . import greatcircle

gdal.UseExceptions()

class EventFigures(object):
    """Plots of alert mag/loc error, MMI (obs vs pred), etc.
    """
    def __init__(self, config, event):
        """
        :type config: ConfigParser
        :param config: Configuration for application.
        """
        self.config = config        

        self.event = event
        return

    def _save(self, figure, label):
        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = analysis_utils.analysis_event_label(self.config, self.event["event_id"])
        filename += "-{}.png".format(label)
        figure.savefig(os.path.join(plotsDir, filename))
        return
    
    def alert_error(self, alerts):
        """Create map with observed MMI with contours.
        """
        originTime = numpy.datetime64(self.event["origin_time"])
        
        alertsMag = [alert["magnitude"] for alert in alerts]
        alertsTime = [numpy.datetime64(alert["timestamp"]) for alert in alerts]
        alertsLon = numpy.array([alert["longitude"] for alert in alerts])
        alertsLat = numpy.array([alert["latitude"] for alert in alerts])
        alertsDepthKm = numpy.array([alert["depth_km"] for alert in alerts])
        if len(alerts) > 0 and "num_stations" in alerts[0].keys():
            alertsNumStations = numpy.array([alert["num_stations"] for alert in alerts])
        else:
            alertsNumStations = None
        
        horizDistKmError = 1.0e-3*greatcircle.distance(self.event["longitude"], self.event["latitude"], alertsLon, alertsLat)

        if len(alertsTime) > 0:
            t = (alertsTime - originTime).astype("timedelta64[us]").astype("float32")/1.0e+6
        else:
            t = []                
        tmax = max(30.0, t[0]+30.0) if len(t) > 0 else 30.0
        
        figure = pyplot.figure(figsize=(10.0, 3.5))
        
        # Magnitude
        rectFactory = matplotlib_extras.axes.RectFactory(figure, ncols=2, margins=((0.7, 0, 5.0), (0.50, 0, 0.25)))
        ax = figure.add_axes(rectFactory.rect())
        ax.plot(t, alertsMag, marker="o", mec="c_red", mfc="c_ltred", lw=0, alpha=0.67)
        ax.set_xlabel("Time after origin time (s)")
        ax.set_ylabel("Magnitude (Mw)")
        ax.set_title("Magnitude")
        ax.set_xlim(0, tmax)
        ax.axhline(self.event["magnitude"], linestyle="--", linewidth=1.0, color="c_ltblue")
        ax.text(ax.get_xlim()[1], self.event["magnitude"], "ANSS", ha="right", va="bottom", color="c_ltblue")

        if not alertsNumStations is None:
            ax2 = ax.twinx()
            ax2.semilogy(t, alertsNumStations, lw=0.5, color="c_green")
            ax2.set_ylabel("# Stations")
            ax2.set_xlim(0, tmax)
        
        # Horizontal location error
        rectFactory = matplotlib_extras.axes.RectFactory(figure, ncols=2, margins=((4.5, 0.8, 0.2), (0.50, 0, 0.25)))
        ax = figure.add_axes(rectFactory.rect(col=1))
        ax.plot(t, horizDistKmError, marker="o", mec="c_red", mfc="c_ltred", lw=0, alpha=0.67)
        ax.axhline(0.0, linestyle="--", linewidth=1.0, color="c_ltblue")
        ax.set_xlabel("Time after origin time (s)")
        ax.set_ylabel("Distance (km)")
        ax.set_xlim(0, tmax)
        ax.set_title("Horiz. Location Error")
        
        # Depth location error
        ax = figure.add_axes(rectFactory.rect(col=2))
        ax.plot(t, self.event["depth_km"]-alertsDepthKm, marker="o", mec="c_red", mfc="c_ltred", lw=0, alpha=0.67)
        ax.axhline(0.0, linestyle="--", linewidth=1.0, color="c_ltblue")
        ax.set_xlabel("Time after origin time (s)")
        ax.set_ylabel("Distance (km)")
        ax.set_xlim(0, tmax)
        ax.set_title("Depth Error (Obs-Pred)")

        self._save(figure, "alert_error")
        pyplot.close(figure)
        return

    def mmi_correlation(self):
        """Plot observed versus predicted MMI.
        """
        cacheDir = self.config.get("files", "analysis_cache_dir")
        filename = "analysis_" + analysis_utils.analysis_event_label(self.config, self.event["event_id"]) + ".tiff"
        rasterData = gdal.Open(os.path.join(cacheDir, filename), gdal.GA_ReadOnly)

        layers = {}
        for iband in range(rasterData.RasterCount):
            band = rasterData.GetRasterBand(1+iband)
            description = band.GetDescription()
            data = numpy.array(band.ReadAsArray())
            data = numpy.ma.masked_values(data, band.GetNoDataValue())
            layers[description] = data


        if numpy.isscalar(layers["mmi_pred"].mask):
            mmiObs = layers["mmi_obs"].ravel().data
            mmiPred = layers["mmi_pred"].ravel().data
        else:
            mask = ~layers["mmi_pred"].ravel().mask
            mmiObs = layers["mmi_obs"].ravel()[mask]
            mmiPred = layers["mmi_pred"].ravel()[mask]

        figure = pyplot.figure(figsize=(4.0, 4.0))
        rectFactory = matplotlib_extras.axes.RectFactory(figure, margins=((0.60, 0, 0.2), (0.45, 0, 0.1)))
        
        # Correlation
        maxMMI = 10.0
        if mmiObs.shape[0] > 0:
            maxMMI = numpy.maximum(5.0, numpy.maximum(numpy.max(mmiPred), numpy.max(mmiObs)))
        ax = figure.add_axes(rectFactory.rect())
        ax.plot(mmiPred, mmiObs, marker="o", ms=2, mec="c_red", mfc="c_ltred", lw=0, alpha=0.67, zorder=1)
        ax.plot([1,maxMMI],[1,maxMMI], "--", color="c_ltblue", zorder=2)

        ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
        ax.xaxis.set_ticks_position("bottom")

        ax.set_xlabel("Predicted MMI")
        ax.set_ylabel("Observed MMI")
        ax.set_xlim(1, maxMMI)
        ax.set_ylim(1, maxMMI)
        ax.set_aspect("equal")

        if mmiObs.shape[0] > 0:
        
            # Mean and std in bins
            bwidth = 0.5
            bins = numpy.arange(1.0-0.5*bwidth, 10.01+0.5*bwidth, bwidth)
            count, bedges = numpy.histogram(mmiPred, bins=bins)
            sum1, bedges = numpy.histogram(mmiPred, bins=bins, weights=mmiObs)
            sum2, bedges = numpy.histogram(mmiPred, bins=bins, weights=mmiObs**2)
            mmiMean = sum1 / count
            mmiStd = numpy.sqrt(sum2/count - mmiMean**2)
            bcenters = 0.5*(bedges[1:] + bedges[:-1])
            ax.errorbar(bcenters, mmiMean, yerr=mmiStd, fmt="none", ecolor="c_ltgreen", elinewidth=2, capthick=1.5, capwidth=6.0, zorder=4)
            ax.plot(bcenters, mmiMean, marker="s", lw=0, ms=6, mec="c_green", mfc="c_ltgreen", zorder=5)
        
            # Inset with residual histogram
            fontsize = 8
            bwidth = 0.25
            bins = numpy.arange(-3.0-0.5*bwidth, +3.001+0.5*bwidth, bwidth)
            residual = mmiObs-mmiPred
            residualMean = numpy.mean(residual)
            residualStd = numpy.std(residual)
            axin = inset_axes(ax, width="33%", height="25%", loc=2, borderpad=1.7)
            axin.hist(residual, bins=bins, density=True, align="mid", color="c_ltred", ec="c_red")
            axin.set_yticks([])
            axin.xaxis.set_ticks_position("bottom")
            for label in axin.xaxis.get_ticklabels():
                label.set_fontsize(fontsize)
            axin.set_xlim(numpy.min(bins), numpy.max(bins))
            axin.text(0.05, 0.95, "mean={m:.2f}\nstd={s:.2f}".format(m=residualMean, s=residualStd), transform=axin.transAxes, va="top", ha="left", fontsize=fontsize)
            axin.set_title("Residual (Obs-Pred)", fontsize=fontsize)

        self._save(figure, "mmi_correlation")
        pyplot.close(figure)
        return
    
    def warning_time_mmi(self):
        """Plot warning time versus observed MMI.
        """
        cacheDir = self.config.get("files", "analysis_cache_dir")
        filename = "analysis_" + analysis_utils.analysis_event_label(self.config, self.event["event_id"]) + ".tiff"
        rasterData = gdal.Open(os.path.join(cacheDir, filename), gdal.GA_ReadOnly)

        layers = {}
        for iband in range(rasterData.RasterCount):
            band = rasterData.GetRasterBand(1+iband)
            description = band.GetDescription()
            data = numpy.array(band.ReadAsArray())
            data = numpy.ma.masked_values(data, band.GetNoDataValue())
            layers[description] = data

        mask = ~layers["warning_time"].ravel().mask
        mmiObs = layers["mmi_obs"].ravel().data[mask]
        warningTime = layers["warning_time"].ravel().data[mask]

        figure = pyplot.figure(figsize=(4.5, 3.5))
        rectFactory = matplotlib_extras.axes.RectFactory(figure, margins=((0.60, 0, 0.2), (0.5, 0, 0.1)))
        
        ax = figure.add_axes(rectFactory.rect())
        mask = warningTime >= 0
        ax.plot(mmiObs[mask], warningTime[mask], marker="o", ms=2, mec="c_red", mfc="c_ltred", lw=0, alpha=0.67, zorder=1)
        ax.plot(mmiObs[~mask], warningTime[~mask], marker="o", ms=2, mec="black", mfc="c_ltgray", lw=0, alpha=0.67, zorder=1)
        ax.set_xlabel("Observed MMI")
        ax.set_xlim(1, 10)
        ax.set_ylabel("Warning Time (s)")

        ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
        
        if mmiObs.shape[0] > 0:
        
            # Mean and std in bins
            bwidth = 0.5
            bins = numpy.arange(1.0-0.5*bwidth, 10.01+0.5*bwidth, bwidth)
            count, bedges = numpy.histogram(mmiObs, bins=bins)
            sum1, bedges = numpy.histogram(mmiObs, bins=bins, weights=warningTime)
            sum2, bedges = numpy.histogram(mmiObs, bins=bins, weights=warningTime**2)
            wtMean = sum1 / count
            wtStd = numpy.sqrt(sum2/count - wtMean**2)
            bcenters = 0.5*(bedges[1:] + bedges[:-1])
            ax.errorbar(bcenters, wtMean, yerr=wtStd, fmt="none", ecolor="c_ltgreen", elinewidth=2, capthick=1.5, capwidth=6.0, zorder=4)
            ax.plot(bcenters, wtMean, marker="s", lw=0, ms=6, mec="c_green", mfc="c_ltgreen", zorder=5)
        
        self._save(figure, "warning_time_mmi")
        pyplot.close(figure)
        return
    


class SummaryFigures(object):
    """Plots of alert mag/loc error, MMI (obs vs pred), etc.
    """
    COLORS = {
        "area_costsavings_eew": ("local:qarea_fc", "local:qarea_ec",),
        "area_costsavings_perfecteew": ("none", "local:qarea_ec",),
        "population_costsavings_eew": ("local:qpop_fc", "local:qpop_ec",),
        "population_costsavings_perfecteew": ("none", "local:qpop_ec",),
    }

    def __init__(self, config, events, db):
        """
        :type config: ConfigParser
        :param config: Configuration for application.
        """
        self.config = config
        self.events = events
        self.db = db
        return

    def _save(self, figure, label):
        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = "eqset_" + analysis_utils.analysis_label(self.config)
        filename += "_{}.png".format(label)
        figure.savefig(os.path.join(plotsDir, filename))
        return
    
    def optimal_mmithresholds(self):
        """Plot Q-area and Q-pop versus magnitude and MMI thresholds to illustrate optimum thresholds.
        """
        server = self.config.get("shakealert.production", "server")
        gmpe = self.config.get("mmi_predicted", "gmpe")
        fragility = self.config.get("fragility_curves", "object").split(".")[-1]
        
        thresholdStart = self.config.getfloat("optimize", "mmi_threshold_min")
        thresholdStop = self.config.getfloat("optimize", "mmi_threshold_max")
        thresholdStep = self.config.getfloat("optimize", "mmi_threshold_step")
        mmiThresholds = numpy.arange(thresholdStart, thresholdStop+0.1*thresholdStep, thresholdStep)

        thresholdStart = self.config.getfloat("optimize", "magnitude_threshold_min")
        thresholdStop = self.config.getfloat("optimize", "magnitude_threshold_max")
        thresholdStep = self.config.getfloat("optimize", "magnitude_threshold_step")
        magThresholds = numpy.arange(thresholdStart, thresholdStop+0.1*thresholdStep, thresholdStep)

        perfs = numpy.array([self.db.performance_stats(eqId, server, gmpe, fragility) for eqId in self.events]).ravel()
        
        dtype = [
            ("area_metric", "float32",),
            ("population_metric", "float32",),
        ]
        metricAllEqs = numpy.zeros((mmiThresholds.shape[0], magThresholds.shape[0]), dtype=dtype)
        for iMag,magThreshold in enumerate(magThresholds):
            maskMag = numpy.ma.masked_values(perfs["magnitude_threshold"], magThreshold).mask
            perfsMag = perfs[maskMag]
            for iMMI,mmiThreshold in enumerate(mmiThresholds):
                maskMMI = numpy.ma.masked_values(perfsMag["mmi_threshold"], mmiThreshold).mask
                perfsMMI = perfsMag[maskMMI]
                if 0 == len(perfsMMI["area_costsavings_eew"]):
                    raise ValueError("Missing performance data for MMI threshold {}".format(mmiThreshold))

                areaMetric = numpy.sum(perfsMMI["area_costsavings_eew"]) / numpy.sum(perfsMMI["area_costsavings_perfecteew"])
                popMetric = numpy.sum(perfsMMI["population_costsavings_eew"]) / numpy.sum(perfsMMI["population_costsavings_perfecteew"])
                
                metricAllEqs[iMMI, iMag] = (areaMetric, popMetric,)
                

        figure = pyplot.figure(figsize=(8.0, 3.5))
        rectFactory = matplotlib_extras.axes.RectFactory(figure, nrows=1, ncols=2, margins=((0.1, 0.5, 0.5), (0.5, 0, 0.2)))
        magOffset = 0.05
        barW = 0.5

        from mpl_toolkits.mplot3d import Axes3D
        # x and y are transposed from metricAllEqs
        x, y = numpy.meshgrid(mmiThresholds-0.25*barW, magThresholds-0.25*barW)
        z = numpy.zeros(x.shape)
        dx = barW*self.config.getfloat("optimize", "mmi_threshold_step")
        dy = barW*self.config.getfloat("optimize", "magnitude_threshold_step")
        c = numpy.zeros(metricAllEqs.shape, dtype=object)
        indices = numpy.indices(metricAllEqs.shape)
        indicesMMI = indices[0,:].ravel()
        indicesMag = indices[1,:].ravel()
        del indices
        
        # Q-area
        metric = "area_metric"
        metricMasked = numpy.ma.masked_less(metricAllEqs[metric], 0.0)
        c[metricAllEqs[metric] > 0.0] = "c_ltred"
        c[metricAllEqs[metric] <= 0.0] = "c_mdgray"

        iMax = numpy.argmax(metricAllEqs[metric].ravel())
        areaMetric = metricAllEqs[metric][indicesMMI[iMax], indicesMag[iMax]]
        areaOptMMI = 0.25*barW + x.T[indicesMMI[iMax], indicesMag[iMax]]
        areaOptMag = 0.25*barW + magOffset + y.T[indicesMMI[iMax], indicesMag[iMax]]
        c[indicesMMI[iMax], indicesMag[iMax]] = "c_ltblue"

        ax = figure.add_axes(rectFactory.rect(row=1, col=1), projection="3d")
        ax.bar3d(x.ravel(), y.ravel(), z.ravel(), dx, dy, metricMasked.ravel("F"), color=c.ravel("F"), zsort="max")
        ax.set_title("Optimal Thresholds for Q-area")
        ax.set_xlabel("MMI Threshold")
        ax.set_ylabel("Magnitude Threshold")
        ax.set_zlabel("Q-area")
        ax.set_zlim(0, 1)
        ax.set_xticks(mmiThresholds)
        ax.set_yticks(magThresholds+magOffset)
        
        # Q-pop
        metric = "population_metric"
        metricMasked = numpy.ma.masked_less(metricAllEqs[metric], 0.0)
        c[metricAllEqs[metric] > 0.0] = "c_ltred"
        c[metricAllEqs[metric] <= 0.0] = "c_mdgray"

        iMax = numpy.argmax(metricAllEqs[metric].ravel())
        popMetric = metricAllEqs[metric][indicesMMI[iMax], indicesMag[iMax]]
        popOptMMI = 0.25*barW + x.T[indicesMMI[iMax], indicesMag[iMax]]
        popOptMag = 0.25*barW + magOffset + y.T[indicesMMI[iMax], indicesMag[iMax]]
        c[indicesMMI[iMax], indicesMag[iMax]] = "c_ltblue"

        ax = figure.add_axes(rectFactory.rect(row=1, col=2), projection="3d")
        ax.bar3d(x.ravel(), y.ravel(), z.ravel(), dx, dy, metricMasked.ravel("F"), color=c.ravel("F"), zsort="max")
        ax.set_title("Optimal Thresholds for Q-pop")
        ax.set_xlabel("MMI Threshold")
        ax.set_ylabel("Magnitude Threshold")
        ax.set_zlabel("Q-pop")
        ax.set_zlim(0, 1)
        ax.set_xticks(mmiThresholds)
        ax.set_yticks(magThresholds+magOffset)

        self._save(figure, "optimal_threshold")
        pyplot.close(figure)

        print("Q-area: {:.2f}, Magnitude threshold: {:.1f}, MMI threshold: {:.1f}".format(areaMetric, areaOptMag, areaOptMMI))
        print("Q-pop: {:.2f}, Magnitude threshold: {:.1f}, MMI threshold: {:.1f}".format(popMetric, popOptMag, popOptMMI))

        print(metricAllEqs)
        return
    

    def costsavings_versus_time(self):
        """Plot costsavings-area and costsavings-pop versus ANSS origin
        time. Cost savings is shown with filled circles relative to
        the ideal case of perfect EEW shown with hollow circles.

        """
        FIG_SIZE = (8.0, 4.5)
        MARGINS = ((1.1, 0, 0.1), (0.5, 0.7, 0.3)) 
        nrows = 2
        ncols = 1

        figure = pyplot.figure(figsize=FIG_SIZE)
        rectFactory = matplotlib_extras.axes.RectFactory(figure, nrows=nrows, ncols=ncols, margins=MARGINS)
        
        server = self.config.get("shakealert.production", "server")
        gmpe = self.config.get("mmi_predicted", "gmpe")
        fragility = self.config.get("fragility_curves", "object").split(".")[-1]
        mmiThreshold = self.config.getfloat("alerts", "mmi_threshold")
        magThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        
        perfs = numpy.array([self.db.performance_stats(eqId, server, gmpe, fragility, magThreshold, mmiThreshold) for eqId in self.events]).ravel()
        
        originTime = numpy.zeros(perfs.shape, dtype="datetime64[s]")
        magnitude = numpy.zeros(perfs.shape, dtype=numpy.float32)
        for i,p in enumerate(perfs):
            event = self.db.comcat_event(p["comcat_id"])
            originTime[i] = numpy.datetime64(dateutil.parser.parse(event["origin_time"]))
            magnitude[i] = event["magnitude"]

        ms = 5.0e-4 * 10**magnitude

        # Q-area
        ax = figure.add_axes(rectFactory.rect(row=1))
        metrics = ["area_costsavings_eew", "area_costsavings_perfecteew"]
        labels = ["ShakeAlert", "Perfect EEW"]
        for metric, label in zip(metrics, labels):
            fc, ec = self.COLORS[metric]
            ax.scatter(originTime.astype(datetime), perfs[metric], s=ms, c=fc, edgecolors=ec, alpha=0.67, label=label)
        ax.set_title("Q-area vs. Origin Time", weight="bold")
        ax.set_ylim(0.0, ax.get_ylim()[1])
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%4.0f"))
        ax.yaxis.set_label_text("Q-area Cost Savings")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        pyplot.legend(handlelength=0.8, borderpad=0.3, labelspacing=0.2, loc="upper left")

        # Q-pop
        ax = figure.add_axes(rectFactory.rect(row=2))
        metrics = ["population_costsavings_eew", "population_costsavings_perfecteew"]
        for metric, label in zip(metrics, labels):
            fc, ec = self.COLORS[metric]
            ax.scatter(originTime.astype(datetime), perfs[metric], s=ms, c=fc, edgecolors=ec, alpha=0.67, label=label)
        ax.set_title("Q-pop vs. Origin Time", weight="bold")
        ax.set_ylim(0.0, ax.get_ylim()[1])
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%7.1e"))
        ax.yaxis.set_label_text("Q-pop Cost Savings")
        ax.xaxis.set_label_text("Origin Time (UTC)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        pyplot.legend(handlelength=0.8, borderpad=0.3, labelspacing=0.2, loc="upper left")

        self._save(figure, "costsavings_time")
        pyplot.close(figure)
        return
    

    def costsavings_versus_magnitude(self):
        """Plot Q-area and Q-pop versus magnitude.
        """
        FIG_SIZE = (8.0, 4.5)
        MARGINS = ((1.1, 0, 0.1), (0.5, 0.7, 0.3)) 
        nrows = 2
        ncols = 1

        figure = pyplot.figure(figsize=FIG_SIZE)
        rectFactory = matplotlib_extras.axes.RectFactory(figure, nrows=nrows, ncols=ncols, margins=MARGINS)
        
        server = self.config.get("shakealert.production", "server")
        gmpe = self.config.get("mmi_predicted", "gmpe")
        fragility = self.config.get("fragility_curves", "object").split(".")[-1]
        mmiThreshold = self.config.getfloat("alerts", "mmi_threshold")
        magThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        
        perfs = numpy.array([self.db.performance_stats(eqId, server, gmpe, fragility, magThreshold, mmiThreshold) for eqId in self.events]).ravel()

        magnitude = numpy.zeros(perfs.shape, dtype=numpy.float32)
        for i,p in enumerate(perfs):
            event = self.db.comcat_event(p["comcat_id"])
            magnitude[i] = event["magnitude"]

        ms = 5.0e-4 * 10**magnitude

        # Q-area
        ax = figure.add_axes(rectFactory.rect(row=1))
        metrics = ["area_costsavings_eew", "area_costsavings_perfecteew"]
        labels = ["ShakeAlert", "Perfect EEW"]
        for metric, label in zip(metrics, labels):
            fc, ec = self.COLORS[metric]
            ax.scatter(magnitude, perfs[metric], s=ms, c=fc, edgecolors=ec, alpha=0.67, label=label)
        ax.set_title("Q-area vs. Earthquake Magnitude", weight="bold")
        ax.set_ylim(0.0, ax.get_ylim()[1])
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%4.0f"))
        ax.yaxis.set_label_text("Q-area Cost Savings")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%3.1f"))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
        pyplot.legend(handlelength=0.8, borderpad=0.3, labelspacing=0.2, loc="upper left")

        # Q-pop
        ax = figure.add_axes(rectFactory.rect(row=2))
        metrics = ["population_costsavings_eew", "population_costsavings_perfecteew"]
        for metric, label in zip(metrics, labels):
            fc, ec = self.COLORS[metric]
            ax.scatter(magnitude, perfs[metric], s=ms, c=fc, edgecolors=ec, alpha=0.67, label=label)
        ax.set_title("Q-pop vs. Earthquake Magnitude", weight="bold")
        ax.set_ylim(0.0, ax.get_ylim()[1])
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%7.1e"))
        ax.yaxis.set_label_text("Q-pop Cost Savings")
        ax.xaxis.set_label_text("Earthquake Magnitude (Mw)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%3.1f"))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
        pyplot.legend(handlelength=0.8, borderpad=0.3, labelspacing=0.2, loc="upper left")

        self._save(figure, "costsavings_magnitude")
        pyplot.close(figure)
        return
    

    def magnitude_versus_time(self):
        """Plot magnitude versus time.
        """
        FIG_SIZE = (8.0, 3.5)
        MARGINS = ((0.6, 0, 0.1), (0.5, 0, 0.3))
        server = self.config.get("shakealert.production", "server")
        gmpe = self.config.get("mmi_predicted", "gmpe")
        fragility = self.config.get("fragility_curves", "object").split(".")[-1]

        numEvents = len(self.events)
        originTime = numpy.zeros(numEvents, dtype="datetime64[s]")
        magnitude = numpy.zeros(numEvents, dtype=numpy.float32)
        for i,eqId in enumerate(self.events):
            event = self.db.comcat_event(eqId)
            originTime[i] = numpy.datetime64(dateutil.parser.parse(event["origin_time"]))
            magnitude[i] = event["magnitude"]

        from matplotlib.dates import YearLocator,date2num,DateFormatter
        figure = pyplot.figure(figsize=FIG_SIZE)
        rectFactory = matplotlib_extras.axes.RectFactory(figure, margins=MARGINS)
        
        ax = figure.add_axes(rectFactory.rect())
        ms = 5.0e-4 * 10**magnitude
        ot = originTime.astype(datetime)
        ax.scatter(ot, magnitude, s=ms, edgecolors="black", c=date2num(ot), cmap="viridis", alpha=0.67)
        ax.set_title("Magnitude versus Earthquake Origin Time")
        ax.set_xlabel("Origin Time (UTC)")
        ax.set_ylabel("Moment Magnitude")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%3.1f"))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))

        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = "eqset_magnitude_time.png"
        figure.savefig(os.path.join(plotsDir, filename))
        pyplot.close(figure)
        
        return
    
# End of file
