# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import numpy

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as pyplot

from osgeo import gdal, osr
from basemap.Figure import Figure

import analysis_utils
import greatcircle

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
    
    def alert_error(self, alerts):
        """Create map with observed MMI with contours.
        """
        originTime = numpy.datetime64(self.event["origin_time"])
        
        alertsMag = [alert["magnitude"] for alert in alerts]
        alertsTime = [numpy.datetime64(alert["timestamp"]) for alert in alerts]
        alertsLon = numpy.array([alert["longitude"] for alert in alerts])
        alertsLat = numpy.array([alert["latitude"] for alert in alerts])
        alertsDepthKm = numpy.array([alert["depth_km"] for alert in alerts])
        
        horizDistKmError = 1.0e-3*greatcircle.distance(self.event["longitude"], self.event["latitude"], alertsLon, alertsLat)

        if len(alertsTime) > 0:
            t = (alertsTime - originTime).astype("timedelta64[us]").astype("float32")/1.0e+6
        else:
            t = []
                
        figure = Figure()
        figure.open(6.0, 2.5, margins=((0.45, 0.65, 0.2), (0.4, 0.8, 0.3)))
        nrows = 1
        ncols = 3
        irow = 1
        icol = 1

        tmax = max(30.0, t[0]+30.0) if len(t) > 0 else 30.0
        
        # Magnitude
        ax = figure.axes(nrows, ncols, irow, icol)
        ax.plot(t, alertsMag, marker="o", mfc="c_ltred", mec="c_fg", lw=0, alpha=0.5)
        ax.set_xlabel("Time after origin time (s)")
        ax.set_ylabel("Magnitude (Mw)")
        ax.set_title("Magnitude")
        ax.set_xlim(0, tmax)
        ax.axhline(self.event["magnitude"], linestyle="--", linewidth=1.0, color="c_blue")
        ax.text(ax.get_xlim()[1], self.event["magnitude"], "ANSS", ha="right", va="bottom", color="c_blue")
        icol += 1
        
        # Horizontal location error
        ax = figure.axes(nrows, ncols, irow, icol)
        ax.plot(t, horizDistKmError, marker="o", mfc="c_ltred", mec="c_fg", lw=0, alpha=0.5)
        ax.axhline(0.0, linestyle="--", linewidth=1.0, color="c_blue")
        ax.set_xlabel("Time after origin time (s)")
        ax.set_ylabel("Distance (km)")
        ax.set_xlim(0, tmax)
        ax.set_title("Horiz. Location Error")
        icol += 1
        
        # Depth location error
        ax = figure.axes(nrows, ncols, irow, icol)
        ax.plot(t, self.event["depth_km"]-alertsDepthKm, marker="o", mfc="c_ltred", mec="c_fg", lw=0, alpha=0.5)
        ax.axhline(0.0, linestyle="--", linewidth=1.0, color="c_blue")
        ax.set_xlabel("Time after origin time (s)")
        ax.set_ylabel("Distance (km)")
        ax.set_xlim(0, tmax)
        ax.set_title("Depth Error (Obs-Pred)")
        icol += 1

        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = analysis_utils.analysis_label(self.config, self.event["event_id"])
        filename += "-alert_error.png"
        figure.figure.savefig(os.path.join(plotsDir, filename))
        figure.close()
        return

    def mmi_correlation(self):
        """Plot observed versus predicted MMI.
        """
        cacheDir = self.config.get("files", "analysis_cache_dir")
        filename = "analysis_" + analysis_utils.analysis_label(self.config, self.event["event_id"]) + ".tiff"
        rasterData = gdal.Open(os.path.join(cacheDir, filename), gdal.GA_ReadOnly)

        layers = {}
        for iband in range(rasterData.RasterCount):
            band = rasterData.GetRasterBand(1+iband)
            description = band.GetDescription()
            data = numpy.array(band.ReadAsArray())
            data = numpy.ma.masked_values(data, band.GetNoDataValue())
            layers[description] = data

        figure = Figure()
        figure.open(4.0, 4.0, margins=((0.45, 0.65, 0.2), (0.4, 0.8, 0.3)))
        nrows = 1
        ncols = 1
        irow = 1
        icol = 1

        if numpy.isscalar(layers["mmi_pred"].mask):
            mmiObs = layers["mmi_obs"].ravel().data
            mmiPred = layers["mmi_pred"].ravel().data
        else:
            mask = ~layers["mmi_pred"].ravel().mask
            mmiObs = layers["mmi_obs"].ravel()[mask]
            mmiPred = layers["mmi_pred"].ravel()[mask]

        # Correlation
        maxMMI = 10.0
        if mmiObs.shape[0] > 0:
            maxMMI = numpy.maximum(5.0, numpy.maximum(numpy.max(mmiPred), numpy.max(mmiObs)))
        ax = figure.axes(nrows, ncols, irow, icol)
        ax.plot(mmiPred, mmiObs, marker="o", ms=2, mfc="c_ltred", mec="c_fg", lw=0, alpha=0.5, zorder=1)
        ax.plot([1,maxMMI],[1,maxMMI], "--", color="c_ltblue", zorder=2)
        ax.set_xlabel("Predicted MMI")
        ax.xaxis.set_ticks_position("bottom")
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
            bwidth = 0.25
            bins = numpy.arange(-3.0-0.5*bwidth, +3.001+0.5*bwidth, bwidth)
            residual = mmiObs-mmiPred
            residualMean = numpy.mean(residual)
            residualStd = numpy.std(residual)
            axin = inset_axes(ax, width="33%", height="25%", loc=2, borderpad=1.7)
            axin.hist(residual, bins=bins, normed=True, align="mid", color="c_ltred")
            axin.set_yticks([])
            axin.xaxis.set_ticks_position("bottom")
            axin.set_xlim(numpy.min(bins), numpy.max(bins))
            axin.text(0.05, 0.95, "mean={m:.2f}\nstd={s:.2f}".format(m=residualMean, s=residualStd), transform=axin.transAxes, va="top", ha="left", fontsize=6)
            axin.set_title("Residual (Obs-Pred)", fontsize=8)

        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = analysis_utils.analysis_label(self.config, self.event["event_id"])
        filename += "-mmi_correlation.png"
        figure.figure.savefig(os.path.join(plotsDir, filename))
        figure.close()
        
        return
    

class SummaryFigures(object):
    """Plots of alert mag/loc error, MMI (obs vs pred), etc.
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

        perfs = None
        for eqId in self.events:
            p = self.db.performance_stats(eqId, server, gmpe, fragility)
            if perfs is None:
                perfs = p
            else:
                perfs = numpy.append(perfs, p)
        
        dtype = [
            ("area_metric", "float32",),
            ("population_metric", "float32",),
            ("area_metric_posonly", "float32",),
            ("population_metric_posonly", "float32",),
        ]
        metricAllEqs = numpy.zeros((mmiThresholds.shape[0], magThresholds.shape[0]), dtype=dtype)
        for imag,magThreshold in enumerate(magThresholds):
            maskMag = numpy.ma.masked_values(perfs["magnitude_threshold"], magThreshold).mask
            perfsMag = perfs[maskMag]
            for immi,mmiThreshold in enumerate(mmiThresholds):
                maskMMI = numpy.ma.masked_values(perfsMag["mmi_threshold"], mmiThreshold).mask
                perfsMMI = perfsMag[maskMMI]

                savingsEEW = numpy.sum(perfsMMI["area_cost_eew"]) - numpy.sum(perfsMMI["area_cost_noeew"])
                savingsPerfectEEW = numpy.sum(perfsMMI["area_cost_perfecteew"]) - numpy.sum(perfsMMI["area_cost_noeew"])
                areaMetric = savingsEEW / savingsPerfectEEW
                numEvents = len(perfsMMI["area_cost_eew"])
                
                savingsEEW = numpy.sum(perfsMMI["population_cost_eew"]) - numpy.sum(perfsMMI["population_cost_noeew"])
                savingsPerfectEEW = numpy.sum(perfsMMI["population_cost_perfecteew"]) - numpy.sum(perfsMMI["population_cost_noeew"])
                popMetric = savingsEEW / savingsPerfectEEW

                # Limit to events with cost savings
                maskQ =  perfsMMI["area_cost_eew"] > perfsMMI["area_cost_noeew"]
                savingsEEW = numpy.sum(perfsMMI["area_cost_eew"][maskQ]) - numpy.sum(perfsMMI["area_cost_noeew"][maskQ])
                savingsPerfectEEW = numpy.sum(perfsMMI["area_cost_perfecteew"][maskQ]) - numpy.sum(perfsMMI["area_cost_noeew"][maskQ])
                areaMetricPosOnly = savingsEEW / savingsPerfectEEW
                numEventsPosOnly = numpy.sum(maskQ)
                
                maskQ =  perfsMMI["population_cost_eew"] > perfsMMI["population_cost_noeew"]
                savingsEEW = numpy.sum(perfsMMI["population_cost_eew"][maskQ]) - numpy.sum(perfsMMI["population_cost_noeew"][maskQ])
                savingsPerfectEEW = numpy.sum(perfsMMI["population_cost_perfecteew"][maskQ]) - numpy.sum(perfsMMI["population_cost_noeew"][maskQ])
                popMetricPosOnly = savingsEEW / savingsPerfectEEW
                
                metricAllEqs[immi, imag] = (areaMetric, popMetric, areaMetricPosOnly, popMetricPosOnly)

        figure = Figure()
        figure.open(6.0, 2.5, margins=((0.45, 0.65, 0.2), (0.4, 0.8, 0.3)))
        nrows = 1
        ncols = 2
        irow = 1
        icol = 1
        magOffset = 0.05
        barW = 0.5
        
        from mpl_toolkits.mplot3d import Axes3D
        x, y = numpy.meshgrid(mmiThresholds-0.25*barW, magThresholds-0.25*barW)
        z = numpy.zeros(x.shape)
        dx = barW*self.config.getfloat("optimize", "mmi_threshold_step")
        dy = barW*self.config.getfloat("optimize", "magnitude_threshold_step")
        c = numpy.zeros(metricAllEqs.shape, dtype=object)

        # Q-area
        index = numpy.argmax(metricAllEqs["area_metric"].ravel())
        areaMetric = metricAllEqs["area_metric"].transpose().ravel()[index]
        areaOptMMI = x.ravel()[index]
        areaOptMag = y.ravel()[index]

        ax = figure.axes(nrows, ncols, irow, icol, projection="3d")
        dz = numpy.ma.masked_less(metricAllEqs["area_metric"], 0.0).transpose()
        c[metricAllEqs["area_metric"] > 0.0] = "c_ltred"
        c[metricAllEqs["area_metric"] <= 0.0] = "c_bg"
        color = c.transpose().ravel()
        color[index] = "c_ltblue"
        ax.bar3d(x.ravel(), y.ravel(), z.ravel(), dx, dy, dz.ravel(), color=color)
        ax.set_xlabel("MMI Threshold")
        ax.set_ylabel("Magnitude Threshold")
        ax.set_zlabel("Q-area")
        ax.set_zlim(0, 1)
        ax.set_xticks(mmiThresholds)
        ax.set_yticks(magThresholds+magOffset)
        icol += 1
        
        # Q-pop
        index = numpy.argmax(metricAllEqs["population_metric"].ravel())
        popMetric = metricAllEqs["population_metric"].ravel()[index]
        popOptMMI = x.ravel()[index]
        popOptMag = y.ravel()[index]

        ax = figure.axes(nrows, ncols, irow, icol, projection="3d")
        dz = numpy.ma.masked_less(metricAllEqs["population_metric"], 0.0).transpose()
        c[metricAllEqs["population_metric"] > 0.0] = "c_ltred"
        c[metricAllEqs["population_metric"] <= 0.0] = "c_bg"
        color = c.transpose().ravel()
        color[index] = "c_ltblue"
        ax.bar3d(x.ravel(), y.ravel(), z.ravel(), dx, dy, dz.ravel(), color=color)
        ax.set_xlabel("MMI Threshold")
        ax.set_ylabel("Magnitude Threshold")
        ax.set_zlabel("Q-pop")
        ax.set_zlim(0, 1)
        ax.set_xticks(mmiThresholds)
        ax.set_yticks(magThresholds+magOffset)
        icol += 1
        
        plotsDir = self.config.get("files", "plots_dir")
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        filename = ""#analysis_utils.analysis_label(self.config, self.event["event_id"])
        filename += "optimal_threshold.pdf"
        figure.figure.savefig(os.path.join(plotsDir, filename))
        figure.close()


        print("# earthquakes: {}".format(numEvents))
        print("Q-area: {:.2f}, Magnitude threshold: {:.1f}, MMI threshold: {:.1f}".format(areaMetric, areaOptMag, areaOptMMI))

        return
    
# End of file
