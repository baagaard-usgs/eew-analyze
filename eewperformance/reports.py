# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import datetime
import dateutil.parser
import numpy

from PIL import Image

from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.platypus import Table
from reportlab.lib.utils import ImageReader

from .analysisdb import AnalysisData
from .analysis_utils import analysis_event_label, analysis_label, timedelta_to_seconds

class AnalysisSummary(object):
    """Summary of analysis.
    """
    HEADER = 0.5*inch
    MARGIN = 0.5*inch
    PAGE_WIDTH, PAGE_HEIGHT = landscape(letter)
    SPACING = 0.125*inch
    YTOP = PAGE_HEIGHT - MARGIN - HEADER
    YBOT = MARGIN
    XLEFT = MARGIN

    CATALOG_MAGNITUDE = "zero-latency-catalog-magnitude"
    CATALOG_MAGNITUDE_BIAS = "zero-latency-catalog-magnitude-bias"
    
    def __init__(self, config, summary_only=True):
        """Constructor.
        """
        self.config = config
        self.canvas = None
        self.summary_only = summary_only
        return

    def generate(self, eqIds):
        """Generate report.
        """
        filename = self.config.get("files", "report")
        self.canvas = canvas.Canvas(
            filename,
            pagesize=landscape(letter),
            bottomup=1,
        )

        eqIds = self._sort_eqs(eqIds)

        self._render_summary(eqIds)

        if not self.summary_only:
            for eqId in eqIds:
                self._render_event(eqId)

        self.canvas.save()
        return

    def _render_event(self, eqId):
        """Generate summary for event.
        """
        db = AnalysisData(self.config.get("files", "analysis_db"))
        event = db.comcat_event(eqId)
        shakemap = db.comcat_shakemap(eqId)
        
        # Get alerts
        server = self.config.get("shakealert.production", "server")
        alerts = db.alerts(eqId, server)
        
        # Get performance information
        #comcat_id, eew_server, gmpe, fragility, magnitude_threshold
        gmpe = self.config.get("mmi_predicted", "gmpe")
        fragility = self.config.get("fragility_curves", "label")
        alertLatency = self.config.getfloat("alerts", "alert_latency_sec")
        perfEEW = numpy.array(db.performance_stats(eqId, server, gmpe, fragility, alertLatency))
        perfTheoryMag = numpy.array(db.performance_stats(eqId, self.CATALOG_MAGNITUDE, gmpe, fragility, alert_latency_sec=0.0))
        perfTheoryMagBias = numpy.array(db.performance_stats(eqId, self.CATALOG_MAGNITUDE_BIAS, gmpe, fragility, alert_latency_sec=0.0))

        # Page 1
        self._render_event_header(event)
        self._render_event_info(event, shakemap, alerts)
        self._render_event_perf_table(perfEEW, perfTheoryMag, perfTheoryMagBias)
        self.canvas.showPage()

        # Page 2
        self._render_event_header(event)
        self._render_event_mmi_maps(event)
        self._render_event_alert_maps(event, server)
        self.canvas.showPage()

        # Page 3
        self._render_event_header(event)
        self._render_event_mmi_correlation(event)
        self._render_event_mmi_warningtime(event)
        self._render_event_error(event)
        self.canvas.showPage()

        return

    def _render_event_header(self, event):
        ot = dateutil.parser.parse(event["origin_time"])
        fragilityFn = self.config.get("fragility_curves", "label")
        magThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        mmiThreshold = self.config.getfloat("alerts", "mmi_threshold")
        header = (
            "{ot:%Y-%m-%d %H:%M:%S} {event[event_id]:s} M{event[magnitude]:.1f} {event[description]:s}".format(ot=ot, event=event),
            "Damage/Action: {fn}, Alert thresholds: M{mag:.1f}, MMI {mmi:.1f}".format(fn=fragilityFn, mag=magThreshold, mmi=mmiThreshold),
        )
        self.canvas.saveState()
        self.canvas.translate(self.MARGIN, self.PAGE_HEIGHT-self.MARGIN)
        text = self.canvas.beginText()
        text.setFont("Helvetica-Bold", 12)
        text.textLines("\n".join(header))
        self.canvas.drawText(text)
        self.canvas.restoreState()

        header = (
            "Brad Aagaard",
            "Version " + datetime.datetime.now().strftime("%Y-%m-%d %I:%M%p"),
        )
        self.canvas.saveState()
        self.canvas.translate(self.PAGE_WIDTH-self.MARGIN, self.PAGE_HEIGHT-self.MARGIN)
        self.canvas.setFont("Helvetica", 8)
        self.canvas.drawRightString(0, 0, header[0])
        self.canvas.drawRightString(0, -10, header[1])
        self.canvas.restoreState()
        return

    def _render_event_mmi_maps(self, event):
        """Maps of MMI observed, predicted, and residual (observed-predicted).
        """
        nrows = 2
        ncols = 3
        row = 1
        
        plots_dir = self.config.get("files", "plots_dir")
        label = analysis_event_label(self.config, event["event_id"])

        (x, y, mapW, mapH) = self._calc_figure_size(nrows, ncols, row=row, col=1)
        filename = os.path.join(plots_dir, label+"-map_mmi_obs.jpg")
        self._render_image(filename, x, y, width=mapW)
    
        (x, y, mapW, mapH) = self._calc_figure_size(nrows, ncols, row=row, col=2)
        filename = os.path.join(plots_dir, label+"-map_mmi_pred.jpg")
        self._render_image(filename, x, y, width=mapW)

        (x, y, mapW, mapH) = self._calc_figure_size(nrows, ncols, row=row, col=3)
        filename = os.path.join(plots_dir, label+"-map_mmi_residual.jpg")
        self._render_image(filename, x, y, width=mapW)

        return

    def _render_event_alert_maps(self, event, server):
        """Maps of alert categories w/warning time for EEW and
        theoretical with catalog magnitude and theoretical with catalog
        magnitude with event bias.
        """
        nrows = 2
        ncols = 3
        row = 2

        plots_dir = self.config.get("files", "plots_dir")
        label = analysis_event_label(self.config, event["event_id"])

        (x, y, mapW, mapH) = self._calc_figure_size(nrows, ncols, row=row, col=1)
        filename = os.path.join(plots_dir, label+"-map_cost_savings.jpg")
        imageWidth, imageHeight = self._render_image(filename, x, y, width=mapW)
        self._figure_label(x, y+imageHeight, "ShakeAlert")
        
        (x, y, mapW, mapH) = self._calc_figure_size(nrows, ncols, row=row, col=2)
        filenameMag = filename.replace(server, self.CATALOG_MAGNITUDE)
        if os.path.isfile(filenameMag):
            imageWidth, imageHeight, self._render_image(filenameMag, x, y, width=mapW)
            self._figure_label(x, y+imageHeight, "Theoretical Ideal: Catalog Mag.")

        (x, y, mapW, mapH) = self._calc_figure_size(nrows, ncols, row=row, col=3)
        filenameMagBias = filename.replace(server, self.CATALOG_MAGNITUDE_BIAS)
        if os.path.isfile(filenameMagBias):
            imageWidth, imageHeight = self._render_image(filenameMagBias, x, y, width=mapW)
            self._figure_label(x, y+imageHeight, "Theoretical Ideal: Catalog Mag. w/Bias")
        return

    def _render_event_error(self, event):
        """Figure of magnitude and location error.
        """
        (x, y, figW, figH) = self._calc_figure_size(nrows=2, ncols=1.3, row=2, col=1)

        plots_dir = self.config.get("files", "plots_dir")
        label = analysis_event_label(self.config, event["event_id"])

        filename = os.path.join(plots_dir, label+"-alert_error.png")
        self._render_image(filename, x, y, width=figW)
        return

    def _render_event_mmi_correlation(self, event):
        """Figure of MMI correlation (observed vs predicted) with inset
        of residual histogram.
        """
        (x, y, figW, figH) = self._calc_figure_size(nrows=2, ncols=2, row=1, col=1)

        plots_dir = self.config.get("files", "plots_dir")
        label = analysis_event_label(self.config, event["event_id"])

        filename = os.path.join(plots_dir, label+"-mmi_correlation.png")
        self._render_image(filename, x, y, height=figH)
        return

    def _render_event_mmi_warningtime(self, event):
        """Figure of MMI correlation (observed vs predicted) with inset
        of residual histogram.
        """
        (x, y, figW, figH) = self._calc_figure_size(nrows=2, ncols=2, row=1, col=2)

        plots_dir = self.config.get("files", "plots_dir")
        label = analysis_event_label(self.config, event["event_id"])

        filename = os.path.join(plots_dir, label+"-warning_time_mmi.png")
        self._render_image(filename, x, y, height=figH)
        return

    def _render_event_info(self, event, shakemap, alerts):
        x = self.XLEFT
        y = self.YTOP
        
        self.canvas.saveState()
        self.canvas.translate(x, y)
        
        ot = dateutil.parser.parse(event["origin_time"])
        text = self.canvas.beginText()
        text.setFont("Courier", 8)
        text.textLine("ANSS")
        text.textLine("   {event[magnitude_type]}{event[magnitude]:.1f}".format(event=event))
        text.textLine("   {ot:%Y-%m-%d %H:%M:%S.%f}".format(ot=ot))
        text.textLine("   {event[depth_km]:.1f} km depth".format(event=event))
        text.textLine("ShakeMap")
        text.textLine("   MMI bias {:6.2f}".format(shakemap["mmi_bias"]))
        text.textLine("")
        
        # First alert
        vs = self.config.getfloat("shaking_time", "vs_kmps")
        vp = self.config.getfloat("shaking_time", "vp_kmps")
        text.textLine("P wave arrival at epicenter: {:.1f}s after OT".format(event["depth_km"]/vp))
        if len(alerts) > 0:
            alert = alerts[0]
            text.textLine("First alert")
            text.textLine("   M{:.1f}".format(alert["magnitude"]))
            at = dateutil.parser.parse(alert["timestamp"])
            dt = timedelta_to_seconds(numpy.timedelta64(at-ot))
            text.textLine("   {at:%Y-%m-%d %H:%M:%S.%f} ({dt:.1f}s after OT)".format(at=at, dt=dt))
            text.textLine("   Alert at epicenter: {:.1f}s after P wave".format(dt-event["depth_km"]/vp))
            ts = dt - event["depth_km"]/vs
            tsLabel = "after" if ts > 0 else "before"
            text.textLine("                       {:.1f}s {} S wave".format(numpy.abs(ts), tsLabel))
            blindDist = ((dt*vs)**2 - event["depth_km"]**2)**0.5 if dt*vs > event["depth_km"] else 0.0
            text.textLine("   Radius of late alert zone: {:.0f} km".format(blindDist))

        # First alert exceeding threshold
        magThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        thresholdAlert = False
        for alert in alerts:
            if alert["magnitude"] > magThreshold:
                if alert != alerts[0]:
                    thresholdAlert = True
                break
        if thresholdAlert:
            text.textLine("First alert exceeding threshold")
            text.textLine("   M{:.1f}".format(alert["magnitude"]))
            at = dateutil.parser.parse(alert["timestamp"])
            dt = timedelta_to_seconds(numpy.timedelta64(at-ot))
            text.textLine("   {at:%Y-%m-%d %H:%M:%S.%f} ({dt:.1f}s after OT)".format(at=at, dt=dt))
            text.textLine("   Time at epicenter: {:.1f}s after P wave".format(dt-event["depth_km"]/vp))
            ts = dt - event["depth_km"]/vs
            tsLabel = "after" if ts > 0 else "before"
            text.textLine("                      {:.1f}s {} S wave".format(ts, tsLabel))
            blindDist = ((dt*vs)**2 - event["depth_km"]**2)**0.5 if dt > 0 else 0.0
            text.textLine("   Radius of late alert zone: {:.0f} km".format(blindDist))

        self.canvas.drawText(text)
        self.canvas.restoreState()
        return

    def _render_event_perf_table(self, perfEEW, perfTheoryMag, perfTheoryMagBias):
        data = [
            ["Thresholds", "", "Alert Region", "", "Cost Savings", "", "", ""],
            ["Mw", "MMI", "Area\n(km^2)", "Population", "ShakeAlert", "", "Catalog Mag.", "", "Catalog Mag. w/Bias", ""],
            [""]*4 + ["Area", "Pop",]*3,
        ]
        numTableHeaderRows = len(data)
        tableCols = (
            ("magnitude_threshold", "{:3.1f}"),
            ("mmi_threshold", "{:3.1f}"),
            ("area_alert", "{:7.1e}"),
            ("population_alert", "{:7.1e}"),
            ("area_costsavings_eew", "{:5.0f}"),
            ("population_costsavings_eew", "{:8.2e}"),
        )
        
        for pEEW in perfEEW:
            row = [s.format(pEEW[v]) for v,s in tableCols]
            mask = numpy.logical_and(perfTheoryMag["magnitude_threshold"] == pEEW["magnitude_threshold"],
                                         perfTheoryMag["mmi_threshold"] == pEEW["mmi_threshold"])
            if numpy.sum(mask) > 0:
                row += [s.format(perfTheoryMag[mask][0][v]) for v,s in tableCols[-2:]]
            else:
                row += [""]*2
            mask = numpy.logical_and(perfTheoryMagBias["magnitude_threshold"] == pEEW["magnitude_threshold"],
                                         perfTheoryMagBias["mmi_threshold"] == pEEW["mmi_threshold"])
            if numpy.sum(mask) > 0:
                row += [s.format(perfTheoryMagBias[mask][0][v]) for v,s in tableCols[-2:]]
            else:
                row += [""]*2
            data.append(row)

        # Perfect EEW in last row
        areaText = "{:5.0f}".format(perfEEW["area_costsavings_perfecteew"][-1])
        popText = "{:8.2e}".format(perfEEW["population_costsavings_perfecteew"][-1])
        row = ["Perfect EEW"] + [""]*3 + [areaText, popText]*3
        data.append(row)
            
        style = [
            ("FONT", (0,0), (-1,-1), "Courier", 7),
            ("LEFTPADDING", (0,0), (-1,-1), 2),
            ("RIGHTPADDING", (0,0), (-1,-1), 2),
            ("BOTTOMPADDING", (0,0), (-1,-1), 2),
            ("TOPPADDING", (0,0), (-1,-1), 2),
            ("GRID", (0,0), (-1,-1), 0.5, (0.5, 0.5, 0.5)),
            ("SPAN", (0,0), (1,0)), # Thresholds
            ("SPAN", (2,0), (3,0)), # Alert Region
            ("SPAN", (4,0), (9,0)), # Cost Savings
            ("SPAN", (0,1), (0,2)), # Mw
            ("SPAN", (1,1), (1,2)), # MMI
            ("SPAN", (2,1), (2,2)), # Area
            ("SPAN", (3,1), (3,2)), # Population
            ("SPAN", (4,1), (5,1)),
            ("SPAN", (6,1), (7,1)),
            ("SPAN", (8,1), (9,1)),
            ("SPAN", (0,-1), (3,-1)),
            ("LINEABOVE", (0,-1), (-1,-1), 1.0, (0,0,0)),
            ("BACKGROUND", (4,-1), (-1,-1), (0.7, 1.0, 1.0)),
            ("ALIGN", (0,0), (-1,-1), "CENTER"),
        ]
        # Highlight row matching alert thresholds
        magThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        mmiThreshold = self.config.getfloat("alerts", "mmi_threshold")
        maskMMI = numpy.ma.masked_values(perfEEW["mmi_threshold"], mmiThreshold).mask
        maskMag = numpy.ma.masked_values(perfEEW["magnitude_threshold"], magThreshold).mask
        row = numTableHeaderRows + numpy.argmax(numpy.logical_and(maskMMI, maskMag))
        style.append(("BACKGROUND", (0,row),(-1,row), (1.0, 1.0, 0.0)))
                
        # Highlight all cases with maximum Q-area and Q-pop
        if perfEEW["area_costsavings_eew"].shape[-1] > 0:
            maxQ = numpy.max(perfEEW["area_costsavings_eew"])
            rowsMaxQ = numTableHeaderRows + numpy.where(perfEEW["area_costsavings_eew"] >= maxQ)[0]
            for row in rowsMaxQ:
                style.append(("BACKGROUND", (4,row),(4,row), (0.7, 1.0, 0.7)))
        if perfEEW["population_costsavings_eew"].shape[-1] > 0:
            maxQ = numpy.max(perfEEW["population_costsavings_eew"])
            rowsMaxQ = numTableHeaderRows + numpy.where(perfEEW["population_costsavings_eew"] >= maxQ)[0]
            for row in rowsMaxQ:
                style.append(("BACKGROUND", (5,row),(5,row), (0.7, 1.0, 0.7)))
                             
        t = Table(data, style=style)
        w, h = t.wrap(0.5*self.PAGE_WIDTH, self.PAGE_HEIGHT-2.0*self.MARGIN)

        self.canvas.saveState()
        self.canvas.translate(0.5*self.PAGE_WIDTH, self.PAGE_HEIGHT-(self.MARGIN+self.HEADER)-h)
        t.drawOn(self.canvas, *self._coord(0.0, 0.0, inch)) #0.0, y-2.0, inch))

        return

    def _render_summary(self, eqIds):
        """Generate summary for event set.
        """
        
        # Page 1
        self._render_summary_header()

        (x, y1, mapW, mapH) = self._calc_figure_size(nrows=2, ncols=2, row=1, col=1)
        
        mapWidth, mapHeight = self._render_events_map(x, y1, mapW, mapH)

        figW = self.PAGE_WIDTH - 2*self.MARGIN - mapWidth - self.SPACING
        figH = self.PAGE_HEIGHT - 2*self.MARGIN - mapHeight - self.SPACING
        x = self.XLEFT + mapWidth + self.SPACING
        self._render_events_timeline(x, y1, figW, figH)

        (x, y, figW, figH) = self._calc_figure_size(nrows=2, ncols=2, row=2, col=1)
        self._render_cost_functions(x, y, figW, figH)

        (x, y, figW, figH) = self._calc_figure_size(nrows=2, ncols=2, row=2, col=2)
        self._render_summary_thresholds(x, y, figW, figH)
        self.canvas.showPage()

        # Page 2
        self._render_summary_header()
        self._render_summary_performance_figures()
        self.canvas.showPage()
        
        # Page 4+
        self._render_summary_header()
        self._render_summary_perf_table(eqIds)
        self.canvas.showPage()

        return

    def _render_summary_header(self):
        fragilityFn = self.config.get("fragility_curves", "label")
        magThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        mmiThreshold = self.config.getfloat("alerts", "mmi_threshold")
        header = (
            "Performance Summary",
            "Damage/Action: {fn}, Alert thresholds: M{mag:.1f}, MMI {mmi:.1f}".format(fn=fragilityFn, mag=magThreshold, mmi=mmiThreshold),
        )
        self.canvas.saveState()
        self.canvas.translate(self.MARGIN, self.PAGE_HEIGHT-self.MARGIN)
        text = self.canvas.beginText()
        text.setFont("Helvetica-Bold", 12)
        text.textLines("\n".join(header))
        self.canvas.drawText(text)
        self.canvas.restoreState()

        header = (
            "Brad Aagaard",
            "Version " + datetime.datetime.now().strftime("%Y-%m-%d %I:%M%p"),
        )
        self.canvas.saveState()
        self.canvas.translate(self.PAGE_WIDTH-self.MARGIN, self.PAGE_HEIGHT-self.MARGIN)
        self.canvas.setFont("Helvetica", 8)
        self.canvas.drawRightString(0, 0, header[0])
        self.canvas.drawRightString(0, -10, header[1])
        self.canvas.restoreState()
        return

    def _render_events_map(self, x, y, figW, figH):
        """Map of events in earthquake set."""
        plotsDir = self.config.get("files", "plots_dir")

        filename = os.path.join(plotsDir, "eqset_map_events.jpg")
        imageWidth, imageHeight = self._render_image(filename, x, y, height=figH)
        self._figure_label(x, y+imageHeight, "Earthquake Set")
        return (imageWidth, imageHeight)

    def _render_events_timeline(self, x, y, figW, figH):
        """Map of earthquakes versus origin time.
        """
        plotsDir = self.config.get("files", "plots_dir")
        
        filename = os.path.join(plotsDir, "eqset_magnitude_time.png")
        imageWidth, imageHeight = self._render_image(filename, x, y, width=figW)
        return (imageWidth, imageHeight)


    def _render_cost_functions(self, x, y, figW, figH):
        """Plot of cost functions.
        """
        plotsDir = self.config.get("files", "plots_dir")
        label = analysis_label(self.config)
        
        filename = os.path.join(plotsDir, "eqset_{}_cost_functions.png".format(label))
        imageWidth, imageHeight = self._render_image(filename, x, y, width=figW)
        self._figure_label(x, y+imageHeight, "Cost Functions")
        return (imageWidth, imageHeight)


    def _render_summary_thresholds(self, x, y, figW, figH):
        """Optimum thresholds.
        """
        plotsDir = self.config.get("files", "plots_dir")
        label = analysis_label(self.config)

        # Figure with optimal thresholds
        filename = os.path.join(plotsDir, "eqset_{}_optimal_threshold.png".format(label))
        imageWidth, imageHeight = self._render_image(filename, x, y, width=figW)
        
        return

    def _render_summary_performance_figures(self):
        """ComCat versus ShakeAlert magnitude.
        """
        plotsDir = self.config.get("files", "plots_dir")
        label = analysis_label(self.config)

        # Maps of Q-area and Q-pop
        (x, yArea, figW, figH) = self._calc_figure_size(nrows=2, ncols=2, row=1, col=1)
        filename = os.path.join(plotsDir, "eqset_{}-map_area_costsavings_eew.jpg".format(label))
        imageWidth, imageHeight = self._render_image(filename, x, yArea, height=figH)
        
        (x, yPop, figW, figH) = self._calc_figure_size(nrows=2, ncols=2, row=2, col=1)
        filename = os.path.join(plotsDir, "eqset_{}-map_population_costsavings_eew.jpg".format(label))
        imageWidth, imageHeight = self._render_image(filename, x, yPop, height=figH)

        # Plot of Q-area and Q-pop versus time
        x += imageWidth + self.SPACING
        plotW = self.PAGE_WIDTH-self.MARGIN-x
        filename = os.path.join(plotsDir, "eqset_{}_costsavings_time.png".format(label))
        imageWidth, imageHeight = self._render_image(filename, x, yArea, width=plotW)
        
        # Plot of Q-area and Q-pop versus earthquake magnitude
        filename = os.path.join(plotsDir, "eqset_{}_costsavings_magnitude.png".format(label))
        imageWidth, imageHeight = self._render_image(filename, x, yPop, width=plotW)
        return

    def _render_summary_perf_table(self, eqIds):
        """Table of performance for all earthquakes in set.
        """
        gmpe = self.config.get("mmi_predicted", "gmpe")
        server = self.config.get("shakealert.production", "server")
        fragility = self.config.get("fragility_curves", "label")
        magThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        mmiThreshold = self.config.getfloat("alerts", "mmi_threshold")
        alertLatency = self.config.getfloat("alerts", "alert_latency_sec")
        
        db = AnalysisData(self.config.get("files", "analysis_db"))
        perfs = numpy.array([db.performance_stats(eqId, server, gmpe, fragility, alertLatency, magThreshold, mmiThreshold) for eqId in eqIds]).ravel()

        from numpy.lib.recfunctions import append_fields
        area_metric = numpy.where(perfs["area_costsavings_perfecteew"] > 0.0, perfs["area_costsavings_eew"] / perfs["area_costsavings_perfecteew"], numpy.nan)
        population_metric = numpy.where(perfs["population_costsavings_perfecteew"] > 0.0, perfs["population_costsavings_eew"] / perfs["population_costsavings_perfecteew"], numpy.nan)
        perfs = append_fields(perfs, ("area_metric_eew", "population_metric_eew"), (area_metric, population_metric), dtypes=("float32", "float32"), usemask=False)

        alertLatency = 0.0
        perfsTheoryMag = numpy.array([db.performance_stats(eqId, self.CATALOG_MAGNITUDE, gmpe, fragility, alertLatency, magThreshold, mmiThreshold) for eqId in eqIds]).ravel()
        perfsTheoryMagBias = numpy.array([db.performance_stats(eqId, self.CATALOG_MAGNITUDE_BIAS, gmpe, fragility, alertLatency, magThreshold, mmiThreshold) for eqId in eqIds]).ravel()

        theader = [
            [
                "Earthquake\nId",
                "Mw",
                "Origin Time\n(UTC)",
                "Description",
                "Alert Region", "",
                "Cost Savings", "", "", "", "", "", "", "", "", "",
            ],
            [""]*4 + ["Area\n(sq. km)", "Population"] + ["ShakeAlert", "", "", "", "Catalog Mag.", "", "Catalog Mag. w/Bias", "", "Perfect EEW", ""],
            [""]*6 + ["Area", "Q-area", "Pop", "Q-pop"] + ["Area", "Pop"]*3,
        ]

        eventCols = (
            "{event[event_id]:s}",
            "{event[magnitude]:4.2f}",
            "{ot:%Y-%m-%d %H:%M:%S}",
            "{event[description]:s}",
        )
        perfCols = (
            ("area_alert", "{:7.1e}"),
            ("population_alert", "{:7.1e}"),
            ("area_costsavings_eew", "{:5.0f}"),
            ("area_metric_eew", "{:5.2f}"),
            ("population_costsavings_eew", "{:8.2e}"),
            ("population_metric_eew", "{:5.2f}"),
        )
        perfTheoryCols = (
            ("area_costsavings_eew", "{:5.0f}"),
            ("population_costsavings_eew", "{:8.2e}"),
        )
        perfPerfectCols = (
            ("area_costsavings_perfecteew", "{:5.0f}"),
            ("population_costsavings_perfecteew", "{:8.2e}"),
        )

        tdata = []
        for iperf, perf in enumerate(perfs):
            event = db.comcat_event(perf["comcat_id"].decode())
            ot = dateutil.parser.parse(event["origin_time"])
            row = [s.format(ot=ot, event=event) for s in eventCols]

            row += [s.format(perf[v]) for v,s in perfCols]

            if len(perfsTheoryMag) > 0:
                assert(perfsTheoryMag[iperf]["comcat_id"] == perf["comcat_id"])
                row += [s.format(perfsTheoryMag[iperf][v]) for v,s in perfTheoryCols]
            else:
                row += [""]*len(perfTheoryCols)

            if len(perfsTheoryMagBias) > 0:
                assert(perfsTheoryMagBias[iperf]["comcat_id"] == perf["comcat_id"])
                row += [s.format(perfsTheoryMagBias[iperf][v]) for v,s in perfTheoryCols]
            else:
                row += [""]*len(perfTheoryCols)

            row += [s.format(perf[v]) for v,s in perfPerfectCols]
            tdata.append(row)

        tfooter = []
        row = ["Total"] + [""]*5
        for p in [perfs, perfsTheoryMag, perfsTheoryMagBias]:
            row += [
                "{:8.2e}".format(numpy.sum(p["area_costsavings_eew"])),
                "{:8.2e}".format(numpy.sum(p["population_costsavings_eew"])),
            ]
        row.insert(7, "")
        row.insert(9, "")
        row += [
                "{:8.2e}".format(numpy.sum(perfs["area_costsavings_perfecteew"])),
                "{:8.2e}".format(numpy.sum(perfs["population_costsavings_perfecteew"])),
        ]
        tfooter.append(row)

        row = ["Q"] + [""]*5
        for p in [perfs, perfsTheoryMag, perfsTheoryMagBias]:
            areaQ = numpy.sum(p["area_costsavings_eew"]) / numpy.sum(p["area_costsavings_perfecteew"])
            popQ = numpy.sum(p["population_costsavings_eew"]) / numpy.sum(p["population_costsavings_perfecteew"])
            row += [
                "{:4.2f}".format(areaQ),
                "{:4.2f}".format(popQ),
            ]
        row.insert(7, "")
        row.insert(9, "")
        row += [""]*2
        tfooter.append(row)

        style = [
            ("FONT", (0,0), (-1,-1), "Courier", 6.5),
            ("LEFTPADDING", (0,0), (-1,-1), 2),
            ("RIGHTPADDING", (0,0), (-1,-1), 2),
            ("BOTTOMPADDING", (0,0), (-1,-1), 2),
            ("TOPPADDING", (0,0), (-1,-1), 2),
            ("GRID", (0,0), (-1,-1), 0.5, (0.5, 0.5, 0.5)),
            ("SPAN", (4,0), (5,0)), # Alert Region
            ("SPAN", (6,0), (15,0)), # Cost Savings
            ("SPAN", (0,0), (0,2)), # Id
            ("SPAN", (1,0), (1,2)), # Mw
            ("SPAN", (2,0), (2,2)), # Origin time
            ("SPAN", (3,0), (3,2)), # Description
            ("SPAN", (4,1), (4,2)), # Area
            ("SPAN", (5,1), (5,2)), # Population
            ("SPAN", (6,1), (9,1)), # ShakeAlert
            ("SPAN", (10,1), (11,1)), # Catalog Mag.
            ("SPAN", (12,1), (13,1)), # Catalog Mag. + Bias
            ("SPAN", (14,1), (15,1)), # Perfect EEW
            ("SPAN", (0,-2), (5,-2)), # Total
            ("SPAN", (0,-1), (5,-1)), # Q
            ("SPAN", (-2,-1), (-1,-1)), # Q perfect
            ("BACKGROUND", (-2,-1), (-1,-1), (0.9, 0.9, 0.9)),
            ("LINEABOVE", (0,-2), (-1,-2), 1.0, (0,0,0)),
            ("ALIGN", (0,0), (-1,-1), "CENTER"),
            ("ALIGN", (0,-2), (5,-1), "RIGHT"),
            ("BACKGROUND", (6,-1),(6,-1), (0.7, 1.0, 0.7)),
            ("BACKGROUND", (7,-2), (7,-1), (0.9, 0.9, 0.9)),
            ("BACKGROUND", (8,-1),(8,-1), (0.7, 1.0, 0.7)),
            ("BACKGROUND", (9,-2), (9,-1), (0.9, 0.9, 0.9)),
        ]

        table = Table(theader + tdata + tfooter, style=style, repeatRows=len(theader))
        self._render_multipage(table, width=self.PAGE_WIDTH-2*self.MARGIN, height=self.PAGE_HEIGHT-2*self.MARGIN-self.HEADER)
        return

    def _figure_label(self, x, y, label):
        """
        """
        self.canvas.saveState()
        self.canvas.translate(x, y)
        text = self.canvas.beginText()
        text.setFont("Helvetica-Bold", 10)
        text.textLines(label)
        self.canvas.drawText(text)
        self.canvas.restoreState()
        return

    def _render_image(self, filename, x, y, width=None, height=None):
        image = Image.open(filename)
        if width is None and height is None:
            raise ValueError("Missing image 'width' or 'height'.")
            
        w = width if width else height*image.width/image.height
        h = height if height else width*image.height/image.width

        self.canvas.drawImage(ImageReader(image), x, y, width=w, height=h)
        return (w,h)

    def _coord(self, x, y, unit=1):
        return x * unit, y * unit

    def _render_multipage(self, flowable, width, height):
        flowables = flowable.split(width, height)
        flowable = flowables[0]
        w, h = flowable.wrap(width, height)
        flowable.drawOn(self.canvas, self.XLEFT, self.YBOT+(height-h))
        if len(flowables) > 1:
            self.canvas.showPage()
            self._render_summary_header()
            self._render_multipage(flowables[1], width, height)
        return

    def _sort_eqs(self, eqIds, key="origin_time"):
        db = AnalysisData(self.config.get("files", "analysis_db"))
        values = numpy.array([dateutil.parser.parse(db.comcat_event(eqId)[key]) for eqId in eqIds])
        return numpy.array(eqIds)[numpy.argsort(values)]

    def _calc_figure_size(self, nrows=1, ncols=1, row=1, col=1):
        figW = (self.PAGE_WIDTH - 2*self.MARGIN - (ncols-1)*self.SPACING) / ncols
        figH = (self.PAGE_HEIGHT - 2*self.MARGIN - self.HEADER - (nrows-1)*self.SPACING) / nrows
        xLeft = self.XLEFT + (col-1)*(figW + self.SPACING)
        yBot = self.YBOT + (nrows-row)*(figH + self.SPACING)
        return (xLeft, yBot, figW, figH)

# End of file
