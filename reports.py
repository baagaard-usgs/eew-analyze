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

from analysisdb import AnalysisData
from analysis_utils import analysis_label

class AnalysisSummary(object):
    """Summary of analysis.
    """
    HEADER = 0.5*inch
    MARGIN = 0.5*inch
    PAGE_WIDTH, PAGE_HEIGHT = landscape(letter)
    
    def __init__(self, config):
        """Constructor.
        """
        self.config = config
        self.canvas = None
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

        # Title page with damage/action figure

        for eqId in eqIds:
            self._render_event(eqId)

        self.canvas.save()
        return

    def _render_event(self, eqId):
        """Generate summary for event.
        """
        db = AnalysisData(self.config.get("files", "analysis_db"))
        event = db.comcat_event(eqId)
        
        self._render_event_header(event)
        self._render_event_plots(event)
        self._render_event_stats(event)

        self.canvas.showPage()
        return

    def _render_event_header(self, event):
        ot = dateutil.parser.parse(event["origin_time"])
        fragilityFn = self.config.get("fragility_curves", "object").split(".")[-1]
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

    def _render_event_plots(self, event):
        MAP_SIZE = 3.25*inch
        SPACING = 0.125*inch
        
        # Maps of MMI observed, predicted, and residual (observed-predicted)
        x = self.MARGIN
        y = self.PAGE_HEIGHT-(self.MARGIN+self.HEADER)-MAP_SIZE

        plots_dir = self.config.get("files", "plots_dir")
        label = analysis_label(self.config, event["event_id"])

        filename = os.path.join(plots_dir, label+"-map_mmi_obs.jpg")
        self._render_image(filename, x, y, width=MAP_SIZE)
    
        x += MAP_SIZE + SPACING
        filename = os.path.join(plots_dir, label+"-map_mmi_pred.jpg")
        self._render_image(filename, x, y, width=MAP_SIZE)

        x += MAP_SIZE + SPACING
        filename = os.path.join(plots_dir, label+"-map_mmi_residual.jpg")
        self._render_image(filename, x, y, width=MAP_SIZE)

        # Map of alert categories w/warning time
        x = self.MARGIN
        y -= MAP_SIZE + SPACING
        filename = os.path.join(plots_dir, label+"-map_alert_category.jpg")
        self._render_image(filename, x, y, width=MAP_SIZE)

        # Figure of magnitude and location error
        x += MAP_SIZE + SPACING
        filename = os.path.join(plots_dir, label+"-alert_error.png")
        self._render_image(filename, x, y, width=MAP_SIZE)

        return

    def _render_event_stats(self, event):
        self.canvas.saveState()
        self.canvas.translate(7.25*inch, 3.5*inch)

        ot = dateutil.parser.parse(event["origin_time"])
        text = self.canvas.beginText()
        text.setFont("Courier", 8)
        text.textLine("ANSS")
        text.textLine("   {event[magnitude_type]}{event[magnitude]:.1f}".format(event=event))
        text.textLine("   {ot:%Y-%m-%d %H:%M:%S:%f}".format(ot=ot))

        text.textLine("First alert")
        text.textLine("   MwX.X")
        text.textLine("   YYYY-MM-DD HH:MM:SS.MSEC (X.Xs after OT)")
        self.canvas.drawText(text)

        self._render_event_stats_table(event)

        self.canvas.restoreState()
        return

    def _render_event_stats_table(self, event):
        cache_dir = self.config.get("files", "analysis_cache_dir")
        label = analysis_label(self.config, event["event_id"])
        filename = os.path.join(cache_dir, "threshold_optimization_"+label+".txt")
        cols = [
            ("mmi_threshold", "float32"),
            ("area_damage", "float32"),
            ("area_alert", "float32"),
            ("area_metric", "float32"),
            ("population_damage", "float32"),
            ("population_alert", "float32"),
            ("population_metric", "float32"),
        ]
        results = numpy.loadtxt(filename, dtype=cols)
        
        data = [
            ["MMI\nThreshold", "Alert", "", "Q", ""],
            ["", "Area (km^2)", "Population", "Area", "Population"],
        ]
        tableCols = (
            ("mmi_threshold", "{:3.1f}"),
            ("area_alert", "{:7.1e}"),
            ("population_alert", "{:7.1e}"),
            ("area_metric", "{:5.2f}"),
            ("population_metric", "{:5.2f}"),
        )
        rowPreferred = 2 + numpy.argmax(results["population_metric"])
        for result in results:
            data.append([s.format(result[v]) for v,s in tableCols])
            
        style = [
            ("FONT", (0,0), (-1,-1), "Courier", 8),
            ("LEFTPADDING", (0,0), (-1,-1), 2),
            ("RIGHTPADDING", (0,0), (-1,-1), 2),
            ("BOTTOMPADDING", (0,0), (-1,-1), 2),
            ("TOPPADDING", (0,0), (-1,-1), 2),
            ("GRID", (0,0), (-1,-1), 0.5, (0.5, 0.5, 0.5)),
            ("SPAN", (0,0), (0,1)),
            ("SPAN", (1,0), (2,0)),
            ("SPAN", (3,0), (4,0)),
            ("ALIGN", (0,0), (-1,-1), "CENTER"),
            ("BACKGROUND", (0,rowPreferred),(-1,rowPreferred), (0.7, 1.0, 0.7)),
        ]
        t = Table(data, style=style)
        t.wrapOn(self.canvas, 3.0*inch, 2.0*inch)
        t.drawOn(self.canvas, *self._coord(0.0, -2.0, inch)) #0.0, y-2.0, inch))
        return

    def _render_image(self, filename, x, y, width=None, height=None):
        image = Image.open(filename)
        if width is None and height is None:
            raise ValueError("Missing image 'width' or 'height'.")
            
        w = width if width else height*image.width/image.height
        h = height if height else width*image.height/image.width

        self.canvas.drawImage(ImageReader(image), x, y, width=w, height=h)
        return
    
    def _coord(self, x, y, unit=1):
        return x * unit, y * unit

        
# End of file
