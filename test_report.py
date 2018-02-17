import os

from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.platypus import Table

HEADER_HEIGHT = 0.5*inch
MARGIN = 0.5*inch

page_width, page_height = landscape(letter)

def drawHeader(canvas):
    header = (
        "YYYY-MM-DD HH:MM:SS ncXXXXXXXX MX.X LOCATION",
        "Damage/Action: PublicAnxiety, Alert thresholds: MX.X, MMI 2.0",
    )
    canvas.saveState()
    canvas.translate(MARGIN, page_height-MARGIN)
    text = canvas.beginText()
    text.setFont("Helvetica-Bold", 12)
    text.textLines("\n".join(header))
    canvas.drawText(text)
    canvas.restoreState()

    header = (
        "Brad Aagaard",
        "Version 2018-02-15 01:02pm"
    )
    canvas.saveState()
    canvas.translate(page_width-MARGIN, page_height-MARGIN)
    canvas.setFont("Helvetica", 8)
    canvas.drawRightString(0, 0, header[0])
    canvas.drawRightString(0, -10, header[1])
    canvas.restoreState()
    return

def coord(x, y, unit=1):
    x, y = x * unit, y * unit
    return x, y


def go():

    c = canvas.Canvas(
        "report.pdf",
        pagesize=landscape(letter),
        bottomup=1,
        )

    drawHeader(c)

    plots_dir = "data/plots"
    filename_prefix = "nc72948801-eew-bk-prod1-ASK2014-PublicAnxiety-M3.0-MMI0.0"
    map_size = 3.25*inch

    spacing = 0.125*inch

    # Maps of MMI observed, predicted, and residual (observed-predicted)
    x = MARGIN
    y = page_height-(MARGIN+HEADER_HEIGHT)-map_size
    
    filename = os.path.join(plots_dir, filename_prefix+"-map_mmi_obs.jpg")
    c.drawImage(filename, x, y, width=map_size, height=map_size)
    
    x += map_size + spacing
    filename = os.path.join(plots_dir, filename_prefix+"-map_mmi_pred.jpg")
    c.drawImage(filename, x, y, width=map_size, height=map_size)

    x += map_size + spacing
    filename = os.path.join(plots_dir, filename_prefix+"-map_mmi_residual.jpg")
    c.drawImage(filename, x, y, width=map_size, height=map_size)

    # Map of alert categories w/warning time
    x = MARGIN
    y -= map_size + spacing
    filename = os.path.join(plots_dir, filename_prefix+"-map_alert_category.jpg")
    c.drawImage(filename, x, y, width=map_size, height=map_size)

    # Figure of magnitude and location error
    x += map_size + spacing
    filename = os.path.join(plots_dir, filename_prefix+"-alert_error.png")
    #c.drawImage(filename, x, y, width=map_size, )    
    
    # Stats
    c.saveState()
    c.translate(7.25*inch, 3.5*inch)

    text = c.beginText()
    text.setFont("Courier", 8)
    text.textLine("ANSS")
    text.textLine("   Mw4.4")
    text.textLine("   YYYY-MM-DD HH:MM:SS.MSEC")
    text.textLine("First alert")
    text.textLine("   MwX.X")
    text.textLine("   YYYY-MM-DD HH:MM:SS.MSEC (X.Xs after OT)")
    c.drawText(text)

    data = [
        ["MMI\nThreshold", "Alert", "", "Q", ""],
        ["", "Area (km^2)", "Population", "Area", "Population"],
        ["2.0", "1.0e+7", "100k", "0.4", "0.5"],
        ["2.5", "1.0e+7", "100k", "0.4", "0.5"],
        ["3.0", "1.0e+7", "100k", "0.4", "0.5"],
    ]
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
        ("BACKGROUND", (0,3),(-1,3), (0.7, 1.0, 0.7)),
    ]
    t = Table(data, style=style)
    t.wrapOn(c, 3.0*inch, 2.0*inch)
    t.drawOn(c, *coord(0.0, -2.0, inch)) #0.0, y-2.0, inch))
    
    c.restoreState()
    
    c.showPage()
    c.save()
    
if __name__ == "__main__":
    go()
