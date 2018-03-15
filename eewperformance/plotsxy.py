# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import numpy

from basemap.Figure import Figure

import analysis_utils
import greatcircle


class Figures(object):
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

        t = (alertsTime - originTime).astype("timedelta64[us]").astype("float32")/1.0e+6
        
        figure = Figure()
        figure.open(6.0, 5.0, margins=((0.45, 0.8, 0.2), (0.4, 0.8, 0.3)))
        nrows = 2
        ncols = 2
        irow = 1
        icol = 1

        tmax = max(30.0, t[0]+30.0)
        
        # Magnitude
        ax = figure.axes(nrows, ncols, irow, icol)
        ax.plot(t, alertsMag, marker="o", mfc="c_ltred", mec="c_fg", lw=0, alpha=0.5)
        ax.set_xlabel("Time after origin time (s)")
        ax.set_ylabel("Magnitude (Mw)")
        ax.set_title("Magnitude")
        ax.set_xlim(0, tmax)

        ax.axhline(self.event["magnitude"], linestyle="--", linewidth=1.0, color="c_blue")
        ax.text(ax.get_xlim()[1], self.event["magnitude"], "ANSS", ha="right", va="bottom", color="c_blue")
        icol = 1
        irow += 1
        
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
        pyplot.close(figure)
        return

# End of file
