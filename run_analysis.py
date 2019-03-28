#!/usr/bin/env python

import shutil
import argparse
import os

import analyzer
import downloader

EQSETS = [
    "eqsets/sanfrancisco.cfg",
    "eqsets/losangeles.cfg",
]

COSTFNS = [
    "fragility_fearavoidance.cfg",
    "fragility_fearavoidance_step.cfg",
    "fragility_fearavoidance_sigmoid.cfg",
    "fragility_fearavoidance_r20.cfg",
    "fragility_injury.cfg",
    "fragility_injury_w2.cfg",
]

class App(object):
    """Application for running EEW cost/benefit analysis for various cost functions.
    """

    def __init__(self):
        """Constructor.
        """
        self.show_progress = False
        self.nthreads = 0
        return
    
    def main(self):
        """Main entry point.
        """
        args = self._parse_command_line()
        self.nthreads = args.nthreads
        self.color_style = args.color_style
        if args.show_progress:
            self.show_progress = True

        if self.show_progress:
            print("Removing analyzer log file...")
            shutil.rmtree("analyzer.log", ignore_errors=True)

        if args.init_perf_db:
            if self.show_progress:
                print("Removing performance results from analysis database...")
            app = downloader.DownloaderApp()
            app.main(**{
                "config": None,
                "fetch_eewalerts": None,
                "fetch_events": None,
                "all": False,
                "db_init": "performance",
                "fetch_shakemaps": None,
                "db_summary": None,
                "db_populate": None,
                "db_replace_rows": False,
                "show_matches": False,
                "show_parameters": True,
                "show_progress": True,
                "debug": True,
            })

        costfns = COSTFNS if args.all else args.costfns.split(",")
        for costfn in costfns:
            self._run_costfn(costfn)
            
        return


    def _run_costfn(self, costfn):
        defaults = {
            "config": None,
            "show_parameters": False,
            "process_events": False,
            "optimize_events": False,
            "plot_alert_maps": None,
            "plot_event_maps": None,
            "plot_event_figures": None,
            "plot_summary_maps": None,
            "plot_summary_figures": None,
            "generate_report": None,
            "nthreads": self.nthreads,
            "all": False,
            "debug": False,
            "show_progress": True,
            "color_style": self.color_style,
        }
            
        # ShakeAlert
        args = defaults.copy()
        args.update({
            "config": ",".join(EQSETS + [costfn]),
            "nthreads": self.nthreads,
            "optimize_events": True,
            "plot_event_maps": "all",
            "plot_event_figures": "all",
        })
        self._run_analyzer(**args)
        
        # Catalog magnitude
        args = defaults.copy()
        args.update({
            "config": ",".join(EQSETS + [costfn, "catalog_magnitude.cfg"]),
            "nthreads": self.nthreads,
            "optimize_events": True,
            "plot_event_maps": "alert",
        })
        self._run_analyzer(**args)

        # Catalog magitude + bias
        args.update({"config": ",".join(EQSETS + [costfn, "catalog_magnitude_bias.cfg"])})
        self._run_analyzer(**args)

        # Summary
        args = defaults.copy()
        args.update({
            "config": ",".join(EQSETS + [costfn]),
            "plot_summary_figures": "all",
            "plot_summary_maps": "all",
            "generate_report": "full",
        })
        self._run_analyzer(**args)
        if os.path.isfile("report.pdf"):
            shutil.move("report.pdf", "report_SF+LA_{}.pdf")
            
        return

    def _run_analyzer(self, **kwargs):
        app = analyzer.EEWAnalyzeApp()
        if self.show_progress:
            print("Running analyzer app w/config: {}".format(kwargs))
        #app.main(**kwargs)
        return
    
    def _parse_command_line(cls):
        """Parse command line arguments.
        """
        parser = argparse.ArgumentParser()
        parser.add_argument("--costfns", action="store", dest="costfns", choices=COSTFNS+[None])
        parser.add_argument("--all", action="store_true", dest="all")

        parser.add_argument("--color-style", action="store", dest="color_style", default="lightbg")
        parser.add_argument("--num-threads", action="store", type=int, dest="nthreads", default=16)
        parser.add_argument("--init-perf-db", action="store_true", dest="init_perf_db", default=False)

        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug", default=True)
        return parser.parse_args()

        
    
# ======================================================================
if __name__ == "__main__":
    App().main()


    
