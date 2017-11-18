#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import logging

DEFAULTS = """
[events]
# Example: nc72923380 = Mw 4.6 Paicines, 2017-11-13

[warning_time]
module = WarningTimeSWave

[gmpe]
module = GMPEUserDisplay

[map]

[shakealert]
server = eew-bk-prod1

[files]
event = event.geojson
shakemap = shakemap_shapefiles.zip
"""

# ----------------------------------------------------------------------


def _config_get_list(list_string):
    """Convert list as string to list.

    :type list_string: list
    :param list_string: List as string.
    :returns: List of strings.
    """
    l = [f.strip() for f in list_string[1:-1].split(",")]
    return l


def _data_filename(params, option, eqId, args=None, makeDir=False):
    """Construct relative path for file.

    :type params: ConfigParser
    :param params: Application parameters.
    :type option: str
    :param option: Option in application parameters for filename.
    :type eqId: str
    :param ComCat event id (e.g., nc72923380).
    :param params: Application parameters.
    :type args: tuple
    :param args: Tuple for arguments for substitution in name of parameter file.
    :type makeDir: bool
    :param makeDir: Create directory for file if True.
    """
    eventDir = os.path.join("data", eqId)
    if makeDir or not os.path.isdir(eventDir):
        os.makedirs(eventDir)

    filename = params.get("files", option) if args is None else params.get("files", option) % args
    return os.path.join(eventDir, filename)
    

# ----------------------------------------------------------------------
class EEWAnalyzeApp(object):
    """
    Analyze ShakeAlert performance using ShakeMap.
    """
    
    def __init__(self):
        """Constructor.
        """
        self.params = None
        self.event = None
        self.shakemap = None
        self.alert = None
        self.warningTime = None
        self.gmpe = None
        return

    def main(self):
        """Main entry point
        """
        # Initialization
        args = self._parseCommandLine()
        logLevel = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=logLevel, filename="analyze_events.log")
        if args.show_progress:
            self.showProgress = True
        self.initialize(args.config)

        # Show parameters
        if args.show_parameters or args.all:
            self.show_parameters()

        if args.fetch or args.process_data or args.plot or args.generate_report or args.all:
            for eqId in self.params.options("events"):
                self.process_event(eqId, args)

        return

    def process_event(self, eqId, args):
        """For given event, fetch data, process data, generate plots, and generate report.
        
        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).
        :type args: argparse.ArgumentParser
        :param args: Command line arguments.
        """
        if self.showProgress:
            print("Processing event %s..." % eqId)
            
        # Fetch data
        if args.fetch == "event" or args.fetch == "all" or args.all:
            self.fetch_event(eqId)
        if args.fetch == "shakemap" or args.fetch == "all" or args.all:
            self._loadEvent(eqId)
            self.fetch_shakemap()
        if args.fetch == "eewalert" or args.fetch == "all" or args.all:
            self._loadEvent(eqId)
            self.fetch_eewalert()

        # Process data
        if args.process_data or args.all:
            self._loadData(eqId)
            self.process_data()

        # Plot data
        if args.plot == "map" or args.plot == "all" or args.all:
            self._loadData(eqId)
            self.plot_map()
        if args.plot == "histograms" or args.plot == "all" or args.all:
            self._loadData(eqId)
            self.plot_histograms()
        if args.plot == "mmierror" or args.plot == "all" or args.all:
            self._loadData(eqId)
            self.plot_mmierror()
            
        # Generate report
        if args.generate_report or args.all:
            self.generate_report(eqId)

        return
    
    def initialize(self, config_filenames):
        """Set parameters from config file and DEFAULTS.

        :type config_filename: str
        :param config_filename: Name of configuration (INI) file with parameters.
        """
        import ConfigParser
        import io
        config = ConfigParser.SafeConfigParser()
        config.readfp(io.BytesIO(DEFAULTS))
        for filename in config_filenames.split(","):
            if self.showProgress:
                print("Fetching parameters from %s..." % filename)
            config.read(filename)

        self.params = config

        return
    
    def show_parameters(self):
        """Write parameters to stdout.
        """
        import sys
        self.params.write(sys.stdout)
        return

    def fetch_event(self, eqId):
        """Fetch geojson event file from USGS ComCat using web services.
        
        :type eqId: string
        :param eqId: ComCat event id (e.g., nc72923380).
        """

        if self.showProgress:
            print("Fetching event file...")

        from comcat import DetailEvent
        self.event = DetailEvent()
        self.event.fetch(eqId, _data_filename(self.params, "event", eqId))
        return
    
    def fetch_shakemap(self):
        """Fetch geojson event file from USGS ComCat using web services.
        """
        if self.showProgress:
            print("Fetching ShakeMap data...")

        shakemap = self.event.get_product("shakemap")
        if len(shakemap) > 1:
            raise ValueError("Expected to get preferred ShakeMap.")
        filename = _data_filename(self.params, "shakemap", self.event.id)
        shakemap[0].fetch("shape.zip", filename)
        return
        
    def fetch_dmlog(self):
        """Fetch DM log from EEW production system.
        """
        if self.showProgress:
            print("Fetching DM log from %s..." % server)
        raise NotImplementedError(":TODO: @brad")
        return
    
    def process_data(self):
        """Analyze ShakeMap and DM log to assess ShakeAlert performance.
        """
        if self.showProgress:
            print("Processing data...")
        self._loadEvent()
        raise NotImplementedError(":TODO: @brad")
        return
    
    
    def plot_map(self):
        """Plot map of MMI with ShakeAlert information.
        """
        if self.showProgress:
            print("Plotting maps...")
        self._loadEvent()
        raise NotImplementedError(":TODO: @brad")
        return
    
    def plot_histograms(self):
        """Plot 3-D and 2-D histograms of area/population versus warning time and MMI.
        """
        if self.showProgress:
            print("Plotting histograms...")
        self._loadEvent()
        raise NotImplementedError(":TODO: @brad")
        return
    
    def plot_mmierror(self):
        """Plot estimated versus observed MMI. Compute median and standard deviation.
        """
        if self.showProgress:
            print("Plotting MMI error...")
        self._loadEvent()
        raise NotImplementedError(":TODO: @brad")
        return
    
    def generate_report(self):
        """Assemble plots, etc into PDF file.
        """
        if self.showProgress:
            print("Generating report...")
        raise NotImplementedError(":TODO: @brad")
        return
        
    def _parseCommandLine(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config", required=True)
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--fetch", action="store", dest="fetch", default=None, choices=[None, "all", "event", "shakemap", "eewalert"])
        parser.add_argument("--process-data", action="store_true", dest="process_data")
        parser.add_argument("--plot", action="store", dest="plot", default=None, choices=[None, "all", "map", "histograms", "mmierror"])
        parser.add_argument("--generate-report", action="store_true", dest="generate_report")
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug")
        return parser.parse_args()

    def _loadData(self, event):
        """Load data for analysis.
        """
        self._loadEvent(event)
        self._loadShakeMap(event)
        self._loadAlert(event)
        return
    
    def _loadEvent(self, eqId):
        """Load GeoJSON event file if not already loaded.

        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).

        """
        if not self.event:
            filename = _data_filename(self.params, "event", eqId)
            if not os.path.isfile(filename):
                raise ValueError("Could not load event file '%s'. Did you forget to --fetch-event?" % filename)
            from comcat import DetailEvent
            self.event = DetailEvent(filename)
        return
        

    def _loadShakeMap(self, eqId):
        """Load event ShakeMap if not already loaded.

        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).
        """
        if not self.shakemap:
            raise NotImplementedError(":TODO: @brad")
        return
        

    def _loadAlert(self, eqId):
        """Load EEW alert if not already loaded.
        
        :type event: str
        :param eqId: ComCat Earthquake id (e.g., nc72923380).
        """
        if not self.alert:
            raise NotImplementedError(":TODO: @brad")
        return
        

# ======================================================================
if __name__ == "__main__":
    EEWAnalyzeApp().main()


# End of file
