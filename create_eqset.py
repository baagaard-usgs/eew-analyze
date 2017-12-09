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
[sanfrancisco]
longitude = -122.419
latitude = 37.775
max_radius_deg = 2.0

[losangeles]
longitude = -117.396
latitude = 33.953
max_radius_deg = 2.5

[fdsn.event]
starttime = 2012-01-27
minmagnitude = 4.0
producttype = shakemap
orderby = magnitude

[fdsn.client]
debug = False

[files]
filename_template = ./eqsets/[DOMAIN].cfg
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


    if makeDir or not os.path.isdir(eventDir):
        os.makedirs(eventDir)

    filename = params.get("files", option) if args is None else params.get("files", option) % args
    return os.path.join(eventDir, filename)
    

# ----------------------------------------------------------------------
class CreateEqSetApp(object):
    """
    Download list of events matching criteria and create earthquake set.
    """
    
    def __init__(self):
        """Constructor.
        """
        self.params = None
        return

    def main(self):
        """Main entry point
        """
        # Initialization
        args = self._parseCommandLine()
        logLevel = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=logLevel, filename="create_eqset.log")
        if args.show_progress:
            self.showProgress = True
        self.initialize(args.config)

        # Show parameters
        if args.show_parameters or args.all:
            self.show_parameters()

        if args.create or args.all:
            self.create(args.domain)

        return

    def create(self, domain):
        """Fetch events from ComCat catalog and create earthquake set .cfg file.

        :type domain: str
        :param domain: Name of domain in parameters to use for search parameters.
        """
        
        if self.showProgress:
            print("Fetching events for earthquake set '%s'..." % domain)

        from obspy.clients.fdsn import Client
        client = Client("USGS", debug=self.params.getboolean("fdsn.client", "debug"))
        
        kwargs = {
            "latitude": self.params.get(domain, "latitude"),
            "longitude": self.params.get(domain, "longitude"),
            "maxradius": self.params.get(domain, "max_radius_deg"),
        }
        kwargs.update(dict(self.params.items("fdsn.event")))
        catalog = client.get_events(**kwargs)

        filename = self.params.get("files", "filename_template").replace("[DOMAIN]", domain)
        path = os.path.split(filename)[0]
        if not os.path.isdir(path):
            os.makedirs(path)
        with open(filename, "w") as fout:
            fout.write("[events]\n")
            for event in catalog:
                eqid = event.preferred_origin_id.id.split("/")[4]
                if "event_descriptions" in event:
                    description = event["event_descriptions"][0].text
                else:
                    description = ""
                info = {
                    "eqid": eqid,
                    "mag": event.preferred_magnitude().mag,
                    "date": str(event.preferred_origin().time.date),
                    "description": description,
                }
                fout.write("%(eqid)s = M%(mag).1f %(description)s, %(date)s\n" % info)

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
        if config_filenames:
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

    def _parseCommandLine(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--config", action="store", dest="config")
        parser.add_argument("--domain", action="store", dest="domain", default="sanfrancisco")
        parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
        parser.add_argument("--create", action="store_true", dest="create")
        parser.add_argument("--all", action="store_true", dest="all")
        parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true", dest="debug")
        return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    CreateEqSetApp().main()


# End of file
