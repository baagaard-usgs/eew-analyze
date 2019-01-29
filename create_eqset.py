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


DEFAULTS = u"""
[sanfrancisco]
longitude = -122.419
latitude = 37.775
max_radius_deg = 2.0

[losangeles]
longitude = -117.396
latitude = 33.953
max_radius_deg = 2.5

[eureka]
longitude = -124.4
latitude = 40.4
max_radius_deg = 1.5

[fdsn.event]
starttime = 2012-01-27
minmagnitude = 3.95
producttype = shakemap
orderby = magnitude

[fdsn.client]
debug = False

[blacklist]
ci11129826 = deleted event
# Ignore events near Nevada
nn00642964 = M4.5 15km WNW of Sandy Valley, Nevada, 2018-07-05
# Ignore events in Mexico
ci11129914 = M4.7 11km W of Alberto Oviedo Mota, B.C., MX, 2012-07-01
ci37373442 = M4.4 1km SE of Delta, B.C., MX, 2018-09-29
ci37359304 = M4.3 3km NNW of Delta, B.C., MX, 2015-04-08
ci15115905 = M4.2 11km WSW of Alberto Oviedo Mota, B.C., MX, 2012-02-29
ci38230144 = M4.2 32km ENE of Ensenada, B.C., MX, 2018-07-25
ci11286130 = M4.1 10km W of Alberto Oviedo Mota, B.C., MX, 2013-04-17
ci11364370 = M4.1 7km NNW of Delta, B.C., MX, 2013-09-14
ci37359312 = M4.0 8km NW of Delta, B.C., MX, 2015-04-08
ci38232616 = M4.0 71km ENE of Maneadero, B.C., MX, 2018-07-29
ci37265488 = M4.0 56km WSW of Rosarito, B.C., MX, 2014-09-07
ci11248258 = M4.0 66km WSW of Rosarito, B.C., MX, 2013-02-23

[files]
filename_template = ./eqsets/[DOMAIN].cfg
"""

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
            print("Fetching events for earthquake set '{}'...".format(domain))

        from obspy.clients.fdsn import Client
        client = Client("USGS", debug=self.params.getboolean("fdsn.client", "debug"))
        
        kwargs = {
            "latitude": self.params.get(domain, "latitude"),
            "longitude": self.params.get(domain, "longitude"),
            "maxradius": self.params.get(domain, "max_radius_deg"),
        }
        kwargs.update(dict(self.params.items("fdsn.event")))
        catalog = client.get_events(**kwargs)


        blacklist = self.params.options("blacklist")
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
                if eqid in blacklist:
                    fout.write("#") # Blacklisted events will be commented out.
                fout.write("{info[eqid]} = M{info[mag]:.1f} {info[description]}, {info[date]}\n".format(info=info))
                    

        return

    def initialize(self, config_filenames):
        """Set parameters from config file and DEFAULTS.

        :type config_filename: str
        :param config_filename: Name of configuration (INI) file with parameters.
        """
        import configparser
        import io
        config = configparser.SafeConfigParser()
        config.readfp(io.StringIO(DEFAULTS))
        if config_filenames:
            for filename in config_filenames.split(","):
                if self.showProgress:
                    print("Fetching parameters from {}...".format(filename))
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
