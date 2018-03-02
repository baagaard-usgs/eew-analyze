# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import sqlite3
import sys
import logging
import datetime
import dateutil.parser
import pytz

import numpy

TABLES = [
    ("eew_alerts", [
        "server TEXT NOT NULL",
        "event_id INTEGER NOT NULL",
        "category TEXT NOT NULL",
        "message_type TEXT NOT NULL",
        "timestamp TEXT NOT NULL",
        "version INTEGER NOT NULL",
        "magnitude REAL NOT NULL",
        "magnitude_type TEXT DEFAULT Mw",
        "latitude REAL NOT NULL",
        "longitude REAL NOT NULL",
        "depth_km REAL NOT NULL",
        "origin_time TEXT NOT NULL",
        "num_stations INTEGER DEFAULT 0",
        "UNIQUE(server, event_id, version, category) ON CONFLICT FAIL",
    ]),
    ("comcat_events", [
        "event_id TEXT NOT NULL PRIMARY KEY",
        "latitude REAL NOT NULL",
        "longitude REAL NOT NULL",
        "depth_km REAL NOT NULL",
        "origin_time TEXT NOT NULL",
        "magnitude REAL NOT NULL",
        "magnitude_type TEXT DEFAULT Mw",
        "description TEXT",
        "UNIQUE(event_id) ON CONFLICT FAIL",
    ]),
    ("comcat_shakemaps", [
        "event_id TEXT NOT NULL PRIMARY KEY",
        "mmi_bias REAL",
        "mmi_max REAL",
        "pga_bias REAL",
        "pga_max REAL",
        "pgv_bias REAL",
        "pgv_max REAL",
        "psa03_bias REAL",
        "psa03_max REAL",
        "psa10_bias REAL",
        "psa10_max REAL",
        "psa30_bias REAL",
        "psa30_max REAL",
        "gmpe TEXT",
        "pgm2mi TEXT",
        "software_version TEXT",
        "UNIQUE(event_id) ON CONFLICT FAIL",
    ]),
    ("performance", [
        "comcat_id TEXT NOT NULL",
        "eew_server TEXT NOT NULL",
        "dm_id INTEGER NOT NULL",
        "dm_timestamp TEXT NOT NULL",
        "gmpe TEXT NOT NULL",
        "fragility TEXT NOT NULL",
        "magnitude_threshold REAL NOT NULL",
        "mmi_threshold REAL NOT NULL",
        "area_damage REAL NOT NULL",
        "area_alert REAL NOT NULL",
        "area_cost_eew REAL NOT NULL",
        "area_cost_noeew REAL NOT NULL",
        "area_cost_perfecteew REAL NOT NULL",
        "area_metric REAL NOT NULL",
        "population_damage REAL NOT NULL",
        "population_alert REAL NOT NULL",
        "population_cost_eew REAL NOT NULL",
        "population_cost_noeew REAL NOT NULL",
        "population_cost_perfecteew REAL NOT NULL",
        "population_metric REAL NOT NULL",
        "UNIQUE(comcat_id, eew_server, dm_id, gmpe, fragility, magnitude_threshold, mmi_threshold) ON CONFLICT FAIL",
    ]),
]


class AnalysisData(object):
    """SQLite database with DM alerts and ComCat events.
    """

    def __init__(self, filename):
        """Constructor with filename.

        :type filename: str
        :param filename: Filename of SQLite database

        """
        self.connection = sqlite3.connect(filename)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()
        return

    def init(self, key):
        """Create database.
        """
        if key == "all":
            for name, columns in TABLES[::-1]:
                self.cursor.execute("DROP TABLE IF EXISTS {}".format(name))
            for name,columns in TABLES:
                self.cursor.execute("CREATE TABLE {name} ({fields})".format(name=name, fields=", ".join(columns)))
        else:
            for name, columns in TABLES:
                if name == key:
                    self.cursor.execute("DROP TABLE IF EXISTS {}".format(name))
                    self.cursor.execute("CREATE TABLE {name} ({fields})".format(name=name, fields=", ".join(columns)))
        self.connection.commit()
        return

    def add_alerts(self, alerts, replace=False):
        """Add alert info to database.
        """
        COLUMNS = (
            "server",
            "event_id",
            "category",
            "message_type",
            "timestamp",
            "version",
            "magnitude",
            "magnitude_type",
            "latitude",
            "longitude",
            "depth_km",
            "origin_time",
            "num_stations",
            )
        insertCols = ", ".join(COLUMNS)
        valueCols = ", ".join([":{}".format(col) for col in COLUMNS])
        cmd = "INSERT"
        if replace:
            cmd += " OR REPLACE"
        try:
            #self.cursor.executemany("INSERT INTO eew_alerts({}) VALUES({})".format(insertCols, valueCols), alerts)
            for alert in alerts:
                self.cursor.execute("{} INTO eew_alerts({}) VALUES({})".format(cmd, insertCols, valueCols), alert)
            self.connection.commit()
        except sqlite3.IntegrityError as ex:
            logging.getLogger(__name__).debug(str(ex))
            #logging.getLogger(__name__).debug(str(alerts))
            logging.getLogger(__name__).debug(str(alert))
        return
    
    def add_event(self, event, replace=False):
        """Add ComCat event to database.

        :type event: Detail event
        :param event: ComCat event to add to database.
        """
        COLUMNS = (
            "event_id",
            "latitude",
            "longitude",
            "depth_km",
            "origin_time",
            "magnitude",
            "magnitude_type",
            "description",
            )

        originTime = event.time if event.time.tzinfo else event.time.replace(tzinfo=pytz.UTC)
        insertCols = ", ".join(COLUMNS)
        eventValues = (event.id, event.latitude, event.longitude, event.depth, originTime, event.magnitude, event["magType"], event.location)
        valueCols = ",".join("?"*len(eventValues))
        cmd = "INSERT"
        if replace:
            cmd += " OR REPLACE"
        try:
            self.cursor.execute("{0} INTO comcat_events({1}) VALUES({2})".format(cmd, insertCols, valueCols), eventValues)
            self.connection.commit()
        except sqlite3.IntegrityError as ex:
            logging.getLogger(__name__).debug(str(ex))
            logging.getLogger(__name__).debug(str(event))
        return

    def add_shakemap_info(self, info, replace=False):
        """Add ComCat ShakeMap info to database.

        :type info: dict
        :param info: ComCat ShakeMap info (from info.json or info.xml).
        """
        COLUMNS = (
            "event_id",
            "mmi_bias",
            "mmi_max",
            "pga_bias",
            "pga_max",
            "pgv_bias",
            "pgv_max",
            "psa03_bias",
            "psa03_max",
            "psa10_bias",
            "psa10_max",
            "psa30_bias",
            "psa30_max",
            "gmpe",
            "pgm2mi",
            "software_version",
            )

        insertCols = ", ".join(COLUMNS)

        gmmod = info["processing"]["ground_motion_modules"]
        infoDict = {
            "event_id": info["event_id"],
            "mmi_bias": 0.0,
            "mmi_max": 0.0,
            "pga_bias": 0.0,
            "pga_max": 0.0,
            "pgv_bias": 0.0,
            "pgv_max": 0.0,
            "psa03_bias": 0.0,
            "psa03_max": 0.0,
            "psa10_bias": 0.0,
            "psa10_max": 0.0,
            "psa30_bias": 0.0,
            "psa30_max": 0.0,
            "gmpe": gmmod["gmpe"]["module"],
            "pgm2mi": gmmod["pgm2mi"]["module"],
            "software_version": info["processing"]["shakemap_versions"]["shakemap_revision"],
        }
        # Update infoDict with available values.
        gm = info["output"]["ground_motions"]
        for key,value in gm.items():
            infoDict[key] = value
        infoValues = tuple([infoDict[col] for col in COLUMNS])
        valueCols = ",".join("?"*len(infoValues))
        cmd = "INSERT"
        if replace:
            cmd += " OR REPLACE"
        try:
            self.cursor.execute("{0} INTO comcat_shakemaps({1}) VALUES({2})".format(cmd, insertCols, valueCols), infoValues)
            self.connection.commit()
        except sqlite3.IntegrityError as ex:
            logging.getLogger(__name__).debug(str(ex))
            logging.getLogger(__name__).debug(str(info))
        return

    def add_performance(self, stats, replace=False):
        """Add performance stats to database.

        :type stats: dict
        :param stats: Performance stats to add to database.
        """
        COLUMNS = (
            "comcat_id",
            "eew_server",
            "dm_id",
            "dm_timestamp",
            "gmpe",
            "fragility",
            "magnitude_threshold",
            "mmi_threshold",
            "area_damage",
            "area_alert",
            "area_cost_eew",
            "area_cost_noeew",
            "area_cost_perfecteew",
            "area_metric",
            "population_damage",
            "population_alert",
            "population_cost_eew",
            "population_cost_noeew",
            "population_cost_perfecteew",
            "population_metric",
            )

        insertCols = ", ".join(COLUMNS)
        perfValues = [stats[col] for col in COLUMNS]
        perfCols = ",".join("?"*len(perfValues))
        cmd = "INSERT"
        if replace:
            cmd += " OR REPLACE"
        try:
            self.cursor.execute("{0} INTO performance({1}) VALUES({2})".format(cmd, insertCols, perfCols), perfValues)
            self.connection.commit()
        except sqlite3.IntegrityError as ex:
            logging.getLogger(__name__).debug(str(ex))
            logging.getLogger(__name__).debug(str(stats))
        return

    def find_match(self, comcatId, server):
        """Find initial alert matching ComCat event.
        """
        import greatcircle
        
        MAX_DISTANCE_DEG = 3.0
        MAX_DISTANCE_KM = 150.0
        MAX_TIME_SECS = 15.0
        VS = 3.0e+3

        self.cursor.execute("SELECT * FROM comcat_events WHERE event_id=?", (comcatId,))
        event = self.cursor.fetchone()

        lat = event["latitude"]
        lon = event["longitude"]
        ot = dateutil.parser.parse(event["origin_time"])
        dt = datetime.timedelta(seconds=MAX_TIME_SECS)
        conditions = [
            "category=?",
            "message_type=?",
            "latitude BETWEEN ? AND ?",
            "longitude BETWEEN ? AND ?",
            "origin_time BETWEEN ? AND ?",
            ]
        values = (
            "live",
            "new",
            lat-MAX_DISTANCE_DEG, lat+MAX_DISTANCE_DEG,
            lon-MAX_DISTANCE_DEG, lon+MAX_DISTANCE_DEG,
            ot-dt, ot+dt,
            )
        self.cursor.execute("SELECT * FROM eew_alerts WHERE " + " AND ".join(conditions), values)
        alerts = self.cursor.fetchall()
        if 0 == len(alerts):
            return None
            
        # Get closest alert, ignoring deleted alerts
        minDist = 1.0e+30
        alertMatch = None
        for alert in alerts:
            # :TODO: Ignore deleted alerts
            # Ignore alerts from non-target servers
            if alert["server"]!="unknown" and alert["server"]!="eew2" and alert["server"]!=server:
                continue
            # Limit to distance range 
            dist = greatcircle.distance(lon, lat, alert["longitude"], alert["latitude"])
            if dist*1e-3 > MAX_DISTANCE_KM:
                continue
            distOT = abs((dateutil.parser.parse(alert["origin_time"])-ot).total_seconds())*VS
            if dist + distOT < minDist:
                alertMatch = alert
                minDist = dist
        return alertMatch

    def alerts(self, comcatId, server):
        """Get ShakeAlert alerts for event matching ComCat id.

        :type comcatId: str
        :param comcatId: ComCat event id.
        """
        alert = self.find_match(comcatId, server)
        if alert is None:
            return []
        
        # Get subsequent alerts matching id  and instance within 10 min
        timestamp = dateutil.parser.parse(alert["timestamp"])
        conditions = [
            "category=?",
            "(message_type=? OR message_type=?)",
            "event_id=?",
            "server=?",
            "timestamp BETWEEN ? AND ?",
            ]
        values = (
            "live",
            "new",
            "update",
            alert["event_id"],
            alert["server"],
            timestamp, timestamp+datetime.timedelta(minutes=10.0),
            )
        self.cursor.execute("SELECT * FROM eew_alerts WHERE " + " AND ".join(conditions), values)
        alerts = self.cursor.fetchall()
        return alerts

    def performance_stats(self, comcatId, server, gmpe, fragility, magnitude_threshold):
        """
        """
        conditions = [
            "comcat_id=?",
            "eew_server=?",
            "gmpe=?",
            "fragility=?",
            "magnitude_threshold=?",
            ]
        values = (
            comcatId,
            server,
            gmpe,
            fragility,
            magnitude_threshold,
            )
        self.cursor.execute("SELECT * FROM performance WHERE " + " AND ".join(conditions) + " ORDER BY mmi_threshold", values)
        # :TODO: :KLUDGE: to get structured array
        # Can we use PRAGMA TABLE_INFO to get column name and type (to map to numpy.dtype)
        dtype = [
            ("comcat_id", "|S32"),
            ("eew_server", "|S32"),
            ("dm_id", "int32"),
            ("dm_timestamp", "|S32"),
            ("gmpe", "|S32"),
            ("fragility", "|S32"),
            ("magnitude_threshold", "float32"),
            ("mmi_threshold", "float32"),
            ("area_damage", "float32"),
            ("area_alert", "float32"),
            ("area_cost_eew", "float32"),
            ("area_cost_noeew", "float32"),
            ("area_cost_perfecteew", "float32"),
            ("area_metric", "float32"),
            ("population_damage", "float32"),
            ("population_alert", "float32"),
            ("population_cost_eew", "float32"),
            ("population_cost_noeew", "float32"),
            ("population_cost_perfecteew", "float32"),
            ("population_metric", "float32"),
            ]
        results = self.cursor.fetchall()
        nrows = len(results)
        stats = numpy.zeros(nrows, dtype=dtype)
        for iresult, result in enumerate(results):
            stats[iresult] = tuple([result[key] for key in result.keys()])
        return stats
    
    def most_recent_alert(self, server):
        """Get most recent alert in database.

        :type server: str
        :param server: Name of EEW server associated with alerts.

        :returns: Most recent alert for server in database.
        """
        self.cursor.execute("SELECT * from eew_alerts ORDER BY date(timestamp) DESC LIMIT 1")
        alert = self.cursor.fetchone()
        return alert
    
    def comcat_event(self, comcatId):
        """Get ComCat event information.

        :type comcatId: str
        :param comcatId: ComCat event id
        """
        self.cursor.execute("SELECT * FROM comcat_events WHERE event_id=?", (comcatId,))
        event = self.cursor.fetchone()
        return event

    
    def summary(self):
        """Returns string with database summary.
        """
        sout = ""
        for name,columns in TABLES:
            self.cursor.execute("PRAGMA TABLE_INFO({})".format(name))
            info = self.cursor.fetchall()
            self.cursor.execute("SELECT COUNT(*) FROM {}".format(name))
            nrows = self.cursor.fetchall()[0][0]

            sout += "Table {}\n".format(name)
            sout += "  Columns\n"
            for column in info:
                sout += "    {name:16} {type:16}\n".format(name=column[1], type=column[2])
            sout += "  Number of rows: {}\n".format(nrows)
        return sout

    def show_matches(self, server):
        """Show matches.
        """
        self.cursor.execute("SELECT * from comcat_events ORDER BY event_id")
        events = self.cursor.fetchall()
        for event in events:
            alert = self.find_match(event["event_id"], server)
            if alert:
                print("COMAT {event[event_id]} M{event[magnitude]:.2f} {event[longitude]:.3f} {event[latitude]:.3f} {event[origin_time]} ALERT {alert[event_id]} {alert[longitude]:.3f} {alert[latitude]:.3f} {alert[origin_time]} {alert[server]}".format(event=event, alert=alert))
            else:
                print("COMAT {event[event_id]} M{event[magnitude]:.2f} {event[longitude]:.3f} {event[latitude]:.3f} {event[origin_time]} ALERT None".format(event=event))
        return

# End of file
