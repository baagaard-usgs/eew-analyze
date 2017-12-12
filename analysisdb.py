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

TABLES = [
    ("dm_alerts", [
        "event_id INTEGER NOT NULL",
        "category TEXT NOT NULL",
        "instance TEXT NOT NULL",
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
        #"UNIQUE(event_id) ON CONFLICT FAIL",
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
    ("matches", [
        "comcat_id TEXT",
        "dm_id INTEGER",
        "dm_timestamp TEXT",
        "UNIQUE(comcat_id) ON CONFLICT FAIL",
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

    def init(self):
        """Create database.
        """
        for name, columns in TABLES[::-1]:
            self.cursor.execute("DROP TABLE IF EXISTS {}".format(name))
        for name,columns in TABLES:
            self.cursor.execute("CREATE TABLE {name} ({fields})".format(name=name, fields=", ".join(columns)))
        self.connection.commit()
        return

    def add_alerts(self, alerts):
        """Add alert info to database.
        """
        COLUMNS = (
            "event_id",
            "category",
            "instance",
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
        try:
            #self.cursor.executemany("INSERT INTO dm_alerts({}) VALUES({})".format(insertCols, valueCols), alerts)
            for alert in alerts:
                self.cursor.execute("INSERT INTO dm_alerts({}) VALUES({})".format(insertCols, valueCols), alert)
            self.connection.commit()
        except sqlite3.IntegrityError as ex:
            logging.getLogger(__name__).debug(str(ex))
            #logging.getLogger(__name__).debug(str(alerts))
            logging.getLogger(__name__).debug(str(alert))
        return
    
    def add_event(self, event):
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
        try:
            self.cursor.execute("INSERT INTO comcat_events({0}) VALUES({1})".format(insertCols, valueCols), eventValues)
            self.connection.commit()
        except sqlite3.IntegrityError as ex:
            logging.getLogger(__name__).debug(str(ex))
            logging.getLogger(__name__).debug(str(event))
        return

    def find_matches(self):
        """Find matches in database between alerts and ComCat events.
        """
        MAX_DISTANCE_DEG = 0.30
        MAX_TIME_SECS = 15.0
        VS = 3.0e+3
        DEG_TO_DIST = 110.0e+3

        self.cursor.execute("select * from comcat_events")
        events = self.cursor.fetchall()
        for event in events:
            lat = event["latitude"]
            lon = event["longitude"]
            ot = dateutil.parser.parse(event["origin_time"])
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
                lat-0.5*MAX_DISTANCE_DEG, lat+0.5*MAX_DISTANCE_DEG,
                lon-0.5*MAX_DISTANCE_DEG, lon+0.5*MAX_DISTANCE_DEG,
                ot-datetime.timedelta(seconds=0.5*MAX_TIME_SECS), ot+datetime.timedelta(seconds=0.5*MAX_TIME_SECS),
                )
            self.cursor.execute("SELECT * FROM dm_alerts WHERE " + " AND ".join(conditions), values)
            alerts = self.cursor.fetchall()
            if 0 == len(alerts):
                self.cursor.execute("INSERT INTO matches(comcat_id, dm_id, dm_timestamp) VALUES(?,?,?)", (event["event_id"], None, None))
                continue
            
            # Get closest alert, ignoring deleted alerts
            minDist = 1.0e+30
            alertMatch = None
            vs = 3.0e+3
            for alert in alerts:
                # :TODO: Ignore deleted alerts
                dist = (((alert["latitude"]-lat)*DEG_TO_DIST)**2
                            + ((alert["longitude"]-lon)*DEG_TO_DIST)**2
                            + 1.0*((dateutil.parser.parse(alert["origin_time"])-ot).total_seconds()*VS)**2)**0.5
                if dist < minDist:
                    alertMatch = alert
                    minDist = dist
            self.cursor.execute("INSERT INTO matches(comcat_id, dm_id, dm_timestamp) VALUES(?,?,?)", (event["event_id"], alertMatch["event_id"], alertMatch["timestamp"]))
        self.connection.commit()
        return

    def get_alerts(self, comcatId):
        """Get ShakeAlert alerts for event matching ComCat id.

        :type comcatId: str
        :param comcatId: ComCat event id.
        """
        import pdb
        pdb.set_trace()
        self.cursor.execute("SELECT * FROM matches WHERE comcat_id=?", (comcatId,))
        match = self.cursor.fetchone()

        alerts = []
        # Get subsequent alerts matching id  and instance within 10 min
        timestamp = dateutil.parser.parse(match["dm_timestamp"])
        conditions = [
            "category=?",
            "(message_type=? OR message_type=?)",
            "event_id=?",
            "timestamp BETWEEN ? AND ?",
            ]
        values = (
            "live",
            "new", "update",
            match["dm_id"],
            timestamp, timestamp+datetime.timedelta(minutes=10.0),
            )
        self.cursor.execute("SELECT * FROM dm_alerts WHERE " + " AND ".join(conditions), values)
        alerts = self.cursor.fetchall()
        return alerts
    
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

    def show_matches(self):
        """Show matches.
        """
        self.cursor.execute("SELECT * FROM matches ORDER BY comcat_id")
        matches = self.cursor.fetchall()
        for match in matches:
            self.cursor.execute("SELECT * FROM comcat_events WHERE event_id=?", (match["comcat_id"],))
            event = self.cursor.fetchone()
            self.cursor.execute("SELECT * FROM dm_alerts WHERE event_id=? AND timestamp=?", (match["dm_id"], match["dm_timestamp"],))
            alert = self.cursor.fetchone()
            if alert:
                print("COMAT {event[event_id]} {event[longitude]:.3f} {event[latitude]:.3f} {event[origin_time]} ALERT {alert[event_id]} {alert[longitude]:.3f} {alert[latitude]:.3f} {alert[origin_time]}".format(event=event, alert=alert))
            else:
                print("COMAT {event[event_id]} {event[longitude]:.3f} {event[latitude]:.3f} {event[origin_time]} ALERT None".format(event=event))
        return

# End of file
