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
        "UNIQUE(event_id, version, instance) ON CONFLICT FAIL",
    ]),
    ("comcat_events", [
        "event_id TEXT NOT NULL",
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
        "id_comcat TEXT",
        "id_dm INTEGER",
        "FOREIGN KEY(id_comcat) REFERENCES comcat_events(event_id)",
        "FOREIGN KEY(id_dm) REFERENCES dm_alerts(event_id)",
        "UNIQUE(id_comcat) ON CONFLICT FAIL",
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
        """Add alert info to database.

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
        insertCols = ", ".join(COLUMNS)
        eventValues = (event.id, event.latitude, event.longitude, event.depth, event.time, event.magnitude, event["magType"], event.location)
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
        MAX_DISTANCE_KM = 100.0
        MAX_TIME_SECS = 15.0
        return
    
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

# End of file
