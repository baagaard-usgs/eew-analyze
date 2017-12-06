import sqlite3

filename = "dmlogs.sqlite"

connection = sqlite3.connect(filename)
cursor = connection.cursor()

columns = [
    "category TEXT",
    "instance TEXT",
    "message_type TEXT",
    "timestamp TEXT",
    "version INTEGER",
    "magnitude REAL",
    "magnitude_type TEXT",
    "latitude REAL",
    "longitude REAL",
    "depth_km REAL",
    "origin_time TEXT",
    "num_stations INTEGER",
    "id_dm INTEGER",
    ]

cursor.execute("CREATE TABLE %(name)s (%(fields)s)" % {"name": "dmalerts", "fields": ", ".join(columns)})
connection.commit()
#connection.close()


cursor.execute("PRAGMA TABLE_INFO(dmalerts)")
info = cursor.fetchall()
print info
