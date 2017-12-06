from lxml import etree
import dateutil.parser
import re

filename = "data/nc72923380/dmevent_20171113.log"

def getChild(el, name):
    """Get child element 'name'. Raises IOError if more than one child element with name is found.

    :type el: XML element
    :param el: XML element

    :type name: str
    :param name: Name of child element.
    """
    elList = el.xpath(name)
    if len(elList) != 1:
        raise IOError("Expect one '%s' child element in XML element. Found %d '%s' elements in %s." % \
                      (name, len(elList), el.tag))
    return elList[0]

with open(filename, "r") as fin:
    bytes = fin.read()
    fin.close()
    pattern = [
        "(?P<timestamp>[0-9]{4}-[0-9]+-[0-9]+T[0-9]{2}:[0-9]{2}:[0-9]{2}.[0-9]+Z)[\s]*",
        "(?P<xml>\<\?xml[\s\S]+?</event_message>)"
        ]
    xmlBlocks = re.findall("".join(pattern), bytes)
    for tstamp,xmlBlock in xmlBlocks:
        elMsg = etree.fromstring(xmlBlock)
        category = elMsg.get("category")
        instance = elMsg.get("instance")
        messageType = elMsg.get("message_type")
        timestamp = dateutil.parser.parse(elMsg.get("timestamp"))
        version = int(elMsg.get("version"))

        elCoreInfo = getChild(elMsg, "core_info")
        id_dm = int(elCoreInfo.get("id"))

        elMag = getChild(elCoreInfo, "mag")
        magnitude = float(elMag.text)
        magnitude_type = elMag.get("units")

        latitude = float(getChild(elCoreInfo, "lat").text)
        longitude = float(getChild(elCoreInfo, "lon").text)
        depthKm = float(getChild(elCoreInfo, "depth").text)
        originTime = dateutil.parser.parse(getChild(elCoreInfo, "orig_time").text)
        likelihood = float(getChild(elCoreInfo, "likelihood").text)
        numStations = int(getChild(elCoreInfo, "num_stations").text)

        print id_dm, magnitude, timestamp
        
