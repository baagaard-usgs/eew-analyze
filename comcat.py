# Stripped down versions of DetailEvent and Product from libcomcat.classes
# https://github.com/usgs/libcomcat/blob/master/libcomcat/classes.py
#
# CHANGES
#
# General
#   * Use urllib2 (Python 2) instead of urllib (Python 3)
#   * Conform to PEP8 style (lowercase and underscores for method names).
#
# DetailEvent
#   * Change constuctor to not automatically download file, but optionally load local file.
#   * Add fetch() to download file.
#   * Remove pandas dependency (rely only on standard modules).
#
# Product
#   * Rename get_content() to fetch().

import json
import urllib2
from datetime import datetime,timedelta
import re
#from enum import Enum

TIMEOUT_SECS = 30 # How many seconds to wait for download
WAIT_SECS = 3

class VersionOption(object):
    LAST = 1
    FIRST = 2
    ALL = 3
    PREFERRED = 4

    

class DetailEvent(object):
    """Wrapper around detailed event as returned by ComCat GeoJSON search results.
    """

    def __init__(self, filename=None):
        """Constructor.
        """
        if filename:
            self.load(filename)
        return

    def fetch(self, event_id, filename):
        """Fetch detailed event GeoJSON object from ComCat.
        
        Documentation for detailed event information is here:
        https://earthquake.usgs.gov/earthquakes/feed/v1.0/geojson_detail.php

        :type event_id: str
        :param event_id: ComCat event id (e.g., nc72923380)
        :type filename: str
        :param filename: Name of file for locally storing GeoJSON event information.
        """
        
        URL_TEMPLATE = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/detail/[EVENTID].geojson"

        url = URL_TEMPLATE.replace("[EVENTID]", event_id)
        
        try:
            fh = urllib2.urlopen(url, timeout=TIMEOUT_SECS)
            data = fh.read().decode('utf-8')
            fh.close()
            self._jdict = json.loads(data)
        except urllib2.HTTPError as htpe:
            time.sleep(WAIT_SECS)
            try:
                fh = urllib.request.urlopen(url,timeout=TIMEOUT_SECS)
                data = fh.read().decode('utf-8')
                fh.close()
                self._jdict = json.loads(data)
            except Exception as msg:
                raise Exception('Could not connect to ComCat server - %s.' % url).with_traceback(msg2.__traceback__)

        with open(filename, "w") as fh:
            fh.write(data)
        return

    def load(self, filename):
        """Load geojson event file.
        """
        with open(filename, "r") as fh:
            self._jdict = json.load(fh)
        return
            
    def __repr__(self):
        tpl = (self.id,str(self.time),self.latitude,self.longitude,self.depth,self.magnitude)
        return '%s %s (%.3f,%.3f) %.1f km M%.1f' % tpl

    @property
    def location(self):
        """Earthquake location string.
        """
        return self._jdict["properties"]["place"]

    @property
    def url(self):
        """ComCat URL.
        """
        return self._jdict["properties"]["url"]
    
    @property
    def latitude(self):
        """Authoritative origin latitude.
        """
        return self._jdict["geometry"]["coordinates"][1]

    @property
    def longitude(self):
        """Authoritative origin longitude.
        """
        return self._jdict["geometry"]["coordinates"][0]

    @property
    def depth(self):
        """Authoritative origin depth.
        """
        return self._jdict["geometry"]["coordinates"][2]

    @property
    def id(self):
        """Authoritative origin ID.
        """
        return self._jdict["id"]

    @property
    def time(self):
        """Authoritative origin time.
        """
        time_in_msec = self._jdict["properties"]["time"]
        time_in_sec = time_in_msec//1000
        msec = time_in_msec - (time_in_sec*1000)
        dtime = datetime.utcfromtimestamp(time_in_sec)
        dt = timedelta(milliseconds=msec)
        dtime = dtime + dt
        return dtime

    @property
    def magnitude(self):
        """Authoritative origin magnitude.
        """
        return self._jdict["properties"]["mag"]
    
    @property
    def properties(self):
        """List of summary event properties (retrievable from object with [] operator).
        """
        return list(self._jdict["properties"].keys())

    def has_product(self, product):
        """Return a boolean indicating whether given product can be extracted from DetailEvent.

        :param product:
          Product to search for.
        :returns:
          Boolean indicating whether that product exists or not.
        """
        if product in self._jdict["properties"]["products"]:
            return True
        return False

    def has_property(self, key):
        """Test to see whether a property with a given key is present in list of properties.
        
        :param key:
          Property to search for.
        :returns:
          Boolean indicating whether that key exists or not.
        """
        if key not in self._jdict["properties"]:
            return False
        return True

    def __getitem__(self, key):
        """Extract DetailEvent property using the [] operator.
        
        :param key:
          Property to extract.
        :returns:
          Desired property.
        """
        if key not in self._jdict["properties"]:
            raise AttributeError('No property %s found for event %s.' % (key,self.id))
        return self._jdict["properties"][key]
        
    
    def get_num_versions(self, productName):
        """Count the number of versions of a product (origin, shakemap, etc.) available for this event.
        
        :param productName:
          Name of product to query.
        :returns:
          Number of versions of a given product.
        """
        if not self.has_product(productName):
            raise AttributeError('Event %s has no product of type %s' % (self.id, productName))
        return len(self._jdict["properties"]["products"][productName])
    
    def get_product(self, productName, source='preferred', version=VersionOption.PREFERRED):
        """Retrieve a Product object from this DetailEvent.

        :param productName: Name of product (origin, shakemap, etc.) to retrieve.
        :param version: An enum value from VersionOption (PREFERRED,FIRST,ALL).
        :param source:
          Any one of: 
            - 'preferred' Get version(s) of products from preferred source.
            - 'all' Get version(s) of products from all sources.
            - Any valid source network for this type of product ('us','ak',etc.)
        :returns: List of Product objects.
        """
        import numpy
        
        if not self.has_product(productName):
            raise AttributeError('Event %s has no product of type %s' % (self.id,productName))

        dtype = [
            ("weight", "int32",),
            ("source", "S64",),
            ("time", "int64",),
            ("index", "int32",),
            ("version", "int32",),
        ]
        nrows = len(self._jdict["properties"]["products"][productName])
        df = numpy.zeros(nrows, dtype=dtype)
        df["weight"] = [product["preferredWeight"] for product in self._jdict["properties"]["products"][productName]]
        df["source"] = [product["source"] for product in self._jdict["properties"]["products"][productName]]
        df["time"] = [product["updateTime"] for product in self._jdict["properties"]["products"][productName]]
        df["index"] = list(range(0,nrows))
        df["version"] = 0
        
        # Add unique version number for each source, ordered by time.
        df.sort(order=["source", "time",])
        psources = []
        pversion = 1
        for idx,row in enumerate(df):
            if row["source"] not in psources:
                psources.append(row["source"])
                pversion = 1
            row["version"] = pversion
            pversion += 1

        if source == 'preferred':
            index = numpy.argmax(df["weight"])
            prefsource = self._jdict["properties"]["products"][productName][index]["source"]
            df = df[df["source"] == prefsource]
            df.sort(order="time")
        elif source == 'all':
            df.sort(order=["source','time"])
        else:
            df = df[df["source"] == source]
            df.sort(order='time')

        if not len(df):
            raise AttributeError('No products found for source "%s".' % source)

        products = []
        usources = set(df["source"])
        if source == 'all': #dataframe includes all sources
            for source in usources:
                df_source = df[df["source"] == source]
                df_source.sort(order="time")
                if version == VersionOption.PREFERRED:
                    df_source.sort(order=["weight","time"])
                    index = df_source.iloc[-1]["index"]
                    pversion = df_source.iloc[-1]["version"]
                    product = Product(productName,pversion,self._jdict["properties"]["products"][productName][index])
                    products.append(product)    
                elif version == VersionOption.LAST:
                    index = df_source.iloc[-1]["index"]
                    pversion = df_source.iloc[-1]["version"]
                    product = Product(productName,pversion,self._jdict["properties"]["products"][productName][index])
                    products.append(product)
                elif version == VersionOption.FIRST:
                    index = df_source.iloc[0]["index"]
                    pversion = df_source.iloc[0]["version"]
                    product = Product(productName,pversion,self._jdict["properties"]["products"][productName][index])
                    products.append(product)
                elif version == VersionOption.ALL:
                    for index,row in df_source.iterrows():
                        index = row["index"]
                        pversion = row["version"]
                        product = Product(productName,pversion,self._jdict["properties"]["products"][productName][index])
                        products.append(product)
                else:
                    raise(AttributeError("No VersionOption defined for %s" % version))
        else: #dataframe only includes one source
            if version == VersionOption.PREFERRED:
                df.sort(order=["weight","time"])
                index = df[-1]["index"]
                pversion = df[-1]["version"]
                product = Product(productName,pversion,self._jdict["properties"]["products"][productName][index])
                products.append(product)    
            elif version == VersionOption.LAST:
                index = df[-1]["index"]
                pversion = df[-1]["version"]
                product = Product(productName,pversion,self._jdict["properties"]["products"][productName][index])
                products.append(product)
            elif version == VersionOption.FIRST:
                index = df[0]["index"]
                pversion = df[0]["version"]
                product = Product(productName,pversion,self._jdict["properties"]["products"][productName][index])
                products.append(product)
            elif version == VersionOption.ALL:
                for index,row in df.iterrows():
                    index = row["index"]
                    pversion = row["version"]
                    product = Product(productName,pversion,self._jdict["properties"]["products"][productName][index])
                    products.append(product)
            else:
                raise(AttributeError("No VersionOption defined for %s" % version))

        return products


class Product(object):
    """Class describing a Product from detailed GeoJSON feed.  Products contain properties and file contents.
    """
    def __init__(self, product_name, version, product):
        """Create a product class from the product found within the detailed event GeoJSON.

        :param product_name:
          Name of Product (origin, shakemap, etc.)
        :param version:
          Best guess as to ordinal version of the product.
        :param product:
          Product data to be copied from DetailEvent.
        """
        self._product_name = product_name
        self._version = version
        self._product = product.copy()
        
    @property
    def preferred_weight(self):
        """The weight assigned to this product by ComCat.
        """
        return self._product["preferredWeight"]

    @property
    def source(self):
        """The contributing source for this product.
        """
        return self._product["source"]

    @property
    def update_time(self):
        """The datetime for when this product was updated.
        """
        time_in_msec = self._product["updateTime"]
        time_in_sec = time_in_msec//1000
        msec = time_in_msec - (time_in_sec*1000)
        dtime = datetime.utcfromtimestamp(time_in_sec)
        dt = timedelta(milliseconds=msec)
        dtime = dtime + dt
        return dtime

    @property
    def version(self):
        """The best guess for the ordinal version number of this product.
        """
        return self._version
    
    @property
    def properties(self):
        """List of product properties (retrievable from object with [] operator).
        """
        return list(self._product["properties"].keys())

    @property
    def contents(self):
        """List of product properties (retrievable from object with getContent() method).
        """
        return list(self._product["contents"].keys())
    
    def __getitem__(self,key):
        """Extract Product property using the [] operator.
        
        :param key:
          Property to extract.
        :returns:
          Desired property.
        """
        if key not in self._product["properties"]:
            raise AttributeError("No property %s found in %s product." % (key,self._product_name))
        return self._product["properties"][key]
    
    def __repr__(self):
        ncontents = len(self._product["contents"])
        tpl = (self._product_name,self.source,self.update_time,ncontents)
        return "Product %s from %s updated %s containing %i content files." % tpl

    def get_contents_matching(self,regexp):
        """Find all contents that match the input regex, ordered by shortest to longest.

        :param regexp:
          Regular expression which should match one of the content files in the Product.
        :returns:
          List of contents matching
        """
        contents = []
        if not len(self._product["contents"]):
            return contents
            
        for contentkey in self._product["contents"].keys():
            url = self._product["contents"][contentkey]["url"]
            parts = urllib2.urlparse.urlparse(url)
            fname = parts.path.split("/")[-1]
            if re.search(regexp+"$",fname):
                contents.append(fname)
        return contents
        
    def get_content_name(self, regexp):
        """Get the shortest filename matching input regular expression.

        For example, if the shakemap product has contents called grid.xml and grid.xml.zip, 
        and the input regexp is grid.xml, then grid.xml will be matched.

        :param regexp:
          Regular expression to use to search for matching contents.
        :returns:
          Shortest file name to match input regexp, or None if no matches found.
        """
        content_name = "a"*1000
        found = False
        for contentkey,content in self._product["contents"].items():
            if re.search(regexp+"$",contentkey) is None:
                continue
            url = content["url"]
            parts = urllib2.urlparse.urlparse(url)
            fname = parts.path.split("/")[-1]
            if len(fname) < len(content_name):
                content_name = fname
                found = True
        if found:
            return content_name
        else:
            return None

    def get_content_url(self,regexp):
        """Get the URL for the shortest filename matching input regular expression.

        For example, if the shakemap product has contents called grid.xml and grid.xml.zip, 
        and the input regexp is grid.xml, then grid.xml will be matched.

        :param regexp:
          Regular expression to use to search for matching contents.
        :returns:
          URL for shortest file name to match input regexp, or None if no matches found.
        """
        content_name = "a"*1000
        found = False
        content_url = ""
        for contentkey,content in self._product["contents"].items():
            if re.search(regexp+"$",contentkey) is None:
                continue
            url = content["url"]
            parts = urllib2.urlparse.urlparse(url)
            fname = parts.path.split("/")[-1]
            if len(fname) < len(content_name):
                content_name = fname
                content_url = url
                found = True
        if found:
            return content_url
        else:
            return None

        
    def fetch(self,regexp, filename):
        """Find and download the shortest file name matching the input regular expression.

        :param regexp:
          Regular expression which should match one of the content files in the Product.
        :param filename:
          Filename to which content should be downloaded.
        :returns:
          The URL from which the content was downloaded.
        :raises:
          Exception if content could not be downloaded from ComCat after two tries.
        """
        content_name = "a"*1000
        content_url = None
        for contentkey,content in self._product["contents"].items():
            if re.search(regexp+"$",contentkey) is None:
                continue
            url = content["url"]
            parts = urllib2.urlparse.urlparse(url)
            fname = parts.path.split("/")[-1]
            if len(fname) < len(content_name):
                content_name = fname
                content_url = url
        if content_url is None:
            raise AttributeError("Could not find any content matching input %s" % regexp)

        try:
            fh = urllib2.urlopen(url, timeout=TIMEOUT_SECS)
            data = fh.read()
            fh.close()
        except urllib2.HTTPError as htpe:
            time.sleep(WAIT_SECS)
            try:
                fh = urllib2.urlopen(url, timeout=TIMEOUT_SECS)
                data = fh.read()
                fh.close()
            except Exception as msg:
                raise Exception("Could not download %s from %s." % (content_name,url))
            
        with open(filename, "w") as f:
            f.write(data)
        return url
    
    def has_property(self,key):
        """Determine if this Product contains a given property.

        :param key:
          Property to search for.
        :returns:
          Boolean indicating whether that key exists or not.
        """
        if key not in self._product["properties"]:
            return False
        return True

# End of file
