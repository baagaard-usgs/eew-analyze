#!/bin/bash

for fileXML in data/*/*info.xml.gz; do
    fileJSON=`echo $fileXML | sed -e s/xml/json/g`
    rm $fileXML $fileJSON
done

exit 0
