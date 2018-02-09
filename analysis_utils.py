# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os

def config_get_list(list_string):
    """Convert list as string to list.

    :type list_string: list
    :param list_string: List as string.
    :returns: List of strings.
    """
    l = [f.strip() for f in list_string[1:-1].split(",")]
    return l


def get_dir(params, name):
    """Get expanded directory name in [files] section.

    :type params: ConfigParser
    :param params: Configuration options

    :type name: str
    :param name: Option in [files] section.
    """
    return os.path.expanduser(params.get("files", name))


def timedelta_to_seconds(value):
    """Convert timedelta to floating point value in seconds.
    
    :type value: numpy.timedelta64
    :param value: Array of time differences.
    """
    return value.astype("timedelta64[us]").astype("float32")/1.0e+6


def analysis_label(params, eqid):
    """Get label for anlysis used in output filenames.
    
    :type params: ConfigParser
    :param params: Configuration options

    :type eqid: str
    :param eqid: ComCat earthquake id.
    """
    server = params.get("shakealert.production", "server")
    gmpe = params.get("mmi_predicted", "gmpe")
    fragility = params.get("fragility_curves", "object").split(".")[-1]
    magThreshold = params.getfloat("alerts", "magnitude_threshold")
    mmiThreshold = params.getfloat("alerts", "mmi_threshold")
    label = "{eqid}-{server}-{gmpe}-{fragility}-M{magThreshold:.1f}-MMI{mmiThreshold:.1f}".format(
        eqid=eqid, server=server, gmpe=gmpe, fragility=fragility, magThreshold=magThreshold, mmiThreshold=mmiThreshold)
    return label
    
