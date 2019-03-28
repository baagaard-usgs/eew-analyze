# ======================================================================
# Brad T. Aagaard
# U.S. Geological Survey
# ======================================================================

"""Add custom applcation specific color names.
"""

from matplotlib.colors import ColorConverter, get_named_colors_mapping

def lightbg():
    """Colors for light background."""
    name_map = get_named_colors_mapping()
    colors = {
        "local:qarea_fc": name_map["c_ltorange"],
        "local:qarea_ec": name_map["c_orange"],

        "local:qpop_fc": name_map["c_ltblue"],
        "local:qpop_ec": name_map["c_blue"],
    }
    ColorConverter.colors.update(colors)
    return


def darkbg():
    """Colors for light background."""
    name_map = get_named_colors_mapping()
    colors = {
        "local:qarea_fc": name_map["c_ltorange"],
        "local:qarea_ec": name_map["c_orange"],

        "local:qpop_fc": name_map["c_ltblue"],
        "local:qpop_ec": name_map["c_blue"],
    }
    ColorConverter.colors.update(colors)
    return


if __name__ == "__main__":
    import matplotlib.pyplot as pyplot
    import matplotlib_extras.colors

    matplotlib_extras.colors.add_general()
    lightbg()
    pyplot.plot([0,1], [1,2], color="local:qarea_fc")
    pyplot.show()
    
