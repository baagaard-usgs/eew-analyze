from lxml import etree

# XML
# <tag name="mi_bias" value="-0.31" desc="magnitude bias for Intensity" />
# <tag name="bias" value="-0.04 -0.14 -0.30 -0.27 -0.45 " desc="magnitude bias (pga pgv psa03 psa10 psa30 )" />
#
# <tag name="pga_max" value="83.83" desc="Max value of grid" />
# <tag name="pgv_max" value="52.53" desc="Max value of grid" />
# <tag name="mi_max" value="7.89" desc="Max value of grid" />
# <tag name="psa03_max" value="112.13" desc="Max value of grid" />
# <tag name="psa10_max" value="54.42" desc="Max value of grid" />
# <tag name="psa30_max" value="5.72" desc="Max value of grid" />
#
# <tag name="GMPE" value="GMPE::BA08" desc="GMPE type" />
# <tag name="pgm2mi" value="GMICE::Wald99 - Wald, et al.; 1999" desc="Intensity Function" />
#
# <tag name="ShakeMap revision" value="3.5.687" desc="ShakeMap source code revision number" />
#
#
# JSON
# output
#    ground_motions
#      intensity
#          bias
#          max
#          units: "intensity"
#      pga
#          bias
#          max
#          units: "%g"
#      pgv
#          bias
#          max
#          units: "cm/s"
#      psa03
#          bias
#          max
#          units: "%g"
#      psa10
#          bias
#          max
#          units: "%g"
#      pds30
#          bias
#          max
#          units: "%g"
# processing
#     ground_motion_modules
#         gmpe
#             module
#             reference
#         pgm2mi
#             module
#             reference
#     shakemap_versions
#         shakemap_revision


with open("info.xml") as fxml:
        bytes = fxml.read()
        elRoot = etree.fromstring(bytes) # info

        # Values
        mmiBias = float(elRoot.xpath("tag[@name='mi_bias']")[0].get("value"))
        pgmBias = map(float, elRoot.xpath("tag[@name='bias']")[0].get("value").split())
        # :TODO: map pgmBias to pgaBias, pgvBias, psa03Bias, etc

        mmiMax = float(elRoot.xpath("tag[@name='mi_max']")[0].get("value"))
        pgvMax = float(elRoot.xpath("tag[@name='pgv_max']")[0].get("value"))
        pgaMax = float(elRoot.xpath("tag[@name='pga_max']")[0].get("value"))
        psa03Max = float(elRoot.xpath("tag[@name='psa03_max']")[0].get("value"))
        psa10Max = float(elRoot.xpath("tag[@name='psa10_max']")[0].get("value"))
        psa30Max = float(elRoot.xpath("tag[@name='psa30_max']")[0].get("value"))
        
        gmpeStr = elRoot.xpath("tag[@name='GMPE']")[0].get("value")
        pgm2miStr = elRoot.xpath("tag[@name='pgm2mi']")[0].get("value")
        shakemapVer = elRoot.xpath("tag[@name='ShakeMap revision']")[0].get("value")

        data = {
                "output": {
                        "ground_motions": {
                                "intensity": {
                                "bias": mmiBias,
                                "max": mmiMax,
                                "units": "intensity",
                        },
                        },
                        },
                        "processing": {
                                "ground_motion_modules": {
                                        "gmpe": {
                                                "module": gmpeStr,
                                                "reference": gmpeStr,
                                                },
                                                "pgm2mi": {
                                                        "module": pgm2miStr,
                                                        "reference": pgm2miStr,
                                                        },
                                },
                                "shakemap_versions": {
                                        "shakemap_revision": shakemapVer,
                                        }
                                        },
                }
        print data
