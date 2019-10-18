# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import logging
from importlib import import_module
import numpy

from . import analysis_utils
from . import gdalraster

class CostSavings(object):
    """Cost savings weighted by area and population.
    """
    
    def __init__(self, config):
        self.config = config
        return

    def compute(self, event, shakemap, alerts, shakingTime, populationDensity, magAlertThreshold, mmiAlertThreshold, plotAlertMaps=False):

        functionPath = self.config.get("mmi_predicted", "function").split(".")
        fn = getattr(import_module(".".join(functionPath[:-1])), functionPath[-1])
            
        shape = shakemap.data["mmi"].shape
        warningTimeZero = numpy.zeros((1,), dtype="timedelta64[us]")
        warningTime = gdalraster.NO_DATA_VALUE * 1.0e+6 * numpy.ones(shape, dtype="timedelta64[us]")
        mmiPred = gdalraster.NO_DATA_VALUE * numpy.ones(shape, numpy.float32)

        gmpe = self.config.get("mmi_predicted", "gmpe")
        gmice = self.config.get("mmi_predicted", "gmice")
        alertLatency = numpy.timedelta64(int(self.config.getfloat("alerts", "alert_latency_sec")*1.0e+3), "ms")
        if gmice == "default":
            gmice = shakemap.gmiceGrid
        
        thresholdReached = False
        for alert in alerts:
            alertTime = numpy.datetime64(alert["timestamp"]) + alertLatency

            if alertTime > numpy.max(shakingTime):
                # Skip alerts with no positive warning times in
                # domain. Changes in estimated earthquake location
                # could result in later alerts having positive warning
                # times.
                logging.getLogger(__name__).debug("Skipping alert version {ver} with no positive warning times.".format(ver=alert["version"]))
                continue
            if not thresholdReached and alert["magnitude"] < magAlertThreshold:
                continue
            else:
                if not thresholdReached:
                    wtime = analysis_utils.timedelta_to_seconds(alertTime - numpy.datetime64(event["origin_time"]))
                    msg = "Alert threshold reached at {tstamp}, {wtime:.1f}s after origin time.".format(tstamp=alertTime, wtime=wtime)
                    logging.getLogger(__name__).info(msg)
                    thresholdReached = True
                
            mmiPredCur = fn(alert, shakemap.data, gmpe, gmice)
            warningTimeCur = shakingTime - alertTime
            
            if plotAlertMaps:
                plotsDir = self.config.get("files", "plots_dir")
                if not os.path.isdir(plotsDir):
                    os.makedirs(plotsDir)
                filename = analysis_utils.analysis_event_label(self.config, self.eqId, magAlertThreshold, mmiAlertThreshold)+"_alert_snapshot.tiff"
                values = [
                    ("mmi_pred", mmiPredCur,),
                    ("warning_time", analysis_utils.timedelta_to_seconds(warningTimeCur),),
                ]
                gdalraster.write(filename, values, shakemap.num_lon(), shakemap.num_lat(), shakemap.spatial_ref(), shakemap.geo_transform())
                mapPanels = maps.MapPanels(self.config)
                mapPanels.load_data(event["event_id"], alert=alert)
                tafterOT = analysis_utils.timedelta_to_seconds(alertTime-numpy.datetime64(event["origin_time"]))
                mapPanels.mmi_warning_time(tafterOT)
            
            # Update alert time if greater than previous
            maskAlert = numpy.bitwise_and(warningTimeCur > warningTime, mmiPredCur >= mmiAlertThreshold)
            warningTime[maskAlert] = warningTimeCur[maskAlert]

            # Update predicted MMI if greater than previous AND
            # positive warning time. Assumes action will be taken if
            # alert threshold is reached (cannot be undone if later
            # updates reduce predicted MMI).
            maskMMI = numpy.bitwise_and(mmiPredCur > mmiPred, warningTimeCur >= warningTimeZero)
            mmiPred[maskMMI] = mmiPredCur[maskMMI]

        filename = "analysis_" + analysis_utils.analysis_event_label(self.config, event["event_id"], magAlertThreshold, mmiAlertThreshold) + ".tiff"
        metrics = self._cost(mmiPred, shakemap, warningTime, populationDensity, mmiAlertThreshold, filename)
        return metrics

    def _cost(self, mmiPred, shakemap, warningTime, populationDensity, mmiAlertThreshold, filename):
        """Compute cost savings metrics.
        """
        mmiObs = shakemap.data["mmi"]
        
        # Compute costNoEEW, costEEW, costPerfectEEW, costSavings
        objectPath = self.config.get("fragility_curves", "object").split(".")
        fragilityOptions = dict(self.config.items("fragility_curves"))
        fragilityOptions.pop("object")
        fragilityOptions.pop("label")
        fragilityOptions = {k: float(v) for k,v in fragilityOptions.items()}
        fragility = getattr(import_module(".".join(objectPath[:-1])), objectPath[-1])(**fragilityOptions)
        
        costDamage = fragility.cost_damage(mmiObs)
        costActionObs = fragility.cost_action(mmiObs)
        costNoEEW = costDamage
        costPerfectEEW = costDamage*(costDamage < costActionObs) + costActionObs*(costDamage >= costActionObs)
        costEEW = fragility.cost_action(mmiPred)*(mmiPred >= mmiAlertThreshold) + costDamage*(mmiPred < mmiAlertThreshold)

        pixelArea = shakemap.pixel_area(self.config.get("shakemap", "projection"))
        areaCostNoEEW = numpy.sum(pixelArea * costNoEEW)
        areaCostPerfectEEW = numpy.sum(pixelArea * costPerfectEEW)
        areaCostEEW = numpy.sum(pixelArea * costEEW)
        areaDamage = numpy.sum(pixelArea * (costDamage > 0.0))
        areaAlert = numpy.sum(pixelArea * (mmiPred >= mmiAlertThreshold))
        areaAlertPerfect = numpy.sum(pixelArea * (costDamage > costActionObs))
        
        popCostNoEEW = numpy.sum(populationDensity * pixelArea * costNoEEW)
        popCostPerfectEEW = numpy.sum(populationDensity * pixelArea * costPerfectEEW)
        popCostEEW = numpy.sum(populationDensity * pixelArea * costEEW)
        popDamage = numpy.sum(pixelArea * populationDensity * (costDamage > 0.0))
        popAlert = numpy.sum(pixelArea * populationDensity * (mmiPred >= mmiAlertThreshold))
        popAlertPerfect = numpy.sum(pixelArea * populationDensity * (costDamage > costActionObs))

        # Alert categories TN(0),FN(1),FP(2),TP(3)
        alertCategory = numpy.zeros(costDamage.shape)
        maskTN = numpy.bitwise_and(mmiPred < mmiAlertThreshold, costDamage < costActionObs)
        maskFN = numpy.bitwise_and(mmiPred < mmiAlertThreshold, costDamage >= costActionObs)
        maskFP = numpy.bitwise_and(mmiPred >= mmiAlertThreshold, costDamage < costActionObs)
        maskTP = numpy.bitwise_and(mmiPred >= mmiAlertThreshold, costDamage >= costActionObs)
        alertCategory = maskTN*0.0 + maskFN*1.0 + maskFP*2.0 + maskTP*3.0

        values = [
            ("mmi_obs", mmiObs,),
            ("mmi_pred", mmiPred,),
            ("warning_time", analysis_utils.timedelta_to_seconds(warningTime),),
            ("population_density", populationDensity,),
            ("cost_no_eew", costNoEEW,),
            ("cost_perfect_eew", costPerfectEEW,),
            ("cost_eew", costEEW,),
            ("alert_category", alertCategory,),
            ("pixel_area", pixelArea,),
            ]
        cacheDir = self.config.get("files", "analysis_cache_dir")
        if not os.path.isdir(cacheDir):
            os.makedirs(cacheDir)
        gdalraster.write(os.path.join(cacheDir, filename), values, shakemap.num_lon(), shakemap.num_lat(), shakemap.spatial_ref(), shakemap.geo_transform())

        metrics = {
            "area_damage": areaDamage,
            "area_alert": areaAlert,
            "area_alert_perfect": areaAlertPerfect,
            "area_costsavings_eew": areaCostNoEEW - areaCostEEW,
            "area_costsavings_perfecteew": areaCostNoEEW - areaCostPerfectEEW,
            "population_damage": popDamage,
            "population_alert": popAlert,
            "population_alert_perfect": popAlertPerfect,
            "population_costsavings_eew": popCostNoEEW - popCostEEW,
            "population_costsavings_perfecteew": popCostNoEEW - popCostPerfectEEW,
            }
        return metrics

# End of file
