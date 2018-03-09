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

import analysis_utils
import gdalraster

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
        warningTimeZero = numpy.zeros((1,), dtype="timedelta64[s]")
        warningTime = gdalraster.NO_DATA_VALUE * numpy.ones(shape, dtype="timedelta64[s]")
        mmiPred = gdalraster.NO_DATA_VALUE * numpy.ones(shape, numpy.float32)
        
        shakingTimeRel = analysis_utils.timedelta_to_seconds(shakingTime - numpy.datetime64(event["origin_time"]))

        thresholdReached = False
        for alert in alerts:
            if numpy.datetime64(alert["timestamp"]) > numpy.max(shakingTime):
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
                    tstamp = numpy.datetime64(alert["timestamp"])
                    wtime = analysis_utils.timedelta_to_seconds(tstamp - numpy.datetime64(event["origin_time"]))
                    msg = "Alert threshold reached at {tstamp}, {wtime:.1f}s after origin time.".format(tstamp=tstamp, wtime=wtime)
                    logging.getLogger(__name__).info(msg)
                    thresholdReached = True

            gmpe = self.config.get("mmi_predicted", "gmpe")
            gmice = self.config.get("mmi_predicted", "gmice")
            mmiPredCur = fn(alert, shakemap.data, gmpe, gmice)
            warningTimeCur = shakingTime - numpy.datetime64(alert["timestamp"])
            
            if plotAlertMaps:
                plotsDir = self.config.get("files", "plots_dir")
                if not os.path.isdir(plotsDir):
                    os.makedirs(plotsDir)
                filename = analysis_utils.analysis_label(self.config, self.eqId, magAlertThreshold, mmiAlertThreshold)+"_alert_snapshot.tiff"
                values = [
                    ("mmi_pred", mmiPredCur,),
                    ("warning_time", analysis_utils.timedelta_to_seconds(warningTimeCur),),
                ]
                gdalraster.write(filename, values, shakemap.num_lon(), shakemap.num_lat(), shakemap.spatial_ref(), shakemap.geo_transform())
                mapPanels = maps.MapPanels(self.config)
                mapPanels.load_data(event["event_id"], alert=alert)
                tafterOT = analysis_utils.timedelta_to_seconds(numpy.datetime64(alert["timestamp"])-numpy.datetime64(event["origin_time"]))
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

        # Compute costNoEEW, costEEW, costPerfectEEW, costSavings
        mmiObs = shakemap.data["mmi"]
        objectPath = self.config.get("fragility_curves", "object").split(".")
        fragilityOptions = dict(self.config.items("fragility_curves"))
        fragilityOptions.pop("object")
        fragilityOptions = {k: float(v) for k,v in fragilityOptions.iteritems()}
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
        areaMetric = (areaCostNoEEW - areaCostEEW) / (areaCostNoEEW - areaCostPerfectEEW)
        areaDamage = numpy.sum(pixelArea * (costDamage > 0.0))
        areaAlert = numpy.sum(pixelArea * (mmiPred >= mmiAlertThreshold))
        
        popCostNoEEW = numpy.sum(populationDensity * pixelArea * costNoEEW)
        popCostPerfectEEW = numpy.sum(populationDensity * pixelArea * costPerfectEEW)
        popCostEEW = numpy.sum(populationDensity * pixelArea * costEEW)
        popMetric = (popCostNoEEW - popCostEEW) / (popCostNoEEW - popCostPerfectEEW)
        popDamage = numpy.sum(pixelArea * populationDensity * (costDamage > 0.0))
        popAlert = numpy.sum(pixelArea * populationDensity * (mmiPred >= mmiAlertThreshold))

        # Alert categories TN(0),FN(1),FP(2),TP(3)
        alertCategory = numpy.zeros(costDamage.shape)
        maskTN = numpy.bitwise_and(mmiPred < mmiAlertThreshold, costDamage < costActionObs)
        maskFN = numpy.bitwise_and(mmiPred < mmiAlertThreshold, costDamage >= costActionObs)
        maskFP = numpy.bitwise_and(mmiPred >= mmiAlertThreshold, costDamage < costActionObs)
        maskTP = numpy.bitwise_and(mmiPred >= mmiAlertThreshold, costDamage >= costActionObs)
        alertCategory = maskTN*0.0 + maskFN*1.0 + maskFP*2.0 + maskTP*3.0
        
        values = [
            ("mmi_obs", shakemap.data["mmi"],),
            ("mmi_pred", mmiPred,),
            ("warning_time", analysis_utils.timedelta_to_seconds(warningTime),),
            ("shaking_time", shakingTimeRel,),
            ("population_density", populationDensity,),
            ("cost_no_eew", costNoEEW,),
            ("cost_perfect_eew", costPerfectEEW,),
            ("cost_eew", costEEW,),
            ("alert_category", alertCategory,),
            ]
        cacheDir = self.config.get("files", "analysis_cache_dir")
        if not os.path.isdir(cacheDir):
            os.makedirs(cacheDir)
        filename = os.path.join(cacheDir, "analysis_" + analysis_utils.analysis_label(self.config, event["event_id"], magAlertThreshold, mmiAlertThreshold) + ".tiff")
        gdalraster.write(filename, values, shakemap.num_lon(), shakemap.num_lat(), shakemap.spatial_ref(), shakemap.geo_transform())

        metrics = {
            "area_damage": areaDamage,
            "area_alert": areaAlert,
            "area_cost_eew": areaCostEEW,
            "area_cost_noeew": areaCostNoEEW,
            "area_cost_perfecteew": areaCostPerfectEEW,
            "area_metric": areaMetric,
            "population_damage": popDamage,
            "population_alert": popAlert,
            "population_cost_eew": popCostEEW,
            "population_cost_noeew": popCostNoEEW,
            "population_cost_perfecteew": popCostPerfectEEW,
            "population_metric": popMetric,
            }
        return metrics

# End of file
