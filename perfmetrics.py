# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import logging
import numpy
from importlib import import_module

import fragilitycurves
import gdalraster

def timedelta_to_seconds(value):
    """Convert timedelta to floating point value in seconds.
    
    :type value: numpy.timedelta64
    :param value: Array of time differences.
    """
    return value.astype("timedelta64[us]").astype("float32")/1.0e+6
    
class CostSavings(object):
    """Cost savings weighted by area and population.
    """
    
    def __init__(self, config, maps):
        self.config = config
        self.maps = maps
        return

    def compute(self, event, shakemap, alerts, shakingTime, populationDensity, mmiAlertThreshold, plotAlertMaps=False):

        functionPath = self.config.get("mmi_predicted", "function").split(".")
        fn = getattr(import_module(".".join(functionPath[:-1])), functionPath[-1])
            
        shape = shakemap.data["mmi"].shape
        warningTimeZero = numpy.zeros((1,), dtype="timedelta64[s]")
        warningTime = gdalraster.NO_DATA_VALUE * numpy.ones(shape, dtype="timedelta64[s]")
        mmiPred = gdalraster.NO_DATA_VALUE * numpy.ones(shape, numpy.float32)
        
        shakingTimeRel = timedelta_to_seconds(shakingTime - numpy.datetime64(event["origin_time"]))

        thresholdReached = False
        magnitudeThreshold = self.config.getfloat("alerts", "magnitude_threshold")
        for alert in alerts:
            if numpy.datetime64(alert["timestamp"]) > numpy.max(shakingTime):
                # Skip alerts with no positive warning times in
                # domain. Changes in estimated earthquake location
                # could result in later alerts having positive warning
                # times.
                logging.getLogger(__name__).debug("Skipping alert version {ver} with no positive warning times.".format(ver=alert["version"]))
                continue
            if not thresholdReached and alert["magnitude"] < magnitudeThreshold:
                continue
            else:
                if not thresholdReached:
                    tstamp = numpy.datetime64(alert["timestamp"])
                    wtime = timedelta_to_seconds(tstamp - numpy.datetime64(event["origin_time"]))
                    msg = "Alert threshold reached at {tstamp}, {wtime:.1f}s after origin time.".format(tstamp=tstamp, wtime=wtime)
                    print msg
                    logging.getLogger(__name__).info(msg)
                thresholdReached = True

            mmiPredCur = fn(alert, shakemap.data, dict(self.config.items("mmi_predicted")))
            warningTimeCur = shakingTime - numpy.datetime64(alert["timestamp"])
            
            if plotAlertMaps:
                filename = self.config.get("files", "analysis_event").replace("[EVENTID]", event["event_id"])
                values = [
                    ("mmi_pred", mmiPredCur,),
                    ("warning_time", timedelta_to_seconds(warningTimeCur),),
                ]
                gdalraster.write(filename, values, shakemap.num_lon(), shakemap.num_lat(), shakemap.spatial_ref(), shakemap.geo_transform())
                self.maps.load_data(event["event_id"], alert=alert)
                tafterOT = timedelta_to_seconds(numpy.datetime64(alert["timestamp"])-numpy.datetime64(event["origin_time"]))
                self.maps.mmi_warning_time(tafterOT)
            
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
        fragility = getattr(import_module(".".join(objectPath[:-1])), objectPath[-1])()
        costDamage = fragility.cost_damage(mmiObs)
        costActionObs = fragility.cost_action(mmiObs)
        costNoEEW = costDamage
        costPerfectEEW = costDamage*(costDamage < costActionObs) + costActionObs*(costDamage >= costActionObs)
        costEEW = fragility.cost_action(mmiPred)*(mmiPred >= mmiAlertThreshold) + costDamage*(mmiPred < mmiAlertThreshold)
        costSavings = costNoEEW - costEEW
        costSavingsPerfect = costNoEEW - costPerfectEEW

        pixelArea = shakemap.pixel_area(self.config.get("shakemap", "projection"))
        costSavingsArea = numpy.sum(pixelArea * costSavings)
        costSavingsPop = numpy.sum(populationDensity * pixelArea * costSavings)
        costSavingsPerfectArea = numpy.sum(pixelArea * costSavingsPerfect)
        costSavingsPerfectPop = numpy.sum(populationDensity * pixelArea * costSavingsPerfect)
        metricArea = costSavingsArea / costSavingsPerfectArea
        metricPop = costSavingsPop / costSavingsPerfectPop
        areaAlert = numpy.sum(pixelArea * (mmiPred >= mmiAlertThreshold))
        areaDamage = numpy.sum(pixelArea * (costDamage > 0.0))
        popAlert = numpy.sum(pixelArea * populationDensity * (mmiPred >= mmiAlertThreshold))
        popDamage = numpy.sum(pixelArea * populationDensity * (costDamage > 0.0))

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
            ("warning_time", timedelta_to_seconds(warningTime),),
            ("shaking_time", shakingTimeRel,),
            ("population_density", populationDensity,),
            ("cost_no_eew", costNoEEW,),
            ("cost_perfect_eew", costPerfectEEW,),
            ("cost_eew", costEEW,),
            ("alert_category", alertCategory,),
            ]
        filename = self.config.get("files", "analysis_event").replace("[EVENTID]", event["event_id"])
        gdalraster.write(filename, values, shakemap.num_lon(), shakemap.num_lat(), shakemap.spatial_ref(), shakemap.geo_transform())

        metrics = {
            "area_alert": areaAlert,
            "population_alert": popAlert,
            "area_damage": areaDamage,
            "population_damage": popDamage,
            "metric_area": metricArea,
            "metric_population": metricPop,
            }
        return metrics

# End of file
