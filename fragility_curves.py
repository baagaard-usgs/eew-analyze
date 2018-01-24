# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import numpy

class Constant(object):

    def __init__(self, y):
        self.y = y
        return

    def cost(self, x):
        return self.y

class LinearRamp(object):
    """
    """
    
    def __init__(self, xLow=0.4, xHigh=0.6, yLow=0.0, yHigh=1.0):
        self.xLow = xLow
        self.xHigh = xHigh
        self.yLow = yLow
        self.yHigh = yHigh
        return

    def cost(self, x):
        maskLow = x <= self.xLow
        maskHigh = x >= self.xHigh
        maskRamp = numpy.bitwise_and(~maskLow, ~maskHigh)
        return maskLow*self.yLow + maskRamp*(x-self.xLow)/(self.xHigh-self.xLow) + maskHigh*self.yHigh


class PublicAnxiety(object):

    def __init__(self, costAction=0.1, mmiLow=2.5, mmiHigh=5.5):
        self.action = Constant(costAction)
        self.damage = LinearRamp(mmiLow, mmiHigh, 0.0, 1.0)
        return

    def cost_action(self, mmiPred):
        return self.action.cost(mmiPred)

    def cost_damage(self, mmiObs):
        return self.damage.cost(mmiObs)
    

class PublicInjury(object):

    def __init__(self, costAction=0.1, mmiLow=.5, mmiHigh=6.5):
        self.action = Constant(costAction)
        self.damage = LinearRamp(mmiLow, mmiHigh, 0.0, 1.0)
        return

    def cost_action(self, mmiPred):
        return self.action.cost(mmiPred)

    def cost_damage(self, mmiObs):
        return self.damage.cost(mmiObs)
    


# End of file
