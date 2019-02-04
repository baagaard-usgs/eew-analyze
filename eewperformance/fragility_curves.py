# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import numpy

class Constant(object):
    """Constant cost function (independent of predictor variable).
    """

    def __init__(self, y):
        """Constructor.

        :type y: float
        :param y: Value of cost function.
        """
        self.y = y
        return

    def cost(self, x):
        """Compute cost.

        :type x: float
        :param x: Predictor variable.

        :returns: Cost at points given by x.
        """
        return self.y


class LinearRamp(object):
    """Linear ramp cost function that is yLow below xLow and then
    increases linearly to yHigh at xHigh.
    """

    def __init__(self, x_low=0.4, x_high=0.6, y_low=0.0, y_high=1.0):
        """Constructor.

        :type x_low: float
        :param x_low: x value at lower end of linear transition region.

        :type x_high: float
        :param x_high: x value at higher end of linear transition region.

        :type y_low: float
        :param y_low: y value at values below lower end of linear transition region.

        :type y_high: float
        :param y_high: y value at values above higher end of linear transition region.
        """
        self.x_low = x_low
        self.x_high = x_high
        self.y_low = y_low
        self.y_high = y_high
        return

    def cost(self, x):
        """Compute cost.

        :type x: float
        :param x: Predictor variable.

        :returns: Cost at points given by x.
        """
        mask_low = x <= self.x_low
        mask_high = x >= self.x_high
        mask_ramp = numpy.bitwise_and(~mask_low, ~mask_high)
        cost = mask_low*self.y_low \
               + mask_ramp*(x-self.x_low)/(self.x_high-self.x_low) \
               + mask_high*self.y_high
        return cost


class Sigmoid(object):
    """Sigmoid function that is y_low at x=-infinity, 0.5*(y_low+y_high) at x=x_middle, and
    y_high at x=+infinity.

    """

    def __init__(self, x_middle=0.4, y_low=0.0, y_high=1.0, slope=2.0):
        """Constructor.

        :type x_middle: float
        :param x_middle: x value at middle value between y_low and y_high.

        :type y_low: float
        :param y_low: y value at x=-infinity.

        :type y_high: float
        :param y_high: y value at x=+infinity.

        :type slope: float
        :param slope: Slope at x=x_middle
        """
        self.x_middle = x_middle
        self.y_low = y_low
        self.y_high = y_high
        self.slope = slope
        return

    def cost(self, x):
        """Compute cost.

        :type x: float
        :param x: Predictor variable.

        :returns: Cost at points given by x.
        """
        cost = 1.0 / (1.0+numpy.exp(-self.slope*(x-self.x_middle)))
        return cost


class CostActionDamage(object):
    """Abstract base class for computing cost of action and damage.
    """

    def __init__(self):
        """Constructor.
        """
        self.action = None
        self.damage = None
        return

    def cost_action(self, mmi_pred):
        """Compute cost of action.

        :type mmi_pred: Numpy array
        :param mmi_pred: Predicted MMI.

        :returns: Numpy array with cost of action at each point.
        """
        return self.action.cost(mmi_pred)

    def cost_damage(self, mmi_obs):
        """Compute cost of damage.

        :type mmi_obs: Numpy array
        :param mmi_obs: Observed MMI.

        :returns: Numpy array with cost of damage at each point.
        """
        return self.damage.cost(mmi_obs)


class StepDamage(CostActionDamage):
    """Damage with uniform cost of action and heavy-side step function for
    cost of damage.
    """

    def __init__(self, cost_action=0.1, damage_mmi=3.5):
        """Constructor.

        :type cost_action: float
        :param cost_action: Cost of action relative to maximum cost of damage.

        :type damage_mmi: float
        :param damage_mmi: MMI associated with onset of damage.
        """
        EPSILON = 1.0e-4

        CostActionDamage.__init__(self)

        self.action = Constant(cost_action)
        self.damage = LinearRamp(damage_mmi - EPSILON, damage_mmi + EPSILON, 0.0, 1.0)
        return


class LinearDamage(CostActionDamage):
    """Damage with uniform cost of action and linear increase in cost of
    damage over some finite interval.
    """

    def __init__(self, cost_action=0.1, damage_low_mmi=2.5, damage_high_mmi=5.5):
        CostActionDamage.__init__(self)

        self.action = Constant(cost_action)
        self.damage = LinearRamp(damage_low_mmi, damage_high_mmi, 0.0, 1.0)
        return


class SigmoidDamage(CostActionDamage):
    """Damage with uniform cost of action and sigmoid function for cost of
    damage.
    """

    def __init__(self, cost_action=0.1, damage_middle_mmi=3.5, damage_slope=3.0):
        """Constructor.

        :type cost_action: float
        :param cost_action: Cost of action relative to maximum cost of damage.

        :type damage_mmi: float
        :param damage_mmi: MMI associated with 50% damage.
        """
        CostActionDamage.__init__(self)

        self.action = Constant(cost_action)
        self.damage = Sigmoid(x_middle=damage_middle_mmi, slope=slope)
        return


# End of file
