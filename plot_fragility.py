#!/usr/bin/env python

import numpy

import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
import matplotlib_extras

import eewperformance.fragility_curves as fragility_curves

MMI_MIDDLE_LOW = 3.5
MMI_MIDDLE_HIGH = 5.5
SIGMOID_SLOPE_LOW = 3.0
SIGMOID_SLOPE_HIGH = 1.5
LINEAR_WIDTH_LOW = 2.0
LINEAR_WIDTH_HIGH = 4.0

def pgaMMI(pga):
    """Use Wald et al. (1999) to convert PGA to MMI.
    
    :type pga: Numpy array
    :param pga: PGA in percent g.
    """
    MIN_FLOAT = 1.0e-20
    G_ACC = 9.80665
    
    mmiPGA = numpy.zeros(pga.shape)
    
    logY = numpy.log10(MIN_FLOAT + pga*G_ACC) # acceleration in cm/s**2
    maskLower = logY <= 1.8193
    mmiPGA = maskLower*(1.00 + 2.1987*logY) + ~maskLower*(-1.6582 + 3.6598*logY)
    
    mmi = numpy.clip(mmiPGA, 1.0, 10.0)
    return mmi
    

pga = numpy.logspace(-1.5, 2.0, 500, base=10.0)
mmi = pgaMMI(pga)

fragilityLowStep = fragility_curves.StepDamage(damage_mmi=MMI_MIDDLE_LOW)
fragilityLowLinear = fragility_curves.PublicFearAvoidance(
    damage_low_mmi=MMI_MIDDLE_LOW-0.5*LINEAR_WIDTH_LOW,
    damage_high_mmi=MMI_MIDDLE_LOW+0.5*LINEAR_WIDTH_LOW)
fragilityLowSigmoid = fragility_curves.SigmoidDamage(damage_middle_mmi=MMI_MIDDLE_LOW, slope=SIGMOID_SLOPE_LOW)

fragilityHighStep = fragility_curves.StepDamage(damage_mmi=MMI_MIDDLE_HIGH)
fragilityHighLinear = fragility_curves.PublicInjury(
    damage_low_mmi=MMI_MIDDLE_HIGH-0.5*LINEAR_WIDTH_HIGH,
    damage_high_mmi=MMI_MIDDLE_HIGH+0.5*LINEAR_WIDTH_HIGH)
fragilityHighSigmoid = fragility_curves.SigmoidDamage(damage_middle_mmi=MMI_MIDDLE_HIGH, slope=SIGMOID_SLOPE_LOW)
fragilityHighSigmoidB = fragility_curves.SigmoidDamage(damage_middle_mmi=MMI_MIDDLE_HIGH, slope=SIGMOID_SLOPE_HIGH)


pyplot.style.use("size-presentation")
pyplot.style.use("color-lightbg")

figure = pyplot.figure(figsize=(6.5, 6.5))
ax_factory = matplotlib_extras.axes.RectFactory(figure, nrows=2, margins=((0.6, 0, 0.1), (0.5, 0.8, 0.3)))

ax = figure.add_axes(ax_factory.rect(row=1))
ax.plot(mmi, fragilityLowStep.cost_damage(mmi), label="Step")
ax.plot(mmi, fragilityLowLinear.cost_damage(mmi), label="Linear")
ax.plot(mmi, fragilityLowSigmoid.cost_damage(mmi), label="Sigmoid")
ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
ax.set_ylabel("Relative Damage Cost")
ax.set_title("Low Damage Threshold (e.g., fear)")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax = figure.add_axes(ax_factory.rect(row=2))
ax.plot(mmi, fragilityHighStep.cost_damage(mmi), label="Step")
ax.plot(mmi, fragilityHighLinear.cost_damage(mmi), label="Linear")
ax.plot(mmi, fragilityHighSigmoid.cost_damage(mmi), label="Sigmoid k={:3.1f}".format(SIGMOID_SLOPE_LOW))
ax.plot(mmi, fragilityHighSigmoidB.cost_damage(mmi), label="Sigmoid k={:3.1f}".format(SIGMOID_SLOPE_HIGH))
ax.set_xlabel("MMI")
ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
ax.set_ylabel("Relative Damage Cost")
ax.set_title("High Damage Threshold (e.g., injuries)")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(loc="upper left")

figure.savefig("damage_functions.pdf")
pyplot.close(figure)
