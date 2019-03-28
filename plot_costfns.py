#!/usr/bin/env python

import configparser
import collections
from importlib import import_module

import numpy

import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
import matplotlib_extras

import eewperformance.fragility_curves as fragility_curves

COSTFNS = {
    "Fear": [
        ("Linear", "fragility_fearavoidance.cfg"),
        ("Sigmoid", "fragility_fearavoidance_sigmoid.cfg"),
        ("Step", "fragility_fearavoidance_step.cfg"),
        #("Linear, r=20", "fragility_fearavoidance_r20.cfg"),
    ],
    "Injury": [
        ("LinearW4", "fragility_injury.cfg"),
        ("LinearW2", "fragility_injury_w2.cfg"),
    ]
}
FIG_SIZE = (7.5, 3.5)
MARGINS = ((0.6, 0.8, 0.1), (0.5, 0.8, 0.3))
   
pyplot.style.use("size-presentation")
pyplot.style.use("color-darkbg")

figure = pyplot.figure(figsize=FIG_SIZE)
nrows = 1
ncols = 2
ax_factory = matplotlib_extras.axes.RectFactory(figure, nrows=nrows, ncols=ncols, margins=MARGINS)

mmi = numpy.arange(1.0, 9.01, 0.01)

col = 1
row = 1
for category in ["Fear", "Injury"]:

    ax = figure.add_axes(ax_factory.rect(row=row, col=col))
    plot_action = True
    
    for label, costfn_config in COSTFNS[category]:
    
        config = configparser.SafeConfigParser()
        config.read(costfn_config)

        objectPath = config.get("fragility_curves", "object").split(".")
        fragilityOptions = dict(config.items("fragility_curves"))
        fragilityOptions.pop("object")
        fragilityOptions.pop("label")
        fragilityOptions = {k: float(v) for k,v in fragilityOptions.items()}
        fragility = getattr(import_module(".".join(objectPath[:-1])), objectPath[-1])(**fragilityOptions)
        
        cost_damage = fragility.cost_damage(mmi)
        cost_action = fragility.cost_action(mmi)

        if plot_action:
            ax.plot(mmi, cost_action, label="Action")
            plot_action = False
        ax.plot(mmi, cost_damage, label="{} Damage".format(label))

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
        ax.set_xlim(1.0, 9.0)
        ax.set_xlabel("MMI")

        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        ax.set_ylim(0.0, 1.02)
        ax.set_ylabel("Relative Cost")

        if row == 1 and col == 1:
            ax.legend(loc="center right")
        
    ax.set_title("{} Reduction".format(category), ha="center", weight="bold")
    col += 1
        
figure.savefig("costfns.pdf")
pyplot.close(figure)
