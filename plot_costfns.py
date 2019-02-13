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
        ("Linear, r=20", "fragility_fearavoidance_r20.cfg"),
    ],
    "Injury": [
        ("LinearW4", "fragility_injury.cfg"),
        ("LinearW2", "fragility_injury_w2.cfg"),
    ]
}
FIG_SIZE = (7.0, 8.5)
MARGINS = ((0.6, 0.8, 0.1), (0.5, 0.8, 0.6))
   

mmi = numpy.arange(1.0, 9.01, 0.01)

pyplot.style.use("size-presentation")
pyplot.style.use("color-lightbg")

figure = pyplot.figure(figsize=FIG_SIZE)
nrows = max(len(COSTFNS["Fear"]), len(COSTFNS["Injury"]))

ax_factory = matplotlib_extras.axes.RectFactory(figure, nrows=nrows, ncols=2, margins=MARGINS)


col = 1
for category in ["Fear", "Injury"]:

    row = 1
    for label, costfn_config in COSTFNS[category]:
    
        ax = figure.add_axes(ax_factory.rect(row=row, col=col))

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

        ax.plot(mmi, cost_action, label="Action")
        ax.plot(mmi, cost_damage, label="Damage")

        ax.set_title(label)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
        ax.set_xlim(1.0, 9.0)
        if row == len(COSTFNS[category]):
            ax.set_xlabel("MMI")

        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        ax.set_ylim(0.0, 1.02)
        ax.set_ylabel("Relative Cost")

        row += 1

    ax.text(0.25+0.5*(col-1), 0.98, "{} Avoidance".format(category),
            transform=figure.transFigure, ha="center", weight="bold")
    col += 1
        
figure.savefig("damage_functions.pdf")
pyplot.close(figure)
