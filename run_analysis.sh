#!/bin/bash

eqset="eqsets/sanfrancisco.cfg,eqsets/losangeles.cfg"

# FearAvoidanceLinear
./analyzer.py --config=${eqset} --num-threads=16 --optimize-events --plot-event-maps=all --plot-event-figures=all
./analyzer.py --config=${eqset},catalog_magnitude.cfg --num-threads=16--optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},catalog_magnitude_bias.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset} --plot-summary-figures=all --plot-summary-maps=all --generate-report && mv report.pdf report_LA+SF-FearAvoidanceLinear.pdf

# FearAvoidanceStep
costfns=fragility_fearavoidance_step.cfg
./analyzer.py --config=${eqset},${costfns} --num-threads=16 --optimize-events --plot-event-maps=all --plot-event-figures=all
./analyzer.py --config=${eqset},${costfns},catalog_magnitude.cfg --num-threads=16--optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns},catalog_magnitude_bias.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns} --plot-summary-figures=all --plot-summary-maps=all --generate-report && mv report.pdf report_LA+SF-FearAvoidanceStep.pdf

# FearAvoidanceSigmoid
costfns=fragility_fearavoidance_sigmoid.cfg
./analyzer.py --config=${eqset},${costfns} --num-threads=16 --optimize-events --plot-event-maps=all --plot-event-figures=all
./analyzer.py --config=${eqset},${costfns},catalog_magnitude.cfg --num-threads=16--optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns},catalog_magnitude_bias.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns} --plot-summary-figures=all --plot-summary-maps=all --generate-report && mv report.pdf report_LA+SF-FearAvoidanceSigmoid.pdf
