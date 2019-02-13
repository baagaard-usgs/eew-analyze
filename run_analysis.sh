#!/bin/bash

./downloader.py --db-init=performance
rm -f analyzer.log

# pass kwargs to main
# if None, call _parse_command_line(), otherwise args = argparse.Namespace(*kwargs)


eqset="eqsets/sanfrancisco.cfg,eqsets/losangeles.cfg"

# FearAvoidanceLinear
costfns="fragility_fearavoidance.cfg"
./analyzer.py --config=${eqset},${costfns} --num-threads=16 --optimize-events --plot-event-maps=all --plot-event-figures=all
./analyzer.py --config=${eqset},${costfns},catalog_magnitude.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns},catalog_magnitude_bias.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns} --plot-summary-figures=all --plot-summary-maps=all --generate-report=full && mv -f report.pdf report_LA+SF-FearAvoidanceLinear.pdf

# FearAvoidanceStep
costfns="fragility_fearavoidance_step.cfg"
./analyzer.py --config=${eqset},${costfns} --num-threads=16 --optimize-events --plot-event-maps=all --plot-event-figures=all
./analyzer.py --config=${eqset},${costfns},catalog_magnitude.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns},catalog_magnitude_bias.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns} --plot-summary-figures=all --plot-summary-maps=all --generate-report=full && mv -f report.pdf report_LA+SF-FearAvoidanceStep.pdf

# FearAvoidanceSigmoid
costfns="fragility_fearavoidance_sigmoid.cfg"
./analyzer.py --config=${eqset},${costfns} --num-threads=16 --optimize-events --plot-event-maps=all --plot-event-figures=all
./analyzer.py --config=${eqset},${costfns},catalog_magnitude.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns},catalog_magnitude_bias.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns} --plot-summary-figures=all --plot-summary-maps=all --generate-report=full && mv -f report.pdf report_LA+SF-FearAvoidanceSigmoid.pdf

# FearAvoidanceLinearR20
costfns="fragility_fearavoidance_r20.cfg"
./analyzer.py --config=${eqset},${costfns} --num-threads=16 --optimize-events --plot-event-maps=all --plot-event-figures=all
./analyzer.py --config=${eqset},${costfns},catalog_magnitude.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns},catalog_magnitude_bias.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns} --plot-summary-figures=all --plot-summary-maps=all --generate-report=full && mv -f report.pdf report_LA+SF-FearAvoidanceLinearR20.pdf

# InjuryLinear
costfns="fragility_injury.cfg"
./analyzer.py --config=${eqset},${costfns} --num-threads=16 --optimize-events --plot-event-maps=all --plot-event-figures=all
./analyzer.py --config=${eqset},${costfns},catalog_magnitude.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns},catalog_magnitude_bias.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns} --plot-summary-figures=all --plot-summary-maps=all --generate-report=full && mv -f report.pdf report_LA+SF-InjuryLinear.pdf

# InjuryLinearW2
costfns="fragility_injury_w2.cfg"
./analyzer.py --config=${eqset},${costfns} --num-threads=16 --optimize-events --plot-event-maps=all --plot-event-figures=all
./analyzer.py --config=${eqset},${costfns},catalog_magnitude.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns},catalog_magnitude_bias.cfg --num-threads=16 --optimize-events --plot-event-maps=alert
./analyzer.py --config=${eqset},${costfns} --plot-summary-figures=all --plot-summary-maps=all --generate-report=full && mv -f report.pdf report_LA+SF-InjuryLinearW2.pdf
