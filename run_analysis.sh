#!/bin/bash

# FearAvoidance
# SF
./analyzer.py --config=security/shakealert_login.cfg,eqsets/sanfrancisco.cfg --num-threads=16 --optimize-events --plot-maps=all --plot-figures=all
./analyzer.py --config=catalog_magnitude.cfg,eqsets/sanfrancisco.cfg --num-threads=16--optimize-events --plot-maps=alert
./analyzer.py --config=catalog_magnitude_bias.cfg,eqsets/sanfrancisco.cfg --num-threads=16 --optimize-events --plot-maps=alert
./analyzer.py --config=security/shakealert_login.cfg,eqsets/sanfrancisco.cfg --generate-report && mv report.pdf report_SF-AvoidFear.pdf

# LA
./analyzer.py --config=security/shakealert_login.cfg,eqsets/losangeles.cfg --num-threads=16 --optimize-events --plot-maps=all --plot-figures=all
./analyzer.py --config=catalog_magnitude.cfg,eqsets/losangeles.cfg --num-threads=16 --optimize-events --plot-maps=alert
./analyzer.py --config=catalog_magnitude_bias.cfg,eqsets/losangeles.cfg --num-threads=16 --optimize-events --plot-maps=alert
./analyzer.py --config=security/shakealert_login.cfg,eqsets/losangeles.cfg --generate-report && mv report.pdf report_LA-AvoidFear.pdf


# Injuries
# SF
./analyzer.py --config=security/shakealert_login.cfg,eqsets/sanfrancisco.cfg,fragility_injury.cfg --num-threads=16 --optimize-events --plot-maps=all --plot-figures=all
./analyzer.py --config=catalog_magnitude.cfg,eqsets/sanfrancisco.cfg,fragility_injury.cfg --num-threads=16 --optimize-events --plot-maps=alert
./analyzer.py --config=catalog_magnitude_bias.cfg,eqsets/sanfrancisco.cfg,fragility_injury.cfg --num-threads=16 --optimize-events --plot-maps=alert
./analyzer.py --config=security/shakealert_login.cfg,eqsets/sanfrancisco.cfg,fragility_injury.cfg --generate-report && mv report.pdf report_SF-AvoidInjury.pdf

# LA
./analyzer.py --config=security/shakealert_login.cfg,eqsets/losangeles.cfg,fragility_injury.cfg --num-threads=16 --optimize-events --plot-maps=all --plot-figures=all
./analyzer.py --config=catalog_magnitude.cfg,eqsets/losangeles.cfg,fragility_injury.cfg --num-threads=16 --optimize-events --plot-maps=alert
./analyzer.py --config=catalog_magnitude_bias.cfg,eqsets/losangeles.cfg,fragility_injury.cfg --num-threads=16 --optimize-events --plot-maps=alert
./analyzer.py --config=security/shakealert_login.cfg,eqsets/losangeles.cfg,fragility_injury.cfg --generate-report && mv report.pdf report_SF-AvoidInjury.pdf
