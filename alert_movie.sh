#!/bin/bash
prefix=nc72948801-current
plots_dir=data/plots

mencoder mf://${plots_dir}/${prefix}-alert???_map_mmi_warning.png -mf fps=8 -ovc x264 -x264encopts subq=4:threads=auto -of lavf -o ${plots_dir}/${prefix}.mp4
