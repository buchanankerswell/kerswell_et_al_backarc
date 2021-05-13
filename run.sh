#!/bin/zsh
echo 'Compiling data ...'
Rscript compile.R
echo '\nPlotting global summary ...'
Rscript global_plots.R
echo '\nDecoding genes ...'
Rscript decode_genes.R
echo '\nInterpolating differences ...'
Rscript interp_diff.R
echo '\nInterpolating differences (GA) ...'
Rscript interp_diff_ga.R
echo '\nComparing variograms ...'
Rscript vgrm_compare.R
echo '\nSummarizing interpolations ...'
Rscript summary.R
echo '\nSummarizing interpolations (GA) ...'
Rscript summary_ga.R
echo '\nVisualizing results ...'
Rscript visualize_diff.R
Rscript visualize_diff_comp.R
echo '\nVisualizing results (GA) ...'
Rscript visualize_diff_ga.R
Rscript visualize_diff_comp_ga.R
echo '\nDone'
echo '\nRun ./knit.sh to render manuscript'