# Workflow

1. `compile.R`: Compile heat flow and subduction zone segments
2. `cluster_ga.R` or `ga.R`: Use genetic algorithm to optimize Kriging parameters
3. `decode_genes.R`: Process genetic algorithm results
4. `interp_diff_ga.R` and `interp_diff_byeye.R`: Compute differences between Kriging and similarity
5. `summary_ga.R` and `summary_byeye.R`: Summarise interpolation differences
6. `visualize_diff_x.R`: Visualize interpolation differences