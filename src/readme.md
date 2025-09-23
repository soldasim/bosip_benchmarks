# src/ readme

## Main scripts

- `main.jl`: Main script for manually running BOSIP. On cluster, the scripts is "main_scripts/" are used instead.
- `generate_starts.jl`: Script for generating experiment starting data.
- `plots.jl`: Scripts for plotting results.
- `calculate_score.jl`: Scripts for re-calculating the scores with different metrics.
- `opt_mmd_ref.jl`: Script for calculating reference "zero" values for the OptMMD metric.

## Other code

All other code is included by `include("include_code.jl")` in each of the main scripts.
