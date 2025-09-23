#!/bin/sh

# ASSUMES pwd == sbi_benchmarks
#
# ARGS[1] = problem: The name of the `AbstractProblem` subtype.
# ARGS[2] = run_name: The name of the whole run-set. The folders are named after this.
# ARGS[3] = run_idx: The index of this particular run. The starting data are selected based on this.
# ARGS[4] = metric: The name of the `DistributionMetric` subtype to calculate.

# something to avoid segmentation faults in juliacall
export PYTHON_JULIACALL_HANDLE_SIGNALS=yes

# activate python env
ml Python/3.10
. venv/bin/activate

# start julia
julia -e "
    include(\"activate.jl\")
    include(\"cluster_scripts/script_score.jl\")
" $1 $2 $3 $4
