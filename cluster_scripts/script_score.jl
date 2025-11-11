@show ARGS

# ARGS are:
# #1: problem name (matches a subtype of `AbstractProblem`)
# #2: run name (describes the used BOSIP setup)
# #3: run index
# #4: metric type name
# #5: estimator function
const problem_name = ARGS[1]
const run_name = ARGS[2]
const run_idx = parse(Int, ARGS[3])
const metric_name = ARGS[4]
const estimator_name = ARGS[5]

### Calculate the scores
include("../src/calculate_score.jl")

const problem = getfield(Main, Symbol(problem_name))()
const metricT = getfield(Main, Symbol(metric_name))
const estimator = getfield(BOSIP, Symbol(estimator_name))

# main(problem; run_name, save_data=true, metric=true, data, run_idx)
calculate_scores(problem, estimator, metricT, run_name, run_idx; save_score=true)
