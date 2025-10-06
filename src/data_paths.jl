
# Data are stored in: "data/{problem}/{run name}_{run index}_{some suffix}.jld2"

# default problem name
get_name(problem::AbstractProblem) = problem |> typeof |> string

data_dir(problem::AbstractProblem) = "data/" * get_name(problem)

data_base_filepath(problem::AbstractProblem, run_name::String, run_idx) = data_dir(problem) * "/" * base_filename(problem, run_name, run_idx)

base_filename(problem::AbstractProblem, run_name::String, run_idx::Int) = run_name * "_" * string(run_idx)
base_filename(problem::AbstractProblem, run_name::String, run_idx::Nothing) = run_name

starts_dir(problem::AbstractProblem) = data_dir(problem) * "/starts"

plot_dir() = "plots"
