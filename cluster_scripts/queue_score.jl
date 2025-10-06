### Queue jobs to recalculate scores for finished runs.
### (based on the "data/{problem}/..._iters.jld2" files)

include("../src/calculate_score.jl")

function queue_score_calc(problem::AbstractProblem, run_name::String, metric::Type{<:DistributionMetric};
    selected_runs = nothing,
)
    pname = get_name(problem)

    iter_files = Glob.glob(data_dir(problem) * "/" * base_filename(problem, run_name, nothing) * "_*_iters.jld2")
    n_runs = isnothing(selected_runs) ? length(iter_files) : length(selected_runs)
    @info "Running $(n_runs) runs of the $(pname) ..."

    for iter_file in iter_files
        m = match(r"_(\d+)_iters\.jld2$", iter_file)
        run_idx = parse(Int, m.captures[1])

        if !isnothing(selected_runs)
            (run_idx in selected_runs) || continue
        end

        # main(; run_name, save_data=true, data, run_idx)
        @info "Queuing calculation of \"$(metric)\" for problem:\"$(pname)\", run_name:\"$(run_name)\", run_idx:\"$(run_idx)\""
        # TODO --mem (the code was failing with the `SimpleProblem` with the default 4G memory)
        job_name = "$(metric)_$(pname)_$(run_name)_$(run_idx)"
        Base.run(`sbatch -p cpulong --mem=16G --job-name=$job_name cluster_scripts/run_score.sh $pname $run_name $run_idx $metric`)
    end

    nothing
end
