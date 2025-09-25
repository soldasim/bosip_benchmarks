# ASSUMES pwd == sbi_benhmarks

include("../src/main.jl")

function queue_jobs(problem::AbstractProblem, run_name::String;
    selected_runs = nothing,
    continued = false,
)
    problem_name = string(typeof(problem))

    start_files = Glob.glob(starts_dir(problem) * "/start_*.jld2")
    n_runs = isnothing(selected_runs) ? length(start_files) : length(selected_runs)
    @info "Running $(n_runs) runs of the $(typeof(problem)) ..."

    for start_file in start_files
        m = match(r"start_(\d+)\.jld2$", start_file)
        run_idx = parse(Int, m.captures[1])
        # data = load(start_file, "data")

        if !isnothing(selected_runs)
            (run_idx in selected_runs) || continue
        end

        # main(; run_name, save_data=true, data, run_idx)
        @info "Queuing run: problem:\"$(problem_name)\", run_name:\"$(run_name)\", run_idx:\"$(run_idx)\""
        # TODO --mem (the code was failing with the `SimpleProblem` with the default 4G memory)
        job_name = "$(problem_name)_$(run_name)_$(run_idx)"
        job_name = continued ? job_name * "_cont" : job_name
        cont = continued ? 1 : 0
        Base.run(`sbatch -p cpulong --mem=16G --job-name=$job_name cluster_scripts/run.sh $problem_name $run_name $run_idx $cont`)
    end

    nothing
end
