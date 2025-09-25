
"""
    SaveCallback(; kwargs...)

Saves the run data after every iteration (by overwriting the data stored in the previous iteration).

# Keywords
- `dir::String`: Directory to save the data to.
- `filename::String`: The base filename (without the file extension).
- `run_idx::Union{Int, Nothing} = nothing`
- `continued::Bool = false`: Whether this is a continued run (i.e., not starting from scratch).
"""
@kwdef struct SaveCallback <: BosipCallback
    dir::String
    filename::String
    run_idx::Union{Int, Nothing} = nothing
    continued::Bool = false
end

function (cb::SaveCallback)(problem::BosipProblem; first, model_fitter, acq_maximizer, term_cond, options)
    mkpath(cb.dir)

    (first && cb.continued) && backup_data(cb)

    # problem & data
    save(cb.dir * "/" * cb.filename * "_problem.jld2", Dict(
        "problem" => problem,
    ))
    save(cb.dir * "/" * cb.filename * "_data.jld2", Dict(
        "run_idx" => cb.run_idx,
        "data" => (problem.problem.data.X, problem.problem.data.Y),
    ))

    # other
    save(cb.dir * "/" * cb.filename * "_extras.jld2", Dict(
        "run_idx" => cb.run_idx,
        "problem" => problem,
        "model_fitter" => model_fitter,
        "acq_maximizer" => acq_maximizer,
        "term_cond" => term_cond,
        "options" => options,
    ))

    # metric
    if !isempty(options.callback.callback.callbacks)
        metric_cb = options.callback.callback.callbacks[1]
        if metric_cb isa MetricCallback
            metricT = Base.typename(typeof(metric_cb.metric)).wrapper
            save(cb.dir * "/" * cb.filename * "_$(metric_fname(metricT)).jld2", Dict(
                "score" => metric_cb.score_history,
                "metric" => metric_cb,
            ))
        end
    end

    # iters data
    if first && cb.continued
        # initialization for continued runs
        # the `problem` is discarded - it should be equal to the last one in the loaded list
        iters = reload_iters(cb)
    elseif first
        iters = [problem]
    else
        iters = load(cb.dir * "/" * cb.filename * "_iters.jld2")["problems"]
        push!(iters, problem)
    end
    save(cb.dir * "/" * cb.filename * "_iters.jld2", Dict("problems" => iters))
end

### for continuing runs ###
function reload_iters(cb::SaveCallback)
    iters_file = cb.dir * "/" * cb.filename * "_iters.jld2"
    @assert isfile(iters_file)

    iters_data = load(iters_file)
    ps_ = iters_data["problems"]
    return ps_
end

function backup_data(cb::SaveCallback)
    ps = reload_iters(cb)
    iters = length(ps) - 1

    # backup previous data files by renaming them with a `.xxx` extension
    # with the number of iterations in the backed up data
    for file in filter(f -> startswith(f, cb.filename), readdir(cb.dir))
        fullpath = joinpath(cb.dir, file)
        if isfile(fullpath)
            cp(fullpath, fullpath * ".$iters"; force=true)
        end
    end
    return
end
