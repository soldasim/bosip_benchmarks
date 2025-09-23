
"""
    SaveCallback(; kwargs...)

Saves the run data after every iteration (by overwriting the data stored in the previous iteration).

# Keywords
- `dir::String`: Directory to save the data to.
- `filename::String`: The base filename (without the file extension).
- `run_idx::Union{Int, Nothing} = nothing`
"""
@kwdef struct SaveCallback <: BosipCallback
    dir::String
    filename::String
    run_idx::Union{Int, Nothing} = nothing
end

function (cb::SaveCallback)(problem::BosipProblem; first, model_fitter, acq_maximizer, term_cond, options)
    mkpath(cb.dir)

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

    if first
        iters = [problem]
    else
        iters = load(cb.dir * "/" * cb.filename * "_iters.jld2")["problems"]
        push!(iters, problem)
    end
    save(cb.dir * "/" * cb.filename * "_iters.jld2", Dict("problems" => iters))
end
