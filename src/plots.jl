include("main.jl")

using ProgressMeter

plot_results(problem::AbstractProblem; kwargs...) = plot_results([problem]; kwargs...)

function plot_results(problems::AbstractVector; save_plot=false, xscale=log10, yscale=log)
    ### setings
    plotted_groups = ["loglike", "loglike-imiqr", "standard", "eiv", "eiig", "nongp"] # TODO
    group_labels = Dict([
        "standard" => "GP - output - MaxVar",
        "loglike" => "GP - loglike - MaxVar",
        "loglike-imiqr" => "GP - loglike - IMIQR",
        "eiv" => "GP - output - EIV",
        "eiig" => "GP - output - IMMD",
        "nongp" => "nonGP - output - MaxVar",
    ]) # TODO
    title = "AB Problem" # TODO
    ylabel = "OptMMD" # TODO
    ###

    # TODO
    scores_by_group = load_stored_scores!(problems)
    # scores_by_group = recalculate_scores!(problems; plotted_groups)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Iteration", ylabel, title, xscale, yscale)

    colors = Makie.wong_colors()  # or use any preferred color palette
    group_names = collect(keys(scores_by_group))
    color_map = Dict(group => colors[i] for (i, group) in enumerate(group_names))

    for group in plotted_groups
        scores = scores_by_group[group]

        color = color_map[group]
        # scores is a Vector of score histories (each is a Vector)
        # Pad with `missing` to equal length if needed
        maxlen = maximum(length.(scores))
        if !allequal(length.(scores))
            maxlen = maximum(length.(scores))
            successful = sum(length.(scores) .== maxlen)
            @warn "Scores for group \"$group\" have different lengths! Padding with `missing`.
            ($successful/$(length(scores)) runs have the max length of $maxlen)"
        end
        padded = [vcat(s, fill(missing, maxlen - length(s))) for s in scores]
        arr = reduce(hcat, padded)

        # avoid zero if xscale if logarithmic
        if xscale(0) |> isinf
            xs = 1:maxlen |> collect
        else
            xs = 0:maxlen-1 |> collect
        end
        
        # Plot median line
        label = get(group_labels, group, group)
        median_scores = mapslices(median∘skipmissing, arr; dims=2)[:]
        lines!(ax, xs, median_scores; label, color=color)
        
        # # Plot quantile band with alpha
        # lq = mapslices(x -> quantile(skipmissing(x), 0.1), arr; dims=2)[:]
        # uq = mapslices(x -> quantile(skipmissing(x), 0.9), arr; dims=2)[:]
        # band!(ax, xs, lq, uq; color=color, alpha=0.6)

        # # Plot min/max dotted lines
        # maxs = mapslices(x -> maximum(skipmissing(x)), arr; dims=2)[:]
        # mins = mapslices(x -> minimum(skipmissing(x)), arr; dims=2)[:]
        # lines!(ax, xs, maxs; color=color, linestyle=:dot, linewidth=1)
        # lines!(ax, xs, mins; color=color, linestyle=:dot, linewidth=1)
    end

    # plot reference opt_mmd values
    p = problems[1]
    ref_optmmd_file = joinpath(data_dir(p), "opt_mmd", "mmd_vals.jld2")
    if isfile(ref_optmmd_file)
        mmd_vals = load(ref_optmmd_file)["mmd_vals"]
        lq = quantile(mmd_vals, 0.1)
        med = median(mmd_vals)
        uq = quantile(mmd_vals, 0.9)

        hlines!(ax, [lq, uq]; color=:black, linestyle=:dot)
        hlines!(ax, [med]; color=:black, linestyle=:dash)
    else
        @warn "Reference OptMMD file not found: $ref_optmmd_file"
    end

    # axislegend(ax)
    axislegend(ax; position=:rc)

    save_plot && save(plot_dir() * "/" * plot_name(problems) * ".png", fig)
    fig
end

plot_name(problems::AbstractVector) = join(plot_name.(problems), "_")
plot_name(problem::AbstractProblem) = string(typeof(problem))

function load_stored_scores!(problems::AbstractVector; kwargs...)
    return merge(load_stored_scores!.(problems; kwargs...)...)
end
function load_stored_scores!(problem::AbstractProblem)
    # TODO dir
    dir = data_dir(problem)
    # dir = "data/archive/data_01/" * string(typeof(problems[1]))
    files = sort(Glob.glob(joinpath(dir, "*.jld2")))
    
    scores_by_group = Dict{String, Vector{Vector{Float64}}}()

    for file in files
        fname = split(basename(file), ".")[1]
        group = split(fname, "_")[1]
        
        # startswith(fname, "start") && continue  # skip start files
        endswith(fname, "metric") || continue     # only consider the metric files

        data = load(file)
        @assert haskey(data, "score")
        if !haskey(scores_by_group, group)
            scores_by_group[group] = Vector{Vector{Float64}}()
        end
        push!(scores_by_group[group], data["score"])
    end

    return scores_by_group
end

function recalculate_scores!(problems::AbstractVector; kwargs...)
    return merge(recalculate_scores!.(problems; kwargs...)...)
end
function recalculate_scores!(problem::AbstractProblem; plotted_groups=nothing)
    dir = data_dir(problem)
    files = sort(Glob.glob(joinpath(dir, "*.jld2")))

    scores_by_group = Dict{String, Vector{Vector{Float64}}}()

    # TODO
    xs = rand(x_prior(problem), 20 * 10^x_dim(problem))
    ws = exp.( (0.) .- logpdf.(Ref(x_prior(problem)), eachcol(xs)) )

    # metric = MMDMetric(;
    #     kernel = with_lengthscale(GaussianKernel(), (bounds[2] .- bounds[1]) ./ 3),
    # )
    # metric = OptMMDMetric(;
    #     kernel = GaussianKernel(),
    #     bounds,
    #     algorithm = BOBYQA(),
    # )
    metric = TVMetric(;
        grid = xs,
        ws = ws,
    )

    # TODO
    sampler = AMISSampler(;
            iters = 10,
            proposal_fitter = BOSIP.AnalyticalFitter(), # re-fit the proposal analytically
            # proposal_fitter = OptimizationFitter(;      # re-fit the proposal by MAP optimization
            #     algorithm = NEWUOA(),
            #     multistart = 24,
            #     parallel = parallel(),
            #     static_schedule = true, # issues with PRIMA.jl
            #     rhoend = 1e-2,
            # ),
            # gauss_mix_options = nothing,                # use Laplace approximation for the 0th iteration
            gauss_mix_options = GaussMixOptions(;       # use Gaussian mixture for the 0th iteration
                algorithm = BOBYQA(),
                multistart = 24,
                parallel = parallel(),
                static_schedule = true, # issues with PRIMA.jl
                cluster_ϵs = nothing,
                rel_min_weight = 1e-8,
                rhoend = 1e-4,
            ),
        )

    for file in files
        fname = split(basename(file), ".")[1]
        group = split(fname, "_")[1]
        run_idx = split(fname, "_")[2]
        
        # only calculate for the groups from the list
        if !(group in plotted_groups)
            scores_by_group[group] = Vector{Float64}[] # just to keep plot colors consistent
            continue
        end
        
        # startswith(fname, "start") && continue  # skip start files
        endswith(fname, "iters") || continue  # only consider the iters files

        iters = load(file)["problems"]

        if !haskey(scores_by_group, group)
            scores_by_group[group] = Vector{Vector{Float64}}()
        end

        # Recalculate the score using the MetricCallback
        scores = zeros(length(iters))
        @showprogress desc="Calculating run_$run_idx" for (idx, itr) in enumerate(iters)
            scores[idx] = calculate_score(metric, problem, itr, sampler)
        end
        push!(scores_by_group[group], scores)
    end

    return scores_by_group
end

function calculate_score(metric::SampleMetric, problem::AbstractProblem, p::BosipProblem, sampler::DistributionSampler)
    # TODO
    sample_count = 2 * 10^x_dim(problem)

    ref = reference(problem)
    if ref isa Function
        true_samples = sample_posterior_pure(sampler, ref, p.problem.domain, sample_count)
    else
        true_samples = ref
    end

    est_logpost = log_posterior_estimate()(p)
    approx_samples = sample_posterior_pure(sampler, est_logpost, p.problem.domain, sample_count)

    score = calculate_metric(metric, true_samples, approx_samples)
    return score
end
function calculate_score(metric::PDFMetric, problem::AbstractProblem, p::BosipProblem, sampler::DistributionSampler)
    ### retrieve the true and approx logpdf
    ref = reference(problem)
    @assert ref isa Function
    
    true_logpdf = ref
    approx_logpdf = log_posterior_estimate()(p)

    ### calculate metric
    score = calculate_metric(metric, true_logpdf, approx_logpdf)
    return score
end
