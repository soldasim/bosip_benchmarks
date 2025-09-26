include("main.jl")

axis_size() = (400, 300)

function get_run_label(abbr::String)
    # TODO
    group_labels = Dict([
        "standard" => "GP - output - MaxVar",
        "loglike" => "GP - loglike - MaxVar",
        "loglike-imiqr" => "GP - loglike - IMIQR",
        "eiv" => "GP - output - EIV",
        "eiig" => "GP - output - IMMD",
        "nongp" => "nonGP - output - MaxVar",
    ])
    return get(group_labels, abbr, abbr)
end

### Create a grid of plots for the individual problems
### (Edit plotted runs in `plotted groups`.)
#
# ## Important Keywords
# - xscale
# - yscale
# - max_iters
function plot_results(; save_plot=false, kwargs...)
    # TODO
    problems = [
        ABProblem()             SimpleProblem()         BananaProblem()         BimodalProblem()
        :legend                 SIRProblem()            DuffingProblem()        "DiffusionProblem()"
    ]

    nrows, ncols = size(problems)
    ax_width, ax_height = axis_size()
    fig = Figure(;
        size = (ax_width * ncols, ax_height * nrows),
    )

    for idx in CartesianIndices(problems)
        if problems[idx] isa AbstractProblem
            ps = AbstractProblem[problems[idx]]
            plot_result_axis!(fig[idx.I...], ps; legend=false, kwargs...)
        else
            val = problems[idx]
            plot_result_special!(fig[idx.I...], fig, val; kwargs...)
        end
    end

    # Set all columns and rows to the same size
    for col in 1:ncols
        colsize!(fig.layout, col, Relative(1 / ncols))
    end
    for row in 1:nrows
        rowsize!(fig.layout, row, Relative(1 / nrows))
    end

    save_plot && save(plot_dir() * "/all_problems.png", fig)
    return fig
end

### Create a single plot with the runs of the given problems
### (Edit plotted runs in `plotted groups`.)
plot_results(problem::AbstractProblem; kwargs...) = plot_results(AbstractProblem[problem]; kwargs...)
function plot_results(problems::AbstractVector; save_plot=false, kwargs...)
    fig = Figure()
    
    plot_result_axis!(fig[1,1], problems; kwargs...)
    
    save_plot && save(plot_dir() * "/" * plot_name(problems) * ".png", fig)
    return fig
end

function plot_result_special!(figpos::GridPosition, fig::Figure, symbol::Symbol; kwargs...)
    (symbol == :nothing) && return

    if symbol == :legend
        Legend(figpos, fig.content[1], "Legend"; titleposition=:left)

    else
        @warn "Unknown special plot symbol: $symbol"
        @assert false
    end
end
function plot_result_special!(figpos::GridPosition, fig::Figure, str::String; kwargs...)
    Label(figpos, str)
end

function plot_result_axis!(figpos::GridPosition, problems::AbstractVector{<:AbstractProblem};
    xscale = log10,
    yscale = log,
    legend = true,
    max_iters = typemax(Int),
)
    ################
    ### SETTINGS ###
    ################

    # TODO metric
    # metric = OptMMDMetric
    metric = TVMetric

    # a list of all groups is needed to keep plot colors consistent
    colors = Makie.wong_colors()
    all_groups = ["loglike-imiqr", "loglike", "standard", "eiv", "eiig", "nongp"]
    color_map = Dict(group => colors[i] for (i, group) in enumerate(all_groups))
    
    # TODO groups
    plotted_groups = ["loglike-imiqr", "loglike", "standard", "eiv", "eiig", "nongp"]
    # plotted_groups = ["loglike-imiqr", "loglike", "standard"]
    # plotted_groups = ["standard", "eiv", "eiig"]
    # plotted_groups = ["standard", "nongp"]

    # include log versions of the problems as well
    add_log_versions!(problems)
    
    title = problems[1] |> typeof |> nameof |> string # TODO
    title = title[1:end-7]  # remove "Problem" suffix
    ylabel = string(metric)[1:end-6]  # remove "Metric" suffix
    ###

    scores_by_group = load_stored_scores!(problems, metric)

    ax = Axis(figpos; xlabel="Iteration", ylabel, title, xscale, yscale)

    for group in plotted_groups
        haskey(scores_by_group, group) || continue
        scores = scores_by_group[group]

        color = color_map[group]
        # scores is a Vector of score histories (each is a Vector)
        # Pad with `missing` to equal length if needed
        maxlen = maximum(length.(scores))
        maxlen = min(maxlen, max_iters)
        padded = pad_with_missing.(scores, Ref(maxlen))
        arr = reduce(hcat, padded)

        if !allequal(length.(scores))
            max_run_len = maximum(length.(scores))
            successful = sum(length.(scores) .== max_run_len)
            failed = findall(length.(scores) .!= max_run_len)
            @warn "Scores for group \"$group\" have different lengths! Padding with `missing`.
            ($successful/$(length(scores)) runs have the max length of $max_run_len,
            runs $failed have shorter lengths.)"
        end

        # avoid zero if xscale if logarithmic
        if xscale(0) |> isinf
            xs = 1:maxlen |> collect
        else
            xs = 0:maxlen-1 |> collect
        end
        
        # Plot median line
        label = get_run_label(group)
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
    if metric == OptMMDMetric
        p = problems[1]
        ref_sample_count = 200
        @warn "CHECK THAT THE OptMMD METRIC HAS BEEN CALCULATED WITH $ref_sample_count SAMPLES!"
        ref_optmmd_file = joinpath(data_dir(p), "opt_mmd", "mmd_vals_$(ref_sample_count).jld2")
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
    end

    if legend
        # axislegend(ax)
        # axislegend(ax; position=:rt)
        # axislegend(ax; position=:rc)
        axislegend(ax; position=:lb)
    end

    return ax
end

plot_name(problems::AbstractVector) = join(plot_name.(problems), "_")
plot_name(problem::AbstractProblem) = string(typeof(problem))

function add_log_versions!(problems::AbstractVector{<:AbstractProblem})
    for p in problems
        pname = p |> typeof |> nameof |> string
        startswith(pname, "Log") && continue
        logpname = "Log" * pname
        if hasproperty(Main, Symbol(logpname))
            logp = getproperty(Main, Symbol(logpname))()
            @warn "Adding $logpname to the list of plotted problems."
            push!(problems, logp)
        end
    end
    return problems
end

function load_stored_scores!(problems::AbstractVector, metricT::Type{<:DistributionMetric}; kwargs...)
    return merge(load_stored_scores!.(problems, Ref(metricT); kwargs...)...)
end
function load_stored_scores!(problem::AbstractProblem, metricT::Type{<:DistributionMetric})
    # TODO dir
    dir = data_dir(problem)
    # dir = "data/archive/data_01/" * string(typeof(problems[1]))
    files = sort(Glob.glob(joinpath(dir, "*.jld2")))
    
    scores_by_group = Dict{String, Vector{Vector{Float64}}}()

    for file in files
        fname, suffixes... = split(basename(file), ".")
        group = split(fname, "_")[1]
        
        (length(suffixes) == 1) || continue                 # ignore backup data (which have an extra suffix)
        endswith(fname, metric_fname(metricT)) || continue  # only consider the metric files

        if !haskey(scores_by_group, group)
            scores_by_group[group] = Vector{Vector{Float64}}()
        end
        push!(scores_by_group[group], load(file, "score"))
    end

    return scores_by_group
end

function pad_with_missing(v::AbstractVector, len::Int)
    if length(v) < len
        return vcat(v, fill(missing, len - length(v)))
    else
        return v[1:len]
    end
end
