###
### Plotting scripts for benchmark results
###

include("main.jl")

axis_size() = (400, 300)

# Set larger font sizes for better readability
function set_theme_fonts!(; base_fontsize=20)
    set_theme!(
        fontsize = base_fontsize,
        Axis = (
            titlesize = base_fontsize + 2,
            xlabelsize = base_fontsize,
            ylabelsize = base_fontsize,
            xticklabelsize = base_fontsize - 2,
            yticklabelsize = base_fontsize - 2,
        ),
        Legend = (
            titlesize = base_fontsize,
            labelsize = base_fontsize - 1,
        ),
        Label = (
            textsize = base_fontsize,
        )
    )
end

function get_run_label(abbr::AbstractString)
    group_labels = Dict([
        "standard" => "GP - output - MaxVar",
        "loglike" => "GP - loglike - MaxVar",
        "loglike-imiqr" => "GP - loglike - IMIQR",
        "eiv" => "GP - output - EIV",
        "eiig" => "GP - output - IMMD",
        "nongp" => "nonGP - output - MaxVar",

        # TODO comment out
        # ### for estimator plot
        # "standard" => "GP - output - MaxVar - exp.",
        # "est" => "GP - output - MaxVar - MAP",
        # "loglike" => "GP - loglike - MaxVar - MAP",
        # "loglike-imiqr" => "GP - loglike - IMIQR - MAP",
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
function plot_results(; save_plot=false, base_fontsize=20, kwargs...)
    # Set theme with larger fonts
    set_theme_fonts!(base_fontsize=base_fontsize)
    
    ### problem grid
    problems = [
        ABProblem()     SimpleProblem()     BananaProblem()     BimodalProblem()
        :legend         SIRProblem()        DuffingProblem()    DiffusionProblem10()
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
        elseif problems[idx] isa Tuple
            ps = AbstractProblem[problems[idx]...]
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

    if save_plot
        save(plot_dir() * "/all_problems.png", fig)
        save(plot_dir() * "/all_problems.pdf", fig)
    end
    return fig
end

### Create a single plot with the runs of the given problems
### (Edit plotted runs in `plotted groups`.)
plot_results(problem::AbstractProblem; kwargs...) = plot_results(AbstractProblem[problem]; kwargs...)
function plot_results(problems::AbstractVector; save_plot=false, base_fontsize=20, kwargs...)
    # Set theme with larger fonts
    set_theme_fonts!(base_fontsize=base_fontsize)
    
    fig = Figure()
    
    plot_result_axis!(fig[1,1], problems; kwargs...)
    
    if save_plot
        save(plot_dir() * "/" * plot_name(problems) * ".png", fig)
        save(plot_dir() * "/" * plot_name(problems) * ".pdf", fig)
    end
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

    ### metric
    # metric = OptMMDMetric
    metric = TVMetric

    ### init data count
    init_data = 3

    ### max plotted iters
    maxiter = nothing
    # maxiter = 100

    # a list of all groups is needed to keep plot colors consistent
    colors = Makie.wong_colors()
    main_groups = ["loglike", "standard", "eiv", "eiig", "nongp", "alt"]
    color_map = Dict(group => colors[i] for (i, group) in enumerate(main_groups))
    
    # TODO groups
    # plotted_groups = ["loglike-imiqr", "loglike", "standard", "eiv", "eiig", "nongp"]
    plotted_groups = ["standard", "est", "loglike", "loglike-imiqr"]
    # plotted_groups = ["standard", "loglike"]
    # plotted_groups = ["standard", "eiv", "eiig"]
    # plotted_groups = ["standard", "nongp"]

    # include log versions of the problems as well
    add_log_variants!(problems)
    
    title = problems[1] |> typeof |> nameof |> string
    title = title[1:end-7]  # remove "Problem" suffix
    ylabel = string(metric)[1:end-6]  # remove "Metric" suffix
    ###

    scores_by_group = load_stored_scores(problems, metric)

    xticks = ([10^x for x in 0.5:0.5:2.0], [L"10^\mathbf{%$x}" for x in 0.5:0.5:2.0])
    # yminorticks = CustomLogMinorTicks(1)

    ax = Axis(figpos; xlabel="simulations",
        ylabel,
        title,
        xscale,
        yscale,
        xticks,
        # ygridvisible = false,
        # yminorticks,
        # yminorticksvisible = true,
        # yminorgridvisible = true,
        # yminorgridcolor = RGBAf(0, 0, 0, 0.12),
    )

    # Super hacky way to order the data same as in `plotted_groups`
    function prepare_data(problem_group, scores)
        pname, group = split(problem_group, "_")
        return pname, group, scores
    end
    plot_data_ = [prepare_data(problem_group, scores) for (problem_group, scores) in scores_by_group]
    plot_data = similar(plot_data_, 0)
    for group in plotted_groups
        for t in plot_data_
            (t[2] == group) && push!(plot_data, t)
        end
    end

    # plotting
    for (pname, group, scores) in plot_data
        label = get_run_label(group)
        color = get(color_map, group, nothing)
        style = :solid
        
        ### proxy variants
        if pname in ["AbsABProblem", "DiffusionProblem2", "ProxySIRProblem"]
            @assert group == "standard"
            color = color_map["alt"]
            # TODO
            label = "GP - alt. proxy - MaxVar"
            # label = "GP - bad proxy - MaxVar"
            # label = "GP - good proxy - MaxVar"
        end
        ### special styles
        if group == "est"
            color = color_map["standard"]
            style = :dash
        end
        if group == "loglike-imiqr"
            color = color_map["loglike"]
            style = :dash
        end

        # scores is a Vector of score histories (each is a Vector)
        # Pad with `missing` to equal length if needed
        maxlen = maximum(length.(scores))
        maxlen = min(maxlen, max_iters)
        
        if !allequal(length.(scores))
            max_run_len = maximum(length.(scores))
            successful = sum(length.(scores) .== max_run_len)
            failed = findall(length.(scores) .!= max_run_len)
            @warn "Scores for group \"$group\" have different lengths!
            ($successful/$(length(scores)) runs have the max length of $max_run_len,
            runs $failed have shorter lengths.)"

            ### pad or align
            # @warn "Padding with `missing`."
            # padded = pad_with_missing.(scores, Ref(maxlen))
            # arr = reduce(hcat, padded)
            @warn "Discarding run ends to include only iterations with full data."
            aligned = align_scores(scores)
            arr = reduce(hcat, aligned)
        else
            arr = reduce(hcat, scores)
        end

        xs = init_data:init_data+maxlen-1
        
        # Plot median line
        median_scores = mapslices(median∘skipmissing, arr; dims=2)[:]
        isnothing(maxiter) || (median_scores = median_scores[1:min(length(median_scores), maxiter)])
        lines!(ax, xs[eachindex(median_scores)], median_scores; label, color=color, linestyle=style, linewidth=2)
        
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

function add_log_variants!(problems::AbstractVector{<:AbstractProblem})
    for p in problems
        pname = p |> typeof |> nameof |> string
        startswith(pname, "Log") && continue
        logpname = "Log" * pname
        if hasproperty(Main, Symbol(logpname))
            logp = getproperty(Main, Symbol(logpname))()
            @warn "Adding $logpname to the list of plotted problems."
            (logp in problems) || push!(problems, logp)
        end
    end
    return problems
end

function load_stored_scores(problems::AbstractVector, metricT::Type{<:DistributionMetric}; kwargs...)
    dicts = load_stored_scores.(problems, Ref(metricT); kwargs...)

    # Find all keys and check for duplicates
    all_keys = reduce(vcat, [collect(keys(d)) for d in dicts])
    key_counts = Dict(k => count(==(k), all_keys) for k in unique(all_keys))
    common_keys = [k for (k, v) in key_counts if v > 1]
    @assert isempty(common_keys)

    return merge(dicts...)
end
function load_stored_scores(problem::AbstractProblem, metricT::Type{<:DistributionMetric})
    ### dir
    dir = data_dir(problem)
    # dir = "data/archive/data_01/" * string(typeof(problems[1]))
    files = sort(Glob.glob(joinpath(dir, "*.jld2")))
    
    scores_by_group = Dict{String, Vector{Vector{Float64}}}()

    for file in files
        fname, suffix = split(basename(file), ".")
        # group = split(fname, "_")[1]
        group = get_name(problem) * "_" * split(fname, "_")[1]
        
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

function align_scores(scores::AbstractVector{<:AbstractVector})
    lens = length.(scores)
    minlen = minimum(lens)
    return [s[1:minlen] for s in scores]
end

function common_prefix(strs::AbstractVector{String})
    isempty(strs) && return ""
    minlen = minimum(length.(strs))
    prefix = strs[1][1:minlen]
    for s in strs[2:end]
        for i in eachindex(prefix)
            if s[i] != prefix[i]
                prefix = prefix[1:i-1]
                break
            end
        end
    end
    return prefix
end

@kwdef struct CustomLogMinorTicks
    logstep::Real = 1
end

function Makie.get_minor_tickvalues(ticks::CustomLogMinorTicks, scale, tickvalues, vmin, vmax)
    return get_minor_ticks(ticks.logstep, vmin, vmax)
end
function get_minor_ticks(logstep, vmin, vmax)
    start = log(vmin) ÷ logstep * logstep
    stop = log(vmax) ÷ logstep * logstep
    ticks = [exp(t) for t in start:logstep:stop]
    # labels = [L"e^{%$t}" for t in start:logstep:stop]
    return ticks
end
