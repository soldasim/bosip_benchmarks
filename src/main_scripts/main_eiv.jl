### The setup for using the EIV acquisition.
# `EIV` instead of `MaxVar`

using BOSS
using BOSIP
using Distributions
using KernelFunctions
using LinearAlgebra
using OptimizationPRIMA
using Bijectors

using JLD2
using Glob
using CairoMakie
using ProgressMeter

using Random
Random.seed!(555)

parallel() = false # PRIMA.jl causes StackOverflow when parallelized on Linux

include(pwd() * "/src/include_code.jl")

### PROBLEMS ###
#
# ABProblem
# SimpleProblem
# SIRProblem
#
# LogABProblem
# LogSimpleProblem
# LogSIRProblem

### START A NEW RUN ###
function main(problem::AbstractProblem; data=nothing, kwargs...)
    ### SETTINGS ###
    init_data_count = 3 # TODO

    ### INIT DATA ###
    if isnothing(data)
        data = get_init_data(problem, init_data_count)
    else
        @assert data isa AbstractMatrix{<:Real}
        sim = simulator(problem)
        X = data
        Y = reduce(hcat, (sim(x) for x in eachcol(X)))[:,:]
        data = BOSS.ExperimentData(X, Y)
    end

    # ### domain mean as only initial point
    # X = hcat(mean(domain(problem).bounds))
    # sim = simulator(problem)
    # Y = reduce(hcat, (sim(x) for x in eachcol(X)))[:,:]
    # data = ExperimentData(X, Y)

    @info "Initial data:"
    for (x, y) in zip(eachcol(data.X), eachcol(data.Y))
        println("  $x -> $y")
    end

    
    ### SURROGATE MODEL ###
    model = GaussianProcess(;
        mean = prior_mean(problem),
        kernel = BOSS.Matern52Kernel(),
        lengthscale_priors = get_lengthscale_priors(problem),
        amplitude_priors = get_amplitude_priors(problem),
        noise_std_priors = get_noise_std_priors(problem),
    )
    # model = NonstationaryGP(;
    #     mean = prior_mean(problem),
    #     lengthscale_model = BOSS.default_lengthscale_model(domain(problem).bounds, y_dim(problem)),
    #     amplitude_model = get_amplitude_priors(problem),
    #     noise_std_model = get_noise_std_priors(problem),
    # )
    
    
    ### ACQUISITION ###
    # acquisition = MaxVar()
    # acquisition = LogMaxVar()
    # acquisition = IMMD(;
    #     y_samples = 20,
    #     x_samples = 2 * 10^x_dim(problem),
    #     x_proposal = x_prior(problem),
    #     y_kernel = BOSS.GaussianKernel(),
    #     p_kernel = BOSS.GaussianKernel(),
    # )
    acquisition = EIV(
        y_samples = 20,
        x_samples = 2 * 10^x_dim(problem),
        x_proposal = x_prior(problem),
    )
    # acquisition = IMIQR(;
    #     p_u = 0.75,
    #     x_samples = 2 * 10^x_dim(problem),
    #     x_proposal = x_prior(problem),
    # )

    
    ### BOSIP PROBLEM ###
    bosip = construct_bosip_problem(;
        problem,
        data,
        acquisition,
        model,
    )

    iters = 100 # TODO
    data_max = size(data.X, 2) + iters

    return main(problem, bosip; data_max, kwargs...)
end

### CONTINUE A RUN ###
function main_continue(problem::AbstractProblem, run_name::String, run_idx::Union{Nothing, Int}; kwargs...)
    # # check the filename just to be sure
    # fname = basename(@__FILE__)
    # fname_split = split(fname, ['.', '_'])
    # @assert length(fname_split) == 3
    # @assert fname_split[2] == run_name

    # load the saved BOSIP problem
    if isnothing(run_idx)
        file = joinpath(data_dir(problem), "$(run_name)_problem.jld2")
    else
        file = joinpath(data_dir(problem), "$(run_name)_$(run_idx)_problem.jld2")
    end
    bosip = load(file)["problem"]
    @assert bosip isa BosipProblem

    # assert iters
    data_count = size(bosip.problem.data.X, 2)
    @assert data_count >= 3 + 100 # TODO
    data_max = 3 + 200 # TODO

    # continue
    return main(problem, bosip; continued=true, run_name, run_idx, data_max, kwargs...)
end

function main(problem::AbstractProblem, bosip::BosipProblem;
    run_name = "test",
    save_data = false,
    metric = false,
    plots = false,
    run_idx = nothing,
    continued = false,
    data_max = 1,
)
    bounds = bosip.problem.domain.bounds

    ### ALGORITHMS ###
    model_fitter = OptimizationMAP(;
        algorithm = NEWUOA(),
        multistart = 24,
        parallel = parallel(),
        rhoend = 1e-4,
    )
    acq_maximizer = OptimizationAM(;
        algorithm = BOBYQA(),
        multistart = 24,
        parallel = parallel(),
        rhoend = 1e-4,
    )

    
    ### TERMINATION CONDITION ###
    term_cond = DataLimit(data_max)


    ### SAMPLER ###
    # sampler = RejectionSampler(;
    #     logpdf_maximizer = LogpdfMaximizer(;
    #         algorithm = BOBYQA(),
    #         multistart = 24,
    #         parallel = parallel(),
    #         static_schedule = true, # issues with PRIMA.jl
    #         rhoend = 1e-4,
    #     ),
    # )
    sampler = AMISSampler(;
        iters = 10,
        proposal_fitter = BOSIP.AnalyticalFitter(), # re-fit the proposal analytically
        # proposal_fitter = OptimizationFitter(;      # re-fit the proposal by MAP optimization
        #     algorithm = NEWUOA(),
        #     multistart = 6,
        #     parallel = parallel(),
        #     static_schedule = true, # issues with PRIMA.jl
        #     rhoend = 1e-2,
        # ),
        # gauss_mix_options = nothing,                # use Laplace approximation for the 0th iteration
        gauss_mix_options = GaussMixOptions(;       # use Gaussian mixture for the 0th iteration
            algorithm = BOBYQA(),
            multistart = 24,
            parallel = parallel(),
            cluster_ϵs = nothing,
            rel_min_weight = 1e-8,
            rhoend = 1e-4,
        ),
    )

    
    ### PERFORMANCE METRIC ###
    xs = rand(bosip.x_prior, 20 * 10^x_dim(problem))
    ws = exp.( (0.) .- logpdf.(Ref(bosip.x_prior), eachcol(xs)) )

    # metric_ = MMDMetric(;
    #     kernel = with_lengthscale(GaussianKernel(), (bounds[2] .- bounds[1]) ./ 3),
    # )
    # metric_ = OptMMDMetric(;
    #     kernel = GaussianKernel(),
    #     bounds,
    #     algorithm = BOBYQA(),
    #     rhoend = 1e-4,
    # )
    metric_ = TVMetric(;
        grid = xs,
        ws = ws,
        true_logpost = true_logpost(problem),
    )

    # Get a reference appropriate for the used metric is available.
    if metric_ isa PDFMetric
        ref = true_logpost(problem)
    else
        ref = reference_samples(problem)
        isnothing(ref) && (ref = true_logpost(problem))
    end
    @assert !isnothing(ref)

    if continued
        metric_cb = reload_metric_cb(metric_, problem, run_name, run_idx)
    else
        metric_cb = MetricCallback(;
            reference = ref,
            logpost_estimator = log_posterior_estimate(problem),
            sampler,
            sample_count = 2 * 10^x_dim(problem),
            metric = metric_,
        )
    end
    # first callback in `callbacks` (this is important for `SaveCallback`)
    callbacks = BosipCallback[]
    metric && push!(callbacks, metric_cb)


    ### PLOTS ###
    plot_cb = PlotModule.PlotCB(;
        problem,
        sampler,
        sample_count = 2 * 10^x_dim(problem),
        resolution = 200,
        plot_each = 10, # TODO
        save_plots = true,
    )
    plots && push!(callbacks, plot_cb)

    
    ### STORING RESULTS ###
    data_cb = SaveCallback(;
        dir = data_dir(problem),
        filename = base_filename(problem, run_name, run_idx),
        continued,
    )
    save_data && push!(callbacks, data_cb)

    options = BosipOptions(;
        callback = CombinedCallback(callbacks...),
    )

    
    ### RUN ###
    bosip!(bosip; model_fitter, acq_maximizer, term_cond, options)
    return bosip
end

### for continuing runs ###
function reload_metric_cb(metric::DistributionMetric, problem::AbstractProblem, run_name::String, run_idx::Union{Nothing, Int})
    dir = data_dir(problem)
    if isnothing(run_idx)
        metric_file = dir * "/$(run_name)_$(metric_fname(Base.typename(typeof(metric)).wrapper)).jld2"
    else
        metric_file = dir * "/$(run_name)_$(run_idx)_$(metric_fname(Base.typename(typeof(metric)).wrapper)).jld2"
    end
    
    metric_data = load(metric_file)
    score_ = metric_data["score"]
    metric_cb_ = metric_data["metric"]

    @assert typeof(metric_cb_.metric) == typeof(metric)
    return metric_cb_
end
