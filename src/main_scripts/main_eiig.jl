### The setup for using the IMMD acquisition.
# `IMMD` instead of `MaxVar`

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

using Random
Random.seed!(555)

# TODO
log_posterior_estimate() = log_posterior_mean
# log_posterior_estimate() = log_approx_posterior

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

function main(problem::AbstractProblem; data=nothing, kwargs...)
    ### SETTINGS ###
    init_data_count = 3

    
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
    #     lengthscale_model = BOSS.default_lengthscale_model(bounds(problem), y_dim(problem)),
    #     amplitude_model = get_amplitude_priors(problem),
    #     noise_std_model = get_noise_std_priors(problem),
    # )
    
    
    ### ACQUISITION ###
    # acquisition = MaxVar()
    # acquisition = LogMaxVar()
    acquisition = IMMD(;
        y_samples = 20,
        x_samples = 2 * 10^x_dim(problem),
        x_proposal = x_prior(problem),
        y_kernel = BOSS.GaussianKernel(),
        p_kernel = BOSS.GaussianKernel(),
    )
    # acquisition = EIV(
    #     y_samples = 20,
    #     x_samples = 2 * 10^x_dim(problem),
    #     x_proposal = x_prior(problem),
    # )
    # acquisition = IMIQR(;
    #     p_u = 0.75,
    #     x_samples = 2 * 10^x_dim(problem),
    #     x_proposal = x_prior(problem),
    # )

    
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

    @info "Initial data:"
    for (x, y) in zip(eachcol(data.X), eachcol(data.Y))
        println("  $x -> $y")
    end

    
    ### BOSIP PROBLEM ###
    bosip = construct_bosip_problem(;
        problem,
        data,
        acquisition,
        model,
    )

    return main(problem, bosip; kwargs...)
end

# for continuing an experiment (mainly for debugging)
function main(problem::AbstractProblem, bosip::BosipProblem; run_name="_test", save_data=false, metric=false, plots=false, run_idx=nothing)
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
    term_cond = IterLimit(100) # TODO


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

    metric_cb = MetricCallback(;
        reference = reference(problem),
        logpost_estimator = log_posterior_estimate(),
        sampler,
        sample_count = 2 * 10^x_dim(problem),
        # metric = MMDMetric(;
        #     kernel = with_lengthscale(GaussianKernel(), (bounds[2] .- bounds[1]) ./ 3),
        # ),
        # metric = OptMMDMetric(;
        #     kernel = GaussianKernel(),
        #     bounds,
        #     algorithm = BOBYQA(),
        #     rhoend = 1e-4,
        # ),
        metric = TVMetric(;
            grid = xs,
            ws = ws,
        ),
    )
    # first callback in `callbacks` (this is important for `SaveCallback`)
    callbacks = BosipCallback[]
    metric && push!(callbacks, metric_cb)


    ### PLOTS ###
    plot_cb = PlotModule.PlotCB(;
        problem,
        sampler,
        sample_count = 2 * 10^x_dim(problem),
        plot_each = 10,
        save_plots = true,
    )
    plots && push!(callbacks, plot_cb)

    
    ### STORING RESULTS ###
    data_cb = SaveCallback(;
        dir = data_dir(problem),
        filename = base_filename(problem, run_name, run_idx),
    )
    save_data && push!(callbacks, data_cb)

    options = BosipOptions(;
        callback = CombinedCallback(callbacks...),
    )

    
    ### RUN ###
    bosip!(bosip; model_fitter, acq_maximizer, term_cond, options)
    return bosip
end

function run(problem::AbstractProblem)
    run_name = "loglike" # the name used for storing data from this run

    start_files = Glob.glob(starts_dir(problem) * "/start_*.jld2")
    @info "Running $(length(start_files)) runs of the $(typeof(problem)) ..."
    
    for start_file in start_files
        m = match(r"start_(\d+)\.jld2$", start_file)
        run_idx = parse(Int, m.captures[1])
        data = load(start_file, "data")
        main(; run_name, save_data=true, data, run_idx)
    end

    nothing
end
