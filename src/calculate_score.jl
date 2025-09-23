include("main.jl")

function calculate_scores(
    problem::AbstractProblem,
    metricT::Type{<:DistributionMetric},
    run_name::String,
    run_idx::Int;
    save_score = false,
)
    dir = data_dir(problem)
    metric = get_metric(metricT, problem)

    # TODO
    sample_count = 2 * 10^x_dim(problem)

    # TODO
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

    # enclose everything in a MetricCallback (to save it compactly later)
    cb = MetricCallback(;
        reference = reference(problem),
        logpost_estimator = log_posterior_estimate(problem),
        sampler,
        sample_count,
        metric,
    )

    # load the stored `BosipProblem`s at each iteration of the run
    iters_file = joinpath(dir, "$(run_name)_$(run_idx)_iters.jld2")
    data = load(iters_file)
    @assert haskey(data, "problems")
    bosip_states = data["problems"]

    # calculate the metric scores
    scores = zeros(length(bosip_states))
    @showprogress "Calculating scores..." for (i, p) in enumerate(bosip_states)
        scores[i] = calculate_score(cb, p)
    end

    # save the scores
    if save_score
        score_file = joinpath(dir, "$(run_name)_$(run_idx)_$(metric_fname(metricT)).jld2")
        @save score_file score=scores metric=cb
    end

    return scores
end

function calculate_score(metric::DistributionMetric, problem::AbstractProblem, sampler::DistributionSampler, p::BosipProblem)
    ref = reference(problem)
    logpost_est = log_posterior_estimate(problem)
    return calculate_score(metric.metric, ref, logpost_est, sampler, p)
end
function calculate_score(cb::MetricCallback, p::BosipProblem)
    return calculate_score(cb.metric, cb.reference, cb.logpost_estimator, cb.sampler, p)
end

function calculate_score(metric::SampleMetric, ref, logpost_est, sampler::DistributionSampler, p::BosipProblem)
    if ref isa Function
        true_samples = sample_posterior_pure(sampler, ref, p.problem.domain, sample_count)
    else
        true_samples = ref
    end

    approx_samples = sample_posterior_pure(sampler, logpost_est, p.problem.domain, sample_count)

    score = calculate_metric(metric, true_samples, approx_samples)
    return score
end
function calculate_score(metric::PDFMetric, ref, logpost_est, sampler::DistributionSampler, p::BosipProblem)
    ### retrieve the true and approx logpdf
    @assert ref isa Function
    true_logpdf = ref
    approx_logpdf = logpost_est(p)

    ### calculate metric
    score = calculate_metric(metric, true_logpdf, approx_logpdf)
    return score
end
