include("main.jl")

# Compute the MMD score of two sample sets, both from the true (reference) posterior,
# multiple times, to determine the sensitivity of the "maximum MMD" metric.
function sample_ref_mmd(problem::AbstractProblem, mmd_samples::Int, x_samples::Int;
    save_data = false,
    parallel = false, # PRIMA.jl causes StackOverflow when parallelized on Linux,
)
    mmd_vals = sample_ref_mmd(reference(problem), domain(problem), mmd_samples, x_samples; parallel)

    if save_data
        dir = data_dir(problem) * "/opt_mmd"
        file = dir * "/mmd_vals.jld2"
        mkpath(dir)
        @save file mmd_vals=mmd_vals
    end
    return mmd_vals
end

function sample_ref_mmd(ref::Function, domain::Domain, mmd_samples::Int, x_samples::Int;
    parallel = false,
)
    ref_logpost = ref
    bounds = domain.bounds

    metric = OptMMDMetric(;
        kernel = GaussianKernel(),
        bounds = bounds,
        algorithm = BOBYQA(),
    )

    sampler = AMISSampler(;
        iters = 10,
        proposal_fitter = BOSIP.AnalyticalFitter(), # re-fit the proposal analytically
        # proposal_fitter = OptimizationFitter(;      # re-fit the proposal by MAP optimization
        #     algorithm = NEWUOA(),
        #     multistart = 6,
        #     parallel,
        #     rhoend = 1e-2,
        # ),
        # gauss_mix_options = nothing,                # use Laplace approximation for the 0th iteration
        gauss_mix_options = GaussMixOptions(;       # use Gaussian mixture for the 0th iteration
            algorithm = BOBYQA(),
            multistart = 24,
            parallel,
            cluster_ϵs = nothing,
            rel_min_weight = 1e-8,
            rhoend = 1e-4,
        ),
    )


    ### the calculation
    mmd_vals = Float64[]
    for i in 1:mmd_samples
        xs_a = sample_posterior_pure(sampler, ref_logpost, domain, x_samples)
        xs_b = sample_posterior_pure(sampler, ref_logpost, domain, x_samples)

        val = calculate_metric(metric, xs_a, xs_b)
        push!(mmd_vals, val)
    end

    return mmd_vals
end

function sample_ref_mmd(ref::AbstractMatrix{<:Real}, domain::Domain, mmd_samples::Int, x_samples::Int;
    parallel = false,
)
    ref_samples = size(ref, 2)
    @assert ref_samples >= 20 * (2 * x_samples) "Not enough reference samples for subsampling."
    bounds = domain.bounds

    metric = OptMMDMetric(;
        kernel = GaussianKernel(),
        bounds = bounds,
        algorithm = BOBYQA(),
    )

    ### the calculation
    mmd_vals = Float64[]
    for i in 1:mmd_samples
        _perm = randperm(ref_samples)
        xs_a = ref[:,_perm[1:x_samples]]
        xs_b = ref[:,_perm[x_samples+1:2*x_samples]]

        val = calculate_metric(metric, xs_a, xs_b)
        push!(mmd_vals, val)
    end

    return mmd_vals
end
