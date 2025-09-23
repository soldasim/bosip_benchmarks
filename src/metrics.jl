
metric_fname(::Type{OptMMDMetric}) = "metric"  # I am too lazy to rename the files.
metric_fname(::Type{TVMetric}) = "TVmetric"

function get_metric(::Type{MMDMetric}, problem::AbstractProblem)
    # TODO
    bounds = domain(problem).bounds
    λs = (bounds[2] .- bounds[1]) ./ 3

    return MMDMetric(;
        kernel = with_lengthscale(GaussianKernel(), λs),
    )
end
function get_metric(::Type{OptMMDMetric}, problem::AbstractProblem)
    return OptMMDMetric(;
        kernel = GaussianKernel(),
        bounds,
        algorithm = BOBYQA(),
    )
end
function get_metric(::Type{TVMetric}, problem::AbstractProblem)
    # TODO
    xs = rand(x_prior(problem), 20 * 10^x_dim(problem))
    ws = exp.( (0.) .- logpdf.(Ref(x_prior(problem)), eachcol(xs)) )

    return TVMetric(;
        grid = xs,
        ws = ws,
    )
end
