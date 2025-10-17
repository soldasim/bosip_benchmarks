
"""
    get_init_data(::AbstractProblem, count::Int) -> ExperimentData

Generate `count` initial data points sampled from the `x_prior` of the given problem.
"""
function get_init_data(problem::AbstractProblem, count::Int)
    prior = x_prior(problem)
    @assert extrema(prior) == domain(problem).bounds
    @show domain(problem).bounds
    sim = simulator(problem)

    if count == 1
        x = mean(domain(problem).bounds)
        y = sim(x)
        return BOSS.ExperimentData(hcat(x), hcat(y))
    else
        X = rand(prior, count)
        Y = reduce(hcat, (sim(x) for x in eachcol(X)))[:,:]
        return BOSS.ExperimentData(X, Y)
    end
end
