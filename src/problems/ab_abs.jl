
"""
    AbsABProblem()

The analytical toy problem of inferring the parameters `a`, `b`
given the observation `z_obs = [1.]`.

The blackbox simulator realizes the function `y = a * b`.
In contrast to the `ABProblem`, here, only the absolute value `|y|` is modeled.
See also the `LogABProblem` for a version of the problem,
where only the log-likelihood is returned by the simulator.

The likelihood is Gaussian.
"""
struct AbsABProblem <: AbstractProblem end


module AbsABProblemModule

import ..AbsABProblem

import ..simulator
import ..domain
import ..y_max
import ..likelihood
import ..prior_mean
import ..x_prior
import ..est_amplitude
import ..est_noise_std
import ..true_f
import ..reference_samples

using BOSS
using BOSIP
using Distributions


# --- API ---

simulator(::AbsABProblem) = ab_simulation

domain(::AbsABProblem) = Domain(;
    bounds = _get_bounds(),
)

# TODO model
# likelihood(::AbsABProblem) = NormalLikelihood(; z_obs, std_obs)
likelihood(::AbsABProblem) = NormalDiffLikelihood(; std_obs)

# TODO model
# prior_mean(::AbsABProblem) = z_obs
prior_mean(::AbsABProblem) = zero(z_obs)

x_prior(::AbsABProblem) = _get_trunc_x_prior()

est_amplitude(::AbsABProblem) = [20.]

# TODO noise
est_noise_std(::AbsABProblem) = nothing

true_f(::AbsABProblem) = x -> ab_simulation(x; noise_std=zero(std_sim))


# --- UTILS ---

const z_obs = [1.]
const std_obs = [0.2]

# TODO noise
# (not using noise in order to compare with loglike modeling more fairly)
const std_sim = [0.]
# const std_sim = [0.1]

# the true blackbox function
f_(x) = [x[1] * x[2]]

function ab_simulation(x; noise_std=std_sim)
    y = f_(x)
    y .+= rand(Normal(0., noise_std[1]))
    
    # TODO model
    # return y
    return abs.(y - z_obs)
end

_get_bounds() = ([-5., -5.], [5., 5.])

function _get_trunc_x_prior()
    prior = _get_x_prior()
    bounds = _get_bounds()
    return truncated(prior; lower=bounds[1], upper=bounds[2])
end
_get_x_prior() = Product(fill(Normal(0., 5/3), 2))

end # module AbsABProblemModule
