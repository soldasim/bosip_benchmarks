"""
    ProxySIRProblem()

The SIR (Susceptible-Infected-Recovered) epidemic model for simulation-based inference.

The SIR model is a system of differential equations:
    dS/dt = -β * S * I / N
    dI/dt = β * S * I / N - γ * I  
    dR/dt = γ * I

where:
- S(t): number of susceptible individuals at time t
- I(t): number of infected individuals at time t
- R(t): number of recovered individuals at time t
- β: contact rate (transmission parameter)
- γ: recovery rate
- N: total population (constant)

The parameters to infer are [β, γ] where:
- β: contact rate (how fast the disease spreads)
- γ: recovery rate (how fast infected individuals recover)

The observations are binomial samples with probability p = I(t)/N at multiple time points,
subsampled every 17 days. The model predicts the probability p, and observations are 
sampled from Binomial(trials=1000, p) at each time point.
"""
struct ProxySIRProblem <: AbstractProblem end

module ProxySIRModule

import ..ProxySIRProblem

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
using DifferentialEquations

### Use the analytical reference posterior from the original SIRProblem
import ..true_logpost
import ..SIRProblem
true_logpost(::ProxySIRProblem) = true_logpost(SIRProblem())


# --- API ---

simulator(::ProxySIRProblem) = model_target

domain(::ProxySIRProblem) = Domain(;
    bounds = _get_bounds(),
)

# TODO proxy
function log_ψ(δ::AbstractVector{<:Real}, x::AbstractVector{<:Real})
    # transform δ back to probabilities
    y = (exp.(δ) .- 1) ./ trials

    y = clamp.(y, 0.0, 1.0)
    return mapreduce((t, p, z) -> logpdf(Binomial(t, p), z), +, trials, y, z_obs)
end

# TODO proxy
# likelihood(::ProxySIRProblem) = BinomialLikelihood(;
#     z_obs = Int64.(z_obs),
#     trials,
#     int_grid_size = 200,
# )
likelihood(::ProxySIRProblem) = CustomLikelihood(;
    log_ψ,
    mc_samples = 1000, # TODO
)

prior_mean(::ProxySIRProblem) = _get_prior_mean()

x_prior(::ProxySIRProblem) = _get_trunc_x_prior()

est_amplitude(::ProxySIRProblem) = _get_est_amplitude()

# TODO noise
est_noise_std(::ProxySIRProblem) = nothing

true_f(::ProxySIRProblem) = model_target


# --- UTILS ---

# Population parameters
const N = 1_000_000.0   # total population
const I0 = 1.0          # initial infected
const R0 = 0.0          # initial recovered
const S0 = N - I0 - R0  # initial susceptible

# Time parameters
const t_span = (0.0, 160.0)  # simulation time span (160 days)
const save_freq = 1.0        # save every day
const subsample_freq = 17    # observe every 17th day

# Observation parameters
const n_obs = 10        # number of time observations (every 17 days from day 0 to day 153)

# Binomial likelihood parameters
const trials = fill(1000, n_obs)  # number of trials for each observation

# Reference parameters: [β, γ]
const x_ref = [0.6148, 0.1917]

"""
    sir_ode!(du, u, p, t)

SIR model ODE system.
u[1] = S (susceptible), u[2] = I (infected), u[3] = R (recovered)
p = [β, γ]: parameters (contact rate, recovery rate)
"""
function sir_ode!(du, u, p, t)
    β, γ = p
    S, I, R = u
    
    # SIR equations
    infection_rate = β * S * I / N
    recovery_rate = γ * I
    
    du[1] = -infection_rate           # dS/dt
    du[2] = infection_rate - recovery_rate  # dI/dt  
    du[3] = recovery_rate             # dR/dt
end

"""
    model_target(x)

The SIR problem simulator together with the mapping to the model target variable.
Returns probabilities p = I/N at subsampled time points.
"""
function model_target(x)
    return sir_simulation(x)
end

"""
    sir_simulation(x)

Simulate the SIR epidemic model with parameters x = [β, γ].
Returns probabilities p = I/N at subsampled time points for binomial likelihood.
"""
function sir_simulation(x)
    β, γ = x
    
    # Initial conditions
    u0 = [S0, I0, R0]
    
    # Set up and solve ODE
    prob = ODEProblem(sir_ode!, u0, t_span, [β, γ])
    sol = solve(prob, Tsit5(), saveat=save_freq)
    
    # Extract infected population I(t) and normalize to get probabilities
    I_series = [sol.u[i][2] for i in 1:length(sol.t)]
    p_series = I_series ./ N
    
    # Subsample every 17 days (indices 1, 18, 35, ..., up to 10 observations)
    indices = 1:subsample_freq:(subsample_freq*n_obs-16)
    p = p_series[indices[1:n_obs]]
    
    # Clamp to [0, 1] to ensure valid probabilities
    p = clamp.(p, 0.0, 1.0)
    
    # TODO proxy
    # return p
    return log.(1 .+ (p .* trials))
end

"""
Generate reference observation data with known parameters.
Returns binomial samples based on the true probabilities.
"""
function _generate_reference_data()
    # Get true probabilities from SIR simulation
    p_true = sir_simulation(x_ref)
    
    # Sample from binomial distribution for each time point
    z_samples = Int64[]
    for i in 1:n_obs
        # Sample from Binomial(trials[i], p_true[i])
        z_sample = rand(Binomial(trials[i], p_true[i]))
        push!(z_samples, z_sample)
    end
    
    return z_samples
end

# Generate synthetic observation data
# const z_obs = _generate_reference_data()
const z_obs = [0,1,311,44,0,1,0,0,0,0]  # from sbibm

"""
Parameter bounds: [β_min, γ_min], [β_max, γ_max]
"""
_get_bounds() = ([0.0, 0.0], [2.0, 0.5])

"""
Prior mean based on binomial probabilities.
"""
# TODO proxy
# _get_prior_mean() = z_obs ./ trials
_get_prior_mean() = log.(1 .+ z_obs)

"""
Estimated amplitude for each observation dimension.
For binomial likelihood, this represents the scale of probability values.
"""
# TODO proxy
# _get_est_amplitude() = fill(1.0, n_obs)
_get_est_amplitude() = log.(1 .+ trials)

"""
Truncated prior distribution for parameters.
"""
function _get_trunc_x_prior()
    prior = _get_x_prior()
    bounds = _get_bounds()
    return truncated(prior; lower=bounds[1], upper=bounds[2])
end

"""
Prior distribution matching Python implementation:
- β ~ LogNormal(log(0.4), 0.5): contact rate
- γ ~ LogNormal(log(0.125), 0.2): recovery rate
"""
function _get_x_prior()
    return Product([
        LogNormal(log(0.4), 0.5),    # β: contact rate
        LogNormal(log(0.125), 0.2),  # γ: recovery rate
    ])
end

end # module ProxySIRModule
