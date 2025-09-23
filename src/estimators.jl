
### Defines which posterior estimator is used.
# "Log" problems use the point-estimate
# because the expected posterior is numerically unstable in their case.
function log_posterior_estimate(p::AbstractProblem)
    pname = p |> typeof |> nameof |> string
    if startswith(pname, "Log")
        return log_approx_posterior
    else
        return log_posterior_mean
    end
end
