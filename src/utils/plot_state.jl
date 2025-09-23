module PlotModule

using BOSIP
using CairoMakie

import ..AbstractProblem
import ..reference
import ..log_posterior_estimate

include("../data_paths.jl")

@kwdef struct PlotCB <: BosipCallback
    problem::AbstractProblem
    sampler::DistributionSampler
    sample_count::Int
    plot_each::Int = 10
    resolution::Int = 500
    save_plots::Bool = false
end

function (cb::PlotCB)(bosip::BosipProblem; term_cond, first, kwargs...)
    # first && return
    (term_cond.iter % cb.plot_each == 0) || return

    plot_state(bosip, cb.problem, cb.sampler, cb.sample_count, term_cond.iter; cb.resolution, cb.save_plots)
end

function plot_state(bosip::BosipProblem, p::AbstractProblem, sampler::DistributionSampler, sample_count::Int, iter::Int; resolution=500, save_plots=false)
    ref = reference(p)
    domain = bosip.problem.domain
    lb, ub = domain.bounds
    X = bosip.problem.data.X
    
    ### log-posteriors
    est_logpost = log_posterior_estimate(p)(bosip)
    if ref isa Function
        ref_logpost = ref
    end

    ### compute grid `Z` for contours & normalization constant `M`
    x = range(lb[1], ub[1], length=resolution)
    y = range(lb[2], ub[2], length=resolution)
    Z_est = [est_logpost([xi, yi]) for xi in x, yi in y]
    M_est = maximum(Z_est)
    Z_est .= exp.(Z_est .- M_est)
    if ref isa Function
        Z_ref = [ref_logpost([xi, yi]) for xi in x, yi in y]
        M_ref = maximum(Z_ref)
        Z_ref .= exp.(Z_ref .- M_ref)
    end

    ### posteriors
    est_post = x -> exp(est_logpost(x) - M_est)
    if ref isa Function
        ref_post = x -> exp(ref_logpost(x) - M_ref)
    end

    ### credible regions
    grid_ = hcat([[t...] for t in Iterators.product(x, y)]...)

    ws_est_ = est_post.(eachcol(grid_))
    (sum(ws_est_) == 0.) && (ws_est_ .= 1.)
    cs_est = find_cutoff.(Ref(est_post), Ref(grid_), Ref(ws_est_), [0.8, 0.95])
    if ref isa Function
        ws_ref_ = ref_post.(eachcol(grid_))
        cs_ref = find_cutoff.(Ref(ref_post), Ref(grid_), Ref(ws_ref_), [0.8, 0.95])
    end

    ### sample posterior
    # TODO
    # if ref isa Function
    #     true_samples = sample_posterior_pure(sampler, ref, domain, sample_count)
    # else
    #     true_samples = ref
    # end
    # approx_samples = sample_posterior_pure(sampler, est_logpost, domain, sample_count) 


    ### ### ### THE FIGURE ### ### ###
    fig = Figure()
    # TODO title
    title = string(log_posterior_estimate(p)) |> exp_title
    # title = "Loglikelihood Model"
    
    ax = Axis(
        fig[1, 1],
        xlabel="x₁", ylabel="x₂",
        title=title,
        xlabelsize=20, ylabelsize=20, titlesize=20,
        xticklabelsize=16, yticklabelsize=16
    )

    ### contour fill
    # contourf!(ax, x, y, Z)

    ### samples
    # TODO
    # scatter!(ax, true_samples[1, :], true_samples[2, :], color=:blue, marker=:x, markersize=4) # label="Reference Samples"
    # scatter!(ax, approx_samples[1, :], approx_samples[2, :], color=:red, marker=:x, markersize=4) # label="Approx. Samples"

    ### data
    scatter!(ax, X[1,:], X[2,:], color=:white, markersize=6)
    scatter!(ax, X[1,:], X[2,:], color=:black, markersize=4, label="Data")
    
    ### contours
    if ref isa Function
        contour!(ax, x, y, Z_ref; levels=cs_ref, color=:blue)
    end
    contour!(ax, x, y, Z_est; levels=cs_est, color=:red)

    ### legend
    # axislegend(ax)

    if save_plots
        dir = plot_dir() * "/state_plots"
        mkpath(dir)
        save(dir * "/" * string(typeof(p)) * "_$iter.png", fig)
    else
        display(fig)
    end
end

function exp_title(est_logpost_name::String)
    @assert startswith(est_logpost_name, "log_")
    est_post_name = est_logpost_name[5:end]
    return est_post_name
end

end # module PlotModule
