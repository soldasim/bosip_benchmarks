# analytical
include("ab.jl")
include("ab_log.jl")

# SBIBM
include("sbibm/sir_det/sir.jl")
include("sbibm/sir_det/sir_log.jl")

# Jarvenpaa & Gutmann
include("jarvenpaa/simple.jl")
include("jarvenpaa/simple_log.jl")
include("jarvenpaa/banana.jl")
include("jarvenpaa/banana_log.jl")
include("jarvenpaa/bimodal.jl")
include("jarvenpaa/bimodal_log.jl")

# Physical systems
include("physical/duffing.jl")
include("physical/duffing_log.jl")
include("physical/diffusion.jl")
