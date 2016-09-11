#

# wrap in individual modules to avoid name conflicts.
module TestAutomaticDifferentiation
include("automatic_differentiation_test.jl")
end

module TestDual
include("dual_n.jl")
end
