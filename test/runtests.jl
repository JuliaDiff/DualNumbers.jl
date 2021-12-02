using DualNumbers
using Test

@test_skip isempty(detect_ambiguities(DualNumbers, Base, Core))

@test checkindex(Bool, 1:3, dual(2))

# wrap in individual modules to avoid name conflicts.
module TestAutomaticDifferentiation
include("automatic_differentiation_test.jl")
end
