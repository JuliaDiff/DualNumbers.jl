using DualNumbers, Base.Test

if VERSION >= v"0.5.0-dev+5429"
    @test checkindex(Bool, 1:3, dual(2))
end

# wrap in individual modules to avoid name conflicts.
module TestAutomaticDifferentiation
include("automatic_differentiation_test.jl")
end
