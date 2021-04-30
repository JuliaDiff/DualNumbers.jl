using DualNumbers
using Test

@test isempty(detect_ambiguities(DualNumbers, Base, Core))

@test checkindex(Bool, 1:3, dual(2))

# wrap in individual modules to avoid name conflicts.
module TestAutomaticDifferentiation
include("automatic_differentiation_test.jl")
end

@testset "complex" begin
    for D in [Dual32, Dual64, Dual128]
        @test typeof(complex(zero(D))) == complex(D)
    end
end
