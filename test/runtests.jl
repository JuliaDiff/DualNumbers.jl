using DualNumbers
using Test

@test isempty(detect_ambiguities(DualNumbers, Base, Core))

@test checkindex(Bool, 1:3, dual(2))

# wrap in individual modules to avoid name conflicts.
module TestAutomaticDifferentiation
include("automatic_differentiation_test.jl")
end

@testset "complex" begin
    for T in [Float16, Float32, Float64, BigFloat, Int8, Int16, Int32, Int64, Int128, BigInt, Bool]
        D = Dual{T}
        @test typeof(complex(zero(D))) == complex(D)
    end
end
