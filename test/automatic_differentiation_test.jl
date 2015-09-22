using DualNumbers
using Base.Test

x = dual(2, 1)
y = x^3

@test_approx_eq real(y) 2.0^3
@test_approx_eq epsilon(y) 3.0*2^2

y = x^3.0

@test_approx_eq real(y) 2.0^3
@test_approx_eq epsilon(y) 3.0*2^2

y = sin(x)+exp(x)
@test_approx_eq real(y) sin(2)+exp(2)
@test_approx_eq epsilon(y) cos(2)+exp(2)

@test x > 1
y = abs(-x)
@test_approx_eq real(y) 2.0
@test_approx_eq epsilon(y) 1.0

@test isequal(1.0,dual(1.0))

y = 1/x
@test_approx_eq real(y) 1/2
@test_approx_eq epsilon(y) -1/2^2

Q = [1.0 0.1; 0.1 1.0]
x = dual([1.0,2.0])
x[1] = dual(1.0,1.0)
y = (1/2)*dot(x,Q*x)
@test_approx_eq real(y) 2.7
@test_approx_eq epsilon(y) 1.2

function squareroot(x)
    it = x
    while abs(it*it - x) > 1e-13
        it = it - (it*it-x)/(2it)
    end
    return it
end

@test_approx_eq epsilon(squareroot(Dual(10000.0,1.0))) 0.005

@test_approx_eq epsilon(exp(1)^Dual(1.0,1.0)) exp(1)
@test_approx_eq epsilon(NaNMath.pow(exp(1),Dual(1.0,1.0))) exp(1)
@test_approx_eq epsilon(NaNMath.sin(Dual(1.0,1.0))) cos(1)

@test Dual(1.0,3) == Dual(1.0,3.0)
x = Dual(1.0,1.0)
@test eps(x) == eps(1.0)
@test eps(Dual{Float64}) == eps(Float64)
@test one(x) == Dual(1.0,0.0)
@test one(Dual{Float64}) == Dual(1.0,0.0)
@test inf(Dual{Float64}) == inf(Float64)
@test isnan(nan(Dual{Float64}))

@test convert(Dual{Float64},Dual(1,2)) == Dual(1.0,2.0)
@test convert(Float64, Dual(10.0,0.0)) == 10.0
@test convert(Dual{Int}, Dual(10.0,0.0)) == Dual(10,0)

# tests for constant du
@test epsilon(1.0 + du) == 1.0
@test epsilon(1.0 + 0.0du) == 0.0
z(x, y) = x^2 + y
@test z(1.0 + du, 1.0) == 2.0 + 2.0du
@test z(1.0, 1.0 + du) == 2.0 + 1.0du

@test real(mod(dual(15.23, 1), 10)) == 5.23
@test epsilon(mod(dual(15.23, 1), 10)) == 1
