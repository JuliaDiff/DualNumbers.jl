using DualNumbers
using Base.Test

x = dual(2, 1)
y = x^3

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
