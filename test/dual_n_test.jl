using DualNumbers
import DualNumbers: DualN, real, epsilon
using Base.Test

function test_dus(z::DualN, f::Function)
    for i=1:neps(z)
        @test_approx_eq epsilon(z,i) f(i)
    end  
end

@test isequal(1.0, DualN(1.0))

x = DualN(2, 1, 2, 3, 4)
@test x - 2 == DualN(0,1,2,3,4)
@test x > 1

y = x^3
@test_approx_eq real(y) 2.0^3
test_dus(y, i -> i*3.0*2^2)

y = x^3.0
@test_approx_eq real(y) 2.0^3
test_dus(y, i -> i*3.0*2^2)

y = NaNMath.pow(x,3)
@test_approx_eq real(y) 2.0^3
test_dus(y, i -> i*3.0*2^2)

y = sin(x)+exp(x)
@test_approx_eq real(y) sin(2)+exp(2)
test_dus(y, i -> i*(cos(2)+exp(2)))

y = x/2
@test_approx_eq real(y) 1.0
test_dus(y, i -> i/2)

y = abs(-x)
@test_approx_eq real(y) 2.0
test_dus(y, i -> i)

y = 1/x
@test_approx_eq real(y) 1/2
test_dus(y, i -> -i/2^2)

function squareroot(x)
    it = x
    while abs(it*it - x) > 1e-13
        it = it - (it*it-x)/(2it)
    end
    return it
end

test_dus(squareroot(DualN(10000.0,1.0,2.0,3.0,4.0)), i -> i*0.005)
test_dus(exp(1)^DualN(1.0,1.0,2.0,3.0,4.0), i -> i*exp(1))

@test_approx_eq epsilon(NaNMath.sin(DualN(1.0,1.0,2.0,3.0,4.0)), 2) 2*cos(1)

test_dus(NaNMath.pow(exp(1),DualN(1.0,1.0,2.0,3.0,4.0)), i -> i*exp(1))