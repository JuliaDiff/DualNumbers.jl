using DualNumbers
import DualNumbers: Dual4, real, epsilon1, epsilon2, epsilon3, epsilon4
import NaNMath
using Base.Test


x = DualNumbers.Dual4(2, 1, 2, 3, 4)
y = x^3

@test_approx_eq real(y) 2.0^3
@test_approx_eq epsilon1(y) 3.0*2^2
@test_approx_eq epsilon2(y) 2*3.0*2^2
@test_approx_eq epsilon3(y) 3*3.0*2^2
@test_approx_eq epsilon4(y) 4*3.0*2^2

y = x^3.0

@test_approx_eq real(y) 2.0^3
@test_approx_eq epsilon1(y) 3.0*2^2
@test_approx_eq epsilon2(y) 2*3.0*2^2
@test_approx_eq epsilon3(y) 3*3.0*2^2
@test_approx_eq epsilon4(y) 4*3.0*2^2

y = NaNMath.pow(x,3)

@test_approx_eq real(y) 2.0^3
@test_approx_eq epsilon1(y) 3.0*2^2
@test_approx_eq epsilon2(y) 2*3.0*2^2
@test_approx_eq epsilon3(y) 3*3.0*2^2
@test_approx_eq epsilon4(y) 4*3.0*2^2

y = sin(x)+exp(x)
@test_approx_eq real(y) sin(2)+exp(2)
@test_approx_eq epsilon1(y) cos(2)+exp(2)
@test_approx_eq epsilon2(y) 2*(cos(2)+exp(2))
@test_approx_eq epsilon3(y) 3*(cos(2)+exp(2))
@test_approx_eq epsilon4(y) 4*(cos(2)+exp(2))

y = x/2
@test_approx_eq real(y) 1.0
@test_approx_eq epsilon1(y) 0.5
@test_approx_eq epsilon2(y) 1.0
@test_approx_eq epsilon3(y) 3/2
@test_approx_eq epsilon4(y) 2.0

@test x - 2 == Dual4(0,1,2,3,4)


@test x > 1
y = abs(-x)
@test_approx_eq real(y) 2.0
@test_approx_eq epsilon1(y) 1.0
@test_approx_eq epsilon2(y) 2.0
@test_approx_eq epsilon3(y) 3.0
@test_approx_eq epsilon4(y) 4.0

@test isequal(1.0,Dual4(1.0))

y = 1/x
@test_approx_eq real(y) 1/2
@test_approx_eq epsilon1(y) -1/2^2
@test_approx_eq epsilon2(y) -2/2^2
@test_approx_eq epsilon3(y) -3/2^2
@test_approx_eq epsilon4(y) -4/2^2

function squareroot(x)
    it = x
    while abs(it*it - x) > 1e-13
        it = it - (it*it-x)/(2it)
    end
    return it
end

sq = squareroot(Dual4(10000.0,1.0,2.0,3.0,4.0))
@test_approx_eq epsilon1(sq) 0.005
@test_approx_eq epsilon2(sq) 2*0.005
@test_approx_eq epsilon3(sq) 3*0.005
@test_approx_eq epsilon4(sq) 4*0.005

@test_approx_eq epsilon1(exp(1)^Dual4(1.0,1.0,2.0,3.0,4.0)) exp(1)
@test_approx_eq epsilon2(exp(1)^Dual4(1.0,1.0,2.0,3.0,4.0)) 2exp(1)
@test_approx_eq epsilon3(exp(1)^Dual4(1.0,1.0,2.0,3.0,4.0)) 3exp(1)
@test_approx_eq epsilon4(exp(1)^Dual4(1.0,1.0,2.0,3.0,4.0)) 4exp(1)
@test_approx_eq epsilon2(NaNMath.sin(Dual4(1.0,1.0,2.0,3.0,4.0))) 2cos(1)

@test_approx_eq epsilon1(NaNMath.pow(exp(1),Dual4(1.0,1.0,2.0,3.0,4.0))) exp(1)
@test_approx_eq epsilon2(NaNMath.pow(exp(1),Dual4(1.0,1.0,2.0,3.0,4.0))) 2exp(1)
@test_approx_eq epsilon3(NaNMath.pow(exp(1),Dual4(1.0,1.0,2.0,3.0,4.0))) 3exp(1)
@test_approx_eq epsilon4(NaNMath.pow(exp(1),Dual4(1.0,1.0,2.0,3.0,4.0))) 4exp(1)
