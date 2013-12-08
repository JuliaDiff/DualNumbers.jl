using DualNumbers
using Base.Test

x = dual(2, 1)
y = x^3

@test_approx_eq real(y) 2.0^3
@test_approx_eq imag(y) 3.0*2^2

y = sin(x)+exp(x)
@test_approx_eq real(y) sin(2)+exp(2)
@test_approx_eq imag(y) cos(2)+exp(2)
