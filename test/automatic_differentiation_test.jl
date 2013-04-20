## Calculation of f(2) and f'(2) of function f:R->R given by f(x)=x^3

using DualNumbers

x = dual(2, 1)
f(x) = x^3
y = f(x)

println("f(x) = x^3")
println("f(2) = ", real(y))
println("f'(2) = ", imag(y))
