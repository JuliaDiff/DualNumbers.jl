using DualNumbers, Base.Test
import DualNumbers: value
import NaNMath
using Compat

x = Dual(2, 1)
y = x^3

@test_approx_eq value(y) 2.0^3
@test_approx_eq epsilon(y) 3.0*2^2

y = x^3.0

@test_approx_eq value(y) 2.0^3
@test_approx_eq epsilon(y) 3.0*2^2

y = sin(x)+exp(x)
@test_approx_eq value(y) sin(2)+exp(2)
@test_approx_eq epsilon(y) cos(2)+exp(2)

@test x > 1
y = abs(-x)
@test_approx_eq value(y) 2.0
@test_approx_eq epsilon(y) 1.0

@test isequal(1.0,Dual(1.0))

y = 1/x
@test_approx_eq value(y) 1/2
@test_approx_eq epsilon(y) -1/2^2

Q = [1.0 0.1; 0.1 1.0]
x = @compat dual.([1.0,2.0])
try 
    x[1] = Dual(1.0,1.0) 
catch e
    @test isa(e, MethodError) # This assignment will fail now because of the extra type parameter
end

y = (1/2)*dot(x,Q*x)
@test_approx_eq value(y) 2.7 
@test isempty(epsilon(y)) # Since there are no partials, this will return empty container

function squareroot(x)
    it = x
    while abs(it*it - x) > 1e-13
        it = (it+x/it)/2
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
#@test eps(Dual{Float64}) == eps(Float64)
@test one(x) == Dual(1.0,0.0)
#@test one(Dual{Float64}) == Dual(1.0,0.0)
#@test convert(Dual{Float64}, Inf) == convert(Float64, Inf)
#@test isnan(convert(Dual{Float64}, NaN))

#@test convert(Dual{Float64},Dual(1,2)) == Dual(1.0,2.0)
#@test convert(Float64, Dual(10.0,0.0)) == 10.0
#@test convert(Dual{Int}, Dual(10.0,0.0)) == Dual(10,0)

x = Dual(1.2,1.0)
@test floor(Int, x) == 1
@test ceil(Int, x)  == 2
@test trunc(Int, x) == 1
@test round(Int, x) == 1

z = Dual(1,1) + im*Dual(0,1)
f = exp(z)
@test abs(f) == exp(abs(z))

#g = sinpi(z)
#@test abs(g) == sinpi(abs(z))

#h = z^4
#@test abs(h) == (abs(z))^4

a = abs2(z)
@test abs(a) == abs2(abs(z))

#l = log(z)
#@test abs(l) == log(abs(z))

s = sign(z)
@test abs(s) == sign(abs(z))

a = angle(z)
@test abs(a) == angle(abs(z))


# Tests limit definition. Let z = a + b ɛ, where a and b ∈ C.
#
# The dual of |z| is lim_{h→0} (|a + bɛh| - |a|)/h
#
# and it depends on the direction (i.e. the complex value of epsilon(z)).
#

z = Dual(1,1) + im * Dual(1,0)
@test abs(z) ≡ sqrt(2) + 1/sqrt(2)*ɛ
z = (1 + 1im) * Dual(1,0) + cis(π/4) * Dual(0,1)
@test abs(z) == sqrt(2) + 2/sqrt(2)^2*ɛ
z = (1.0 + 1.0im) * Dual(1,0) + cis(π/2) * Dual(0,1)
@test abs(z) == sqrt(2) + 1/sqrt(2)*ɛ

# tests vectorized methods
const zv = @compat dual.(collect(1.0:10.0), ones(10))

f = @compat exp.(zv)
@test all(@compat value.(f) .== exp.(value.(zv)))
@test all(@compat epsilon.(f) .== epsilon.(zv) .* exp.(value.(zv)))

# tests norms and inequalities
@test norm(f,Inf) ≤ norm(f) ≤ norm(f,1)

# tests for constant ɛ
@test epsilon(1.0 + ɛ) == [1.0] # Returns container
@test epsilon(1.0 + 0.0ɛ) == [0.0] # Returns container 
test(x, y) = x^2 + y
@test test(1.0 + ɛ, 1.0) == 2.0 + 2.0ɛ
@test test(1.0, 1.0 + ɛ) == 2.0 + 1.0ɛ

@test ɛ*im == Dual(false,false) + Complex(false,true)*Dual(false,true)

@test value(mod(Dual(15.23, 1), 10)) == 5.23
@test epsilon(mod(Dual(15.23, 1), 10)) == [1] # Returns container

@test epsilon(Dual(-2.0,1.0)^2.0) == [-4] # Returns container
@test epsilon(Dual(-2.0,1.0)^Dual(2.0,0.0)) == [-4] # Returns container
