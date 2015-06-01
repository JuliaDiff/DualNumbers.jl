immutable DualN{N,T<:Real} <: Number
    re::T
    dus::NTuple{N,T}
end

DualN{T<:Real}(re::T) = DualN{0,T}(re, tuple())
DualN{T<:Real}(re::T, dus::T...) = DualN(re, dus)
DualN(x::Real, dus::Real...) = DualN(promote(x, dus...)...)

real(z::DualN) = z.re
epsilon(z::DualN) = z.dus
epsilon(z::DualN, i) = z.dus[i]

neps{N}(::DualN{N}) = N
neps{N,T}(::Type{DualN{N,T}}) = N

eps(z::DualN) = eps(real(z))
eps{N,T}(::Type{DualN{N,T}}) = eps(T)

# Type stable method for generating
# zero-filled  NTuples of a given
# length and type.
@generated function zero_tup{N,T}(::Type{NTuple{N,T}})
    z = zero(T)
    ex = "tuple(" * repeat("zero(T),", N) * ")"
    return parse(ex)
end

zero(z::DualN) = DualN(zero(real(z)), zero_tup(typeof(epsilon(z))))
zero{N,T}(::Type{DualN{N,T}}) = DualN(zero(T), zero_tup(NTuple{N,T}))
one(z::DualN) = DualN(one(real(z)), zero_tup(typeof(epsilon(z))))
one{N,T}(::Type{DualN{N,T}}) = DualN(one(T), zero_tup(NTuple{N,T}))

inf(z::DualN) = DualN(inf(real(z)))
nan(z::DualN) = DualN(nan(real(z)))

isnan(z::DualN) = isnan(real(z))
isdual(::DualN) = true

isreal(z::DualN) = any(x -> x == 0, epsilon(z))
isreal(z::DualN{0}) = true
isfinite(z::DualN) = isfinite(real(z))

convert{N,T<:Real}(::Type{DualN{N,T}}, x::Real) = DualN(convert(T, x))
convert{N,T<:Real}(::Type{DualN{N,T}}, z::DualN{N,T}) = z
convert{N,T<:Real}(::Type{DualN{N,T}}, z::DualN{N}) = DualN{N,T}(real(z), epsilon(z))
convert{T<:Real}(::Type{T}, z::DualN{0}) = convert(T, real(z))
convert{T<:Real}(::Type{T}, z::DualN) = isreal(z) ? convert(T, real(z)) : throw(InexactError())

promote_rule{N, A<:Real, B<:Real}(::Type{DualN{N,A}}, ::Type{DualN{N,B}}) = DualN{N,promote_type(A, B)}
promote_rule{N, T<:Real}(::Type{DualN{N,T}}, ::Type{T}) = DualN{N,T}
promote_rule{N, A<:Real, B<:Real}(::Type{DualN{N,A}}, ::Type{B}) = DualN{N,promote_type(A, B)}

#######################################
## Generic functions of dual numbers ##

convert(::Type{DualN}, z::Dual) = z
convert(::Type{DualN}, x::Real) = DualN(x)

==(z::DualN, w::DualN) = real(z) == real(w) && epsilon(z) == epsilon(w)
==(z::DualN, x::Real) = real(z) == x
==(x::Real, z::DualN) = z == x

isequal(z::DualN, w::DualN) = isequal(real(z), real(w)) && isequal(epsilon(z), epsilon(w))
isequal(z::DualN, x::Real) = isreal(z) && isequal(real(z), x)
isequal(x::Real, z::DualN) = isequal(z, x)

isless(z::DualN, w::DualN) = real(z) < real(w)
isless(z::Real, w::DualN) = z < real(w)
isless(z::DualN, w::Real) = real(z) < w

hash(z::DualN) = isreal(z) ? hash(real(z)) : hash(real(z), hash(epsilon(z)))

# we don't support Dual{Complex}, so conj is a noop
conj(z::DualN) = z

# Why is abs defined in this manner?
abs(z::DualN)  = (real(z) >= 0) ? z : -z
abs2(z::DualN) = z*z

# How should these be defined?
#conjdual(z::Dual) = Dual(real(z),-epsilon(z))
#absdual(z::Dual) = abs(real(z))
#abs2dual(z::Dual) = abs2(real(z))

scale_tuple(x::Number, tup::Tuple) = map(i -> i * x, tup)

# What about the mixed N case? See comment above near the convert method definitions
+{N}(z::DualN{N}, w::DualN{N}) = DualN(real(z)+real(w), map(+, epsilon(z), epsilon(w)))
+(x::Real, z::DualN) = DualN(z+real(z), epsilon(z))
+(z::DualN, x::Real) = w+z

-(z::DualN) = DualN(-real(z), map(-, epsilon(z)))
-{N}(z::DualN{N}, w::DualN{N}) = DualN(real(z)-real(w), map(-, epsilon(z), epsilon(w)))
-(x::Real, z::DualN) = DualN(x-real(z), map(-, epsilon(z)))
-(z::DualN, x::Real) = DualN(real(z)-x, epsilon(z))

# avoid ambiguous definition with Bool*Number
# Should we use ifelse rather than standard ternary operations?
*(x::Bool, z::DualN) = ifelse(x, z, ifelse(signbit(real(z))==0, zero(z), -zero(z)))
*(z::DualN, x::Bool) = x*z

function *{N}(z::DualN{N}, w::DualN{N})
    z_r, w_r = real(z), real(w)
    dus =  map((zdu, wdu)->zdu*w_r+z_r*wdu, epsilon(z), epsilon(w))
    return DualN(real(z)*real(w), dus)
end

*(x::Real, z::DualN) = DualN(x*real(z), scale_tuple(x, epsilon(z)))
*(z::DualN, x::Real) = x*z

function /(x::Real, z::DualN)
    df = -x/real(z)^2
    return DualN(x/real(z), scale_tuple(df, epsilon(z)))
end

/(z::DualN, x::Real) = DualN(real(z)/x, scale_tuple(1/x, epsilon(z)))

function /{N}(z::DualN{N}, w::DualN{N})
    z_r, w_r = real(z), real(w)    
    denom = 1/w_r^2
    dus = map((zdu, wdu) -> (zdu*w_r - z_r*wdu) * denom, epsilon(z), epsilon(w))
    return DualN(real(z)/real(w), dus)
end

for f in (:^, :(NaNMath.pow))
    @eval function ($f)(z::DualN, w::DualN)
        re = $f(real(z),real(w))
        powval = real(w)*(($f)(real(z),real(w)-1))
        logval = ($f)(real(z),real(w))*log(real(z))
        dus = map((zdu, wdu) -> zdu*powval + logval*wdu, epsilon(z), epsilon(w))
        return DualN(re, dus)
    end
end

# generate redundant definitions to resolve ambiguity warnings
for T in (:Integer, :Rational, :Real)
    @eval function ^(z::DualN, n::($T))
        zn1 = n*real(z)^(n-1)
        return DualN(real(z)^n, scale_tuple(zn1, epsilon(z)))
    end
end

function NaNMath.pow(z::DualN, x::Real)
    powval = x*NaNMath.pow(real(z),x-1)
    return DualN(NaNMath.pow(real(z),x), scale_tuple(powval, epsilon(z)))
end

function NaNMath.pow(x::Real, z::DualN)
    logval = NaNMath.pow(x,real(z))*log(x) 
    return DualN(NaNMath.pow(x,real(z)), scale_tuple(logval, epsilon(z)))
end

for (funsym, exp) in Calculus.symbolic_derivatives_1arg()
    funsym == :exp && continue
    
    @eval function $(funsym)(z::DualN)
        x = real(z)
        df = $exp
        return DualN($(funsym)(x), scale_tuple(df, epsilon(z)))
    end
    
    # extend corresponding NaNMath methods
    if funsym in (:sin, :cos, :tan, 
                  :asin, :acos, :acosh, 
                  :atanh, :log, :log2, 
                  :log10, :lgamma, :log1p)

        funsym = Expr(:.,:NaNMath,Base.Meta.quot(funsym))
        
        @eval function $(funsym)(z::DualN)
            x = real(z)
            df = $(to_nanmath(exp))
            return DualN($(funsym)(x), scale_tuple(df, epsilon(z)))
        end

    end
end

# only need to compute exp once
function exp(z::DualN)
    df = exp(real(z))
    return DualN(df, scale_tuple(df, epsilon(z)))
end