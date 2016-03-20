# Dual number with multiple epsilon components.
# Fast inlined storage (will switch to tuples when they're fast)
# This is useful because we can avoid recomputing expensive functions and share computations across the epsilon components.
# Interface will improve, so not exported yet.

immutable Dual4{T<:Real} <: Number
    re::T
    du1::T
    du2::T
    du3::T
    du4::T
end
Dual4(x::Real, y1::Real, y2::Real, y3::Real, y4::Real) = Dual4(promote(x,y1,y2,y3,y4)...)
Dual4(x::Real) = Dual4(x, zero(x), zero(x), zero(x), zero(x))

real(z::Dual4) = z.re
epsilon1(z::Dual4) = z.du1
epsilon2(z::Dual4) = z.du2
epsilon3(z::Dual4) = z.du3
epsilon4(z::Dual4) = z.du4

eps(z::Dual4) = eps(real(z))
eps{T}(::Type{Dual4{T}}) = eps(T)
one(z::Dual4) = Dual4(one(real(z)))
one{T}(::Type{Dual4{T}}) = Dual4(one(T))
@deprecate inf{T}(::Type{Dual4{T}}) convert(Dual4{T}, Inf)
@deprecate nan{T}(::Type{Dual4{T}}) convert(Dual4{T}, NaN)
isnan(z::Dual4) = isnan(real(z))

convert{T<:Real}(::Type{Dual4{T}}, x::Real) =
  Dual4{T}(convert(T, x), zero(T), zero(T), zero(T), zero(T))
convert{T<:Real}(::Type{Dual4{T}}, z::Dual4{T}) = z
convert{T<:Real}(::Type{Dual4{T}}, z::Dual4) =
  Dual4{T}(real(z), epsilon1(z), epsilon2(z), epsilon3(z), epsilon4(z))

convert{T<:Real}(::Type{T}, z::Dual4) =
  (epsilon1(z)==0 && epsilon2(z) == 0 && epsilon3(z) == 0 && epsilon4(z) == 0 ? convert(T, real(z)) : throw(InexactError()))

promote_rule{T<:Real, S<:Real}(::Type{Dual4{T}}, ::Type{Dual4{S}}) =
    Dual4{promote_type(T, S)}
# these promotion rules shouldn't be used for scalar operations -- they're slow
promote_rule{T<:Real}(::Type{Dual4{T}}, ::Type{T}) = Dual4{T}
promote_rule{T<:Real, S<:Real}(::Type{Dual4{T}}, ::Type{S}) =
  Dual4{promote_type(T, S)}

isdual(x::Dual4) = true

real_valued{T<:Real}(z::Dual4{T}) = epsilon1(z) == 0 && epsilon2(z) == 0 && epsilon3(z) == 0 && epsilon4(z) == 0

isfinite(z::Dual4) = isfinite(real(z))

## Generic functions of dual numbers ##

convert(::Type{Dual4}, z::Dual) = z
convert(::Type{Dual4}, x::Real) = Dual(x)

==(z::Dual4, w::Dual4) = real(z) == real(w)
==(z::Dual4, x::Real) = real(z) == x
==(x::Real, z::Dual4) = real(z) == x

isequal(z::Dual4, w::Dual4) =
  isequal(real(z),real(w)) && isequal(epsilon1(z), epsilon1(w)) &&
  isequal(epsilon2(z), epsilon2(w)) && isequal(epsilon3(z),epsilon3(w)) &&
  isequal(epsilon4(z), epsilon4(w))
isequal(z::Dual4, x::Real) = real_valued(z) && isequal(real(z), x)
isequal(x::Real, z::Dual4) = real_valued(z) && isequal(real(z), x)

isless(z::Dual4,w::Dual4) = real(z) < real(w)
isless(z::Real,w::Dual4) = z < real(w)
isless(z::Dual4,w::Real) = real(z) < w

#hash(z::Dual) =
#  (x = hash(real(z)); real_valued(z) ? x : bitmix(x,hash(epsilon(z))))

# we don't support Dual{Complex}, so conj is a noop
conj(z::Dual4) = z
abs(z::Dual4)  = (real(z) >= 0) ? z : -z
abs2(z::Dual4) = z*z

# algebraic definitions
#conjdual(z::Dual) = Dual(real(z),-epsilon(z))
#absdual(z::Dual) = abs(real(z))
#abs2dual(z::Dual) = abs2(real(z))

+(z::Dual4, w::Dual4) = Dual4(real(z)+real(w), epsilon1(z)+epsilon1(w), epsilon2(z) + epsilon2(w), epsilon3(z) + epsilon3(w), epsilon4(z) + epsilon4(w))
+(z::Real, w::Dual4) = Dual4(z+real(w), epsilon1(w), epsilon2(w), epsilon3(w), epsilon4(w))
+(z::Dual4, w::Real) = w+z

-(z::Dual4) = Dual4(-real(z), -epsilon1(z), -epsilon2(z), -epsilon3(z), -epsilon4(z))
-(z::Dual4, w::Dual4) = Dual4(real(z)-real(w), epsilon1(z)-epsilon1(w), epsilon2(z)-epsilon2(w), epsilon3(z)-epsilon3(w), epsilon4(z)-epsilon4(w))
-(z::Real, w::Dual4) = Dual4(z-real(w), -epsilon1(w), -epsilon2(w), -epsilon3(w), -epsilon4(w))
-(z::Dual4, w::Real) = Dual4(real(z)-w, epsilon1(z), epsilon2(z), epsilon3(z), epsilon4(z))

# avoid ambiguous definition with Bool*Number
*(x::Bool, z::Dual4) = ifelse(x, z, ifelse(signbit(real(z))==0, zero(z), -zero(z)))
*(x::Dual4, z::Bool) = z*x

*(z::Dual4, w::Dual4) = Dual4(real(z)*real(w), epsilon1(z)*real(w)+real(z)*epsilon1(w),
epsilon2(z)*real(w)+real(z)*epsilon2(w), epsilon3(z)*real(w)+real(z)*epsilon3(w),
epsilon4(z)*real(w)+real(z)*epsilon4(w))
*(x::Real, z::Dual4) = Dual4(x*real(z), x*epsilon1(z), x*epsilon2(z), x*epsilon3(z), x*epsilon4(z))
*(z::Dual4, x::Real) = x*z

function /(z::Real, w::Dual4)
    df = -z/real(w)^2
    Dual4(z/real(w), epsilon1(w)*df, epsilon2(w)*df, epsilon3(w)*df, epsilon4(w)*df)
end
/(z::Dual4, x::Real) = Dual4(real(z)/x, epsilon1(z)/x, epsilon2(z)/x, epsilon3(z)/x, epsilon4(z)/x)
function /(z::Dual4, w::Dual4)
    denom = 1/(real(w)*real(w))
    Dual4(real(z)/real(w), (epsilon1(z)*real(w)-real(z)*epsilon1(w))*denom,
        (epsilon2(z)*real(w)-real(z)*epsilon2(w))*denom,
        (epsilon3(z)*real(w)-real(z)*epsilon3(w))*denom,
        (epsilon4(z)*real(w)-real(z)*epsilon4(w))*denom)
end

for f in [:^, :(NaNMath.pow)]
    @eval function ($f)(z::Dual4, w::Dual4)
        re = $f(real(z),real(w))
        powval = real(w)*(($f)(real(z),real(w)-1))
        logval = ($f)(real(z),real(w))*log(real(z))
        du1 = epsilon1(z)*powval+epsilon1(w)*logval
        du2 = epsilon2(z)*powval+epsilon2(w)*logval
        du3 = epsilon3(z)*powval+epsilon3(w)*logval
        du4 = epsilon4(z)*powval+epsilon4(w)*logval
        Dual4(re, du1, du2, du3, du4)
    end
end

# these two definitions are needed to fix ambiguity warnings
function ^(z::Dual4, n::Integer)
    zn1 = n*real(z)^(n-1)
    Dual4(real(z)^n, epsilon1(z)*zn1, epsilon2(z)*zn1, epsilon3(z)*zn1, epsilon4(z)*zn1)
end
^(z::Dual4, n::Rational) = invoke(^, (Dual4,Real), z,n)

function ^(z::Dual4, n::Real)
    zn1 = n*real(z)^(n-1)
    Dual4(real(z)^n, epsilon1(z)*zn1, epsilon2(z)*zn1, epsilon3(z)*zn1, epsilon4(z)*zn1)
end

function NaNMath.pow(z::Dual4, n::Real)
    powval = n*NaNMath.pow(real(z),n-1)
    Dual4(NaNMath.pow(real(z),n), epsilon1(z)*powval, epsilon2(z)*powval, epsilon3(z)*powval, epsilon4(z)*powval)
end
function NaNMath.pow(z::Real, w::Dual4)
    logval = NaNMath.pow(z,real(w))*log(z) 
    Dual4(NaNMath.pow(z,real(w)), epsilon1(w)*logval, epsilon2(w)*logval, epsilon3(w)*logval, epsilon4(w)*logval)
end

for (funsym, exp) in Calculus.symbolic_derivatives_1arg()
    funsym == :exp && continue
    funsym == :abs2 && continue
    @eval function $(funsym)(z::Dual4)
        x = real(z)
        df = $exp
        Dual4($(funsym)(x),epsilon1(z)*df,epsilon2(z)*df,epsilon3(z)*df,epsilon4(z)*df)
    end
    # extend corresponding NaNMath methods
    if funsym in (:sin, :cos, :tan, :asin, :acos, :acosh, :atanh, :log, :log2, :log10,
          :lgamma, :log1p)
        funsym = Expr(:.,:NaNMath,Base.Meta.quot(funsym))
        @eval function $(funsym)(z::Dual4)
            x = real(z)
            df = $(to_nanmath(exp))
            Dual4($(funsym)(x),epsilon1(z)*df,epsilon2(z)*df,epsilon3(z)*df,epsilon4(z)*df)
        end
    end
end

# only need to compute exp once
function exp(z::Dual4)
    df = exp(real(z))
    return Dual4(df, epsilon1(z)*df, epsilon2(z)*df, epsilon3(z)*df, epsilon4(z)*df)
end
