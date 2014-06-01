immutable Dual{T<:Real} <: Number
    re::T
    du::T
end
Dual(x::Real, y::Real) = Dual(promote(x,y)...)
Dual(x::Real) = Dual(x, zero(x))

typealias Dual128 Dual{Float64}
typealias Dual64 Dual{Float32}
typealias DualPair Dual

real(z::Dual) = z.re
epsilon(z::Dual) = z.du

eps(z::Dual) = eps(real(z))
eps{T}(::Type{Dual{T}}) = eps(T)
one(z::Dual) = dual(one(real(z)))
one{T}(::Type{Dual{T}}) = dual(one(T))
inf{T}(::Type{Dual{T}}) = dual(inf(T))
nan{T}(::Type{Dual{T}}) = nan(T)
isnan(z::Dual) = isnan(real(z))

convert{T<:Real}(::Type{Dual{T}}, x::Real) =
  Dual{T}(convert(T, x), convert(T, 0))
convert{T<:Real}(::Type{Dual{T}}, z::Dual{T}) = z
convert{T<:Real}(::Type{Dual{T}}, z::Dual) =
  Dual{T}(convert(T, real(z)), convert(T, epsilon(z)))

convert{T<:Real}(::Type{T}, z::Dual) =
  (epsilon(z)==0 ? convert(T, real(z)) : throw(InexactError()))

promote_rule{T<:Real, S<:Real}(::Type{Dual{T}}, ::Type{Dual{S}}) =
    Dual{promote_type(T, S)}
# these promotion rules shouldn't be used for scalar operations -- they're slow
promote_rule{T<:Real}(::Type{Dual{T}}, ::Type{T}) = Dual{T}
promote_rule{T<:Real, S<:Real}(::Type{Dual{T}}, ::Type{S}) =
  Dual{promote_type(T, S)}

dual(x, y) = Dual(x, y)
dual(x) = Dual(x)

@vectorize_1arg Real dual

dual128(x::Float64, y::Float64) = Dual{Float64}(x, y)
dual128(x::Real, y::Real) = dual128(float64(x), float64(y))
dual128(z) = dual128(real(z), epsilon(z))
dual64(x::Float32, y::Float32) = Dual{Float32}(x, y)
dual64(x::Real, y::Real) = dual64(float32(x), float32(y))
dual64(z) = dual64(real(z), epsilon(z))

isdual(x::Dual) = true
isdual(x::Number) = false

real_valued{T<:Real}(z::Dual{T}) = epsilon(z) == 0
integer_valued(z::Dual) = real_valued(z) && integer_valued(real(z))

isfinite(z::Dual) = isfinite(real(z))
reim(z::Dual) = (real(z), epsilon(z))

function dual_show(io::IO, z::Dual, compact::Bool)
    x, y = reim(z)
    if isnan(x) || isfinite(y)
        compact ? showcompact(io,x) : show(io,x)
        if signbit(y)==1 && !isnan(y)
            y = -y
            print(io, compact ? "-" : " - ")
        else
            print(io, compact ? "+" : " + ")
        end
        compact ? showcompact(io, y) : show(io, y)
        if !(isa(y,Integer) || isa(y,Rational) ||
             isa(y,FloatingPoint) && isfinite(y))
            print(io, "*")
        end
        print(io, "du")
    else
        print(io, "dual(", x, ",", y, ")")
    end
end
show(io::IO, z::Dual) = dual_show(io, z, false)
showcompact(io::IO, z::Dual) = dual_show(io, z, true)

function read{T<:Real}(s::IO, ::Type{Dual{T}})
    x = read(s, T)
    y = read(s, T)
    Dual{T}(x, y)
end
function write(s::IO, z::Dual)
    write(s, real(z))
    write(s, epsilon(z))
end


## Generic functions of dual numbers ##

convert(::Type{Dual}, z::Dual) = z
convert(::Type{Dual}, x::Real) = dual(x)

==(z::Dual, w::Dual) = real(z) == real(w)
==(z::Dual, x::Real) = real(z) == x
==(x::Real, z::Dual) = real(z) == x

isequal(z::Dual, w::Dual) =
  isequal(real(z),real(w)) && isequal(epsilon(z), epsilon(w))
isequal(z::Dual, x::Real) = real_valued(z) && isequal(real(z), x)
isequal(x::Real, z::Dual) = real_valued(z) && isequal(real(z), x)

isless(z::Dual,w::Dual) = real(z) < real(w)
isless(z::Number,w::Dual) = z < real(w)
isless(z::Dual,w::Number) = real(z) < w

hash(z::Dual) =
  (x = hash(real(z)); real_valued(z) ? x : bitmix(x,hash(epsilon(z))))

# we don't support Dual{Complex}, so conj is a noop
conj(z::Dual) = z
abs(z::Dual)  = (real(z) >= 0) ? z : -z
abs2(z::Dual) = z*z

# algebraic definitions
conjdual(z::Dual) = Dual(real(z),-epsilon(z))
absdual(z::Dual) = abs(real(z))
abs2dual(z::Dual) = abs2(real(z))

+(z::Dual, w::Dual) = dual(real(z)+real(w), epsilon(z)+epsilon(w))
+(z::Number, w::Dual) = dual(z+real(w), epsilon(w))
+(z::Dual, w::Number) = dual(real(z)+w, epsilon(z))

-(z::Dual) = dual(-real(z), -epsilon(z))
-(z::Dual, w::Dual) = dual(real(z)-real(w), epsilon(z)-epsilon(w))
-(z::Number, w::Dual) = dual(z-real(w), -epsilon(w))
-(z::Dual, w::Number) = dual(real(z)-w, epsilon(z))

# avoid ambiguous definition with Bool*Number
*(x::Bool, z::Dual) = ifelse(x, z, ifelse(signbit(real(z))==0, zero(z), -zero(z)))
*(x::Dual, z::Bool) = z*x

*(z::Dual, w::Dual) = dual(real(z)*real(w), epsilon(z)*real(w)+real(z)*epsilon(w))
*(x::Real, z::Dual) = dual(x*real(z), x*epsilon(z))
*(z::Dual, x::Real) = dual(x*real(z), x*epsilon(z))

/(z::Real, w::Dual) = dual(z/real(w), -z*epsilon(w)/real(w)^2)
/(z::Dual, x::Real) = dual(real(z)/x, epsilon(z)/x)
/(z::Dual, w::Dual) =
  dual(real(z)/real(w), (epsilon(z)*real(w)-real(z)*epsilon(w))/(real(w)*real(w)))

function ^(z::Dual, w::Dual)
  re = real(z)^real(w)
  
  du =
    epsilon(z)*real(w)*(real(z)^(real(w)-1))+epsilon(w)*(real(z)^real(w))*log(real(z))
    
  dual(re, du)
end

# these two definitions are needed to fix ambiguity warnings
^(z::Dual, n::Integer) = dual(real(z)^n, epsilon(z)*n*real(z)^(n-1))
^(z::Dual, n::Rational) = dual(real(z)^n, epsilon(z)*n*real(z)^(n-1))

^(z::Dual, n::Real) = dual(real(z)^n, epsilon(z)*n*real(z)^(n-1))

for (funsym, exp) in Calculus.derivative_rules
    @eval function $(funsym)(z::Dual)
        xp = epsilon(z)
        x = real(z)
        Dual($(funsym)(x),$exp)
    end
end

