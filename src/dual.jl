const ReComp = Union{Real,Complex}

immutable Dual{T<:ReComp} <: Number
    value::T
    epsilon::T
end
Dual{S<:ReComp,T<:ReComp}(x::S, y::T) = Dual(promote(x,y)...)
Dual(x::ReComp) = Dual(x, zero(x))

const ɛ = Dual(false, true)
const imɛ = Dual(Complex(false, false), Complex(false, true))

const Dual128 = Dual{Float64}
const Dual64  = Dual{Float32}
const Dual32  = Dual{Float16}
const DualComplex256 = Dual{Complex128}
const DualComplex128 = Dual{Complex64}
const DualComplex64  = Dual{Complex32}

convert{T<:ReComp}(::Type{Dual{T}}, z::Dual{T}) = z
convert{T<:ReComp}(::Type{Dual{T}}, z::Dual) = Dual{T}(convert(T, value(z)), convert(T, epsilon(z)))
convert{T<:ReComp}(::Type{Dual{T}}, x::Number) = Dual{T}(convert(T, x), convert(T, 0))
convert{T<:ReComp}(::Type{T}, z::Dual) = (epsilon(z)==0 ? convert(T, value(z)) : throw(InexactError()))

promote_rule{T<:ReComp, S<:ReComp}(::Type{Dual{T}}, ::Type{Dual{S}}) = Dual{promote_type(T, S)}
promote_rule{T<:ReComp, S<:ReComp}(::Type{Dual{T}}, ::Type{S}) = Dual{promote_type(T, S)}
promote_rule{T<:ReComp}(::Type{Dual{T}}, ::Type{T}) = Dual{T}

widen{T}(::Type{Dual{T}}) = Dual{widen(T)}

value(z::Dual) = z.value
epsilon(z::Dual) = z.epsilon

dual(x::ReComp, y::ReComp) = Dual(x, y)
dual(x::ReComp) = Dual(x)
dual(z::Dual) = z

Compat.@dep_vectorize_1arg ReComp dual
Compat.@dep_vectorize_2arg ReComp dual
Compat.@dep_vectorize_1arg Dual dual
Compat.@dep_vectorize_1arg Dual value
Compat.@dep_vectorize_1arg Dual epsilon

const realpart = value
const dualpart = epsilon

isnan(z::Dual) = isnan(value(z))
isinf(z::Dual) = isinf(value(z))
isfinite(z::Dual) = isfinite(value(z))
isdual(x::Dual) = true
isdual(x::Number) = false
eps(z::Dual) = eps(value(z))
eps{T}(::Type{Dual{T}}) = eps(T)

function dual_show{T<:Real}(io::IO, z::Dual{T}, compact::Bool)
    x, y = value(z), epsilon(z)
    if isnan(x) || isfinite(y)
        compact ? showcompact(io,x) : show(io,x)
        if signbit(y)==1 && !isnan(y)
            y = -y
            print(io, compact ? "-" : " - ")
        else
            print(io, compact ? "+" : " + ")
        end
        compact ? showcompact(io, y) : show(io, y)
        printtimes(io, y)
        print(io, "ɛ")
    else
        print(io, "Dual{",T,"}(", x, ",", y, ")")
    end
end

function dual_show{T<:Complex}(io::IO, z::Dual{T}, compact::Bool)
    x, y = value(z), epsilon(z)
    xr, xi = reim(x)
    yr, yi = reim(y)
    if isnan(x) || isfinite(y)
        compact ? showcompact(io,x) : show(io,x)
        if signbit(yr)==1 && !isnan(y)
            yr = -yr
            print(io, " - ")
        else
            print(io, " + ")
        end
        if compact
            if signbit(yi)==1 && !isnan(y)
                yi = -yi
                showcompact(io, yr)
                printtimes(io, yr)
                print(io, "ɛ-")
                showcompact(io, yi)
            else
                showcompact(io, yr)
                print(io, "ɛ+")
                showcompact(io, yi)
            end
        else
            if signbit(yi)==1 && !isnan(y)
                yi = -yi
                show(io, yr)
                printtimes(io, yr)
                print(io, "ɛ - ")
                show(io, yi)
            else
                show(io, yr)
                print(io, "ɛ + ")
                show(io, yi)
            end
        end
        printtimes(io, yi)
        print(io, "imɛ")
    else
        print(io, "Dual{",T,"}(", x, ",", y, ")")
    end
end

function dual_show{T<:Bool}(io::IO, z::Dual{T}, compact::Bool)
    x, y = value(z), epsilon(z)
    if !value(z) && epsilon(z)
        print(io, "ɛ")
    else
        print(io, "Dual{",T,"}(", x, ",", y, ")")
    end
end

function dual_show{T<:Bool}(io::IO, z::Dual{Complex{T}}, compact::Bool)
    x, y = value(z), epsilon(z)
    xr, xi = reim(x)
    yr, yi = reim(y)
    if !xr
        if xi*!yr*!yi
            print(io, "im")
        elseif !xi*yr*!yi
            print(io, "ɛ")
        elseif !xi*!yr*yi
            print(io, "imɛ")
        end
    else
        print(io, "Dual{",T,"}(", x, ",", y, ")")
    end
end

function printtimes(io::IO, x::Real)
    if !(isa(x,Integer) || isa(x,Rational) ||
         isa(x,AbstractFloat) && isfinite(x))
        print(io, "*")
    end
end

show(io::IO, z::Dual) = dual_show(io, z, false)
showcompact(io::IO, z::Dual) = dual_show(io, z, true)

function read{T<:ReComp}(s::IO, ::Type{Dual{T}})
    x = read(s, T)
    y = read(s, T)
    Dual{T}(x, y)
end
function write(s::IO, z::Dual)
    write(s, value(z))
    write(s, epsilon(z))
end


## Generic functions of dual numbers ##

convert(::Type{Dual}, z::Dual) = z
convert(::Type{Dual}, x::Number) = Dual(x)

==(z::Dual, w::Dual) = value(z) == value(w)
==(z::Dual, x::Number) = value(z) == x
==(x::Number, z::Dual) = value(z) == x

isequal(z::Dual, w::Dual) = isequal(value(z),value(w)) && isequal(epsilon(z), epsilon(w))
isequal(z::Dual, x::Number) = isequal(value(z), x) && isequal(epsilon(z), zero(x))
isequal(x::Number, z::Dual) = isequal(z, x)

isless{T<:Real}(z::Dual{T},w::Dual{T}) = value(z) < value(w)
isless{T<:Real}(z::Real,w::Dual{T}) = z < value(w)
isless{T<:Real}(z::Dual{T},w::Real) = value(z) < w

hash(z::Dual) = (x = hash(value(z)); epsilon(z)==0 ? x : bitmix(x,hash(epsilon(z))))

float{T<:AbstractFloat}(z::Union{Dual{T},Dual{Complex{T}}})=z
complex{T<:Real}(z::Dual{Complex{T}})=z

floor{T<:Real}(::Type{T}, z::Dual) = floor(T, value(z))
ceil{ T<:Real}(::Type{T}, z::Dual) = ceil( T, value(z))
trunc{T<:Real}(::Type{T}, z::Dual) = trunc(T, value(z))
round{T<:Real}(::Type{T}, z::Dual) = round(T, value(z))

for op in (:real,:imag,:conj,:float,:complex)
    @eval begin
        $op(z::Dual) = Dual($op(value(z)),$op(epsilon(z)))
    end
end

abs(z::Dual) = sqrt(abs2(z))
abs2(z::Dual) = real(conj(z)*z)

real{T<:Real}(z::Dual{T}) = z
abs{T<:Real}(z::Dual{T}) = z ≥ 0 ? z : -z

angle{T<:Real}(z::Dual{T}) = z ≥ 0 ? zero(z) : one(z)*π
angle{T<:Real}(z::Dual{Complex{T}}) = z == 0 ? (imag(epsilon(z)) == 0 ? Dual(zero(T),zero(T)) : Dual(zero(T),convert(T, Inf))) : real(log(sign(z))/im)

# algebraic definitions
conjdual(z::Dual) = Dual(value(z),-epsilon(z))
absdual(z::Dual) = abs(value(z))
abs2dual(z::Dual) = abs2(value(z))

# algebra

+(z::Dual, w::Dual) = Dual(value(z)+value(w), epsilon(z)+epsilon(w))
+(z::Number, w::Dual) = Dual(z+value(w), epsilon(w))
+(z::Dual, w::Number) = Dual(value(z)+w, epsilon(z))

-(z::Dual) = Dual(-value(z), -epsilon(z))
-(z::Dual, w::Dual) = Dual(value(z)-value(w), epsilon(z)-epsilon(w))
-(z::Number, w::Dual) = Dual(z-value(w), -epsilon(w))
-(z::Dual, w::Number) = Dual(value(z)-w, epsilon(z))

# avoid ambiguous definition with Bool*Number
*(x::Bool, z::Dual) = ifelse(x, z, ifelse(signbit(real(value(z)))==0, zero(z), -zero(z)))
*(x::Dual, z::Bool) = z*x

*(z::Dual, w::Dual) = Dual(value(z)*value(w), epsilon(z)*value(w)+value(z)*epsilon(w))
*(x::Number, z::Dual) = Dual(x*value(z), x*epsilon(z))
*(z::Dual, x::Number) = Dual(x*value(z), x*epsilon(z))

/(z::Dual, w::Dual) = Dual(value(z)/value(w), (epsilon(z)*value(w)-value(z)*epsilon(w))/(value(w)*value(w)))
/(z::Number, w::Dual) = Dual(z/value(w), -z*epsilon(w)/value(w)^2)
/(z::Dual, x::Number) = Dual(value(z)/x, epsilon(z)/x)

for f in [:^, :(NaNMath.pow)]
    @eval function ($f)(z::Dual, w::Dual)
        if epsilon(w) == 0.0
            return $f(z,value(w))
        end
        val = $f(value(z),value(w))

        du =
        epsilon(z)*value(w)*(($f)(value(z),value(w)-1))+epsilon(w)*($f)(value(z),value(w))*log(value(z))

        Dual(val, du)
    end
end

mod(z::Dual, n::Number) = Dual(mod(value(z), n), epsilon(z))

# these two definitions are needed to fix ambiguity warnings
^(z::Dual, n::Integer) = Dual(value(z)^n, epsilon(z)*n*value(z)^(n-1))
^(z::Dual, n::Rational) = Dual(value(z)^n, epsilon(z)*n*value(z)^(n-1))

^(z::Dual, n::Number) = Dual(value(z)^n, epsilon(z)*n*value(z)^(n-1))
NaNMath.pow(z::Dual, n::Number) = Dual(NaNMath.pow(value(z),n), epsilon(z)*n*NaNMath.pow(value(z),n-1))
NaNMath.pow(z::Number, w::Dual) = Dual(NaNMath.pow(z,value(w)), epsilon(w)*NaNMath.pow(z,value(w))*log(z))

# force use of NaNMath functions in derivative calculations
function to_nanmath(x::Expr)
    if x.head == :call
        funsym = Expr(:.,:NaNMath,Base.Meta.quot(x.args[1]))
        return Expr(:call,funsym,[to_nanmath(z) for z in x.args[2:end]]...)
    else
        return Expr(:call,[to_nanmath(z) for z in x.args]...)
    end
end
to_nanmath(x) = x

for (funsym, exp) in Calculus.symbolic_derivatives_1arg()
    funsym == :exp && continue
    funsym == :abs2 && continue
    @eval function $(funsym)(z::Dual)
        x = value(z)
        xp = epsilon(z)
        Dual($(funsym)(x),xp*$exp)
    end
    # extend corresponding NaNMath methods
    if funsym in (:sin, :cos, :tan, :asin, :acos, :acosh, :atanh, :log, :log2, :log10,
          :lgamma, :log1p)
        funsym = Expr(:.,:NaNMath,Base.Meta.quot(funsym))
        @eval function $(funsym)(z::Dual)
            x = value(z)
            xp = epsilon(z)
            Dual($(funsym)(x),xp*$(to_nanmath(exp)))
        end
    end
end

# only need to compute exp/cis once
exp(z::Dual) = (expval = exp(value(z)); Dual(expval, epsilon(z)*expval))
cis(z::Dual) = (cisval = cis(value(z)); Dual(cisval, im*epsilon(z)*cisval))

## TODO: should be generated in Calculus
sinpi(z::Dual) = Dual(sinpi(value(z)),epsilon(z)*cospi(value(z))*π)
cospi(z::Dual) = Dual(cospi(value(z)),-epsilon(z)*sinpi(value(z))*π)

if VERSION >= v"0.5.0-dev+5429"
    Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, i::Dual) = checkindex(Bool, inds, value(i))
end
