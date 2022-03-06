const ReComp = Union{Real,Complex}

struct Dual{T<:ReComp} <: Number
    value::T
    epsilon::T
end
Dual(x::S, y::T) where {S<:ReComp,T<:ReComp} = Dual(promote(x,y)...)
Dual(x::ReComp) = Dual(x, zero(x))
Dual{T}(x::ReComp) where T<:ReComp = Dual{T}(T(x), zero(T))

const ɛ = Dual(false, true)
const imɛ = Dual(Complex(false, false), Complex(false, true))

const Dual128 = Dual{Float64}
const Dual64  = Dual{Float32}
const Dual32  = Dual{Float16}
const DualComplex256 = Dual{ComplexF64}
const DualComplex128 = Dual{ComplexF32}
const DualComplex64  = Dual{ComplexF16}

Base.convert(::Type{Dual{T}}, z::Dual{T}) where {T<:ReComp} = z
Base.convert(::Type{Dual{T}}, z::Dual) where {T<:ReComp} = Dual{T}(convert(T, value(z)), convert(T, epsilon(z)))
Base.convert(::Type{Dual{T}}, x::Number) where {T<:ReComp} = Dual{T}(convert(T, x), convert(T, 0))
Base.convert(::Type{T}, z::Dual) where {T<:ReComp} = (epsilon(z)==0 ? convert(T, value(z)) : throw(InexactError()))

Base.promote_rule(::Type{Dual{T}}, ::Type{Dual{S}}) where {T<:ReComp,S<:ReComp} = Dual{promote_type(T, S)}
Base.promote_rule(::Type{Dual{T}}, ::Type{S}) where {T<:ReComp,S<:ReComp} = Dual{promote_type(T, S)}
Base.promote_rule(::Type{Dual{T}}, ::Type{T}) where {T<:ReComp} = Dual{T}

Base.widen(::Type{Dual{T}}) where {T} = Dual{widen(T)}

value(z::Dual) = z.value
epsilon(z::Dual) = z.epsilon

value(x::Number) = x
epsilon(x::Number) = zero(typeof(x))

dual(x::ReComp, y::ReComp) = Dual(x, y)
dual(x::ReComp) = Dual(x)
dual(z::Dual) = z

function Base.complex(x::Dual, y::Dual)
    dual(complex(value(x), value(y)), complex(epsilon(x), epsilon(y)))
end
Base.complex(x::Real, y::Dual) = complex(dual(x), y)
Base.complex(x::Dual, y::Real) = complex(x, dual(y))
Base.complex(::Type{Dual{T}}) where {T} = Dual{complex(T)}

const realpart = value
const dualpart = epsilon

Base.isnan(z::Dual) = isnan(value(z))
Base.isinf(z::Dual) = isinf(value(z))
Base.isfinite(z::Dual) = isfinite(value(z))
isdual(x::Dual) = true
isdual(x::Number) = false
Base.eps(z::Dual) = eps(value(z))
Base.eps(::Type{Dual{T}}) where {T} = eps(T)

function dual_show(io::IO, z::Dual{T}, compact::Bool) where T<:Real
    x, y = value(z), epsilon(z)
    if isnan(x) || isfinite(y)
        compact ? show(IOContext(io, :compact=>true), x) : show(io, x)
        if signbit(y)==1 && !isnan(y)
            y = -y
            print(io, compact ? "-" : " - ")
        else
            print(io, compact ? "+" : " + ")
        end
        compact ? show(IOContext(io, :compact=>true), y) : show(io, y)
        printtimes(io, y)
        print(io, "ɛ")
    else
        print(io, "Dual{",T,"}(", x, ",", y, ")")
    end
end

function dual_show(io::IO, z::Dual{T}, compact::Bool) where T<:Complex
    x, y = value(z), epsilon(z)
    xr, xi = reim(x)
    yr, yi = reim(y)
    if isnan(x) || isfinite(y)
        compact ? show(IOContext(io, :compact=>true), x) : show(io, x)
        if signbit(yr)==1 && !isnan(y)
            yr = -yr
            print(io, " - ")
        else
            print(io, " + ")
        end
        if compact
            if signbit(yi)==1 && !isnan(y)
                yi = -yi
                show(IOContext(io, :compact=>true), yr)
                printtimes(io, yr)
                print(io, "ɛ-")
                show(IOContext(io, :compact=>true), yi)
            else
                show(IOContext(io, :compact=>true), yr)
                print(io, "ɛ+")
                show(IOContext(io, :compact=>true), yi)
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

function dual_show(io::IO, z::Dual{T}, compact::Bool) where T<:Bool
    x, y = value(z), epsilon(z)
    if !value(z) && epsilon(z)
        print(io, "ɛ")
    else
        print(io, "Dual{",T,"}(", x, ",", y, ")")
    end
end

function dual_show(io::IO, z::Dual{Complex{T}}, compact::Bool) where T<:Bool
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

Base.show(io::IO, z::Dual) = dual_show(io, z, get(IOContext(io), :compact, false))

function Base.read(s::IO, ::Type{Dual{T}}) where T<:ReComp
    x = read(s, T)
    y = read(s, T)
    Dual{T}(x, y)
end
function Base.write(s::IO, z::Dual)
    write(s, value(z))
    write(s, epsilon(z))
end


## Generic functions of dual numbers ##

Base.convert(::Type{Dual}, z::Dual) = z
Base.convert(::Type{Dual}, x::Number) = Dual(x)

Base.:(==)(z::Dual, w::Dual) = value(z) == value(w)
Base.:(==)(z::Dual, x::Number) = value(z) == x
Base.:(==)(x::Number, z::Dual) = value(z) == x

Base.isequal(z::Dual, w::Dual) = isequal(value(z),value(w)) && isequal(epsilon(z), epsilon(w))
Base.isequal(z::Dual, x::Number) = isequal(value(z), x) && isequal(epsilon(z), zero(x))
Base.isequal(x::Number, z::Dual) = isequal(z, x)

Base.isless(z::Dual{<:Real},w::Dual{<:Real}) = value(z) < value(w)
Base.isless(z::Real,w::Dual{<:Real}) = z < value(w)
Base.isless(z::Dual{<:Real},w::Real) = value(z) < w

Base.hash(z::Dual) = (x = hash(value(z)); epsilon(z)==0 ? x : bitmix(x,hash(epsilon(z))))

Base.float(z::Union{Dual{T}, Dual{Complex{T}}}) where {T<:AbstractFloat} = z
Base.complex(z::Dual{<:Complex}) = z

Base.floor(z::Dual) = floor(value(z))
Base.ceil(z::Dual)  = ceil(value(z))
Base.trunc(z::Dual) = trunc(value(z))
Base.round(z::Dual) = round(value(z))
Base.floor(::Type{T}, z::Dual) where {T<:Real} = floor(T, value(z))
Base.ceil( ::Type{T}, z::Dual) where {T<:Real} = ceil( T, value(z))
Base.trunc(::Type{T}, z::Dual) where {T<:Real} = trunc(T, value(z))
Base.round(::Type{T}, z::Dual) where {T<:Real} = round(T, value(z))

for op in (:real, :imag, :conj, :float, :complex)
    @eval Base.$op(z::Dual) = Dual($op(value(z)), $op(epsilon(z)))
end

Base.abs(z::Dual) = sqrt(abs2(z))
Base.abs2(z::Dual) = real(conj(z)*z)

Base.real(z::Dual{<:Real}) = z
Base.abs(z::Dual{<:Real}) = z ≥ 0 ? z : -z

Base.angle(z::Dual{<:Real}) = z ≥ 0 ? zero(z) : one(z)*π
function Base.angle(z::Dual{Complex{T}}) where T<:Real
    if z == 0
        if imag(epsilon(z)) == 0
            Dual(zero(T), zero(T))
        else
            Dual(zero(T), convert(T, Inf))
        end
    else
        real(log(sign(z)) / im)
    end
end

Base.flipsign(x::Dual,y::Dual) = y == 0 ? flipsign(x, epsilon(y)) : flipsign(x, value(y))
Base.flipsign(x, y::Dual) = y == 0 ? flipsign(x, epsilon(y)) : flipsign(x, value(y))
Base.flipsign(x::Dual, y) = dual(flipsign(value(x), y), flipsign(epsilon(x), y))

# algebraic definitions
conjdual(z::Dual) = Dual(value(z),-epsilon(z))
absdual(z::Dual) = abs(value(z))
abs2dual(z::Dual) = abs2(value(z))

# algebra

Base.:+(z::Dual, w::Dual) = Dual(value(z)+value(w), epsilon(z)+epsilon(w))
Base.:+(z::Number, w::Dual) = Dual(z+value(w), epsilon(w))
Base.:+(z::Dual, w::Number) = Dual(value(z)+w, epsilon(z))

Base.:-(z::Dual) = Dual(-value(z), -epsilon(z))
Base.:-(z::Dual, w::Dual) = Dual(value(z)-value(w), epsilon(z)-epsilon(w))
Base.:-(z::Number, w::Dual) = Dual(z-value(w), -epsilon(w))
Base.:-(z::Dual, w::Number) = Dual(value(z)-w, epsilon(z))

# avoid ambiguous definition with Bool*Number
Base.:*(x::Bool, z::Dual) = ifelse(x, z, ifelse(signbit(real(value(z)))==0, zero(z), -zero(z)))
Base.:*(x::Dual, z::Bool) = z*x

Base.:*(z::Dual, w::Dual) = Dual(value(z)*value(w), epsilon(z)*value(w)+value(z)*epsilon(w))
Base.:*(x::Number, z::Dual) = Dual(x*value(z), x*epsilon(z))
Base.:*(z::Dual, x::Number) = Dual(x*value(z), x*epsilon(z))

Base.:/(z::Dual, w::Dual) = Dual(value(z)/value(w), (epsilon(z)*value(w)-value(z)*epsilon(w))/(value(w)*value(w)))
Base.:/(z::Number, w::Dual) = Dual(z/value(w), -z*epsilon(w)/value(w)^2)
Base.:/(z::Dual, x::Number) = Dual(value(z)/x, epsilon(z)/x)

for f in [:(Base.:^), :(NaNMath.pow)]
    @eval function ($f)(z::Dual{T1}, w::Dual{T2}) where {T1, T2}
        T = promote_type(T1, T2) # for type stability in ? : statements
        val = $f(value(z), value(w))

        ezvw = epsilon(z) * value(w) # for using in ? : statement
        du1 = iszero(ezvw) ? zero(T) : ezvw * $f(value(z), value(w) - 1)
        ew = epsilon(w) # for using in ? : statement
        # the float is for type stability because log promotes to floats
        du2 = iszero(ew) ? zero(float(T)) : ew * val * log(value(z))
        du = du1 + du2

        Dual(val, du)
    end
end

Base.mod(z::Dual, n::Number) = Dual(mod(value(z), n), epsilon(z))

# introduce a boolean !iszero(n) for hard zero behaviour to combat NaNs
function pow(z::Dual, n::AbstractFloat)
    return Dual(value(z)^n, !iszero(n) * (epsilon(z) * n * value(z)^(n - 1)))
end
function pow(z::Dual{T}, n::Integer) where T
    iszero(n) && return Dual(one(T), zero(T)) # avoid DomainError Int^(negative Int)
    isone(z) && return Dual(one(T), epsilon(z) * n)
    return Dual(value(z)^n, epsilon(z) * n * value(z)^(n - 1))
end
# these first two definitions are needed to fix ambiguity warnings
for T1 ∈ (:Integer, :Rational, :Number)
    @eval Base.:^(z::Dual{T}, n::$T1) where T = pow(z, n)
end


NaNMath.pow(z::Dual{T}, n::Number) where T = Dual(NaNMath.pow(value(z),n), epsilon(z)*n*NaNMath.pow(value(z),n-1))
NaNMath.pow(z::Number, w::Dual{T}) where T = Dual(NaNMath.pow(z,value(w)), epsilon(w)*NaNMath.pow(z,value(w))*log(z))

Base.inv(z::Dual) = dual(inv(value(z)),-epsilon(z)/value(z)^2)

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
    funsym == :inv && continue
    if isdefined(SpecialFunctions, funsym)
        @eval function SpecialFunctions.$(funsym)(z::Dual)
            x = value(z)
            xp = epsilon(z)
            Dual($(funsym)(x),xp*$exp)
        end
    elseif isdefined(Base, funsym)
        @eval function Base.$(funsym)(z::Dual)
            x = value(z)
            xp = epsilon(z)
            Dual($(funsym)(x),xp*$exp)
        end
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
Base.exp(z::Dual) = (expval = exp(value(z)); Dual(expval, epsilon(z)*expval))
Base.cis(z::Dual) = (cisval = cis(value(z)); Dual(cisval, im*epsilon(z)*cisval))

Base.exp10(x::Dual) = (y = exp10(value(x)); Dual(y, y * log(10) * epsilon(x)))

## TODO: should be generated in Calculus
Base.sinpi(z::Dual) = Dual(sinpi(value(z)),epsilon(z)*cospi(value(z))*π)
Base.cospi(z::Dual) = Dual(cospi(value(z)),-epsilon(z)*sinpi(value(z))*π)

Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, i::Dual) = checkindex(Bool, inds, value(i))
