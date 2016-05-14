__precompile__()

module DualNumbers
  importall Base

  import NaNMath
  import Calculus

  include("dual.jl")
  include("dual_n.jl")

  Base.@deprecate_binding du ɛ
  @deprecate inf{T}(::Type{Dual{T}}) convert(Dual{T}, Inf)
  @deprecate nan{T}(::Type{Dual{T}}) convert(Dual{T}, NaN)

  export
    Dual,
    Dual128,
    Dual64,
    Dual32,
    DualComplex256,
    DualComplex128,
    DualComplex64,
    dual,
    epsilon,
    realpart,
    dualpart,
    isdual,
    dual_show,
    conjdual,
    absdual,
    abs2dual,
    ɛ,
    imɛ
end
