__precompile__()

module DualNumbers
  importall Base

  import NaNMath
  import Calculus
  using Compat

  const NANSAFE_MODE_ENABLED = false                                              
                                                                                  
  const IS_MULTITHREADED_JULIA = VERSION >= v"0.5.0-dev+923" 
  const AUTO_DEFINED_UNARY_FUNCS = map(first, Calculus.symbolic_derivatives_1arg())

  const NANMATH_FUNCS = (:sin, :cos, :tan, :asin, :acos, :acosh,                  
                         :atanh, :log, :log2, :log10, :lgamma, :log1p)

  include("partials.jl")
  include("dual.jl")
  include("deprecate.jl")

  export Dual,
  value,
  partials

  export
      Dual128,
      Dual64,
      Dual32,
      epsilon
      conjdual,
      absdual,
      abs2dual
      #=DualComplex256,
      DualComplex128,
      DualComplex64,
      dual,
      realpart,
      dualpart,
      isdual,
      dual_show=#
end
