__precompile__()

module DualNumbers
  importall Base

  import NaNMath
  import Calculus

  include("dual.jl")
  include("dual_n.jl")

  export
    Dual,
    Dual128,
    Dual64,
    DualPair,
    isdual,
    dual_show,
    epsilon,
    conjdual,
    absdual,
    abs2dual,
    du
end
