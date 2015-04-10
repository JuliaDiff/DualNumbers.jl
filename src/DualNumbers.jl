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
    dual,
    dual128,
    dual64,
    isdual,
    dual_show,
    epsilon,
    conjdual,
    absdual,
    abs2dual
end
