module DualNumbers
  importall Base

  import NaNMath
  import Calculus

  include("dual.jl")
  include("dual_4.jl")
  include("dual_n.jl")

  export
    Dual,
    DualN,
    Dual128,
    Dual64,
    DualPair,
    dual,
    dual128,
    dual64,
    isdual,
    dual_show,
    epsilon,
    neps,
    conjdual,
    absdual,
    abs2dual,
    du
end
