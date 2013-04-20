module DualNumbers
  importall Base

  include("dual.jl")

  export
    Dual,
    Dual128,
    Dual64,
    DualPair,
    dual,
    dual128,
    dual64,
    isdual,
    dual_show
end
