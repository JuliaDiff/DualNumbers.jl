import Base.depwarn, Base.@deprecate_binding

@deprecate epsilon(d) partials(d)
@deprecate Dual128(a,b) Dual(a,b) 
@deprecate Dual64(a,b) Dual(a,b) 
@deprecate Dual32(a,b) Dual(a,b) 
@deprecate_binding ɛ Dual(false, true)
@deprecate_binding imɛ Dual(0,1)*im
