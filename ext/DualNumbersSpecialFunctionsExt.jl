module DualNumbersSpecialFunctionsExt

using DualNumbers
using DualNumbers: value, epsilon, to_nanmath
using NaNMath
using Calculus
using SpecialFunctions

for (funsym, expr) in Calculus.symbolic_derivatives_1arg()
    if isdefined(SpecialFunctions, funsym) && !isdefined(Base, funsym)
        @eval function SpecialFunctions.$(funsym)(z::Dual)
            x = value(z)
            xp = epsilon(z)
            Dual($(funsym)(x),xp*$expr)
        end
    end
    # extend corresponding NaNMath methods
	if funsym in (:lgamma,)
	    funsym = Expr(:.,:NaNMath,Base.Meta.quot(funsym))
	    @eval function $(funsym)(z::Dual)
	        x = value(z)
	        xp = epsilon(z)
	        Dual($(funsym)(x),xp*$(to_nanmath(expr)))
	    end
	end
end


end
