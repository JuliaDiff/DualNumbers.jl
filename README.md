### Scope of DualNumbers.jl

The `DualNumbers` package defines the `Dual` type to represent dual numbers and 
supports the standard mathematical operations on them. Conversions and 
promotions are defined to allow performing operations on combinations of 
dual numbers with predefined Julia numeric types.

Dual numbers extend the real numbers, similar to complex numbers. They adjoin a 
new element `du` such that `du*du=0`, in a similar way that complex numbers 
adjoin the imaginary unit `i` with the property `i*i=-1`. So the typical 
representation of a dual number takes the form `x+y*du`, where `x` and `y` are 
real numbers.

Apart from their mathematical role in algebraic and differential geometry (they 
are mainly interpreted as angles between lines), they also find applications in 
physics (the real part of a dual represents the bosonic direction, while the 
imaginary part represents the fermionic direction), in screw theory, in motor 
and spatial vector algebra, and in computer science due to its relation with the 
forward mode of automatic differentiation.

## List of functions for numbers of `Dual` types

All the typical mathematical operations are supported. These include
* `real`,
* `imag`,
* `real_valued`,
* `integer_valued`,
* `isfinite`,
* `reim`,
* `show`,
* `showcompact`,
* `read`,
* `write`,
* `==`,
* `isequal`,
* `hash`,
* `conj`,
* `abs`,
* `abs2`,
* `inv`,
* `+`,
* `-`,
* `*`,
* `/`,
* `sqrt`,
* `cbrt`,
* `^`,
* `exp`,
* `log`,
* `log2`,
* `log10`,
* `sin`,
* `cos`,
* `tan`,
* `asin`,
* `acos`,
* `atan`,
* `sinh`,
* `cosh`,
* `tanh`,
* `asinh`,
* `acosh`,
* `atanh`.

The following functions are specific to dual numbers:
* `dual`,
* `dual128`,
* `dual64`,
* `isdual`,
* `dual_show`.

### A walk-through example

The example below demonstrates basic usage of dual numbers by using them to 
perform automatic differentiation. The code for this example can be found in 
`test/automatic_differentiation_test.jl`.

First install the package by using the Julia package manager:
    Pkg.update()
    Pkg.add("DualNumbers")
    
Then make the package available via
    using DualNumbers

Use the `dual()` function to define the dual number `2+1*du`:
    x = dual(2., 1.)

Define a function that will be differentiated, say
    f(x) = x^3

Perform automatic differentiation by simply passing the dual number `x` as 
argument to `f`    
    y = f(x)

Use the functions `real()` and `imag()` to get the real and imaginary (dual) 
parts of `x`, respectively:  
    println("f(x) = x^3")
    println("f(2) = ", real(y))
    println("f'(2) = ", imag(y))
