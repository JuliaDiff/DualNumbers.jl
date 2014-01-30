### Scope of DualNumbers.jl

The `DualNumbers` package defines the `Dual` type to represent dual numbers and 
supports standard mathematical operations on them. Conversions and promotions 
are defined to allow performing operations on combinations of dual numbers with 
predefined Julia numeric types.

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

## Supported functions

We aim for complete support for `Dual` types for numerical functions within Julia's 
`Base`. Currently, basic mathematical operations and trigonometric functions are
supported.


The following functions are specific to dual numbers:
* `dual`,
* `dual128`,
* `dual64`,
* `epsilon`,
* `isdual`,
* `dual_show`,
* `conjdual`,
* `absdual`,
* `abs2dual`.

In some cases the mathematical definition of functions of ``Dual`` numbers
is in conflict with their use as a drop-in replacement for calculating
numerical derivatives, for example, ``conj``, ``abs`` and ``abs2``. In these
cases, we choose to follow the rule ``f(x::Dual) = Dual(f(real(x)),epsilon(x)*f'(real(x)))``,
where ``f'`` is the derivative of ``f``. The mathematical definitions are
available using the functions with the suffix ``dual``.


### A walk-through example

The example below demonstrates basic usage of dual numbers by employing them to 
perform automatic differentiation. The code for this example can be found in 
`test/automatic_differentiation_test.jl`.

First install the package by using the Julia package manager:

    Pkg.update()
    Pkg.add("DualNumbers")
    
Then make the package available via

    using DualNumbers

Use the `dual()` function to define the dual number `2+1*du`:

    x = dual(2, 1)

Define a function that will be differentiated, say

    f(x) = x^3

Perform automatic differentiation by passing the dual number `x` as argument to 
`f`:

    y = f(x)

Use the functions `real()` and `epsilon()` to get the real and imaginary (dual) 
parts of `x`, respectively:

    println("f(x) = x^3")
    println("f(2) = ", real(y))
    println("f'(2) = ", epsilon(y))

[![Build Status](https://travis-ci.org/scidom/DualNumbers.jl.png)](https://travis-ci.org/scidom/DualNumbers.jl)
