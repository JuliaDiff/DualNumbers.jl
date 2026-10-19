[![Build Status](https://github.com/JuliaDiff/DualNumbers.jl/workflows/CI/badge.svg)](https://github.com/JuliaDiff/DualNumbers.jl/actions)
[![codecov](https://codecov.io/gh/JuliaDiff/DualNumbers.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaDiff/DualNumbers.jl)


### DualNumbers.jl 的定位

`DualNumbers` 是一个 [Julia](http://julialang.org) 包，定义了表示[对偶数](https://en.wikipedia.org/wiki/Dual_number)的 `Dual` 类型，并支持对偶数的标准数学运算。包中还定义了转换和提升规则，使对偶数可以与 Julia 内置数值类型混合运算。

对偶数是对实数的扩展，方式类似复数。它添加一个新元素 `ɛ`，满足 `ɛ*ɛ=0`，就像复数添加虚数单位 `i` 并满足 `i*i=-1` 一样。对偶数的典型表示形式是 `x+y*ɛ`，其中 `x` 和 `y` 是实数。对偶数还可以进一步扩展复数：在实部和虚部各添加一个这样的新元素。

除了在代数和微分几何中的数学意义（主要用来解释直线之间的夹角）之外，对偶数在物理（对偶数的实部表示玻色方向，epsilon 部分表示费米方向）、螺旋理论、运动与空间向量代数，以及计算机科学中都有应用，后者源于它与自动微分前向模式的关系。

就前向模式自动微分而言，本包已被 [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) 取代。

`DualNumbers` 目前没有活跃的维护者。

## 支持的函数

我们的目标是为 Julia `Base` 中的数值函数提供完整的 `Dual` 类型支持。目前支持基本数学运算和三角函数。

以下函数是对偶数专用的：
* `dual`
* `realpart`
* `dualpart`
* `epsilon`
* `isdual`
* `dual_show`
* `conjdual`
* `absdual`
* `abs2dual`

对偶数函数值 `f(a+bɛ)` 由如下极限定义：

    f(a+bɛ) := f(a) + lim_{h→0} (f(a + bɛh) - f(a))/h .

对复可微函数，这等价于求导：

    f(a+bɛ) := f(a) + b f'(a) ɛ.

对不是复可微的函数，对偶部分返回该极限值，可以看作 `R²` 中的一个方向导数。

有时 `Dual` 数的数学定义会与它作为数值求导替代品的用法冲突，例如 `conj`、`abs` 和 `abs2`。数学定义通过带 `dual` 后缀的函数提供。类似地，比较运算符 `<`、`>`、`==` 被重载为只比较数值部分。

### 示例

下面的例子通过自动微分展示对偶数的基本用法。示例代码见 `test/automatic_differentiation_test.jl`。

首先用 Julia 的包管理器安装本包：

    Pkg.update()
    Pkg.add("DualNumbers")

然后加载包：

    using DualNumbers

用 `Dual()` 构造对偶数 `2+1*ɛ`：

    x = Dual(2, 1)

定义一个要求导的函数，例如：

    f(x) = x^3

把对偶数 `x` 作为参数传给 `f`，就完成了自动微分：

    y = f(x)

用 `realpart` 和 `dualpart` 分别取出 `x` 的数值部分和对偶部分：

    println("f(x) = x^3")
    println("f(2) = ", realpart(y))
    println("f'(2) = ", dualpart(y))
