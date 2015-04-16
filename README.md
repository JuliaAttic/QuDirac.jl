[![Build Status](https://travis-ci.org/JuliaQuantum/QuDirac.jl.svg?branch=master)](https://travis-ci.org/JuliaQuantum/QuDirac.jl)
# QuDirac.jl

QuDirac.jl is a [Julia](http://julialang.org/) library for using Dirac notation to perform 
quantum mechanics computations. 

## Features

Below are some toy examples; actual documentation is coming soon! 

#### State types (`Ket`,`Bra`) and Operator types (`GenericOp`,`OuterProduct`)

```julia
julia> bell = 1/√2 * (ket(0,0) + ket(1,1))
Ket{KroneckerDelta,2,Float64} with 2 state(s):
  0.7071067811865475 | 0,0 ⟩
  0.7071067811865475 | 1,1 ⟩

julia> bell'
Bra{KroneckerDelta,2,Float64} with 2 state(s):
  0.7071067811865475 ⟨ 0,0 |
  0.7071067811865475 ⟨ 1,1 |

julia> op = bell * bell'
OuterProduct with 4 operator(s); Ket{KroneckerDelta,2,Float64} * Bra{KroneckerDelta,2,Float64}:  
  0.4999999999999999 | 0,0 ⟩⟨ 0,0 |
  0.4999999999999999 | 0,0 ⟩⟨ 1,1 |
  0.4999999999999999 | 1,1 ⟩⟨ 0,0 |
  0.4999999999999999 | 1,1 ⟩⟨ 1,1 |

julia> ptrace(op, 1)
GenericOp{KroneckerDelta,1,Float64} with 2 operator(s):
  0.4999999999999999 | 0 ⟩⟨ 0 |
  0.4999999999999999 | 1 ⟩⟨ 1 |
```

#### Support for abstract inner products

```julia
# tells QuDirac to use the rule for undefined inner products
julia> default_inner(UndefinedInner());
 
julia> k = 1/√2 * (ket('a') + ket('b'))
Ket{UndefinedInner,1,Float64} with 2 state(s):
  0.7071067811865475 | 'b' ⟩
  0.7071067811865475 | 'a' ⟩

julia> bra('a') * k
((0.7071067811865475 * ⟨ 'a' | 'b' ⟩) + (0.7071067811865475 * ⟨ 'a' | 'a' ⟩))
```

#### Custom inner product rules

```julia
julia> immutable MyInner <: AbstractInner end

julia> QuDirac.inner_rule(::MyInner, ktlabel, brlabel) = sqrt(ktlabel[1]+brlabel[1])
inner_rule (generic function with 3 methods)

julia> default_inner(MyInner());

julia> bra(π) * ket(e) # eval ⟨ π | e ⟩ with MyInner rule -> sqrt(π + e)
2.420717761749361
```

#### Generate operator representations from user-defined functions

For example, here's the functional generation of a lowering operator: 

```julia
# State in a number basis
julia> k = normalize(sum(ket, 0:4))
Ket{KroneckerDelta,1,Float64} with 5 state(s):
  0.4472135954999579 | 4 ⟩
  0.4472135954999579 | 3 ⟩
  0.4472135954999579 | 2 ⟩
  0.4472135954999579 | 0 ⟩
  0.4472135954999579 | 1 ⟩

julia> f(label) = (sqrt(label), label-1)
f (generic function with 1 method)

# generate â in the basis of k
julia> â = func_op(f, k, 1)
GenericOp{KroneckerDelta,1,Float64} with 4 operator(s):
  2.0 | 3 ⟩⟨ 4 |
  1.0 | 0 ⟩⟨ 1 |
  1.4142135623730951 | 1 ⟩⟨ 2 |
  1.7320508075688772 | 2 ⟩⟨ 3 |

julia> â*k
Ket{KroneckerDelta,1,Float64} with 4 state(s):
  0.8944271909999159 | 3 ⟩
  0.7745966692414833 | 2 ⟩
  0.4472135954999579 | 0 ⟩
  0.6324555320336759 | 1 ⟩

julia> â*ket(4)
Ket{KroneckerDelta,1,Float64} with 1 state(s):
  2.0 | 3 ⟩
```

#### `d" ... "` syntax for natural Dirac notation input

```julia
julia> d" < 0,0 | *  (| 0,0 > + | 1,1 >)/√2 "
0.7071067811865475

julia> default_inner(UndefinedInner());

julia> d" < 'a','b' | *  (| 0,0 > + | 1,1 >)/√2 "
((⟨ 'a','b' | 0,0 ⟩ + ⟨ 'a','b' | 1,1 ⟩) / 1.4142135623730951)

julia> d"""
       ψ = 1/√2 * (| 0,0 > + | 1,1 >)
       a = purity(ptrace(ψ*ψ', 2))
       ϕ = normalize!( 1/5 * < 0 | + 4/5 * < 1 | )
       result = normalize!(a * act_on(ϕ, ψ, 2))
       """

julia> result
Ket{KroneckerDelta,1,Float64} with 2 state(s):
  0.24253562503633297 | 0 ⟩
  0.9701425001453319 | 1 ⟩
```

#### ...and other stuff (examples and documentation coming soon)

- `xsubspace` allows easy selection of excitation subspaces of states and operators
- `permute` and `switch` allows generic permutation of factor labels for states
- `filter`/`filter!` are supported on both the labels and coefficients of operators/states
- Arbitrary mapping functions (`map`/`maplabels`/`mapcoeffs`) are also provided for applying functions to labels and coefficients
