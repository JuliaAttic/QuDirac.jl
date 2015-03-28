# QuDirac.jl

QuDirac.jl is a [Julia](http://julialang.org/) library for using Dirac notation to perform 
quantum mechanics computations. 

## Features

Below are some toy examples; actual documentation is coming soon! 

#### State types (`Ket`,`Bra`) and Operator types (`GenericOp`,`Projector`)

```julia
julia> bell = 1/√2 * (ket(0,0) + ket(1,1))
Ket{Orthonormal,2} with 2 state(s):
  0.7071067811865475 | 1,1 ⟩
  0.7071067811865475 | 0,0 ⟩

julia> bell'
Bra{Orthonormal,2} with 2 state(s):
  0.7071067811865475 ⟨ 1,1 |
  0.7071067811865475 ⟨ 0,0 |

julia> op = bell * bell'
Projector{Orthonormal,2} with 4 operator(s):
  0.4999999999999999 | 1,1 ⟩⟨ 1,1 |
  0.4999999999999999 | 1,1 ⟩⟨ 0,0 |
  0.4999999999999999 | 0,0 ⟩⟨ 1,1 |
  0.4999999999999999 | 0,0 ⟩⟨ 0,0 |

julia> ptrace(op, 1)
GenericOp{Orthonormal,1} with 2 operator(s):
  0.4999999999999999 | 0 ⟩⟨ 0 |
  0.4999999999999999 | 1 ⟩⟨ 1 |
```

#### User-extensible inner product rules + support for abstract inner products

Abstract inner product example:

```julia
julia> QuDirac.set_default_inner(AbstractInner);

julia> k = 1/√2 * (ket('a') + ket('b'))
Ket{AbstractInner,1} with 2 state(s):
  0.7071067811865475 | 'b' ⟩
  0.7071067811865475 | 'a' ⟩

julia> bra('a') * k
((0.7071067811865475 * ⟨ 'a' | 'b' ⟩) + (0.7071067811865475 * ⟨ 'a' | 'a' ⟩))
```

Custom inner product example:

```julia
julia> abstract MyInner <: AbstractInner

julia> QuDirac.inner_rule{T<:MyInner}(::Type{T}, ktlabel, brlabel) = sqrt(ktlabel[1]+brlabel[1])
inner_rule (generic function with 3 methods)

julia> QuDirac.set_default_inner(MyInner);

julia> bra(π) * ket(e) # eval ⟨ π | e ⟩ with MyInner rule -> sqrt(π + e)
2.420717761749361
```

#### Generate operator representations from user-defined functions

Functional generation of a lowering operator: 

```julia
julia> k = normalize!(sum(ket, 0:4))
Ket{Orthonormal,1} with 5 state(s):
  0.4472135954999579 | 0 ⟩
  0.4472135954999579 | 2 ⟩
  0.4472135954999579 | 3 ⟩
  0.4472135954999579 | 4 ⟩
  0.4472135954999579 | 1 ⟩

julia> a = GenericOp(label->(sqrt(label[1]),{label[1]-1}), k)
GenericOp{Orthonormal,1} with 4 operator(s):
  1.0 | 0 ⟩⟨ 1 |
  2.0 | 3 ⟩⟨ 4 |
  1.7320508075688772 | 2 ⟩⟨ 3 |
  1.4142135623730951 | 1 ⟩⟨ 2 |

julia> a*k
Ket{Orthonormal,1} with 4 state(s):
  0.4472135954999579 | 0 ⟩
  0.8944271909999159 | 3 ⟩
  0.7745966692414833 | 2 ⟩
  0.6324555320336759 | 1 ⟩

julia> a*ket(4)
Ket{Orthonormal,1} with 4 state(s):
  0.0 | 0 ⟩
  2.0 | 3 ⟩
  0.0 | 2 ⟩
  0.0 | 1 ⟩
```

#### `@drc_str` for inputting expressions in natural Dirac notation

```julia
julia> drc" < 0,0 | *  (| 0,0 > + | 1,1 >)/√2 "
0.7071067811865475

julia> QuDirac.set_default_inner(AbstractInner)
AbstractInner

julia> drc" < 'a','b' | *  (| 0,0 > + | 1,1 >)/√2 "
((⟨ 'a','b' | 1,1 ⟩ + ⟨ 'a','b' | 0,0 ⟩) / 1.4142135623730951)
```

#### ...and other stuff (examples and documentation coming soon)

- `xsubspace` allows easy selection of excitation subspaces of states and operators
- `permute!`/`permute` and `switch!`/`switch` allows generic permutation of factor labels for states
- `filter`/`filter!` are supported on both the labels and coefficients of operators/states
- Arbitrary `map` functions (`map`/`maplabels`/`mapcoeffs`) are also provided for applying functions to labels and coefficients
- Wave function generation from QuDirac states by using label functionals
