[![Build Status](https://travis-ci.org/JuliaQuantum/QuDirac.jl.svg?branch=master)](https://travis-ci.org/JuliaQuantum/QuDirac.jl)
# QuDirac.jl

QuDirac.jl is a [Julia](http://julialang.org/) library for using Dirac notation to perform 
quantum mechanics computations. 

Documentation for the current release version (v0.1) can be found [**here**](http://qudiracjl.readthedocs.org/en/release-0.1/).

## Installation

To install QuDirac.jl, you should have [a working build of Julia v0.3](https://github.com/JuliaLang/julia#source-download-and-compilation). Then, you can grab QuDirac.jl via the package manager:

```julia
julia> Pkg.add("QuDirac")
```

## Features

These are toy examples for demoing features. See [below for more involved examples](https://github.com/JuliaQuantum/QuDirac.jl#examples).

#### Ket, Bra, and Operator types

```julia
julia> bell = d" 1/√2 * (| 0,0 > + | 1,1 >) "
Ket{KronDelta,2,Float64} with 2 state(s):
  0.7071067811865475 | 0,0 ⟩
  0.7071067811865475 | 1,1 ⟩

julia> bell'
Bra{KronDelta,2,Float64} with 2 state(s):
  0.7071067811865475 ⟨ 0,0 |
  0.7071067811865475 ⟨ 1,1 |

julia> ptrace(bell * bell', 1)
OpSum{KronDelta,1,Float64} with 2 operator(s):
  0.4999999999999999 | 0 ⟩⟨ 0 |
  0.4999999999999999 | 1 ⟩⟨ 1 |
```

#### Support for undefined inner products

```julia
# tells QuDirac to use the rule for undefined inner products
julia> default_inner(UndefinedInner)
INFO: QuDirac's default inner product type is currently UndefinedInner.

julia> d" < 0,0 | *  (| 0,0 > + | 1,1 >)/√2 "
((⟨ 0,0 | 0,0 ⟩ + ⟨ 0,0 | 1,1 ⟩) / 1.4142135623730951)

julia> s = d" (e^( < 1,2 | 3,4 > ) + < 5,6 | 7,8 > * im)^4 "
(((exp(⟨ 1,2 | 3,4 ⟩)) + (⟨ 5,6 | 7,8 ⟩ * im))^4)

julia> inner_eval((b, k) -> sum(k) - sum(b), s)
8.600194553751864e6 + 2.5900995362955774e6im
```

#### Custom inner product rules

```julia
julia> @def_inner MyInner Float64
INFO: MyInner is now defined as an inner product type.
INFO: Inner products using the MyInner type should return values of type Float64.

julia> MyInner(a::Float64, b::Float64) = sqrt(a+b)
MyInner (constructor with 2 methods)

julia> MyInner(a, b) = MyInner(float(a), float(b))
MyInner (constructor with 3 methods)

julia> default_inner(MyInner)
INFO: QuDirac's default inner product type is currently MyInner.

julia> d" < π | e > " # sqrt(π + e)
2.420717761749361

julia> d" < 4.4,5 | 2,3.42 > " # sqrt(4.4 + 2) * sqrt(5 + 3.42)
7.340844638050856
```

#### Functional operator construction

```julia
# define a₂ on Ket
julia> @def_op " a₂ | x, y, z > = √y * | x, y - 1, z > "
a₂ (generic function with 1 method)

# define a₂ on Bra
julia> @def_op " < x, y, z | a₂ = √(y + 1) * < x, y + 1, z | "
a₂ (generic function with 2 methods)

julia> d" a₂ * | 3,5,5 > "
Ket{KronDelta,3,Float64} with 1 state(s):
  2.23606797749979 | 3,4,5 ⟩

julia> d" a₂' * | 3,4,5 > "
Ket{KronDelta,3,Float64} with 1 state(s):
  2.23606797749979 | 3,5,5 ⟩

julia> d" < 3,4,5 | * a₂ * | 3,5,5 > "
2.23606797749979

# Hadamard operator
julia> @rep_op " H | n > = 1/√2 * ( | 0 > + (-1)^n *| 1 > ) " 0:1;

julia> H
OpSum{KronDelta,1,Float64} with 4 operator(s):
  0.7071067811865475 | 1 ⟩⟨ 0 |
  0.7071067811865475 | 0 ⟩⟨ 0 |
  0.7071067811865475 | 0 ⟩⟨ 1 |
  -0.7071067811865475 | 1 ⟩⟨ 1 |
```

#### ...and other stuff

- Implementation of common operations like partial trace (`ptrace`) and partial transpose (`ptranspose`)
- Treat states and operators as map-like data structures, enabling label-based analysis for spectroscopy purposes
- `xsubspace` allows easy selection of excitation subspaces of states and operators
- `permute` and `switch` allows generic permutation of factor labels for states
- `filter`/`filter!` for the filtering out component states/operators via predicate functions
- Arbitrary mapping functions (`map`/`maplabels`/`mapcoeffs`) for applying functions to labels and coefficients

## Examples

There are currently two example files, `qho.jl` and `randwalk.jl`. The former implements methods for plotting quantum harmonic oscillator wave functions using Plotly. The latter is a simple implementation of a quantum random 
walk.

To run the examples, one can do the following (using `qho.jl` as an example):

```julia
julia> cd(Pkg.dir("QuDirac"))
julia> include("examples/qho.jl")
```
