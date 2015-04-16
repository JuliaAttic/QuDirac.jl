# QuDirac.jl

QuDirac.jl is a [Julia](http://julialang.org/) library for using Dirac notation to perform 
quantum mechanics computations. 

Current Features:

- Implementations of state types (`Ket`,`Bra`), and a variety of operator types (`GenericOp`,`OuterProduct`)
- Support for abstract/undefined inner products
- Users can define custom inner product rules
- Subspace selection/transformation via functions on state labels and coefficients:
    - `xsubspace` allows easy selection of excitation subspaces of states and operators
    - `permute` and `switch` allows generic permutation of factor labels for states
    - `filter`/`filter!` are supported on both the labels and coefficients of operators/states
    - Arbitrary mapping functions (`map`/`maplabels`/`mapcoeffs`) are also provided for applying functions to labels and coefficients
- Functional generation of operator representations
- `d" ... "` syntax for natural Dirac notation input