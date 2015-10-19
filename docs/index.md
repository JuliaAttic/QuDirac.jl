# DiracNotation.jl

DiracNotation.jl is a [Julia](http://julialang.org/) library for using Dirac notation to perform 
quantum mechanics computations. 

Current Features:

- Implementations of state types (`Ket`,`Bra`), and a variety of operator types (`OuterSum`,`OuterProduct`)
- Treat states and operators as map-like data structures, enabling label-based analysis for spectroscopy purposes
- Implementation of common operations like partial trace (`ptrace`) and partial transpose (`ptranspose`)
- Support for abstract/undefined inner products
- User-definable custom inner product rules
- Subspace selection/transformation via functions on state labels and coefficients:
    - `xsubspace` allows easy selection of excitation subspaces of states and operators
    - `permute` and `switch` allows generic permutation of factor labels for states
    - `filter`/`filter!` are supported on both the labels and coefficients of operators/states
    - Mapping functions (`map`/`maplabels`/`mapcoeffs`) for applying arbitrary functions to labels and coefficients
- Functional generation of operators using `@defop` and `@rep_op`
- `d" ... "` literals for natural Dirac notation input syntax