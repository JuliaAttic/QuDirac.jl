# QuDirac.jl

QuDirac.jl is a [Julia](http://julialang.org/) library for using Dirac notation to perform 
quantum mechanics computations. 

Documentation is coming soon!

Features:

- Implementations of state types (`Ket`,`Bra`), and a variety of operator types (`GenericOp`,`Projector`)
- Extensible inner product rules overloaded on types `T<:AbstractInner` (undefined inner products are lazy evaluated)
- Subspace selection and transformation via user-defined functions on labels and coefficients
- Generate operator representations from user-defined functions
- Wave function generation from QuDirac states by using label functionals