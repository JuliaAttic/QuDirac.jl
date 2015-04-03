## Outer Product of Two States
---

The simplest way to construct a QuDirac operator is to take the outer product of two states:

```
julia> k = 1/√2 * (ket(0,0) - ket(1,1))
Ket{Orthonormal,2} with 2 state(s):
  -0.7071067811865475 | 1,1 ⟩
  0.7071067811865475 | 0,0 ⟩

julia> k*k' # | k ⟩⟨ k |
Projector{Orthonormal,2} with 4 operator(s):
  0.4999999999999999 | 1,1 ⟩⟨ 1,1 |
  -0.4999999999999999 | 1,1 ⟩⟨ 0,0 |
  -0.4999999999999999 | 0,0 ⟩⟨ 1,1 |
  0.4999999999999999 | 0,0 ⟩⟨ 0,0 |
```

Specifically, an outer product of two states will yield an instance of the `Projector` type, as can be seen above. 

The `Projector` type is a lazy representation of the outer product of two states - it simply stores a reference to 
the two factor states, and uses the state information to behave like an operator. Thus, the `Projector` type acts as 
a *view* onto the factor states (see [Views vs. Copies in QuDirac](view_copy.md)). 

Outer products are not the only way to construct operators - certain kinds of operators can be constructed simply by 
defining their action on basis Kets. See [Functionally Defining Operators](func_op_def.md) for details.

---
## Scalar Multiplication
---

Like states, operators can be multiplied by a scalar: 

```
julia> op = ket('a') * bra('b')
Projector{Orthonormal,1} with 1 operator(s):
  1 | 'a' ⟩⟨ 'b' |

julia> im * op
Projector{Orthonormal,1} with 1 operator(s):
  0 + 1im | 'a' ⟩⟨ 'b' |

julia> op/2
Projector{Orthonormal,1} with 1 operator(s):
  0.5 | 'a' ⟩⟨ 'b' |
```
---
## Addition and Subtraction
---

Operators can be added and subtracted just like states:

```
julia> op + op
GenericOp{Orthonormal,1} with 1 operator(s):
  2 | 'a' ⟩⟨ 'b' |

julia> gop = 1/√2 * (op - (ket(0) * bra(1)))
GenericOp{Orthonormal,1} with 2 operator(s):
  0.7071067811865475 | 'a' ⟩⟨ 'b' |
  -0.7071067811865475 | 0 ⟩⟨ 1 |
```

Note that addition and subtraction of `Projector`s results in `GenericOp`s. While the former
simply represents an outer product, the latter can more generally represent a *sum* of operators, 
and is no longer a view. 

---
## Normalization
---

Similarly to states, the Hilbert–Schmidt norm can be computed on an operator using the `norm` function:

```
julia> norm(gop) # using gop from the previous example
0.9999999999999999
```

Of course, the `normalize` and `normalize!` functions are provided for operators:

```
julia> normalize!(sum(i -> ket(i) * bra(i^2), 1:5))
GenericOp{Orthonormal,1} with 5 operator(s):
  0.4472135954999579 | 1 ⟩⟨ 1 |
  0.4472135954999579 | 4 ⟩⟨ 16 |
  0.4472135954999579 | 5 ⟩⟨ 25 |
  0.4472135954999579 | 2 ⟩⟨ 4 |
  0.4472135954999579 | 3 ⟩⟨ 9 |
```

---
## Getting an Operator's Adjoint
---

Take the conjugate transpose of an operator, simply use the `ctranspose` function:

```
julia> gop #from the earlier example
GenericOp{Orthonormal,1} with 2 operator(s):
  0.7071067811865475 | 'a' ⟩⟨ 'b' |
  -0.7071067811865475 | 0 ⟩⟨ 1 |

julia> gop'
DualOp{Orthonormal,1} with 2 operator(s):
  0.7071067811865475 | 'b' ⟩⟨ 'a' |
  -0.7071067811865475 | 1 ⟩⟨ 0 |

julia> gop'' == gop
true
```

---
## Tensor Product
---

---
## Inner Product
---

---
## Expectation values
---

---
## Trace and Partial Trace
---

