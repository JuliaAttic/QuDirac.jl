# Outer Product of Two States
---

The simplest way to construct a QuDirac operator is to take the outer product of two states:

```julia
julia> k = d" 1/√2 * (| 0,0 > - | 1,1 >) "
Ket{KroneckerDelta,2,Float64} with 2 state(s):
  0.7071067811865475 | 0,0 ⟩
  -0.7071067811865475 | 1,1 ⟩

julia> k*k'
OuterProduct with 4 operator(s); Ket{KroneckerDelta,2,Float64} * Bra{KroneckerDelta,2,Float64}:
  0.4999999999999999 | 0,0 ⟩⟨ 0,0 |
  -0.4999999999999999 | 0,0 ⟩⟨ 1,1 |
  -0.4999999999999999 | 1,1 ⟩⟨ 0,0 |
  0.4999999999999999 | 1,1 ⟩⟨ 1,1 |
```

Specifically, an outer product of two states will yield an instance of the `OuterProduct` type, as can be seen above. 

The `OuterProduct` type is a lazy representation of the outer product of two states - it simply stores a reference to 
the two factor states, and uses the state information to behave like an operator. Thus, the `OuterProduct` type acts as 
a *view* onto the factor states.

With the exception of scaling functions, most mutating functions are not defined on `OuterProduct`. The non-mutating versions of these functions will work, however, by converting the operator into the more flexible `OpSum` type. The `OpSum` type represents a sum of operators rather than an outer product of states, and is *not* a view. 

---
# Scalar Multiplication
---

Like states, operators can be multiplied by a scalar: 

```julia
julia> op = d" | 'a' >< 'b' | "
OuterProduct with 1 operator(s); Ket{KroneckerDelta,1,Int64} * Bra{KroneckerDelta,1,Int64}:
  1 | 'a' ⟩⟨ 'b' |

julia> im * op
OuterProduct with 1 operator(s); Ket{KroneckerDelta,1,Int64} * Bra{KroneckerDelta,1,Int64}:
  0 + 1im | 'a' ⟩⟨ 'b' |

julia> op/2
OuterProduct with 1 operator(s); Ket{KroneckerDelta,1,Int64} * Bra{KroneckerDelta,1,Int64}:
  0.5 | 'a' ⟩⟨ 'b' |
```

---
# Addition and Subtraction
---

Operators can be added and subtracted just like states:

```julia
julia> op + op
OpSum{KroneckerDelta,1,Int64} with 1 operator(s):
  2 | 'a' ⟩⟨ 'b' |

julia> d" 1/√2 * (op - | 0 >< 1 |) "
OpSum{KroneckerDelta,1,Float64} with 2 operator(s):
  0.7071067811865475 | 'a' ⟩⟨ 'b' |
  -0.7071067811865475 | 0 ⟩⟨ 1 |
```

As you can see, generic sums of operators are represented using the `OpSum` type rather than the `OuterProduct` type.

---
# Normalization
---

Similarly to states, one can normalize operators using the `normalize` and `normalize!` functions:

```julia
julia> op = normalize(sum(i -> d"| i >< i^2 |", 1:5))
OpSum{KroneckerDelta,1,Float64} with 5 operator(s):
  0.4472135954999579 | 5 ⟩⟨ 25 |
  0.4472135954999579 | 3 ⟩⟨ 9 |
  0.4472135954999579 | 4 ⟩⟨ 16 |
  0.4472135954999579 | 1 ⟩⟨ 1 |
  0.4472135954999579 | 2 ⟩⟨ 4 |

julia> norm(op)
0.9999999999999999
```

For QuDirac operators, the `norm` function specifically computes the [Frobenius norm](http://en.wikipedia.org/wiki/Matrix_norm#Frobenius_norm).

---
# Getting an Operator's Adjoint
---

Take the conjugate transpose of an operator, simply call `ctranspose` on it:

```julia
julia> A = normalize(d"(1+3im)| 1 >< 2 | + (5-2im)| 3 >< 4 |")
OpSum{KroneckerDelta,1,Complex{Float64}} with 2 operator(s):
  0.8006407690254357 - 0.32025630761017426im | 3 ⟩⟨ 4 |
  0.16012815380508713 + 0.48038446141526137im | 1 ⟩⟨ 2 |

julia> A'
DualOpSum{KroneckerDelta,1,Complex{Float64}} with 2 operator(s):
  0.8006407690254357 + 0.32025630761017426im | 4 ⟩⟨ 3 |
  0.16012815380508713 - 0.48038446141526137im | 2 ⟩⟨ 1 |
```

Like the conjugate transpose of a Ket is a Bra, the conjugate transpose of an `OpSum` is a `DualOpSum`.
The `DualOpSum` type is a *view* on the original operator, so mutating a `DualOpSum` will mutate the underlying `OpSum` (one can explicitly make a copy via the `copy` function). 

The dual of an `OuterProduct` is simply an `OuterProduct`, and is still a view on the original factor states.

---
# Inner Product
---

Use the `*` function to take the inner product of states/operators: 

```julia
julia> k = d" 1/√2 * (| 0,0 > - | 1,1 >) "; P = k*k'
OuterProduct with 4 operator(s); Ket{KroneckerDelta,2,Float64} * Bra{KroneckerDelta,2,Float64}:
  0.4999999999999999 | 0,0 ⟩⟨ 0,0 |
  -0.4999999999999999 | 0,0 ⟩⟨ 1,1 |
  -0.4999999999999999 | 1,1 ⟩⟨ 0,0 |
  0.4999999999999999 | 1,1 ⟩⟨ 1,1 |

julia> d" P * | 1,1 > "
Ket{KroneckerDelta,2,Float64} with 2 state(s):
  -0.4999999999999999 | 0,0 ⟩
  0.4999999999999999 | 1,1 ⟩

julia> d" < 0,0 | * P * | 1,1 > "
-0.4999999999999999

julia> d" P * (| 1,1 >< 0,0 |) "
OuterProduct with 2 operator(s); Ket{KroneckerDelta,2,Float64} * Bra{KroneckerDelta,2,Int64}:
  -0.4999999999999999 | 0,0 ⟩⟨ 0,0 |
  0.4999999999999999 | 1,1 ⟩⟨ 0,0 |
```

Thus, expectation values are naturally obtained in this manner:

```julia
julia> d" < 1,1 | * P * | 1,1 > "
0.4999999999999999

julia> k' * P * k
0.9999999999999997
```

Like states, operator inner products have support for both lazy and custom evaluation rules.
See the [Working with Inner Products](inner_products.md) section for information
regarding these features.

---
# Acting an operator on a specific Ket factor
---

One can use the `act_on` function to apply an operator to a specific factor of a Ket:

```julia
julia> a = sum(i-> d"(√i)| i-1 >< i |", 1:5) # lowering operator
OpSum{KroneckerDelta,1,Float64} with 5 operator(s):
  2.23606797749979 | 4 ⟩⟨ 5 |
  2.0 | 3 ⟩⟨ 4 |
  1.0 | 0 ⟩⟨ 1 |
  1.4142135623730951 | 1 ⟩⟨ 2 |
  1.7320508075688772 | 2 ⟩⟨ 3 |

julia> k = d" | 1,2,3 > + 2| 3,5,1> "
Ket{KroneckerDelta,3,Int64} with 2 state(s):
  1 | 1,2,3 ⟩
  2 | 3,5,1 ⟩

julia> act_on(a, k, 2)
Ket{KroneckerDelta,3,Float64} with 2 state(s):
  4.47213595499958 | 3,4,1 ⟩
  1.4142135623730951 | 1,1,3 ⟩
```

---
# Tensor Product
---

Unlike states, one does not take the tensor product of operators using `*`; that function is
already used for inner products. Thus, one must use the `tensor` function:

```julia
julia> op = d" 1/√2 * (| 'a' >< 'b' | + | 'c' >< 'd' |)"
OpSum{KroneckerDelta,1,Float64} with 2 operator(s):
  0.7071067811865475 | 'a' ⟩⟨ 'b' |
  0.7071067811865475 | 'c' ⟩⟨ 'd' |

julia> tensor(op,op,op)
OpSum{KroneckerDelta,3,Float64} with 8 operator(s):
  0.3535533905932737 | 'c','c','a' ⟩⟨ 'd','d','b' |
  0.3535533905932737 | 'c','a','a' ⟩⟨ 'd','b','b' |
  0.3535533905932737 | 'a','c','a' ⟩⟨ 'b','d','b' |
  0.3535533905932737 | 'a','a','c' ⟩⟨ 'b','b','d' |
  0.3535533905932737 | 'a','c','c' ⟩⟨ 'b','d','d' |
  0.3535533905932737 | 'c','c','c' ⟩⟨ 'd','d','d' |
  0.3535533905932737 | 'a','a','a' ⟩⟨ 'b','b','b' |
  0.3535533905932737 | 'c','a','c' ⟩⟨ 'd','b','d' |
```

---
# Trace and Partial Trace
---

To take the trace of an operator, simply use the `trace` function:

```julia
julia> k = normalize(sum(ket, 0:5))
Ket{KroneckerDelta,1,Float64} with 6 state(s):
  0.4082482904638631 | 0 ⟩
  0.4082482904638631 | 2 ⟩
  0.4082482904638631 | 3 ⟩
  0.4082482904638631 | 5 ⟩
  0.4082482904638631 | 4 ⟩
  0.4082482904638631 | 1 ⟩

julia> trace(k*k')
1.0000000000000002
```

Note that the trace calculation is defined to sum the coefficients for which the Ket label and Bra label are equal (analogous to summing over the diagonal of a matrix representation).

The partial trace of an operator can be taken using the `ptrace` function:

```julia
julia> bell = d" 1/√2 * (| 'b','a' > + | 'a','b' >) "
Ket{KroneckerDelta,2,Float64} with 2 state(s):
  0.7071067811865475 | 'a','b' ⟩
  0.7071067811865475 | 'b','a' ⟩

julia> dense = bell * bell'
OuterProduct with 4 operator(s); Ket{KroneckerDelta,2,Float64} * Bra{KroneckerDelta,2,Float64}:
  0.4999999999999999 | 'b','a' ⟩⟨ 'b','a' |
  0.4999999999999999 | 'b','a' ⟩⟨ 'a','b' |
  0.4999999999999999 | 'a','b' ⟩⟨ 'b','a' |
  0.4999999999999999 | 'a','b' ⟩⟨ 'a','b' |

julia> ptrace(dense,1) # trace over the 1st subsystem
OpSum{KroneckerDelta,1,Float64} with 2 operator(s):
  0.4999999999999999 | 'b' ⟩⟨ 'b' |
  0.4999999999999999 | 'a' ⟩⟨ 'a' |

julia> purity(ans) # get the purity of the previous result
0.4999999999999998
```
