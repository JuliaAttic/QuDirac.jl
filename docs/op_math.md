# Outer Product of Two States
---

The simplest way to construct a QuDirac operator is to take the outer product of two states:

```
julia> k = 1/√2 * (ket(0,0) - ket(1,1))
Ket{KroneckerDelta,2,Float64} with 2 state(s):
  -0.7071067811865475 | 1,1 ⟩
  0.7071067811865475 | 0,0 ⟩

julia> k*k'
OuterProduct{KroneckerDelta,2,Float64} with 4 operator(s):
  0.4999999999999999 | 1,1 ⟩⟨ 1,1 |
  -0.4999999999999999 | 1,1 ⟩⟨ 0,0 |
  -0.4999999999999999 | 0,0 ⟩⟨ 1,1 |
  0.4999999999999999 | 0,0 ⟩⟨ 0,0 |
```

Specifically, an outer product of two states will yield an instance of the `OuterProduct` type, as can be seen above. 

The `OuterProduct` type is a lazy representation of the outer product of two states - it simply stores a reference to 
the two factor states, and uses the state information to behave like an operator. Thus, the `OuterProduct` type acts as 
a *view* onto the factor states. This allows quick, memory-efficient construction of density operators and the like.

With the exception of scaling functions, most mutating functions are not defined on `OuterProduct`. The non-mutating versions of these functions will work, however, by converting the operator into the more flexible `GenericOp` type. This type represents a sum of operators rather than an outer product of states, and is *not* a view. 

---
# Scalar Multiplication
---

Like states, operators can be multiplied by a scalar: 

```
julia> op = ket('a') * bra('b')
OuterProduct{KroneckerDelta,1} with 1 operator(s):
  1 | 'a' ⟩⟨ 'b' |

julia> im * op
OuterProduct{KroneckerDelta,1} with 1 operator(s):
  0 + 1im | 'a' ⟩⟨ 'b' |

julia> op/2
OuterProduct{KroneckerDelta,1} with 1 operator(s):
  0.5 | 'a' ⟩⟨ 'b' |
```
---
# Addition and Subtraction
---

Operators can be added and subtracted just like states:

```
julia> op + op
GenericOp{KroneckerDelta,1,Int64} with 1 operator(s):
  2 | 'a' ⟩⟨ 'b' |

julia> gop = 1/√2 * (op - (ket(0) * bra(1)))
GenericOp{KroneckerDelta,1,Float64} with 2 operator(s):
  0.7071067811865475 | 'a' ⟩⟨ 'b' |
  -0.7071067811865475 | 0 ⟩⟨ 1 |
```

As you can see, generic sums of operators are represented using the `GenericOp` type rather than the `OuterProduct` type.

---
# Normalization
---

Similarly to states, one can normalize operators using the `normalize` and `normalize!` functions:

```
julia> op = normalize(sum(i -> ket(i) * bra(i^2), 1:5))
GenericOp{KroneckerDelta,1,Float64} with 5 operator(s):
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

```
julia> gop #from the earlier example
GenericOp{KroneckerDelta,1,Float64} with 2 operator(s):
  0.7071067811865475 | 'a' ⟩⟨ 'b' |
  -0.7071067811865475 | 0 ⟩⟨ 1 |

julia> gop'
DualOp{KroneckerDelta,1,Float64} with 2 operator(s):
  0.7071067811865475 | 'b' ⟩⟨ 'a' |
  -0.7071067811865475 | 1 ⟩⟨ 0 |
```

Like the conjugate transpose of a Ket is a Bra, the conjugate transpose of a `GenericOp` is a `DualOp`.
The `DualOp` type is a *view* on the original operator, so mutating a `DualOp` will mutate the underlying `GenericOp` (one can explicitly make a copy via the `copy` function). 

The dual of a `OuterProduct` is simply a `OuterProduct`, and is still a view on the original factor states.

---
# Inner Product
---

Use the `*` function to take the inner product of states/operators: 

```
julia> k = 1/√2 * (ket(0,0) - ket(1,1)); op = k*k'
OuterProduct{KroneckerDelta,2} with 4 operator(s):
  0.4999999999999999 | 0,0 ⟩⟨ 0,0 |
  -0.4999999999999999 | 0,0 ⟩⟨ 1,1 |
  -0.4999999999999999 | 1,1 ⟩⟨ 0,0 |
  0.4999999999999999 | 1,1 ⟩⟨ 1,1 |

julia> op * ket(1,1)
Ket{KroneckerDelta,2,Float64} with 2 state(s):
  -0.4999999999999999 | 0,0 ⟩
  0.4999999999999999 | 1,1 ⟩

julia> bra(0,0) * op * ket(1,1)
-0.4999999999999999

julia> op * (ket(1,1) * bra(0,0))
OuterProduct{KroneckerDelta,2} with 2 operator(s):
  -0.4999999999999999 | 0,0 ⟩⟨ 0,0 |
  0.4999999999999999 | 1,1 ⟩⟨ 0,0 |
```

Thus, expectation values are naturally obtained in this manner:

```
julia> bra(1,1)*op*ket(1,1)
0.4999999999999999

julia> k'*op*k
0.9999999999999997
```

Like states, operator inner products have support for both lazy and custom evaluation rules.
See the [Working with Inner Products](inner_products.md) section for information
regarding these features.

---
# Acting an operator on a specific Ket factor
---

One can use the `act_on` function to apply an operator to a specific factor of a Ket:

```
julia> a = sum(i->sqrt(i) * ket(i-1) * bra(i), 1:5) # lowering operator
GenericOp{KroneckerDelta,1,Float64} with 5 operator(s):
  2.23606797749979 | 4 ⟩⟨ 5 |
  2.0 | 3 ⟩⟨ 4 |
  1.0 | 0 ⟩⟨ 1 |
  1.4142135623730951 | 1 ⟩⟨ 2 |
  1.7320508075688772 | 2 ⟩⟨ 3 |

julia> k = ket(1,2,3) + 2*ket(3,5,1)
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

```
julia> tensor(op,op,op,op) # op from previous example
OuterProduct{KroneckerDelta,8} with 256 operator(s):
  0.06249999999999996 | 1,1,1,1,0,0,0,0 ⟩⟨ 1,1,1,1,0,0,0,0 |
  -0.06249999999999996 | 1,1,1,1,0,0,0,0 ⟩⟨ 0,0,0,0,1,1,0,0 |
  -0.06249999999999996 | 1,1,1,1,0,0,0,0 ⟩⟨ 1,1,1,1,1,1,0,0 |
  0.06249999999999996 | 1,1,1,1,0,0,0,0 ⟩⟨ 1,1,0,0,0,0,1,1 |
  -0.06249999999999996 | 0,0,0,0,1,1,0,0 ⟩⟨ 1,1,1,1,0,0,0,0 |
  0.06249999999999996 | 0,0,0,0,1,1,0,0 ⟩⟨ 0,0,0,0,1,1,0,0 |
  ⁞

julia> tensor(op, ket('a')*bra('b'))
OuterProduct{KroneckerDelta,3} with 4 operator(s):
  0.4999999999999999 | 0,0,'a' ⟩⟨ 0,0,'b' |
  -0.4999999999999999 | 0,0,'a' ⟩⟨ 1,1,'b' |
  -0.4999999999999999 | 1,1,'a' ⟩⟨ 0,0,'b' |
  0.4999999999999999 | 1,1,'a' ⟩⟨ 1,1,'b' |
```

---
# Trace and Partial Trace
---

To take the trace of an operator, simply use the `trace` function:

```
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

The partial trace of an operator can be taken using the `ptrace` function:

```
julia> bell = 1/√2 * (ket('b','a') + ket('a','b'))
Ket{KroneckerDelta,2,Float64} with 2 state(s):
  0.7071067811865475 | 'a','b' ⟩
  0.7071067811865475 | 'b','a' ⟩

julia> dense = bell * bell'
OuterProduct{KroneckerDelta,2,Float64} with 4 operator(s):
  0.4999999999999999 | 'a','b' ⟩⟨ 'a','b' |
  0.4999999999999999 | 'a','b' ⟩⟨ 'b','a' |
  0.4999999999999999 | 'b','a' ⟩⟨ 'a','b' |
  0.4999999999999999 | 'b','a' ⟩⟨ 'b','a' |

julia> ptrace(dense,1) # trace over the 1st subsystem
GenericOp{KroneckerDelta,1,Float64} with 2 operator(s):
  0.4999999999999999 | 'b' ⟩⟨ 'b' |
  0.4999999999999999 | 'a' ⟩⟨ 'a' |

julia> purity(ans) # get the purity of the previous result
0.4999999999999998
```

