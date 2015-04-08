# QuDirac objects as data structures
---

Under the hood, QuDirac's Kets, Bras, and operator types use `Dict`s to map labels to coefficients.

There are few important things to keep in mind when working with these structures:

- All coefficients are restricted to be a subtype of `Number`.
- All state labels are of type `StateLabel`.
- All operator labels are of type `OuterLabel`, a composite type that holds two `StateLabel`s (one for the Ket, and one for the Bra)
- Because the label-to-coefficient map is stored as a `Dict`, the components of a QuDirac object are "unordered". In other words, iteration through a QuDirac object's labels and coefficients is not guaranteed to be in any particular order.

---
#  Accessing/Assigning coefficients
---

**States**

---

Julia's `getindex` function has been overloaded to allow coefficients of a state to be accessed by the basis states' labels:

```
julia> k = normalize!(sum(i->i*ket(i), 1:3))
Ket{Orthonormal,1} with 3 state(s):
  0.5345224838248488 | 2 ⟩
  0.8017837257372732 | 3 ⟩
  0.2672612419124244 | 1 ⟩

julia> k[3]
0.8017837257372732

```

The coefficients of product states are accessed in the same fashion:

```
julia> k4 = k^4
Ket{Orthonormal,4} with 81 state(s):
  0.06122448979591838 | 3,2,2,1 ⟩
  0.03061224489795919 | 1,2,1,3 ⟩
  0.18367346938775514 | 3,3,2,2 ⟩
  0.045918367346938785 | 1,3,1,3 ⟩
  0.09183673469387757 | 2,1,3,3 ⟩
  0.06122448979591838 | 3,2,1,2 ⟩
  ⁞

julia> k4[3,2,3,3]
0.2755102040816327

```

One can also peform `setindex!` operations on a state: 

```
julia> k = 0*ket(0)
Ket{Orthonormal,1} with 1 state(s):
  0 | 0 ⟩

julia> k[0] = 1; k[1] = 1/2; k[2] = 1/4; normalize!(k)
Ket{Orthonormal,1} with 3 state(s):
  0.8728715609439696 | 0 ⟩
  0.2182178902359924 | 2 ⟩
  0.4364357804719848 | 1 ⟩
```

Performing a `setindex!` operation on a state is faster than constructing new states 
for the sake of addition/subtraction, but is a mutation of the orginal state:

```
julia> k + ket(0) # this requires constructing a new instance of | 0 ⟩
Ket{Orthonormal,1} with 3 state(s):
  1.8728715609439694 | 0 ⟩
  0.2182178902359924 | 2 ⟩
  0.4364357804719848 | 1 ⟩

julia> k[0] += 1 # faster than the above, but mutates k
1.8728715609439694
```

---
**Operators**

---

Operator coefficients are accessed in the same manner as state coefficients, 
except two labels are required; one for the basis Ket, and another for the basis Bra:

```
julia> k = sum(i->i*ket(i), 1:10); op = k*k'
OuterProduct{Orthonormal,1} with 100 operator(s):
  36 | 6 ⟩⟨ 6 |
  12 | 6 ⟩⟨ 2 |
  18 | 6 ⟩⟨ 3 |
  48 | 6 ⟩⟨ 8 |
  12 | 2 ⟩⟨ 6 |
  4 | 2 ⟩⟨ 2 |
  ⁞

julia> op[6,8]
48

julia> op = tensor(op,op,op)
OuterProduct{Orthonormal,3} with 1000000 operator(s):
  6400 | 8,1,10 ⟩⟨ 8,1,10 |
  40320 | 8,1,10 ⟩⟨ 8,9,7 |
  320 | 8,1,10 ⟩⟨ 2,1,2 |
  3200 | 8,1,10 ⟩⟨ 2,5,4 |
  40320 | 8,9,7 ⟩⟨ 8,1,10 |
  254016 | 8,9,7 ⟩⟨ 8,9,7 |
  ⁞

julia> op[(8,1,10),(2,1,2)]
320
```

The above obviously works with `GenericOp`s as well as `OuterProduct`s.

Assigning coefficients, however, only works with `GenericOp`s (due to the 
`OuterProduct` type simply being a view on its underlying state factors):

```
julia> op[(8,1,10),(2,1,2)] = 1
ERROR: `setindex!` has no method matching setindex!(::OuterProduct{Orthonormal,3}, ::Int64, ::Array{Int64,1}, ::Array{Int64,1})

julia> gop = convert(GenericOp, op)
GenericOp{Orthonormal,3} with 1000000 operator(s):
  33600 | 4,7,4 ⟩⟨ 10,5,6 |
  9720 | 9,9,1 ⟩⟨ 5,3,8 |
  5292 | 4,3,1 ⟩⟨ 9,7,7 |
  47628 | 7,3,7 ⟩⟨ 6,9,6 |
  1701 | 3,3,7 ⟩⟨ 3,1,9 |
  80640 | 7,8,6 ⟩⟨ 4,6,10 |
  ⁞

julia> gop[(8,1,10),(2,1,2)] = 1
1

julia> gop[(8,1,10),(2,1,2)]
1
```

There are a few functions, however, that are exclusive to the `OuterProduct` 
type: `getket` and `getbra`. 

These functions are equivalent to taking the inner product of an operator
and a specific basis Ket/Bra, but is generally :

```
julia> getket(op, (2,5,4))
Ket{Orthonormal,3} with 1000 state(s):
  3200 | 8,1,10 ⟩
  20160 | 8,9,7 ⟩
  160 | 2,1,2 ⟩
  1600 | 2,5,4 ⟩
  2880 | 6,4,3 ⟩
  1600 | 1,10,4 ⟩
  ⁞

julia> getket(op, (2,5,4)) == op*ket(2,5,4)
true
```

Likewise with `getbra`:

```
julia> getbra(op, (2,5,4))
Bra{Orthonormal,3} with 1000 state(s):
  3200 ⟨ 8,1,10 |
  20160 ⟨ 8,9,7 |
  160 ⟨ 2,1,2 |
  1600 ⟨ 2,5,4 |
  2880 ⟨ 6,4,3 |
  1600 ⟨ 1,10,4 |
  ⁞

julia> getbra(op, (2,5,4)) == bra(2,5,4)*op
true
```

Why have these methods at all if we can already take inner products? Well, these methods are 
generally more efficient than calculating the inner product if the result you're looking
for is the action on a basis state.

---
#  Using the `get` function
---

If a label is not explictly present in a QuDirac object, then calling `getindex` with that label results in an error: 

```
julia> k = sum(i-> i * ket(i), 1:3)^3
Ket{Orthonormal,3} with 27 state(s):
  12 | 2,2,3 ⟩
  27 | 3,3,3 ⟩
  4 | 1,2,2 ⟩
  6 | 2,1,3 ⟩
  3 | 1,3,1 ⟩
  18 | 2,3,3 ⟩
  ⁞

julia> k['a', 'b', 'c']
ERROR: key not found: StateLabel{3}(('a','b','c'),0x504c9a1142755a02)
 in getindex at /Users/jarrettrevels/.julia/QuDirac/src/dirac/state.jl:58
 in getindex at /Users/jarrettrevels/.julia/QuDirac/src/dirac/state.jl:61
```

QuDirac overloads Julia's `get` method in order to provide more flexibility in the above scenario. By default, `get` is defined to return a `0` instead of an error if the label it's given isn't found:

```
julia> get(k, ('a', 'b', 'c'))
0

julia> get(k, (2,3,3))
18

julia> get(k, (2,3,3)) == k[2,3,3]
true
```

One can pass in an extra argument to `get` that replaces `0` as the default return value: 

```
julia> get(k, ('a', 'b', 'c'), "Not here")
"Not here"
```

---
#  Label and coefficient iterators
---

One can iterate through the labels or coefficients of a QuDirac object using the `labels` and `coeffs` functions:

---
#  Mapping and filtering functions
---


