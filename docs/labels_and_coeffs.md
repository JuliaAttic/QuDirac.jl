# QuDirac objects as data structures
---

Under the hood, QuDirac's Kets, Bras, and operator types use dictionaries to map labels to coefficients.

There are few important things to keep in mind when working with these structures:

- All coefficients are restricted to be a subtype of `Number`.
- All state labels are of type `Vector`.
- All operator labels are of type `OpLabel`, a composite type that holds a pair of state labels
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

Note that performing a `setindex!` operation on a state is faster than constructing new states for the sake of addition/subtraction,
but is a mutation of the orginal state:

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

Operator coefficients are accessed in the same manner as state coefficients, except two labels are required; one for the basis Ket, and another for the basis Bra:

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

julia> op[[8,1,10],[2,1,2]]
320
```

The above obviously works with `GenericOp`s as well as `OuterProduct`s.

Assigning coefficients, however, only works with `GenericOp`s (due to the 
`OuterProduct` type simply being a view on its underlying state factors):

```
julia> op[[8,1,10],[2,1,2]] = 1
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

julia> gop[[8,1,10],[2,1,2]] = 1
1

julia> gop[[8,1,10],[2,1,2]]
1
```

There are a few functions, however, that are exclusive to the `OuterProduct` 
type: `getket` and `getbra`. 

These functions are equivalent to taking the inner product of an operator
and a specific basis Ket/Bra, but is generally :

```
julia> getket(op, {2,5,4})
Ket{Orthonormal,3} with 1000 state(s):
  3200 | 8,1,10 ⟩
  20160 | 8,9,7 ⟩
  160 | 2,1,2 ⟩
  1600 | 2,5,4 ⟩
  2880 | 6,4,3 ⟩
  1600 | 1,10,4 ⟩
  ⁞

julia> getket(op, [2,5,4]) == op*ket(2,5,4)
true
```

Likewise with `getbra`:

```
julia> getbra(op, [2,5,4])
Bra{Orthonormal,3} with 1000 state(s):
  3200 ⟨ 8,1,10 |
  20160 ⟨ 8,9,7 |
  160 ⟨ 2,1,2 |
  1600 ⟨ 2,5,4 |
  2880 ⟨ 6,4,3 |
  1600 ⟨ 1,10,4 |
  ⁞

julia> getbra(op, [2,5,4]) == bra(2,5,4)*op
true
```

Why have these methods at all if we can already take inner products? Well, these methods are 
generally more efficient than simply taking the inner product: 

```
julia> @time op * ket(2,5,4);
elapsed time: 0.001857005 seconds (622152 bytes allocated)

julia> @time getket(op, [2,5,4]);
elapsed time: 0.000877927 seconds (331232 bytes allocated)
```

---
#  Using the `get` function
---

Note that trying to access a coefficient using a label not present in the Ket/Bra/operator will 
result in an error:

```
julia> k
Ket{Orthonormal,3} with 125 state(s):
  4 | 1,4,1 ⟩
  48 | 3,4,4 ⟩
  20 | 1,5,4 ⟩
  45 | 5,3,3 ⟩
  6 | 2,3,1 ⟩
  25 | 1,5,5 ⟩
  ⁞

julia> k['a', 'b', 'c']
ERROR: key not found: Char[a,b,c]
 in getindex at /Users/jarrettrevels/.julia/QuDirac/src/dirac/state.jl:50
 in getindex at /Users/jarrettrevels/.julia/QuDirac/src/dirac/state.jl:52
```

In some sense, if a basis state is not "present" in the Ket, the coefficient for that state could be said to be zero.

As such, QuDirac overloads base Julia's `get` method for safely accessing states that may or may not be present in the specified QuDirac object. 

By default, `get` will return a `0` if the passed in label is not found, and the result of `getindex` if it is found:

```
julia> get(k, ['a', 'b', 'c'])
0

julia> get(k, [2, 3, 4])
24

julia> get(k, [2,3,4]) == k[2,3,4]
true
```

One can pass in an extra argument to `get` that replaces `0` as the default return value: 

```
julia> get(k, ['a', 'b', 'c'], "Not here")
"Not here"
```

---
#  Label and coefficient iterators
---

One can iterate through the labels or coefficients of a QuDirac object using the `labels` and `coeffs` functions:



