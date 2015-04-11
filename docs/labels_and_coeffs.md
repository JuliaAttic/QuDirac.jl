# QuDirac objects as data structures
---

Under the hood, QuDirac's Kets, Bras, and operator types use `Dict`s to map labels to coefficients.

There are few important things to keep in mind when working with these structures:

- All state labels are of type `StateLabel`.
- All operator labels are of type `OuterLabel`, a composite type that holds two `StateLabel`s (one for the Ket, and one for the Bra).
- Because the label-to-coefficient map is stored as a `Dict`, the components of a QuDirac object are unordered. In other words, iteration through a QuDirac object's labels/coefficients is not guaranteed to be in any particular order.

---
#  Accessing/Assigning coefficients
---

**States**

---

Julia's `getindex` function has been overloaded to allow coefficients of a state to be accessed by the basis states' labels:

```
julia> k = normalize(sum(i->i*ket(i), 1:3))
Ket{KroneckerDelta,1,Float64} with 3 state(s):
  0.8017837257372732 | 3 ⟩
  0.5345224838248488 | 2 ⟩
  0.2672612419124244 | 1 ⟩

julia> k[3]
0.8017837257372732

```

The coefficients of product states are accessed in the same fashion:

```
julia> k4 = k^4
Ket{KroneckerDelta,4,Float64} with 81 state(s):
  0.03061224489795919 | 2,1,3,1 ⟩
  0.03061224489795919 | 1,3,2,1 ⟩
  0.09183673469387757 | 3,2,3,1 ⟩
  0.18367346938775514 | 2,2,3,3 ⟩
  0.18367346938775514 | 3,2,2,3 ⟩
  0.06122448979591838 | 1,2,3,2 ⟩
  ⁞

julia> k4[3,2,3,3]
0.2755102040816327

```

One can also use `setindex!` on states: 

```
julia> k = 0.0*ket(0)
Ket{KroneckerDelta,1,Float64} with 1 state(s):
  0.0 | 0 ⟩

julia> k[0] = 1; k[1] = 1/2; k[2] = 1/4; normalize!(k)
Ket{KroneckerDelta,1,Float64} with 3 state(s):
  0.2182178902359924 | 2 ⟩
  0.8728715609439696 | 0 ⟩
  0.4364357804719848 | 1 ⟩
```

Performing a `setindex!` operation on a state is faster than constructing new states 
for the sake of addition/subtraction, but is a mutation of the orginal state:

```
julia> k + ket(0) # this requires constructing a new instance of | 0 ⟩
Ket{KroneckerDelta,1,Float64} with 3 state(s):
  0.2182178902359924 | 2 ⟩
  1.8728715609439694 | 0 ⟩
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
OuterProduct{KroneckerDelta,1,Int64,Int64,Int64} with 100 operator(s):
  81 | 9 ⟩⟨ 9 |
  36 | 9 ⟩⟨ 4 |
  72 | 9 ⟩⟨ 8 |
  27 | 9 ⟩⟨ 3 |
  36 | 4 ⟩⟨ 9 |
  16 | 4 ⟩⟨ 4 |
  ⁞

julia> op[6,8]
48

julia> op = tensor(op,op,op)
OuterProduct{KroneckerDelta,3,Int64,Int64,Int64} with 1000000 operator(s):
  1600 | 4,5,2 ⟩⟨ 4,5,2 |
  1920 | 4,5,2 ⟩⟨ 2,8,3 |
  8640 | 4,5,2 ⟩⟨ 6,6,6 |
  4000 | 4,5,2 ⟩⟨ 10,2,5 |
  1920 | 2,8,3 ⟩⟨ 4,5,2 |
  2304 | 2,8,3 ⟩⟨ 2,8,3 |
  ⁞

julia> op[(8,1,10),(2,1,2)]
320
```

The above obviously works with `GenericOp`s as well as `OuterProduct`s.

Assigning coefficients, however, only works with `GenericOp`s (due to the 
`OuterProduct` type simply being a view on its underlying state factors):

```
julia> op[(8,1,10),(2,1,2)] = 1
ERROR: `setindex!` has no method matching setindex!(::OuterProduct{KroneckerDelta,3,Int64,Int64,Int64}, ::Int64, ::(Int64,Int64,Int64), ::(Int64,Int64,Int64))

julia> gop = convert(GenericOp, op)
GenericOp{KroneckerDelta,3,Int64} with 1000000 operator(s):
  168 | 7,4,1 ⟩⟨ 1,6,1 |
  90 | 1,1,9 ⟩⟨ 5,1,2 |
  2520 | 5,6,2 ⟩⟨ 2,3,7 |
  51840 | 6,4,5 ⟩⟨ 6,9,8 |
  2800 | 10,5,7 ⟩⟨ 8,1,1 |
  26460 | 9,3,7 ⟩⟨ 7,10,2 |
  ⁞

julia> gop[(8,1,10),(2,1,2)] = 1
1

julia> gop[(8,1,10),(2,1,2)]
1
```

---
#  Using the `get` function
---

If a label is not explictly present in a QuDirac object, then calling `getindex` with that label results in an error: 

```
julia> k = sum(i-> i * ket(i), 1:3)^3
Ket{KroneckerDelta,3,Int64}} with 27 state(s):
  12 | 2,2,3 ⟩
  27 | 3,3,3 ⟩
  4 | 1,2,2 ⟩
  6 | 2,1,3 ⟩
  3 | 1,3,1 ⟩
  18 | 2,3,3 ⟩
  ⁞

julia> k['a', 'b', 'c']
ERROR: key not found: StateLabel{3}('a','b','c')
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
#  Mapping and filtering functions
---


