## QuDirac objects as data structures
---

Under the hood, QuDirac's Kets, Bras, and operator types use dictionaries to map labels to coefficients.

There are few important things to keep in mind when working with these structures:

- All coefficients are restricted to be a subtype of `Number`.
- All state labels are of type `Vector{Any}`.
- All operator labels are of type `OpLabel`, a composite type that holds a pair of state labels
- Because the label-to-coefficient map is stored as a `Dict`, the components of a QuDirac object are "unordered". In other words, iteration through a QuDirac object's labels is not guaranteed to be in any particular order.

---
##  Accessing/Assigning coefficients
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
  0.040816326530612256 | 2,1,2,2 ⟩
  0.2755102040816327 | 3,2,3,3 ⟩
  0.13775510204081634 | 3,3,1,3 ⟩
  0.020408163265306128 | 1,2,2,1 ⟩
  0.13775510204081634 | 3,3,3,1 ⟩
  0.09183673469387757 | 3,3,1,2 ⟩
  0.09183673469387757 | 3,3,2,1 ⟩
  0.03061224489795919 | 3,1,1,2 ⟩
  0.03061224489795919 | 1,3,2,1 ⟩
  0.12244897959183676 | 2,3,2,2 ⟩
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

TODO

---
##  Label and coefficient iterators
---

TODO

