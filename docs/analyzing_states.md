## Unordered States, and a Note on Performance
---

Under the hood, QuDirac's Kets and Bras use the `Dict` type to map state labels to coefficients. Because the label --> coefficient map is stored as a `Dict`, state labels are unordered. This is, for the most part, not a problem - ordered bases are mainly useful for mapping basis states to array indices, but since states in QuDirac are uniquely indexed by their labels, we have no need for to establish an order to keep track of the coefficients.

While supporting a lot of `Dict`-like behavior opens up some interesting avenues of analysis (some of which will be demonstrated in the following sections), operations on QuDirac states and operators are generally not as performant as similar operations on pure arrays. Hopefully the tradeoffs will make more sense as the reader makes their way through this documentation. If you find that a given operation is too slow for your purposes, feel free to check if the operation can be performed using [QuBase.jl](https://github.com/JuliaQuantum/QuBase.jl), another Julia library that focuses on quantum operations with `Array`-like structures.

---
## Accessing coefficients
---

Coefficients are accessed by the labels of the basis states:

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

