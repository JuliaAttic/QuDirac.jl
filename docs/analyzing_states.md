## States are basically `Dict`s
---

Under the hood, QuDirac's Kets and Bras use the `Dict` type to map state labels to coefficients. 

There are few important things to keep in mind when working with states:

- State labels are of type `Vector{Any}` and coefficients are restricted to be a subtype of `Number`.
- Because the label-to-coefficient map is stored as a `Dict`, the state labels are "unordered", i.e. iteration through a set of labels is not guaranteed to be in any particular order.

Lastly, a note on performance:

While supporting a lot of `Dict`-like behavior opens up some interesting avenues of analysis (some of which will be demonstrated in the following sections), operations on QuDirac states and operators are generally not as performant as similar operations on pure arrays. Hopefully the tradeoffs will make more sense as the reader makes their way through this documentation. If you find that a given operation is too slow for your purposes, feel free to check if the operation can be performed using [QuBase.jl](https://github.com/JuliaQuantum/QuBase.jl), another Julia library that focuses on quantum operations with `Array`-like structures.

---
## Accessing/Assigning coefficients with labels
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

Tip: Performing a `setindex!` operation on a state is faster than constructing new states for the sake of addition/subtraction:

```
julia> k + ket(0) # this requires constructing a new instance of | 0 ⟩
Ket{Orthonormal,1} with 3 state(s):
  1.8728715609439694 | 0 ⟩
  0.2182178902359924 | 2 ⟩
  0.4364357804719848 | 1 ⟩

julia> k[0] += 1 # faster than the above, but of course mutates k
1.8728715609439694
```

---
## Filtering States
---

QuDirac provides a variety of functions for filtering states. The most generic of these functions is `filter(f::Function, s::AbstractState)` (and its in-place counterpart, `filter!`). This `filter` function acts exactly like a Julia's built-in filtering functions for `Dict`s:

```
julia> s = normalize!(sum(ket, 0:4)^3);

julia> filter((label, c)->label[2]==2, s) # extract labels where the second factor is labeled "2" 
Ket{Orthonormal,3} with 25 state(s):
  0.08944271909999159 | 4,2,0 ⟩
  0.08944271909999159 | 1,2,2 ⟩
  0.08944271909999159 | 3,2,3 ⟩
  0.08944271909999159 | 0,2,1 ⟩
  0.08944271909999159 | 0,2,3 ⟩
  0.08944271909999159 | 0,2,4 ⟩
  0.08944271909999159 | 0,2,0 ⟩
  0.08944271909999159 | 3,2,2 ⟩
  0.08944271909999159 | 1,2,4 ⟩
  0.08944271909999159 | 2,2,3 ⟩
  0.08944271909999159 | 2,2,4 ⟩
  0.08944271909999159 | 2,2,2 ⟩
  0.08944271909999159 | 4,2,1 ⟩
  0.08944271909999159 | 2,2,0 ⟩
  0.08944271909999159 | 2,2,1 ⟩
  0.08944271909999159 | 0,2,2 ⟩
  ⁞
```

For convience, a function `xsubspace(s::AbstractState, x)` is provided which extracts the states whose labels sum to `x`:

```julia
julia> xsubspace(s, 10)
Ket{Orthonormal,3} with 6 state(s):
  0.08944271909999159 | 4,4,2 ⟩
  0.08944271909999159 | 4,3,3 ⟩
  0.08944271909999159 | 3,4,3 ⟩
  0.08944271909999159 | 3,3,4 ⟩
  0.08944271909999159 | 2,4,4 ⟩
  0.08944271909999159 | 4,2,4 ⟩
```

For example's sake, this is equivalent to calling:

```
julia> filter((label,c)->sum(label)==10, s)
Ket{Orthonormal,3} with 6 state(s):
  0.08944271909999159 | 4,4,2 ⟩
  0.08944271909999159 | 4,3,3 ⟩
  0.08944271909999159 | 3,4,3 ⟩
  0.08944271909999159 | 3,3,4 ⟩
  0.08944271909999159 | 2,4,4 ⟩
  0.08944271909999159 | 4,2,4 ⟩
```

Our last filtering function is `filternz(s::AbstractState)`, which filters out the zero components of a state:

```
julia> k = sum(ket, 0:4); k[2] = 0; normalize!(k)
Ket{Orthonormal,1} with 5 state(s):
  0.5 | 0 ⟩
  0.0 | 2 ⟩
  0.5 | 3 ⟩
  0.5 | 4 ⟩
  0.5 | 1 ⟩

julia> filternz(k^3)
Ket{Orthonormal,3} with 64 state(s):
  0.125 | 1,4,1 ⟩
  0.125 | 3,4,4 ⟩
  0.125 | 3,3,1 ⟩
  0.125 | 0,3,1 ⟩
  0.125 | 1,1,0 ⟩
  0.125 | 1,3,3 ⟩
  0.125 | 3,4,0 ⟩
  0.125 | 3,1,4 ⟩
  0.125 | 0,1,3 ⟩
  0.125 | 4,3,1 ⟩
  0.125 | 0,4,1 ⟩
  0.125 | 4,3,3 ⟩
  0.125 | 0,3,4 ⟩
  0.125 | 1,0,4 ⟩
  0.125 | 4,1,3 ⟩
  0.125 | 0,1,1 ⟩
  ⁞
```


