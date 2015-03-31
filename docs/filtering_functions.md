QuDirac provides a variety of functions for filtering out components of states.

### filter/filter!

Methods: `filter(obj::AbstractDirac)`, `filter!(obj::AbstractDirac)` 
Description: This function acts exactly like a Julia's built-in filtering functions for `Dict`s

Example:

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

### xsubspace

Methods: `xsubspace(obj::AbstractDirac, x)`
Description: Extracts the elements of `obj` whose labels sum to `x`.

Example:

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

### filternz/filternz!

Methods: `filter(obj::AbstractDirac)`, `filter!(obj::AbstractDirac)` 
Description: Removes the zero components of `obj`

Example:

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
