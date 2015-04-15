## Constructor functions
---

_**ket**_

_**bra**_

_**func_op**_

---
## Math Functions
---

_**purity**_

_**norm**_

_**normalize/normalize!**_

_**scale/scale!**_

_**act_on**_

_**tensor**_

_**trace**_

_**ptrace**_

_**commmutator**_

_**anticommmutator**_

_**switch**_

_**permute**_

---
## `Dict`-like Functions
---

_**length**_

_**nfactors**_

_**get**_

_**getket/getbra**_

_**haskey**_

_**getindex**_

_**setindex!**_

_**delete!**_

---
## Mapping Functions
---

_**map(f::Function, obj::AbstractDirac)**_

Maps `f` onto the `(label, coefficient)` pairs of `obj`. For states, `f` is called as `f(::StateLabel,::T)` where `T` is the coefficient type of `obj`. For operators, `f` is called as `f(:OpLabel,::T)`.

_**maplabels(f::Function, obj::AbstractDirac)**_

Maps `f` onto the labels of `obj`. For states, `f` is called as `f(::StateLabel)`. For operators, `f` is called as `f(::OpLabel)`.

_**mapcoeffs(f::Function, obj::AbstractDirac)**_

Maps `f` onto the coefficients of `obj`. An in-place version, `mapcoeffs!`, is also provided. The function `f` will be called as `f(::T)` where `T` is the coefficient type of `obj`.

---
## Filtering Functions
---

_**filter(f::Function, obj::AbstractDirac)**_

This function acts exactly like Julia's built-in filtering function for `Dict`s. An in-place version, `filter!`, is also provided.

Example:

```julia
julia> s = normalize!(sum(ket, 0:4)^3);

julia> filter((label, c)->label[2]==2, s) # extract labels where the second factor is labeled "2" 
Ket{KroneckerDelta,3} with 25 state(s):
  0.08944271909999159 | 4,2,0 ⟩
  0.08944271909999159 | 1,2,2 ⟩
  0.08944271909999159 | 3,2,3 ⟩
  0.08944271909999159 | 0,2,1 ⟩
  0.08944271909999159 | 0,2,3 ⟩
  ⁞
```

_**xsubspace(obj::AbstractDirac, x)**_

Extracts the elements of `obj` whose labels sum to `x`

Example:

```julia
julia> xsubspace(s, 10)
Ket{KroneckerDelta,3} with 6 state(s):
  0.08944271909999159 | 4,4,2 ⟩
  0.08944271909999159 | 4,3,3 ⟩
  0.08944271909999159 | 3,4,3 ⟩
  0.08944271909999159 | 3,3,4 ⟩
  0.08944271909999159 | 2,4,4 ⟩
  0.08944271909999159 | 4,2,4 ⟩
```

_**filternz(obj::AbstractDirac)**_

Removes the zero-valued components of `obj`. An in-place version, `filternz!`, is also provided.

Example:

```julia
julia> k = sum(ket, 0:4); k[2] = 0; normalize!(k)
Ket{KroneckerDelta,1} with 5 state(s):
  0.5 | 0 ⟩
  0.0 | 2 ⟩
  0.5 | 3 ⟩
  0.5 | 4 ⟩
  0.5 | 1 ⟩

julia> filternz(k^3)
Ket{KroneckerDelta,3} with 64 state(s):
  0.125 | 1,4,1 ⟩
  0.125 | 3,4,4 ⟩
  0.125 | 3,3,1 ⟩
  0.125 | 0,3,1 ⟩
  0.125 | 1,1,0 ⟩
  0.125 | 1,3,3 ⟩
  ⁞
```

---
## Inner Product Evaluation
---

_**inner_eval**_

_**default_inner**_

_**QuDirac.inner_rule**_

_**QuDirac.inner_rettype**_
