# Constructor functions
---

*ket(labels...)*

Construct a single (i.e. non-superposed) Ket with as many factors as there are `labels`. 
See [Constructing Single Kets](constructing_states/#constructing-single-kets) section for more.

---
*bra(labels...)*

Construct a single (i.e. non-superposed) Bra with as many factors as there are `labels`. 
See the [Constructing Single Bras](constructing_states/#constructing-single-bras) section for more.

---
*func_op(f::Function, k::Ket[, i::Int])*

Generate the representation for the operator specified by the function `f` acting on `k`'s basis. 
One can optionally pass in `i` to specify the operator's action on a specific factor of the basis.
See the [Functionally Defining Operators](func_op_def.md) section for more.

---
# Math Functions
---

*purity(op::DiracOp)*

Calculate `Tr(op^2)`.

---
*norm(obj::AbstractDirac)*

Calculate the Frobenius norm of `obj`.

---
*scale(obj::AbstractDirac, c::Number),* 

*scale(c::Number, obj::AbstractDirac)*

Multiply `obj` by the scalar `c`. The in-place version, `scale!`, is also defined.

---
*normalize(obj::AbstractDirac)*

Normalize `obj`, i.e. scale it by `1/norm(obj)`. The in-place version, `normalize!`, is also defined.

---
*act_on{P}(br::Bra{P,1}, kt::Ket{P}, i::Int)*

Act `br` on the `i`th factor of `kt`. This method is discussed in detail [here](state_math/#acting-a-bra-on-a-specific-ket-factor).

*act_on{P}(op::DiracOp{P,1}, kt::Ket{P}, i::Int)*

Act `op` on the `i`th factor of `kt`. This method is discussed in detail [here](op_math/#acting-an-operator-on-a-specific-ket-factor).

---
*tensor{P}(ops::DiracOp{P}...),* 

*tensor{P}(states::DiracState{P}...)*

Take the tensor product of the given arguments.

---
*trace(op::DiracOp)*

Take the trace of `op`. See [here](op_math/#trace-and-partial-trace) for details.

---
*ptrace(op::DiracOp, i::Int)*

Take the partial trace of `op` over the `i`th subsystem. See [here](op_math/#trace-and-partial-trace) for details.

---
*commmutator(a::DiracOp, b::DiracOp)*

Calculate the commutator `a*b - b*a`. 

---
*anticommmutator(a::DiracOp, b::DiracOp)*

Calculate the anticommutator `a*b + b*a`. 

---
*switch(obj::AbstractDirac, i, j)*

Switch the `i`th and `j`th factors in the labels of `obj`.

Example:

```
julia> k = sum(i -> i * ket(i, i+1, i+2), 1:3)
Ket{KroneckerDelta,3,Int64} with 3 state(s):
  1 | 1,2,3 ⟩
  2 | 2,3,4 ⟩
  3 | 3,4,5 ⟩

julia> switch(k, 2, 3)
Ket{KroneckerDelta,3,Int64} with 3 state(s):
  2 | 2,4,3 ⟩
  1 | 1,3,2 ⟩
  3 | 3,5,4 ⟩
```

---
*permute(obj::AbstractDirac, perm::Vector)*

Apply the given permutation, `perm`, to the labels of `obj`

Example:

```
julia> k = sum(i -> i * ket(i, i+1, i+2), 1:3)
Ket{KroneckerDelta,3,Int64} with 3 state(s):
  1 | 1,2,3 ⟩
  2 | 2,3,4 ⟩
  3 | 3,4,5 ⟩

julia> permute(k, [2, 3, 1])
Ket{KroneckerDelta,3,Int64} with 3 state(s):
  3 | 4,5,3 ⟩
  1 | 2,3,1 ⟩
  2 | 3,4,2 ⟩
```

---
# Dict-like Functions
---

*length*

---
*nfactors*

---
*get*

---
*haskey*

---
*getindex*

---
*setindex!*

---
*delete!*

---
# Mapping Functions
---

*map(f::Function, obj::AbstractDirac)*

Maps `f` onto the `(label, coefficient)` pairs of `obj`. For states, `f` is called as `f(::StateLabel,::T)` where `T` is the coefficient type of `obj`. For operators, `f` is called as `f(:OpLabel,::T)`.

---
*maplabels(f::Function, obj::AbstractDirac)*

Maps `f` onto the labels of `obj`. For states, `f` is called as `f(::StateLabel)`. For operators, `f` is called as `f(::OpLabel)`.

---
*mapcoeffs(f::Function, obj::AbstractDirac)*

Maps `f` onto the coefficients of `obj`. An in-place version, `mapcoeffs!`, is also provided. The function `f` will be called as `f(::T)` where `T` is the coefficient type of `obj`.

---
# Filtering Functions
---

*filter(f::Function, obj::AbstractDirac)*

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

---
*xsubspace(obj::AbstractDirac, x)*

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

---
*filternz(obj::AbstractDirac)*

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
# Inner Product Evaluation
---

*inner_eval*

---
*default_inner*

---
*QuDirac.inner_rule*

---
*QuDirac.inner_rettype*
