# QuDirac's Type Hiearchy
---

This API refers to various types to describe the domain of the included functions. For clarity's sake, here is the type hiearchy for QuDirac objects:

```
abstract AbstractDirac{P<:AbstractInner,N}

abstract DiracOp{P,N} <: AbstractDirac{P,N}
abstract DiracState{P,N} <: AbstractDirac{P,N}

type Ket{P,N,T} <: DiracState{P,N}
type Bra{P,N,T} <: DiracState{P,N}

type OpSum{P,N,T} <: DiracOp{P,N}
type DualOpSum{P,N,T} <: DiracOp{P,N}
type OuterProduct{P,N,S,K,B} <: DiracOp{P,N}
```

---
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
*@def_op(str)*

Using the definition given by `str`, generate a function that acts on Kets via the `*` operator. 

See the [Functionally Defining Operators](adv_cons/#functionally-defining-operators) section for detailed examples.

---
*@repr_op(str, basis_labels)*

Generate an operator representation by applying the definition given by `str` to the labels provided by `basis_labels`.

See the [Functionally Representing Operators](adv_cons/#functionally-representing-operators) section for detailed examples.


---
# Math Functions
---

*nfactors(obj::AbstractDirac)*

Returns the number of factor systems of `obj`. This information is also parameterized in the
type of `obj`, e.g. an instance of type `Ket{KroneckerDelta,3}` has 3 factors.

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

*length(obj::AbstractDirac)*

Returns the number of (label, coefficient) pairs stored in `obj`. This is the same as the 
number of basis states/operators present in `obj`.

---
*collect(obj::AbstractDirac)*

Returns a `Vector` of the (label, coefficient) pairs stored in `obj`. The ordering of the returned pairs is not guaranteed.

---
*get(state::DiracState, label[, default]),*

*get(op::DiracOp, ktlabel, brlabel[, default])*

Get the coefficient specified by the provided label(s), return a default value of 0 if the labels could not be found. This default value can be overwridden by passing in the desired value instead. See [here](labels_and_coeffs/#using-the-get-function) for details.

---
*haskey(state::DiracState, label),*

*haskey(op::DiracOp, ktlabel, brlabel)*

Return `true` if the labels are found in the provided objects, and `false` otherwise. 

---
*getindex(state::DiracState, label...),*

*getindex(op::DiracOp, ktlabel, brlabel)*

Return the coefficient for the basis state/operator with the given labels, erroring if the labels could not be found.
See [here](labels_and_coeffs/#accessing-and-assigning-coefficients) for details.

---
*setindex!(state::DiracState, c, label...),*

*setindex!(op::DiracOp, c, ktlabel, brlabel)*

Set the coefficient `c` for the basis state/operator with the given labels, mutating the first argument.
See [here](labels_and_coeffs/#accessing-and-assigning-coefficients) for details.

This function is not implemented on `OuterProduct`s.

---
*delete!(state::DiracState, label),*

*delete!(op::DiracOp, ktlabel, brlabel)*

Delete the specified label->coefficient mapping from the first argument.

This function is not implemented on `OuterProduct`s.

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
julia> k = normalize(sum(ket, 0:4)^3);

 # extract labels where the second factor is 2
julia> filter((label, c)->label[2]==2, k)
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
julia> k = normalize(sum(ket, 0:4)^3);

julia> xsubspace(k, 10)
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
julia> k0 = sum(ket, 0:10); 

julia> k1 = map((label, v) -> iseven(label[1]) ? (label, v) : (label, 0), k0)
Ket{KroneckerDelta,1,Int64} with 11 state(s):
  1 | 2 ⟩
  0 | 5 ⟩
  0 | 1 ⟩
  1 | 6 ⟩
  0 | 9 ⟩
  1 | 8 ⟩
  0 | 7 ⟩
  1 | 4 ⟩
  0 | 3 ⟩
  1 | 0 ⟩
  1 | 10 ⟩

julia> filternz(k1)
Ket{KroneckerDelta,1,Int64} with 6 state(s):
  1 | 2 ⟩
  1 | 6 ⟩
  1 | 8 ⟩
  1 | 4 ⟩
  1 | 0 ⟩
  1 | 10 ⟩
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
