# QuDirac's Type Hiearchy
---

This API refers to various types to describe the domain of the included functions. For clarity's sake, here is the type hiearchy for QuDirac objects:

```julia
abstract AbstractDirac{P<:AbstractInner,N}

abstract DiracOp{P,N} <: AbstractDirac{P,N}
abstract DiracState{P,N} <: AbstractDirac{P,N}

type Ket{P,N,T} <: DiracState{P,N}
type Bra{P,N,T} <: DiracState{P,N}

type OuterProduct{P,N,S,K,B} <: DiracOp{P,N}

abstract AbsOpSum{P,N,T} <: DiracOp{P,N}
type OpSum{P,N,T} <: AbsOpSum{P,N,T}
type DualOpSum{P,N,T} <: AbsOpSum{P,N,T}
```

---
# Constructor functions
---

**Ket(dict::Dict{StateLabel{N}, T}),**

**Ket(ptype::AbstractInner, dict::Dict{StateLabel{N}, T})**

Construct a Ket from the provided `dict`.

---
**ket(labels...)**

Construct a single (i.e. non-superposed) Ket with as many factors as there are `labels`. 
See [Constructing Single Kets](constructing_states/#constructing-single-kets) section for more.

---
**bra(labels...)**

Construct a single (i.e. non-superposed) Bra with as many factors as there are `labels`. 
See the [Constructing Single Bras](constructing_states/#constructing-single-bras) section for more.

---
**OpSum(dict::Dict{OpLabel{N}, T}),**

**OpSum(ptype::AbstractInner, dict::Dict{OpLabel{N}, T})**

Construct an operator from the provided `dict`.

---
**@def_op(str)**

Using the definition given by `str`, generate a function that acts on states via the `*` operator. 
See the [Defining Operators as Functions](func_op_def/#defining-operators-as-functions) section for detailed examples.

---
**@rep_op(str, basis_labels)**

Generate an operator representation by applying the definition given by `str` to the labels given by the iterable `basis_labels`.
See the [Generating Operator Representations](func_op_def/#generating-operator-representations) section for detailed examples.

---
**StateLabel(labels::Tuple),**

**StateLabel(labels...)**

Construct a `StateLabel` with the given factor labels.
`StateLabel`s are iterable, indexable, and mappable.

---
**OpLabel(ketlabel::StateLabel, bralabel::StateLabel),**

**OpLabel(ketlabel, bralabel)**

Construct an `OpLabel` with the given Ket and Bra labels.

To retrieve an `OpLabel`'s Ket label, call `klabel(::OpLabel)`.
To retrieve an `OpLabel`'s Bra label, call `blabel(::OpLabel)`.

---
# Math Functions
---

**nfactors(obj::AbstractDirac)**

Returns the number of factors of `obj`. This information is also parameterized in the
type of `obj`, e.g. an instance of type `Ket{KronDelta,3}` has 3 factors.

---
**purity(op::DiracOp)**

Calculate `Tr(op^2)`.

---
**norm(obj::AbstractDirac)**

Calculate the Frobenius norm of `obj`.

---
**scale(obj::AbstractDirac, c::Number),** 

**scale(c::Number, obj::AbstractDirac)**

Multiply `obj` by the scalar `c`. The in-place version, `scale!`, is also defined.

---
**normalize(obj::AbstractDirac)**

Normalize `obj`, i.e. scale it by `1/norm(obj)`. The in-place version, `normalize!`, is also defined.

---
**act_on(a::DiracState, b::DiracState, i::Int),**

**act_on(a::DiracOp, b::DiracState, i::Int)**

Act `a` on the `i`th subsystem of `b`. For states acting on states, this method is discussed in detail [here](state_math/#acting-a-bra-on-a-specific-ket-factor-and-vice-versa). For an operator acting on a state, details can be found [here](op_math/#acting-an-operator-on-a-specific-state-factor).

---
**tensor{P}(ops::DiracOp{P}...),**

**tensor{P}(states::DiracState{P}...)**

Take the tensor product of the given arguments.

---
**trace(op::DiracOp)**

Take the trace of `op`. See [here](op_math/#trace-and-partial-trace) for details.

---
**ptrace(op::DiracOp, i::Int)**

Take the partial trace of `op` over the `i`th subsystem. See [here](op_math/#trace-and-partial-trace) for details.

---
**ptranspose(op::DiracOp, i::Int)**

Take the partial transpose of `op` over the `i`th subsystem. See [here](op_math/#partial-transpose) for details.

---
**commmute(a::DiracOp, b::DiracOp)**

Calculate the commutator `a*b - b*a`. 

---
**anticommmute(a::DiracOp, b::DiracOp)**

Calculate the anticommutator `a*b + b*a`. 

---
**switch(obj::AbstractDirac, i, j)**

Switch the `i`th and `j`th factors in the labels of `obj`.

Example:

```
julia> k = sum(i ->d" i * | i,i+1,i+2 > ", 1:3)
Ket{KronDelta,3,Int64} with 3 state(s):
  1 | 1,2,3 ⟩
  2 | 2,3,4 ⟩
  3 | 3,4,5 ⟩

julia> switch(k, 2, 3)
Ket{KronDelta,3,Int64} with 3 state(s):
  2 | 2,4,3 ⟩
  1 | 1,3,2 ⟩
  3 | 3,5,4 ⟩
```

---
**permute(obj::AbstractDirac, perm::Vector)**

Apply the given permutation, `perm`, to the labels of `obj`

Example:

```
julia> k = sum(i -> d" i * | i,i+1,i+2 > ", 1:3)
Ket{KronDelta,3,Int64} with 3 state(s):
  1 | 1,2,3 ⟩
  2 | 2,3,4 ⟩
  3 | 3,4,5 ⟩

julia> permute(k, [2, 3, 1])
Ket{KronDelta,3,Int64} with 3 state(s):
  3 | 4,5,3 ⟩
  1 | 2,3,1 ⟩
  2 | 3,4,2 ⟩
```
---
**raise(state::DiracState)**

Calculate the action of the raising operator on `state`.
This is generally faster than actually constructing and
applying a raising operator.

---
**lower(state::DiracState)**

Calculate the action of the lowering operator on `state`.
This is generally faster than actually constructing and
applying a lowering operator.

---
**represent(op::DiracOp, labels...)**

Compute the matrix representation of `op` by calculating `∑ᵢⱼ ⟨ i | op | j ⟩ for i,j ∈ product(labels...)` where `product` denotes the catesian product of iterables.

For example:

```julia
julia> @rep_op " H | n > = 1/√2 * ( | 0 > + (-1)^n *| 1 > ) " 0:1
OpSum{KronDelta,1,Float64} with 4 operator(s):
  0.7071067811865475 | 1 ⟩⟨ 0 |
  0.7071067811865475 | 0 ⟩⟨ 0 |
  0.7071067811865475 | 0 ⟩⟨ 1 |
  -0.7071067811865475 | 1 ⟩⟨ 1 |

julia> represent(H, 0:1)
2x2 Array{Float64,2}:
 0.707107   0.707107
 0.707107  -0.707107

julia> represent(tensor(H, H), 0:1, 0:1)
4x4 Array{Float64,2}:
 0.5   0.5   0.5   0.5
 0.5  -0.5   0.5  -0.5
 0.5   0.5  -0.5  -0.5
 0.5  -0.5  -0.5   0.5
```

**represent(state::DiracState, labels...)**

Compute the vector representation of a Ket by calculating `∑ᵢ ⟨ i | state ⟩ for i ∈ product(labels...)` where `product` denotes the catesian product of iterables. If a Bra is passed in, compute the same thing, just conjugate transposed.

For example:

```julia
julia> state = normalize(sum(i -> d" (int(i) + int(i) * im) * | i > ", 'a':'c'))
Ket{KronDelta,1,Complex{Float64}} with 3 state(s):
  0.40823412181754837 + 0.40823412181754837im | 'b' ⟩
  0.41239977612180906 + 0.41239977612180906im | 'c' ⟩
  0.4040684675132877 + 0.4040684675132877im | 'a' ⟩

julia> represent(state, 'a':'c')
3-element Array{Complex{Float64},1}:
 0.404068+0.404068im
 0.408234+0.408234im
   0.4124+0.4124im

julia> represent(state', 'a':'c')
1x3 Array{Complex{Float64},2}:
 0.404068-0.404068im  0.408234-0.408234im  0.4124-0.4124im
```

---
# Inner Product Evaluation
---

**@def_inner(type_name, T)**

Define a new inner product type with name `type_name`. Assume that inner products taken using this type 
will return values of type `V` where `V<:T` (this return type assumption is used for operations that 
require pre-allocation). 

See the [Custom Inner Product Types](inner_products/#custom-inner-product-types) section for details.

---
**inner_eval(f::Function, obj::AbstractDirac),**

**inner_eval(f::Function, i::InnerExpr),**

**inner_eval{P<:AbstractInner}(::Type{P}, obj::AbstractDirac),**

**inner_eval{P<:AbstractInner}(::Type{P}, i::InnerExpr)**

Evaluate unresolved `InnerExpr`s using the provided function, whose signature should be `f(bralabel::StateLabel, ketlabel::StateLabel)`. If a product type is provided instead of a function, use that type for evaluation.

See the [Delayed Inner Product Evaluation](inner_products/#delayed-inner-product-evaluation) section for details.

---
**default_inner{P<:AbstractInner}(::Type{P})**

Set QuDirac's default inner product type to `P`. See [here](inner_products/#assigning-inner-product-types-to-qudirac-objects) for details.

---
# Dict-like Functions
---

**length(obj::AbstractDirac)**

Returns the number of (label, coefficient) pairs stored in `obj`. This is the same as the 
number of basis states/operators present in `obj`.

---
**collect(obj::AbstractDirac)**

Returns a `Vector` of the (label, coefficient) pairs stored in `obj`. The ordering of the returned pairs is not guaranteed.

---
**get(state::DiracState, label[, default]),**

**get(op::DiracOp, ktlabel, brlabel[, default])**

Get the coefficient specified by the provided label(s), return a default value of 0 if the labels could not be found. This default value can be overwridden by passing in the desired value instead. See [here](labels_and_coeffs/#using-the-get-function) for details.

---
**haskey(state::DiracState, label),**

**haskey(op::DiracOp, ktlabel, brlabel)**

Return `true` if the labels are found in the provided objects, and `false` otherwise. 

---
**getindex(state::DiracState, label...),**

**getindex(op::DiracOp, ktlabel, brlabel)**

Return the coefficient for the basis state/operator with the given labels, erroring if the labels could not be found.
See [here](labels_and_coeffs/#accessing-and-assigning-coefficients) for details.

---
**setindex!(state::DiracState, c, label...),**

**setindex!(op::DiracOp, c, ktlabel, brlabel)**

Set the coefficient `c` for the basis state/operator with the given labels, mutating the first argument.
See [here](labels_and_coeffs/#accessing-and-assigning-coefficients) for details.

This function is not implemented on `OuterProduct`s.

---
**delete!(state::DiracState, label),**

**delete!(op::DiracOp, ktlabel, brlabel)**

Delete the specified label->coefficient mapping from the first argument.

This function is not implemented on `OuterProduct`s.

---
# Mapping Functions
---

**map(f::Function, obj::AbstractDirac)**

Maps `f` onto the `(label, coefficient)` pairs of `obj`. For states, `f` is called as `f(::StateLabel,::T)` where `T` is the coefficient type of `obj`. For operators, `f` is called as `f(:OpLabel,::T)`.

Example:

```julia
julia> k0 = sum(ket, 0:10); 

julia> k1 = map((label, v) -> iseven(label[1]) ? (label, v) : (label, 0), k0)
Ket{KronDelta,1,Int64} with 11 state(s):
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
```
---
**maplabels(f::Function, obj::AbstractDirac)**

Maps `f` onto the labels of `obj`. For states, `f` is called as `f(::StateLabel)`. For operators, `f` is called as `f(::OpLabel)`.

Example:

```julia
julia> k = sum(ket, 0:1)*sum(ket, -1:1)
Ket{KronDelta,2,Int64} with 6 state(s):
  1 | 0,0 ⟩
  1 | 0,-1 ⟩
  1 | 1,0 ⟩
  1 | 0,1 ⟩
  1 | 1,1 ⟩
  1 | 1,-1 ⟩

# Performs a shift operation such that
# | 0, j ⟩ -> | 0, j - 1 ⟩
# | 1, j ⟩ -> | 1, j + 1 ⟩
julia> shift_label(s) = StateLabel(s[1], s[2] - (-1)^s[1])
shift_label (generic function with 1 method)

julia> maplabels(shift_label, k)
Ket{KronDelta,2,Int64} with 6 state(s):
  1 | 0,-2 ⟩
  1 | 0,0 ⟩
  1 | 0,-1 ⟩
  1 | 1,2 ⟩
  1 | 1,0 ⟩
  1 | 1,1 ⟩
```
---
**mapcoeffs(f::Function, obj::AbstractDirac)**

Maps `f` onto the coefficients of `obj`. An in-place version, `mapcoeffs!`, is also provided. The function `f` will be called as `f(::T)` where `T` is the coefficient type of `obj`.

Example:

```julia
julia> k = sum(i->d" i * | i > ", 1:6)
Ket{KronDelta,1,Int64} with 6 state(s):
  4 | 4 ⟩
  3 | 3 ⟩
  6 | 6 ⟩
  2 | 2 ⟩
  5 | 5 ⟩
  1 | 1 ⟩

julia> mapcoeffs(i->i^2,k)
Ket{KronDelta,1,Int64} with 6 state(s):
  16 | 4 ⟩
  9 | 3 ⟩
  36 | 6 ⟩
  4 | 2 ⟩
  25 | 5 ⟩
  1 | 1 ⟩
```

---
# Filtering Functions
---

**filter(f::Function, obj::AbstractDirac)**

This function acts exactly like Julia's built-in filtering function for `Dict`s. An in-place version, `filter!`, is also provided.

Example:

```julia
julia> k = normalize(sum(ket, 0:4)^3);

 # extract labels where the second factor is 2
julia> filter((label, c)->label[2]==2, k)
Ket{KronDelta,3} with 25 state(s):
  0.08944271909999159 | 4,2,0 ⟩
  0.08944271909999159 | 1,2,2 ⟩
  0.08944271909999159 | 3,2,3 ⟩
  0.08944271909999159 | 0,2,1 ⟩
  0.08944271909999159 | 0,2,3 ⟩
  ⁞
```

---
**xsubspace(obj::AbstractDirac, x)**

Extracts the elements of `obj` whose labels sum to `x`.

Example:

```julia
julia> k = normalize(sum(ket, 0:4)^3);

julia> xsubspace(k, 10)
Ket{KronDelta,3} with 6 state(s):
  0.08944271909999159 | 4,4,2 ⟩
  0.08944271909999159 | 4,3,3 ⟩
  0.08944271909999159 | 3,4,3 ⟩
  0.08944271909999159 | 3,3,4 ⟩
  0.08944271909999159 | 2,4,4 ⟩
  0.08944271909999159 | 4,2,4 ⟩
```

---
**filternz(obj::AbstractDirac)**

Removes the zero-valued components of `obj`. An in-place version, `filternz!`, is also provided.

Example:

```julia
julia> filternz(d" 0*| 'a' > + 2*| 'b' > + 0*| 'c' > + 1.2*| 'd' >")
Ket{KronDelta,1,Float64} with 2 state(s):
  2.0 | 'b' ⟩
  1.2 | 'd' ⟩
```