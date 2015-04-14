# Intro to Inner Products in QuDirac
---

Every QuDirac state or operator has a type parameter `P<:AbstractInner` which describes the inner product behavior for the object. QuDirac comes with two such inner product types:

```
abstract AbstractInner
immutable KroneckerDelta <: AbstractInner end
immutable UndefinedInner <: AbstractInner end
```
Inner products are then evaluated differently based on these types by referring to the `inner_rule` function. This function evaluates the inner product of two basis states given a product type and the states' labels. For example, the definitions of `inner_rule` for the two types above are similar to the following:

```
inner_rule(p::UndefinedInner, b, k) = ScalarExpr(InnerProduct(p, b, k)) # lazy evaluation of inner product
inner_rule(::KroneckerDelta, b, k) = b == k ? 1 : 0
```

The arguments to `inner_rule` are, in order:

1. `p` -> The product type for which this method is defined
2. `b` -> The Bra label of the inner product being evaluated.
3. `k` -> The Ket label of the inner product being evaluated.

The type of `b`/`k` can be a `StateLabel` or a single factor element of a `StateLabel`.

Most examples in this documentation show the `KroneckerDelta` rule in action. To contrast, the following example illustrates the inner product rule for `UndefinedInner`:

```
julia> k = sum(i->i*ket(i), 1:5)
Ket{UndefinedInner,1,Int64} with 5 state(s):
  4 | 4 ⟩
  3 | 3 ⟩
  2 | 2 ⟩
  5 | 5 ⟩
  1 | 1 ⟩

julia> k'*k
(((((((((((((((((((((((((4 * ⟨ 2 | 2 ⟩) + (6 * ⟨ 2 | 3 ⟩)) + 
(10 * ⟨ 2 | 5 ⟩)) + (8 * ⟨ 2 | 4 ⟩)) + (2 * ⟨ 2 | 1 ⟩)) + 
(6 * ⟨ 3 | 2 ⟩)) + (9 * ⟨ 3 | 3 ⟩)) + (15 * ⟨ 3 | 5 ⟩)) + 
(12 * ⟨ 3 | 4 ⟩)) + (3 * ⟨ 3 | 1 ⟩)) + (10 * ⟨ 5 | 2 ⟩)) + 
(15 * ⟨ 5 | 3 ⟩)) + (25 * ⟨ 5 | 5 ⟩)) + (20 * ⟨ 5 | 4 ⟩)) + 
(5 * ⟨ 5 | 1 ⟩)) + (8 * ⟨ 4 | 2 ⟩)) + (12 * ⟨ 4 | 3 ⟩)) + 
(20 * ⟨ 4 | 5 ⟩)) + (16 * ⟨ 4 | 4 ⟩)) + (4 * ⟨ 4 | 1 ⟩)) + 
(2 * ⟨ 1 | 2 ⟩)) + (3 * ⟨ 1 | 3 ⟩)) + (5 * ⟨ 1 | 5 ⟩)) + 
(4 * ⟨ 1 | 4 ⟩)) + ⟨ 1 | 1 ⟩)
```

The parens in the above result may look ugly, but are currently necessary for disambiguating the order of evaluation for more complicated expressions (see [below](#scalar-expressions-and-inner_eval) for details). A more elegant pretty-printing strategy would be preferable, of course - feel free to open a pull request on the QuDirac repo if you have one!

---
# Assigning inner product types to QuDirac objects
---

To create an instance of a Ket or Bra with a specific inner product type, you can simply 
pass an instance of the type to the `ket` or `bra` function:

```
julia> ket(UndefinedInner(), 1, 2)
Ket{UndefinedInner,2,Int64} with 1 state(s):
  1 | 1,2 ⟩
```

If you don't explicitly select a type, a default inner product type is used. The user can set the default inner product type for the current session with the `default_inner` function:

```
julia> default_inner(UndefinedInner());

julia> ket(1,2)
Ket{UndefinedInner,2,Int64} with 1 state(s):
  1 | 1,2 ⟩
```

As you've probably noticed from previous examples, the out-of-the-box default inner product type is `KroneckerDelta`.

---
# Custom inner product types
---

**inner_rule**

---

One can easily define their own inner product type in order to overload QuDirac's built-in behavior:

```
julia> immutable MyInner <: AbstractInner end

# k and b are factor labels
julia> QuDirac.inner_rule(::MyInner, k, b) = sqrt(k+b)
inner_rule (generic function with 4 methods)

julia> default_inner(MyInner());

julia> bra(1)*ket(1)
1.4142135623730951

julia> bra(1,2,3)*ket(1,2,3)
6.928203230275509
```

---
**inner_labels**

---

Overloading `inner_rule` defines the behavior between each factor's Bra and Ket labels. 

There is another method, `inner_labels`, which defines inner product behavior over the *entire* Bra and Ket labels, 
rather than individual factors. By default, this method is defined such that the following is true:

```
⟨ b₁,b₂,b₃... | k₁,k₂,k₃... ⟩ --> ⟨ b₁ | k₁ ⟩ * ⟨ b₂ | k₂ ⟩ * ⟨ b₃ | k₃ ⟩ * ...
```

...where the inner products between the individual factors (denoted by the subscripts above) are evaluated using 
`inner_rule`. 

If one so wishes, `inner_labels` can be overloaded for specific behaviors just like `inner_rule`:

```
julia> immutable SumInner <: AbstractInner end

# b and k are factor labels
julia> QuDirac.inner_rule(::SumInner, b, k) = b + k
inner_rule (generic function with 5 methods)

# sum factor-wise inner products instead of the default behavior (multiplying them)
julia> QuDirac.inner_labels(p::SumInner, b::StateLabel, k::StateLabel) = sum(map((x, y) -> QuDirac.inner_rule(p, x, y), b, k))
inner_labels (generic function with 4 methods)

julia> bra(1,2,3)*ket(4,5,6)
21
```

---
# Scalar expressions and `inner_eval`
---


