# Intro to Inner Products in QuDirac
---

Every QuDirac state or operator has a type parameter `P<:AbstractInner` (e.g., the inner product type of `Ket{AbstractInner,1}` is `AbstractInner`).

QuDirac comes with two such types:

```
abstract AbstractInner
abstract Orthonormal <: AbstractInner
```

---
**The `inner_rule` function**

---

For each type `P<:AbstractInner`, a method of the `inner_rule` function is defined which specifies how
QuDirac should evaluate inner products of basis states with that type. 

For example, the definitions of `inner_rule` for the two types above look similar to the following:

```
inner_rule{P<:AbstractInner}(::Type{P}, b, k) = ScalarExpr(InnerProduct{P}(b, k))
inner_rule{O<:Orthonormal}(::Type{O}, b, k) = b == k ? 1 : 0
```
The arguments to these functions are, in order:

1. `::Type{P}` -> The product type for which this method is defined
2. `b` -> The label of the basis Bra of the inner product.
3. `k` -> The label of the basis Ket of the inner product.

The state labels passed in can be an entire label (i.e., a `Vector{Any}`) or simply a single factor (i.e. an element of a `Vector{Any}`).

---
**Inner product rules in action**

---

The `inner_rule` function is called for each inner product operation that occurs between basis Bras and basis Kets. This is best seen by example.

If we have a Ket with `P<:Orthonormal`, and we take an inner product with itself, we get the following:

```
julia> k
Ket{Orthonormal,1} with 5 state(s):
  2 | 2 ⟩
  3 | 3 ⟩
  5 | 5 ⟩
  4 | 4 ⟩
  1 | 1 ⟩

julia> k'*k
55
```

If that same Ket had an undefined inner product instead (represented by the type `AbstractInner`), the above operation would look like this:

```
julia> k
Ket{AbstractInner,1} with 5 state(s):
  2 | 2 ⟩
  3 | 3 ⟩
  5 | 5 ⟩
  4 | 4 ⟩
  1 | 1 ⟩

julia> k'*k
(((((((((((((((((((((((((4 * ⟨ 2 | 2 ⟩) + (6 * ⟨ 2 | 3 ⟩)) + (10 * ⟨ 2 | 5 ⟩)) + (8 * ⟨ 2 | 4 ⟩)) + (2 * ⟨ 2 | 1 ⟩)) + (6 * ⟨ 3 | 2 ⟩)) + (9 * ⟨ 3 | 3 ⟩)) + (15 * ⟨ 3 | 5 ⟩)) + (12 * ⟨ 3 | 4 ⟩)) + (3 * ⟨ 3 | 1 ⟩)) + (10 * ⟨ 5 | 2 ⟩)) + (15 * ⟨ 5 | 3 ⟩)) + (25 * ⟨ 5 | 5 ⟩)) + (20 * ⟨ 5 | 4 ⟩)) + (5 * ⟨ 5 | 1 ⟩)) + (8 * ⟨ 4 | 2 ⟩)) + (12 * ⟨ 4 | 3 ⟩)) + (20 * ⟨ 4 | 5 ⟩)) + (16 * ⟨ 4 | 4 ⟩)) + (4 * ⟨ 4 | 1 ⟩)) + (2 * ⟨ 1 | 2 ⟩)) + (3 * ⟨ 1 | 3 ⟩)) + (5 * ⟨ 1 | 5 ⟩)) + (4 * ⟨ 1 | 4 ⟩)) + ⟨ 1 | 1 ⟩)

julia> typeof(ans)
ScalarExpr (constructor with 4 methods)
```

Note: The excessive parens in the above result look ugly, but are currently necessary for disambiguating the order of evaluation for more complicated expressions (see [below](#scalar-expressions-and-inner_eval) for details). A more elegant pretty-printing strategy would be preferable, of course - feel free to open a pull request on the QuDirac repo if you have one!

---
# Constructing objects with specific inner product types
---

To create an instance of a Ket or Bra with a specific inner product type, you can simply 
pass that type as the first argument to the `ket` and `bra` functions:

```
julia> ket(AbstractInner, 1, 2)
Ket{AbstractInner,2} with 1 state(s):
  1 | 1,2 ⟩
```

---
# Custom inner product types
---


---
# Scalar expressions and `inner_eval`
---


