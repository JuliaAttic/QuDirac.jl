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
inner_rule(p::UndefinedInner, b::StateLabel, k::StateLabel) = InnerExpr(InnerProduct(p, b, k)) # lazy evaluation of inner product
inner_rule(::KroneckerDelta, b::StateLabel, k::StateLabel) = b == k ? 1 : 0
```

The arguments to `inner_rule` are, in order:

1. `p::AbstractInner` -> An instance of the product type for which this method is defined
2. `b::StateLabel` -> The Bra label of the inner product being evaluated.
3. `k::StateLabel` -> The Ket label of the inner product being evaluated.

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
julia> immutable SumInner <: AbstractInner end

julia> QuDirac.inner_rule(::SumInner, k, b) = sum(k) + sum(b)
inner_rule (generic function with 4 methods)

julia> default_inner(SumInner());

julia> bra(1)*ket(1)
2

julia> bra(1,2,3)*ket(4,-5,6)
11
```

---
**inner_rettype**

---

For certain operations, QuDirac tries to predict the return type of `inner_rule` to allow for more exact coefficient type inferencing. You can help this prediction by defining a method for the `inner_rettype` function that specificies what the return type will be for your custom inner product evaluations. For example, for the `KroneckerDelta` and `UndefinedInner` types, the following are defined:

```
inner_rettype(::UndefinedInner) = InnerExpr
inner_rettype(::KroneckerDelta) = Int64
```

Another example: If my use of `SumInner` will be restricted to cases where the state labels hold `Int64`s, it's easy to see that evaluations using `SumInner` will always produce `Int64`s. I can then share this assumption with QuDirac by defining the following:

```
julia> QuDirac.inner_rettype(::SumInner) = Int64
inner_rettype (generic function with 4 methods)
```

---
# Delayed Inner Product Evaluation
---

Evaluation using the `UndefinedInner` type yields objects of type `InnerExpr`. These objects are merely representations of unevaluated inner products, and can be treated like numbers in most respects:

```
julia> default_inner(UndefinedInner());

julia> act_on(bra('x'), ket('a','b','c'), 2)
Ket{UndefinedInner,2,Number} with 1 state(s):
  ⟨ 'x' | 'b' ⟩ | 'a','c' ⟩

julia> s = √((bra(1) * ket(3))^2 + 1)
(sqrt(((⟨ 1 | 3 ⟩^2) + 1)))

julia> typeof(s)
InnerExpr (constructor with 4 methods)

julia> s * ket(1,2,3)
Ket{UndefinedInner,3,Number} with 1 state(s):
  (sqrt(((⟨ 1 | 3 ⟩^2) + 1))) | 1,2,3 ⟩
```

The `inner_eval` function can be used to re-evaluate `InnerExpr`s by mapping a function to each unresolved inner product:

```
julia> s = (e^(bra(1,2) * ket(3,4)) + (bra(5,6)*ket(7,8))im)^4
(((exp(⟨ 1,2 | 3,4 ⟩)) + (⟨ 5,6 | 7,8 ⟩ * im))^4)

julia> f(b::StateLabel, k::StateLabel) = sum(k) - sum(b)
f (generic function with 1 method)

julia> inner_eval(f, s)
8.600194553751864e6 + 2.5900995362955774e6im

# expanding the function ourselves...
julia> inner_eval(f, s) == ((e^((3+4) - (1+2))) + (((7+8) - (5+6)) * im))^4
true
```

One can pass in inner product type instances instead of functions, which will evaluate the `InnerExpr` using the inner product type's `inner_rule` method:

```
# evaluate s using KroneckerDelta inner_rule
inner_eval(KroneckerDelta(), s)
1.0 + 0.0im
```

Finally, `inner_eval` can be called on states and operators to perform the evaluation on their coefficients:

```
julia> s = act_on(bra('x') + bra('y'), ket('a','b','c') + ket('d', 'e', 'f'), 2)
Ket{UndefinedInner,2,Number} with 2 state(s):
  (⟨ 'x' | 'e' ⟩ + ⟨ 'y' | 'e' ⟩) | 'd','f' ⟩
  (⟨ 'x' | 'b' ⟩ + ⟨ 'y' | 'b' ⟩) | 'a','c' ⟩

julia> inner_eval((b,k) -> int(b[1]) * int(k[1]), s)
Ket{UndefinedInner,2,Int64} with 2 state(s):
  24341 | 'd','f' ⟩
  23618 | 'a','c' ⟩

julia> inner_eval(KroneckerDelta(), s)
Ket{UndefinedInner,2,Int64} with 2 state(s):
  0 | 'd','f' ⟩
  0 | 'a','c' ⟩
```

---
# Functions supported by `InnerExpr`
---

The following is a list of the functions supported for use with `InnerExpr` (keep in mind that `InnerExpr <: Number`):

```
one(::InnerExpr)
zero(::InnerExpr)

abs(::InnerExpr)
abs2(::InnerExpr)

conj(::InnerExpr)
ctranspose(::InnerExpr)

+(::InnerExpr, ::Number)
+(::Number, ::InnerExpr)

-(::InnerExpr, ::Number)
-(::Number, ::InnerExpr)

/(::InnerExpr, ::Number)
/(::Number, ::InnerExpr)

*(::InnerExpr, ::Number)
*(::Number, ::InnerExpr)

^(::InnerExpr, ::Number)
^(::Number, ::InnerExpr)

exp(::InnerExpr)
exp2(::InnerExpr)

sqrt(::InnerExpr)

log(::InnerExpr, ::Number)
log(::Number, ::InnerExpr)

log(::InnerExpr)
log2(::InnerExpr)
```

If you would like support for a function not in the above list, feel free to open an issue or pull request on the QuDirac repo.
