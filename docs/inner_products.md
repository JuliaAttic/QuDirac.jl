# Intro to Inner Products in QuDirac
---

Every QuDirac state or operator has a type parameter `P<:AbstractInner` which denotes the inner product behavior for the object. QuDirac comes with two such inner product types:

```julia
abstract AbstractInner
immutable KroneckerDelta <: AbstractInner end
immutable UndefinedInner <: AbstractInner end
```

For each inner product type, a method of the `inner_rule` function is defined which evaluates the inner product of two basis states given the states' labels. For example, the definitions of `inner_rule` for the two types above are similar to the following:

```julia
inner_rule(p::UndefinedInner, b::StateLabel, k::StateLabel) = InnerExpr(InnerProduct(p, b, k)) # lazy evaluation of inner product
inner_rule(::KroneckerDelta, b::StateLabel, k::StateLabel) = b == k ? 1 : 0
```

The arguments to `inner_rule` are, in order:

1. `p::AbstractInner` -> An instance of the product type for which this method is defined
2. `b::StateLabel` -> The Bra label of the inner product being evaluated.
3. `k::StateLabel` -> The Ket label of the inner product being evaluated.

Later, we'll examine how to overload this function for new inner product types. There's one thing to learn before we get there, however: how to switch between existing inner product types.

---
# Assigning inner product types to QuDirac objects
---

To create an instance of a Ket or Bra with a specific inner product type, you can simply 
pass an instance of the type to the `ket` or `bra` function:

```julia
julia> ket(UndefinedInner(), 1, 2)
Ket{UndefinedInner,2,Int64} with 1 state(s):
  1 | 1,2 ⟩
```

If you don't explicitly select a type, a default inner product type is used. The user can set the default inner product type for the current session with the `default_inner` function:

```julia
julia> default_inner(UndefinedInner())
INFO: QuDirac default inner product type is currently UndefinedInner()

julia> ket(1,2)
Ket{UndefinedInner,2,Int64} with 1 state(s):
  1 | 1,2 ⟩
```

As you've probably noticed from previous examples, the out-of-the-box default inner product type is `KroneckerDelta`.

---
# Custom inner product types
---

Defining a new inner product type is as easy as creating the type and its `inner_rule`:

```julia
julia> immutable SumInner <: AbstractInner end

julia> QuDirac.inner_rule(::SumInner, k, b) = int(sum(k) + sum(b))
inner_rule (generic function with 4 methods)

julia> default_inner(SumInner());
INFO: QuDirac's default inner product type is currently SumInner()

julia> d" < 1 | 1 > "
2

julia> d" < 1,2,3 | 4,-5,6 > "
11
```

---
**inner_rettype**

---

For certain operations, QuDirac tries to predict the return type of `inner_rule` to allow for more exact coefficient type inferencing. You can help this prediction by defining a method for the `inner_rettype` function that specificies what the return type will be for your custom inner product evaluations. For example, for the `KroneckerDelta` and `UndefinedInner` types, the following are defined:

```julia
inner_rettype(::UndefinedInner) = InnerExpr
inner_rettype(::KroneckerDelta) = Int64
```

Another example: It's easy to see that evaluations of `inner_rule` with `SumInner` will always produce `Int64`s (or `Int32`s, if you're on a 32 bit machine), yet Julia isn't necessarily able to inference this (the inferenced return type can be verified using the `code_typed` macro).

Thus, I can manually share the assumption about the return type with QuDirac by defining the following:

```julia
julia> QuDirac.inner_rettype(::SumInner) = Int
inner_rettype (generic function with 4 methods)
```

---
# Delayed Inner Product Evaluation
---

Taking inner products with the `UndefinedInner` type will yield `InnerExpr`s, QuDirac's representations of unevaluated inner products. Instances of `InnerExpr` can be treated like numbers in most respects:

```julia
julia> default_inner(UndefinedInner());
INFO: QuDirac default inner product type is currently UndefinedInner()

julia> d" act_on(< 'x' |, | 'a','b','c' >, 2) "
Ket{UndefinedInner,2,Number} with 1 state(s):
  ⟨ 'x' | 'b' ⟩ | 'a','c' ⟩

julia> s = d" √(< 1 | 3 >)^2 + 1 "
(sqrt(((⟨ 1 | 3 ⟩^2) + 1)))

julia> typeof(s)
InnerExpr (constructor with 4 methods)

julia> d" s * | 1,2,3 > "
Ket{UndefinedInner,3,Number} with 1 state(s):
  (sqrt(((⟨ 1 | 3 ⟩^2) + 1))) | 1,2,3 ⟩
```

The `inner_eval` function can be used to re-evaluate `InnerExpr`s by mapping a function to each unresolved inner product:

```julia
julia> s = d" (e^( < 1,2 | 3,4 > ) + < 5,6 | 7,8 > * im)^4 "
(((exp(⟨ 1,2 | 3,4 ⟩)) + (⟨ 5,6 | 7,8 ⟩ * im))^4)

julia> f(b::StateLabel, k::StateLabel) = sum(k) - sum(b);

julia> inner_eval(f, s)
8.600194553751864e6 + 2.5900995362955774e6im

# expanding the function ourselves...
julia> inner_eval(f, s) == ((e^((3+4) - (1+2))) + (((7+8) - (5+6)) * im))^4
true
```

One can pass in inner product type instances instead of functions, which will evaluate the `InnerExpr` using the inner product type's `inner_rule` method:

```julia
# evaluate s using KroneckerDelta inner_rule
inner_eval(KroneckerDelta(), s)
1.0 + 0.0im
```

Finally, `inner_eval` can be called on states and operators to perform the evaluation on their coefficients:

```julia
julia> s = d" act_on( < 'x' | + < 'y' |, | 'a','b','c' > + | 'd','e','f' >, 2) "
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

```julia
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
