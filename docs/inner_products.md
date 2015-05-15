# Intro to Inner Products in QuDirac
---

Every QuDirac state or operator has a type parameter `P<:AbstractInner` which denotes the inner product 
behavior for the object. QuDirac comes with two such inner product types:

```julia
abstract AbstractInner
immutable KronDelta <: AbstractInner end
immutable UndefinedInner <: AbstractInner end
```

For each inner product type, a constructor is defined that evaluates the inner product of two basis states 
given the states' labels. For example, the constructors for the two types above are similar to the following:

```julia
UndefinedInner{N}(b::StateLabel{N}, k::StateLabel{N}) = InnerExpr(InnerLabel(b, k)) # lazy evaluation of inner product
KronDelta{N}(b::StateLabel{N}, k::StateLabel{N}) = b == k ? 1 : 0
```

---
# Assigning inner product types to QuDirac objects
---

To create an instance of a Ket or Bra with a specific inner product type, you can simply 
pass the desired type to the `ket` or `bra` function:

```julia
julia> ket(UndefinedInner, 1, 2)
Ket{UndefinedInner,2,Int64} with 1 state(s):
  1 | 1,2 ⟩
```

If you don't explicitly select a type, a default inner product type is used. The user can set the default inner product type for the current session with the `default_inner` function:

```julia
julia> default_inner(UndefinedInner)
INFO: QuDirac default inner product type is currently UndefinedInner.

julia> ket(1,2)
Ket{UndefinedInner,2,Int64} with 1 state(s):
  1 | 1,2 ⟩
```

As you've probably noticed from previous examples, the out-of-the-box default inner product type is `KronDelta`.

---
# Custom Inner Product Types
---

Defining a new inner product type is done using the `@def_inner` macro. This macro takes in 
the name of the type to be defined, and the return type of its inner product evaluation (this
return type assumption is utilized by QuDirac to perform efficient pre-allocation operations):

```julia
julia> @def_inner SumInner Int
INFO: SumInner is now defined as an inner product type.
INFO: Inner products using the SumInner type should return values of type Int64.
```

The type `SumInner` is now defined, and is automatically given a constructor that behaves like this:

```julia
SumInner{N}(b::StateLabel{N}, k::StateLabel{N}) = SumInner(b[1], k[1]) * SumInner(b[2], k[2]) * ... * SumInner(b[N], k[N])
```

Thus, we need to give the `SumInner` type a constructor for handling the individual factor labels:

```julia
julia> SumInner(b::Int, k::Int) = b + k
SumInner (constructor with 2 methods)
```

Here's our freshly defined `SumInner` type in action:

```julia
julia> default_inner(SumInner)
INFO: QuDirac's default inner product type is currently SumInner.

julia> d" < 1 | 1 > "
2

julia> d" < 1,2,3 | 4,-5,6 > "
-135

julia> ans == d" < 1 | 4 > * < 2 | -5 > * < 3 | 6 > "
true
```

---
# Delayed Inner Product Evaluation
---

Taking inner products with the `UndefinedInner` type will yield `InnerExpr`s, QuDirac's representations of unevaluated inner products. Instances of `InnerExpr` can be treated like numbers in most respects:

```julia
julia> default_inner(UndefinedInner)
INFO: QuDirac default inner product type is currently UndefinedInner.

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

One can pass in an inner product type instead of a function to evaluate the `InnerExpr`:

```julia
# evaluate s using KronDelta
julia> inner_eval(KronDelta, s)
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

julia> inner_eval(KronDelta, s)
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
