# Natural Dirac notation input syntax
---

QuDirac supports a "natural" input format by implementing `d"..."` literals. Seeing this syntax in action:

```julia
julia> d" | 0 > "
Ket{KronDelta,1,Int64} with 1 state(s):
  1 | 0 ⟩
```

Using `d"..."` calls the macro `@d_str`. This macro performs three transformations on the 
given string:

1. Replace `" < ... | ... > "` with `"(bra(...)*ket(...))"`
2. Replace `" < ... | "` with `"bra(...)"`
2. Replace `" | ... > "` with `"ket(...)"`

...then parses the resulting string into an `Expr`, which is returned back to the calling environment and evaluated.

For example, this statement:

```julia
d" < :a | :b > * | 1 > + 3 * | 2 > "
```

...translates to this code:

```julia
(bra(:a) * ket(:b)) * ket(1) + 3 * ket(2)
```

When using `d"..."` literals, the `|`, `>` and `<` characters *cannot be used for anything other than as symbols for Kets, Bras, and inner products*. Other than that, any other Julia syntax should evaluate properly. For example, assignments and function calls work as expected:

```julia
julia> d" A = tensor(| 0 >< 1 |, | 1 >< 0 |) ";

julia> A
OuterProduct with 1 operator(s); Ket{UndefinedInner,2,Int64} * Bra{UndefinedInner,2,Int64}:
  1 | 0,1 ⟩⟨ 1,0 |
```

---
# Multi-line support
---

With multi-line support, one can write entire chunks of code using the above notation. Just use triple quotes (`"""`) instead of single quotes:

```
julia> d"""
       ψ = 1/√2 * (| 0,0 > + | 1,1 >)
       a = purity(ptrace(ψ*ψ', 2))
       ϕ = normalize!( 1/5 * < 0 | + 4/5 * < 1 | )
       result = normalize!(a * act_on(ϕ, ψ, 2))
       """

julia> result
Ket{KronDelta,1,Float64} with 2 state(s):
  0.24253562503633297 | 0 ⟩
  0.9701425001453319 | 1 ⟩
```

