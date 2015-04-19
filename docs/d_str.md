# Natural Dirac notation input syntax
---

QuDirac supports a "natural" input format by implementing `d"..."` syntax. Seeing this syntax in action:

```julia
julia> d" | 0 > "
Ket{KroneckerDelta,1,Int64} with 1 state(s):
  1 | 0 ⟩
```

Using `d"..."` calls the macro `@d_str`. This macro parses the given string for the `|`, `>` and `<` characters, replacing them where appropriate with the `ket` and `bra` functions. The result is then parsed into an `Expr`, which is returned back to the calling environment and evaluated. For more on the behavior of Julia's non-standard string literals, see [the official docs](http://julia.readthedocs.org/en/release-0.3/manual/metaprogramming/#non-standard-string-literals).

Note that when using `d"..."` syntax, the `|`, `>` and `<` characters *cannot be used for anything other than as indicators of Kets and Bras*. Assignments and function calls work as expected, though:

```julia
julia> d" A = tensor( | 0 >< 1 |, | 1 >< 0 | ) ";

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
Ket{KroneckerDelta,1,Float64} with 2 state(s):
  0.24253562503633297 | 0 ⟩
  0.9701425001453319 | 1 ⟩
```

