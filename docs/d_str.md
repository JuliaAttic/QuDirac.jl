# Natural Dirac notation input syntax
---

QuDirac supports an input format that is even more natural than using the `ket` or `bra` functions by implementing `d" ... "` and `d""" ... """` syntax. The former is for single-line input, while the latter is for multi-line input.

Seeing this syntax in action:

```julia
julia> d" | 0 > "
Ket{KroneckerDelta,1,Int64} with 1 state(s):
  1 | 0 ⟩
```

The macro parses the string for the `|`, `>` and `<` characters, replacing them where appropriate with the `ket` and `bra` functions. Thus, when using this syntax, the `|`, `>` and `<` characters *cannot be used for anything other than as indicators of Kets and Bras*.

Assignments and function calls work as expected, though:

```julia
julia> @d_str " ψ = 1/√2 * (| 0,0 > + | 1,1 >); purity(ptrace(ψ*ψ', 2)) "
0.4999999999999998
```

---
# Multi-line syntax
---

With multi-line support, one can write entire chunks of code using the above notation. Just wrap the code in `d""" ... """`:

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

