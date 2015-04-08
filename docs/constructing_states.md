*Note: Throughout this documentation, I'll refer to `ket` (lowercase) and Ket (uppercase). "Ket" refers to a type, while `ket` refers to a convenience constructor for that type. The same goes for `bra` and Bra.*

---
# Constructing Single Kets
---

To begin, let's make a single state using the `ket` function:


```
julia> ket(0)
Ket{KroneckerDelta,1} with 1 state(s):
  1 | 0 ⟩
```

As you can see, the `ket` function takes labels (in this case, a single zero) as arguments. In QuDirac, anything can be used as a Ket label - primitives, strings, composite types, etc. Simply pass the desired label in as we did `0` above:

```
julia> ket(":)")
Ket{KroneckerDelta,1} with 1 state(s):
  1 | ":)" ⟩

julia> ket(:a)
Ket{KroneckerDelta,1} with 1 state(s):
  1 | :a ⟩

julia> ket([1,2,3])
Ket{KroneckerDelta,1} with 1 state(s):
  1 | [1,2,3] ⟩
```

---
# Constructing Single Bras
---

Bras can be constructed the same way using the `bra` function:

```
julia> bra(0)
Bra{KroneckerDelta,1} with 1 state(s):
  1 ⟨ 0 |
```

Just like Kets, Bra labels can be anything.  


---
# Constructing Single Product States
---

The number of labels passed to the `ket`/`bra` functions determines how many factors there are in the basis of the resulting state. For example, to construct `| 0 ⟩ ⊗| 0 ⟩ ⊗ | 0 ⟩` we can simply do the following:


```
julia> k = ket(0,0,0)
Ket{KroneckerDelta,3} with 1 state(s):
  1 | 0,0,0 ⟩

julia> nfactors(k)
3
```

Just as with single labels, one is free to use labels of any type for multi-factor states:

```
julia> k = ket(0, ":)", :a, [1,2,3])
Ket{KroneckerDelta,4} with 1 state(s):
  1 | 0,":)",:a,[1,2,3] ⟩

julia> nfactors(k)
4
```
