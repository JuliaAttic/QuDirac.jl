## Constructing Single Kets
---

To begin, let's make a single state using the `ket` function:


```
julia> ket(0)
Ket{Orthonormal,1} with 1 state(s):
  1 | 0 ⟩
```

As you can see, the `ket` function takes labels (in this case, a single zero) as arguments. In QuDirac.jl, anything can be used as a Ket label - primitives, strings, composite types, and even other QuDirac objects. Simply pass the desired label in as we did `0` above:

```
julia> ket(":)")
Ket{Orthonormal,1} with 1 state(s):
  1 | ":)" ⟩

julia> ket(:a)
Ket{Orthonormal,1} with 1 state(s):
  1 | :a ⟩

julia> ket([1,2,3])
Ket{Orthonormal,1} with 1 state(s):
  1 | [1,2,3] ⟩
```

---
## Constructing Single Bras
---

Bras can be constructed the same way using the `bra` function:

```
julia> bra(0)
Bra{Orthonormal,1} with 1 state(s):
  1 ⟨ 0 |
```

Just like Kets, Bra labels can be anything.  


---
## Constructing Single Product States
---

The number of labels passed to the `ket`/`bra` functions determines how many factors there are in the basis of the resulting state. For example, to construct `| 0 ⟩ ⊗| 0 ⟩ ⊗ | 0 ⟩` we can simply do the following:


```
julia> k = ket(0,0,0)
Ket{Orthonormal,3} with 1 state(s):
  1 | 0,0,0 ⟩

julia> nfactors(k)
3
```

The number of factors is encoded in the type information of a state (e.g. the `3` in `Ket{Orthonormal, 3}` above) and can be retrieved using the `nfactors` function.

Just as with single labels, one is free to use labels of any type for multi-factor states:

```
julia> k = ket(0, ":)", :a, [1,2,3])
Ket{Orthonormal,4} with 1 state(s):
  1 | 0,":)",:a,[1,2,3] ⟩

julia> nfactors(k)
4
```
