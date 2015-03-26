*Note: All arithmetic on Kets also works on Bras, although we may not explicitly give examples here.*

---
## Addition, Subtraction, and Scalar multiplication
---

Multiplying a state by a scalar modifies the coefficient appropriately:

```
julia> (1+3.4im) * ket(0)
Ket{Orthonormal,1} with 1 state(s):
  1.0 + 3.4im | 0 ⟩
```

States can be also be added and subtracted:

```
julia> ket(0) + ket(0) == 2 * ket(0)
true

julia> ket(0) - ket(0) == 0 * ket(0)
true

julia> 1/√3 * (ket(0) + ket(1) - ket(2))
Ket{Orthonormal,1} with 3 state(s):
  0.5773502691896258 | 0 ⟩
  -0.5773502691896258 | 2 ⟩
  0.5773502691896258 | 1 ⟩
```

Two key observations can be made here: 

1. The basis states of a Ket are unordered. See the [Analyzing States](analyzing_states.md) section for more info about Bras and Kets as data structures.
2. States do not automatically normalize themselves under operations like addition, which leads us to the next section...

---
## Normalization
---

We can normalize a state in-place by using the `normalize!` function:

```
julia> k = ket(0) + ket(1)
Ket{Orthonormal,1} with 2 state(s):
  1 | 0 ⟩
  1 | 1 ⟩

julia> normalize!(k); k
Ket{Orthonormal,1} with 2 state(s):
  0.7071067811865476 | 0 ⟩
  0.7071067811865476 | 1 ⟩
```

The `normalize` function (without the trailing "!") is also provided, which normalizes a copy of the state instead of modifying the original:

```
julia> k = ket(0) + ket(1)
Ket{Orthonormal,1} with 2 state(s):
  1 | 0 ⟩
  1 | 1 ⟩

julia> normalize(k)
Ket{Orthonormal,1} with 2 state(s):
  0.7071067811865475 | 0 ⟩
  0.7071067811865475 | 1 ⟩

julia> k
Ket{Orthonormal,1} with 2 state(s):
  1 | 0 ⟩
  1 | 1 ⟩
```

You can check the norm of a state using the `norm` function:

```
julia> norm(1/√2 * (ket(0) + ket(1)))
0.9999999999999999
```

---
## Getting a State's Dual
---

One can use the `ctranspose` function to construct the dual of a given state:

```
julia> k = im * ket(0)
Ket{Orthonormal,1} with 1 state(s):
  0 + 1im | 0 ⟩

julia> k'
Bra{Orthonormal,1} with 1 state(s):
  0 - 1im ⟨ 0 |
  
julia> k'' == k
true

```

---
## Tensor Product
---

One can take a tensor product of states simply by multiplying them:

```
julia> ket(0) * ket(0)
Ket{Orthonormal,2} with 1 state(s):
  1 | 0,0 ⟩
```

As you might notice, a tensor product of Kets is itself a Ket, and the result
of the above is the same as if we input `ket(0,0)`. Taking the tensor product of 
more complicated states illustrates the tensor product's cartesian properties:

```
julia> normalize!(sum(i->i^2*ket(i), 0:3) * sum(i->i/2*ket(i), -3:3))
Ket{Orthonormal,2} with 28 state(s):
  -0.5154323951168185 | 3,-3 ⟩
  0.0 | 1,0 ⟩
  0.3436215967445456 | 3,2 ⟩
  -0.019090088708030313 | 1,-1 ⟩
  0.22908106449636376 | 2,3 ⟩
  -0.07636035483212125 | 2,-1 ⟩
  0.019090088708030313 | 1,1 ⟩
  0.0 | 2,0 ⟩
  -0.22908106449636376 | 2,-3 ⟩
  0.1718107983722728 | 3,1 ⟩
  -0.0 | 0,-1 ⟩
  -0.05727026612409094 | 1,-3 ⟩
  0.0 | 0,2 ⟩
  0.05727026612409094 | 1,3 ⟩
  0.0 | 0,0 ⟩
  -0.1718107983722728 | 3,-1 ⟩
  ⁞
```

Note that in the above, we used Julia's [`sum`](http://julia.readthedocs.org/en/latest/stdlib/collections/?highlight=sum#Base.sum) function to quickly construct superpositions of states where the coefficients were a function of the labels (and scaled proportionally in the final normalized product). 

---
## Inner Product
---

Similarly to the tensor product, the inner product can be taken simply by multiplying Bras with Kets:

```
julia> bra(0)*ket(1)
0

julia> k = 1/√2 * (ket(0,0) + ket(1,1)); k' * k
0.9999999999999998

julia> bra(0,0) * k
0.7071067811865475
```

It is sometimes useful to take the inner product between a Bra and a *specific factor* of a Ket. A math example might be:

```
| ψ ⟩ =   c₁ | 0, 1 ⟩ + c₂ | 1, 0 ⟩

⟨ 0₂ | ψ ⟩ =  c₁ ⟨ 0₂ | 0, 1 ⟩ + c₂ ⟨ 0₂ | 1, 0 ⟩

⟨ 0₂ | ψ ⟩ =  c₁ ⟨ 0 | 1 ⟩| 0 ⟩ + c₂ ⟨ 0 | 0 ⟩| 1 ⟩
```

If these states are orthonormal, our final result is

```
⟨ 0₂ | ψ ⟩ = 0 | 0 ⟩ + c₂ | 1 ⟩ 
```

QuDirac supports this operation through the use of the `inner` function:

```
julia> ψ = normalize!(ket(0,1) + 2*ket(1,0))
Ket{Orthonormal,2} with 2 state(s):
  0.4472135954999579 | 0,1 ⟩
  0.8944271909999159 | 1,0 ⟩

julia> inner(bra(0), ψ, 2) # ⟨ 0₂ | ψ ⟩
Ket{Orthonormal,1} with 2 state(s):
  0.0 | 0 ⟩
  0.8944271909999159 | 1 ⟩

julia> inner(bra(0), ψ, 1) # ⟨ 0₁ | ψ ⟩
Ket{Orthonormal,1} with 2 state(s):
  0.0 | 0 ⟩
  0.4472135954999579 | 1 ⟩
```

This does, of course, work even when the Bra is a superposition of states:

```
julia> ϕ = 1/√2 * (bra(0) + bra(1))
Bra{Orthonormal,1} with 2 state(s):
  0.7071067811865475 ⟨ 0 |
  0.7071067811865475 ⟨ 1 |

julia> inner(ϕ, ψ, 2)
Ket{Orthonormal,1} with 2 state(s):
  0.3162277660168379 | 0 ⟩
  0.6324555320336758 | 1 ⟩
```

As you can see, the above calculations assume an *orthonormal* inner product for the involved states. This behavior is stored in the state's type information (e.g. `Orthonormal` in `Ket{Orthonormal,1}`), and you may notice that the `bra`/`ket` functions construct states with product type `Orthonormal` by default. QuDirac also has support for arbitrary, lazily evaluated inner products when states aren't simply orthonormal. To learn more, see the [Custom Inner Products](custom_inner_products.md) section.

---
## Outer Product
---

TODO
