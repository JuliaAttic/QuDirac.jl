# Bras and Kets

This section provides an introduction to working with Bras and Kets, the building blocks of Dirac notation. 

*Note: Throughout this documentation, I'll refer to `ket`s (lowercase) and Kets (uppercase). "Ket" refers to a type, while `ket` refers to the convenience constructor for that type. The same goes for `bra` and "Bra".*

## 1. Constructing Simple Bras and Kets

### 1.1 `ket` Constructor
To begin, let's make a single state using the `ket` function:

```
julia> using QuDirac

julia> ket(0)
Ket{Orthonormal,1} with 1 state(s):
  1 | 0 ⟩
```

As you can see, the `ket` function takes labels (in this case, a single zero) as arguments. 

Using QuDirac.jl, *ANYTHING can be used as a Ket label* - primitives, `String`s, composite types, and even other QuDirac objects. Simply pass the label in as we did `0` above:

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

### 1.2 `bra` Constructor


Bras can be constructed the same way using the `bra` function:

```
julia> bra(0)
Bra{Orthonormal,1} with 1 state(s):
  1 ⟨ 0 |
```

Just like Kets, Bra labels can be anything.  


### 1.3 Multi-factor states

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

## 2. Math with States

*Note: All arithmetic on Kets also works on Bras, though we may not explicitly give examples here.*

### 2.1 Addition, Subtraction, and Scalar multiplication

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

1. The basis states of a Ket are unordered. See the [States as Data Structures](States_as_Data_Structures.md) section below.
2. States do not automatically normalize themselves under operations like addition, which leads us to...


### 2.2 Normalization

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

A `normalize` function (without the trailing "!") is also provided, which normalizes a copy of the state instead of modifying the original:

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

### 2.3 Getting a State's Dual

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

### 2.4 Tensor Product

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

Note that in the above, we used Julia's [`sum`](http://julia.readthedocs.org/en/latest/stdlib/collections/?highlight=sum#Base.sum) function to quickly construct superpositions of states where the coefficients were a function of the labels (and thus scaled proportionally in the final normalized product). 

### 2.5 Inner Product

Similarly to the tensor product, the inner product can be taken simply by multiplying Bras with Kets:

```
julia> bra(0)*ket(1)
0

julia> k = 1/√2 * (ket(0,0) + ket(1,1)); k' * k
0.9999999999999998

julia> bra(0,0) * k
0.7071067811865475
```

As you can see, the above calculations assume an *orthonormal* inner product. This behavior is stored in the state's type information (e.g. `Ket{Orthonormal,1}`), and you may notice that the `bra`/`ket` functions construct states with product type `P<:Orthonormal` by default. 

QuDirac.jl has support for arbitrary, lazily evaluated inner products as well. To learn more, see the [Custom Inner Products](Custom_Inner_Products.md) section.

### 2.6 Outer Product




