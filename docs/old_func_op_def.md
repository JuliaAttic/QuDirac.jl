# Constructing operators by defining their action on a basis
---

Mathematically, one can represent an operator `Ô` with the following definition:

```
Ô | i ⟩ = ∑ⱼ cᵢⱼ | j ⟩
```

QuDirac supports the construction of operators in the above manner using the `@repr_op` macro.

This usage is easiest to understand by example, so let's construct a Hadamard operator:

```
julia> flip(label::StateLabel) = 1/√2 * (ket(0) + (-1)^label[1] * ket(1))
flip (generic function with 1 method)

julia> op = repr_op(flip, d"| 0 > + | 1 >")

julia> op * ket(1)
Ket{KroneckerDelta,1,Float64} with 2 state(s):
  0.7071067811865475 | 0 ⟩
  -0.7071067811865475 | 1 ⟩

julia> op * ket(0)
Ket{KroneckerDelta,1,Float64} with 2 state(s):
  0.7071067811865475 | 0 ⟩
  0.7071067811865475 | 1 ⟩

julia> op * (ket(0) + ket(1))
Ket{KroneckerDelta,1,Float64} with 2 state(s):
  1.414213562373095 | 0 ⟩
  0.0 | 1 ⟩
```

As you can see, the `repr_op` function takes two arguments:

1. A function that takes in a `StateLabel` and returns a Ket.
2. A Ket, in whose basis the operator will be represented.

Thus, this expression:

```
Ô = repr_op(f, labels)
```
...is mathematically equivalent to this definition:

```
Ô | l ∈ labels ⟩ = f(l) 
```

...where `f` returns a Ket.

---
# Permutative Operator Representation
---

If the operator you wish to construct can be represented as a generalized permutation matrix, one can use 
the more efficient `repr_permop` function instead of the `repr_op` function. "Generalized permutation", in 
this case, means that the operator representation can be defined as:

```
Ô | i ⟩ = cᵢⱼ | j ⟩
```

...where both `| i ⟩` and `| j ⟩` are basis states, i.e. *not* superpositional states. Common examples of 
operators that adhere to this definition are the identity operator and the ladder operators for the 
quantum harmonic oscillator. 

The use of `repr_permop` is a little different from the more general `repr_op`, so let's 
construct the lowering operator as an example. The lowering operator `â` is defined as 

```
â | n ⟩ = √n | n - 1 ⟩
``` 

Let's say we wish to define `â` on the basis of the following state:

```julia
julia> k = normalize(sum(ket, 1:5))
Ket{KroneckerDelta,1,Float64} with 5 state(s):
  0.4472135954999579 | 4 ⟩
  0.4472135954999579 | 3 ⟩
  0.4472135954999579 | 2 ⟩
  0.4472135954999579 | 0 ⟩
  0.4472135954999579 | 1 ⟩
```

To do so, we construct a function that applies to each basis label, and pass it
to `repr_permop` along with `k`:

```julia
julia> lower(n::StateLabel) = (sqrt(n[1]), StateLabel(n[1]-1))
lower (generic function with 1 method)

julia> â = repr_permop(lower, k)
OpSum{KroneckerDelta,1,Float64} with 5 operator(s):
  2.23606797749979 | 4 ⟩⟨ 5 |
  2.0 | 3 ⟩⟨ 4 |
  1.0 | 0 ⟩⟨ 1 |
  1.4142135623730951 | 1 ⟩⟨ 2 |
  1.7320508075688772 | 2 ⟩⟨ 3 |
```

As you can see, the output of `lower` is structured differently than the output of functions 
passed to `repr_op`. A function passed to `repr_permop` takes in a `StateLabel` and return 
a (coefficient, label) pair of type `(T, StateLabel)`, where `T` is the coefficient type 
of the resulting operator. 

Thus, mixing Julia syntax and math, our functional construction of `â` could be written out like this:

```
for n ∈ k, 

construct â such that

â | n ⟩ = lower(n)[1] | lower(n)[2] ⟩
``` 

Our construction of `â` works as expected:

```julia
julia> â*k
Ket{KroneckerDelta,1,Float64} with 5 state(s):
  1.0 | 4 ⟩
  0.8944271909999159 | 3 ⟩
  0.7745966692414833 | 2 ⟩
  0.4472135954999579 | 0 ⟩
  0.6324555320336759 | 1 ⟩

julia> k' * â
Bra{KroneckerDelta,1,Float64} with 4 state(s):
  0.8944271909999159 ⟨ 4 |
  0.7745966692414833 ⟨ 3 |
  0.6324555320336759 ⟨ 2 |
  1.0 ⟨ 5 |

julia> â*ket(4)
Ket{KroneckerDelta,1,Float64} with 1 state(s):
  2.0 | 3 ⟩

julia> â'*ket(4)
Ket{KroneckerDelta,1,Float64} with 1 state(s):
  2.23606797749979 | 5 ⟩

julia> bra(4) * â
Bra{KroneckerDelta,1,Float64} with 1 state(s):
  2.23606797749979 ⟨ 5 |

julia> bra(4) * â'
Bra{KroneckerDelta,1,Float64} with 1 state(s):
  2.0 ⟨ 3 |
```

Another example is the factor switching operator `Â₁₃₂`, defined as: 

```
Â₁₃₂ | a, b, c ⟩  = | a, c, b ⟩
```

We can construct this operator thusly:

```julia
julia> k = sum(i -> i * ket(i,i+1,i+2), 1:10)
Ket{KroneckerDelta,3,Int64} with 10 state(s):
  4 | 4,5,6 ⟩
  1 | 1,2,3 ⟩
  2 | 2,3,4 ⟩
  5 | 5,6,7 ⟩
  6 | 6,7,8 ⟩
  7 | 7,8,9 ⟩
  3 | 3,4,5 ⟩
  9 | 9,10,11 ⟩
  10 | 10,11,12 ⟩
  8 | 8,9,10 ⟩

julia> Â₁₃₂ = repr_permop(label -> (1, StateLabel(label[[1,3,2]])), k)
OpSum{KroneckerDelta,3,Int64} with 10 operator(s):
  1 | 2,4,3 ⟩⟨ 2,3,4 |
  1 | 4,6,5 ⟩⟨ 4,5,6 |
  1 | 7,9,8 ⟩⟨ 7,8,9 |
  1 | 10,12,11 ⟩⟨ 10,11,12 |
  1 | 1,3,2 ⟩⟨ 1,2,3 |
  1 | 3,5,4 ⟩⟨ 3,4,5 |
  1 | 5,7,6 ⟩⟨ 5,6,7 |
  1 | 6,8,7 ⟩⟨ 6,7,8 |
  1 | 9,11,10 ⟩⟨ 9,10,11 |
  1 | 8,10,9 ⟩⟨ 8,9,10 |
  
julia> Â₁₃₂ * k
Ket{KroneckerDelta,3,Int64} with 10 state(s):
  2 | 2,4,3 ⟩
  9 | 9,11,10 ⟩
  7 | 7,9,8 ⟩
  5 | 5,7,6 ⟩
  10 | 10,12,11 ⟩
  6 | 6,8,7 ⟩
  4 | 4,6,5 ⟩
  1 | 1,3,2 ⟩
  3 | 3,5,4 ⟩
  8 | 8,10,9 ⟩
```

Note that, in practice, QuDirac provides the `switch` and `permute` functions for states and operators, 
which can imitate the above behavior in functional form (and is not restricted to a basis)

```
julia> Â₁₃₂ * k == switch(k, 2, 3) == permute(k, [1, 3, 2])
true
```