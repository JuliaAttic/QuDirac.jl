QuDirac supports the functional construction of operators that can be defined in the form:

```
Ô | a ⟩ = c_ab | b ⟩
```

Where `| a ⟩` and `| b ⟩` are single basis states (i.e. not superpositions).

---
# Constructing operators by defining their action on a basis
---

Operators are often defined by their actions on individual Kets. The quintessential example is the lowering operator, defined
on Kets as:

```
â | n ⟩ = √n | n - 1 ⟩
``` 

To define an operator in this manner using QuDirac, we can use the `func_op` method, 
which takes a function, a Ket, and a factor index as arguments. 

Let's say we have a simple Ket `k`:

```
julia> k = normalize!(sum(ket, 0:4))
Ket{Orthonormal,1} with 5 state(s):
  0.4472135954999579 | 0 ⟩
  0.4472135954999579 | 2 ⟩
  0.4472135954999579 | 3 ⟩
  0.4472135954999579 | 4 ⟩
  0.4472135954999579 | 1 ⟩
```

Now we can define our lowering operator in `k`'s basis:

```
julia> f(label) = (sqrt(label), label-1)

julia> â = func_op(f, k, 1)
GenericOp{Orthonormal,1} with 4 operator(s):
  1.0 | 0 ⟩⟨ 1 |
  2.0 | 3 ⟩⟨ 4 |
  1.7320508075688772 | 2 ⟩⟨ 3 |
  1.4142135623730951 | 1 ⟩⟨ 2 |
```

Let's examine the arguments we passed to the `GenericOp` constructor above. 

First, we have a function `f` defined as `label->(sqrt(label),label-1)`. As you can see, this
function takes in a single argument (a label) and returns two things: a transformation coefficient, 
and the resultant label. The labels passed to this function are simply the labels of the basis
states of `k`. We passed in `1` to denote that we wanted this function to act on the
first (in this case, only) factor of the labels of `k`.

Thus, mixing Julia syntax and math, our functional construction of `â` can be written out like this:

```
for n ∈ k, 

where f(label) = (sqrt(label),label-1):

construct â such that

â | n ⟩ = f(n)[1] | f(n)[2] ⟩
``` 

As we can see, this operator we defined works as expected:

```
julia> â*k
Ket{Orthonormal,1} with 4 state(s):
  0.4472135954999579 | 0 ⟩
  0.8944271909999159 | 3 ⟩
  0.7745966692414833 | 2 ⟩
  0.6324555320336759 | 1 ⟩

julia> â*ket(4)
Ket{Orthonormal,1} with 1 state(s):
  2.0 | 3 ⟩
```

---
# Functionally defining operators on single factors
---

In the previous example, we passed in `1` to the `GenericOp` constructor to indicate that our function `f`
acted on the first factor of the labels in `k`'s basis. This is rather boring, since `k` only had one factor 
to begin with!

To demonstrate factor-specific operator defintions in a more practical context, let's construct a lowering 
operator on the second factor of a three-factor basis.

We'll just reuse `k` from the previous example to construct our basis labels:

```
julia> k3 = k^3
Ket{Orthonormal,3} with 125 state(s):
  0.08944271909999157 | 1,4,1 ⟩
  0.08944271909999157 | 3,4,4 ⟩
  0.08944271909999157 | 2,3,1 ⟩
  0.08944271909999157 | 4,2,0 ⟩
  0.08944271909999157 | 0,1,2 ⟩
  0.08944271909999157 | 2,3,2 ⟩
  ⁞
```

...and we can reuse `f` as well, but define the operator on the second factor by passing in `2`:

```
julia> â₂ = func_op(f, k3, 2)
GenericOp{Orthonormal,3} with 100 operator(s):
  1.4142135623730951 | 3,1,0 ⟩⟨ 3,2,0 |
  1.0 | 3,0,3 ⟩⟨ 3,1,3 |
  1.4142135623730951 | 4,1,2 ⟩⟨ 4,2,2 |
  1.7320508075688772 | 3,2,4 ⟩⟨ 3,3,4 |
  1.7320508075688772 | 4,2,0 ⟩⟨ 4,3,0 |
  1.0 | 1,0,2 ⟩⟨ 1,1,2 |
  ⁞
```

Testing it out:

```
julia> â₂ * k3
Ket{Orthonormal,3} with 100 state(s):
  0.17888543819998315 | 2,3,1 ⟩
  0.15491933384829665 | 4,2,0 ⟩
  0.12649110640673517 | 0,1,2 ⟩
  0.17888543819998315 | 3,3,1 ⟩
  0.17888543819998315 | 2,3,2 ⟩
  0.17888543819998315 | 2,3,4 ⟩
  ⁞

julia> â₂ * ket(1,2,3)
Ket{Orthonormal,3} with 1 state(s):
  1.4142135623730951 | 1,1,3 ⟩
```

---
# Functionally defining operators on all factors
---

If you wish to functionally define an operator on all factors of the basis, you may easily do so by omitting the factor index as an argument:

```
julia> k = sum(i -> i * ket(i,i+1,i+2), 1:10)
Ket{Orthonormal,3} with 10 state(s):
  1 | 1,2,3 ⟩
  5 | 5,6,7 ⟩
  3 | 3,4,5 ⟩
  8 | 8,9,10 ⟩
  4 | 4,5,6 ⟩
  10 | 10,11,12 ⟩
  9 | 9,10,11 ⟩
  6 | 6,7,8 ⟩
  2 | 2,3,4 ⟩
  7 | 7,8,9 ⟩

julia> p132 = func_op(label -> (1, switch(label, 2, 3)), k)
GenericOp{Orthonormal,3} with 10 operator(s):
  1 | 6,8,7 ⟩⟨ 6,7,8 |
  1 | 4,6,5 ⟩⟨ 4,5,6 |
  1 | 2,4,3 ⟩⟨ 2,3,4 |
  1 | 8,10,9 ⟩⟨ 8,9,10 |
  1 | 1,3,2 ⟩⟨ 1,2,3 |
  1 | 9,11,10 ⟩⟨ 9,10,11 |
  1 | 5,7,6 ⟩⟨ 5,6,7 |
  1 | 10,12,11 ⟩⟨ 10,11,12 |
  1 | 3,5,4 ⟩⟨ 3,4,5 |
  1 | 7,9,8 ⟩⟨ 7,8,9 |
```

*Note: The `label` argument above is a `Vector`. See the [Labels and coefficients](labels_and_coeffs.md) section for details.*

As you can see, the above operator simply switches the second and third factors of 
the basis labels (the `switch` function is provided by QuDirac, and can be found in
the API section).

Acting it on `k`:

```
julia> p132 * k
Ket{Orthonormal,3} with 10 state(s):
  4 | 4,6,5 ⟩
  9 | 9,11,10 ⟩
  6 | 6,8,7 ⟩
  7 | 7,9,8 ⟩
  5 | 5,7,6 ⟩
  2 | 2,4,3 ⟩
  1 | 1,3,2 ⟩
  8 | 8,10,9 ⟩
  10 | 10,12,11 ⟩
  3 | 3,5,4 ⟩
```

