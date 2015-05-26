# Defining Operators as Functions
---

Mathematically, one can represent an operator `Ô` with the following definitions:

```
Action of Ô on Ket:
Ô | i ⟩ = ∑ⱼ cᵢⱼ | j ⟩ 

Action of Ô on Bra:
⟨ u | Ô  = ∑ᵥ cᵤᵥ ⟨ v |
```

QuDirac allows us to define normal Julia functions that act like `Ô` using 
the `@def_op` macro:

```julia
# define "a" on Kets
julia> @def_op " a | n > = √n * | n-1 > "
a (generic function with 1 methods)

julia> d" a * | 42 > "
Ket{KronDelta,1,Float64} with 1 state(s):
  6.48074069840786 | 41 ⟩

julia> d" a * (| 1 > + | 2 > - | 3 >) "
Ket{KronDelta,1,Float64} with 3 state(s):
  -1.7320508075688772 | 2 ⟩
  1.0 | 0 ⟩
  1.4142135623730951 | 1 ⟩
```

Action of the operator's dual on Bras is well-defined by the definition of 
the operator on Kets, since ` < 1 | * Ô' == ( Ô * | 1 > )' `:

```
julia> d" < 2 | * a' "
Bra{KronDelta,1,Float64} with 1 state(s):
  1.4142135623730951 ⟨ 1 |
```

To act the operator on a Bra (or its dual on a Ket), we need to define its action on the Bra:

```
# define "a" on Bras
julia> @def_op " < n | a = √(n+1) * < n+1 | "
a (generic function with 2 methods)

julia> d" < 5 | * a "
Bra{KronDelta,1,Float64} with 1 state(s):
  2.449489742783178 ⟨ 6 |

julia> d" a' * (| 5 > + | 2 >)"
Ket{KronDelta,1,Float64} with 2 state(s):
  1.7320508075688772 | 3 ⟩
  2.449489742783178 | 6 ⟩

julia> d" < 5 | * a * | 6 > "
2.449489742783178
```

The `@def_op` macro works for product bases as well:

```
julia> @def_op " a₂ | x,y,z >  = √y * | x,y-1,z > "
a₂ (generic function with 2 methods)

julia> d" a₂ * (| 0,10,8 > - | 12,31,838 >) "
Ket{KronDelta,3,Float64} with 2 state(s):
  3.1622776601683795 | 0,9,8 ⟩
  -5.5677643628300215 | 12,30,838 ⟩
```

For an example of an operator that throws its basis Kets into superpositions, 
here's a function emulating a Hadamard operator:

```julia
julia> @def_op " h | n > = 1/√2 * ( | 0 > + (-1)^n *| 1 > )"
h (generic function with 1 methods)

julia> d" h * | 0 > "
Ket{KronDelta,1,Float64} with 2 state(s):
  0.7071067811865475 | 0 ⟩
  0.7071067811865475 | 1 ⟩

julia> d" h * | 1 > "
Ket{KronDelta,1,Float64} with 2 state(s):
  0.7071067811865475 | 0 ⟩
  -0.7071067811865475 | 1 ⟩
```

---
# Grammar for the Definition String
---

The grammar of the string passed to `@def_op` is:

1. Defining action on Kets:
  
        @def_op " $op_name | $label_args > = f($label_args...) "

      where `f` is an arbitrary expanded function that takes in the `$label_args` and
      returns a Ket.

2. Defining action on Bras:

        @def_op " < $label_args | $op_name  = f($label_args...) "
      
      where `f` is an arbitrary expanded function that takes in the `$label_args` and
      returns a Bra.

Allowable syntax for the right-hand side of the equation
is [exactly the same syntax allowed by `d"..."`](d_str.md).

---
# Generating Operator Representations
---

The operator-functions described in the previous example are no doubt useful, 
but they are just normal Julia functions, and so are quite limited when it comes 
to mimicking the behavior of *actual* quantum operators. For example, they can't 
be added or be factors of a tensor product (this functionality may indeed be implemented
in the future, however).

For those capabilities, we'll need to generate an `OuterSum` representation in a basis. To do so, 
we can use the `@rep_op` macro. Here's a familiar example:

```julia
julia> @rep_op " a | n > = √n * | n-1 > " 1:10;

julia> a
OuterSum{KronDelta,1,Float64} with 10 operator(s):
  2.449489742783178 | 5 ⟩⟨ 6 |
  3.1622776601683795 | 9 ⟩⟨ 10 |
  2.23606797749979 | 4 ⟩⟨ 5 |
  2.8284271247461903 | 7 ⟩⟨ 8 |
  2.6457513110645907 | 6 ⟩⟨ 7 |
  2.0 | 3 ⟩⟨ 4 |
  1.0 | 0 ⟩⟨ 1 |
  3.0 | 8 ⟩⟨ 9 |
  1.4142135623730951 | 1 ⟩⟨ 2 |
  1.7320508075688772 | 2 ⟩⟨ 3 |
```

The `@rep_op` macro takes in a definition string, and an iterable of items to be used as basis labels.
The grammar and allowable syntax of the definition string is *exactly* that of the definition string passed 
to `@def_op`. The only difference between the two is that the `@rep_op` macro feeds in the given basis labels
to produce an `OuterSum`.

To generate a representation on a product basis, one can provide multiple iterables to `@rep_op`.
Their cartesian product will then be used as the basis for the representation: 

```julia
# define P₁₃₂ by it's action on a Bra
julia> @rep_op " < i,j,k | P₁₃₂  = < i,k,j | "  1:10  'a':'f'  -(1:10)
OuterSum{KronDelta,3,Int64} with 600 operator(s):
  1 | 10,'e',-2 ⟩⟨ 10,-2,'e' |
  1 | 4,'c',-5 ⟩⟨ 4,-5,'c' |
  1 | 3,'c',-10 ⟩⟨ 3,-10,'c' |
  1 | 2,'e',-3 ⟩⟨ 2,-3,'e' |
  1 | 7,'c',-9 ⟩⟨ 7,-9,'c' |
  1 | 1,'f',-9 ⟩⟨ 1,-9,'f' |
  1 | 3,'f',-6 ⟩⟨ 3,-6,'f' |
  1 | 9,'b',-2 ⟩⟨ 9,-2,'b' |
  1 | 1,'c',-3 ⟩⟨ 1,-3,'c' |
  1 | 4,'d',-8 ⟩⟨ 4,-8,'d' |

julia> d" P₁₃₂ * | 10,-2,'e' > "
Ket{KronDelta,3,Int64} with 1 state(s):
  1 | 10,'e',-2 ⟩
```

Let's say I already have an operator-function defined, and want
to represent it in a basis. For example, take the Hadamard operator-function 
`h`, constructed in the previous section as:

```julia
julia> @def_op " h | n > = 1/√2 * ( | 0 > + (-1)^n *| 1 > )"
h (generic function with 1 methods)
```

I can easily generate a representation for this function by using `@rep_op` and 
calling `h` on the right-hand side:

```julia
julia> @rep_op " H | n > = h * | n > " 0:1
OuterSum{KronDelta,1,Float64} with 4 operator(s):
  0.7071067811865475 | 1 ⟩⟨ 0 |
  0.7071067811865475 | 0 ⟩⟨ 0 |
  0.7071067811865475 | 0 ⟩⟨ 1 |
  -0.7071067811865475 | 1 ⟩⟨ 1 |
```

The above strategy works to represent any operator-function. Just be aware that
the function and the actual representation need to have unique names:

```julia
julia> @rep_op " h | n > = h * | n > " 0:1
ERROR: invalid redefinition of constant h
 in anonymous at no file:70
```
