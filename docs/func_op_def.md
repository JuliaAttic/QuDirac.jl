# Functionally Defining Operators
---

Mathematically, one can represent an operator `Ô` with the following definition:

```
Ô | i ⟩ = ∑ⱼ cᵢⱼ | j ⟩
```

QuDirac supports the construction of operators in the above manner using the `@repr_op` macro.
