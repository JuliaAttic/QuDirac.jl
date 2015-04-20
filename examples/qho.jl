using QuDirac
##################
# Lower Operator #
##################
# This function can be applied to Kets like "a * k"
@def_op " a | n > = √n * | n - 1 > "

#################################
# Hermite Polynomial Evaluation #
#################################
# This is a very naive implementation, 
# likely extremely prone to numerical error. 
# It should suffice for this example, however.
# Implements the expression found here: 
# http://en.wikipedia.org/wiki/Hermite_polynomials#Explicit_expression

hermite_term(m,n,x) = ((-1)^m / (factorial(m) * factorial(n - 2m))) * (2x)^(n-2m)
hermite(n::BigInt, x::Float64) = factorial(n) * sum(m -> hermite_term(BigInt(m),n,x), 0:floor(n/2))

############
# QHOInner #
############
# Custom inner 
immutable QHOInner <: AbstractInner end

immutable PosX
    val::Float64
end

# using natural units
qho_inner(x::PosX, k::BigInt) = e^((-x.val^2)/2) * hermite(k, x.val) * 1/√(2^k * factorial(k) * √π)
qho_inner(x::PosX, k::Int) = qho_inner(x, BigInt(k))
qho_inner(n, m) = inner_rule(KroneckerDelta(), n, m)

QuDirac.inner_rule(::QHOInner, b, k) = reduce(*, map(qho_inner, b, k))
QuDirac.inner_rettype(::QHOInner) = BigFloat

default_inner(QHOInner())
