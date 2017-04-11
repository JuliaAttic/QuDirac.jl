##############
# InnerLabel #
##############
immutable InnerProduct{N} <: Number
    b::StateLabel{N}
    k::StateLabel{N}
end

blabel(i::InnerProduct) = i.b
klabel(i::InnerProduct) = i.k

Base.repr(i::InnerProduct) = brstr(blabel(i))*ktstr(klabel(i))[2:end]
Base.show(io::IO, i::InnerProduct) = print(io, repr(i))

Base.:(==)(a::InnerProduct, b::InnerProduct) = blabel(a) == blabel(b) && klabel(a) == klabel(b)
Base.:(==)(::InnerProduct, ::Number) = false
Base.:(==)(::Number, ::InnerProduct) = false

Base.hash(i::InnerProduct) = hash(blabel(i), hash(klabel(i)))
Base.hash(i::InnerProduct, h::UInt64) = hash(hash(i), h)

Base.conj(i::InnerProduct) = InnerProduct(klabel(i), blabel(i))

##############
# inner_rule #
##############

# we can cheat here to avoid tempory StateLabel construction
# to evaluate KroneckerDelta inner products
eval_inner_rule(p::KroneckerDelta, b, k) = inner_rule(p, b, k)
eval_inner_rule(p::KroneckerDelta, b::StateLabel, k::StateLabel) = inner_rule(p, b, k)

eval_inner_rule(p::AbstractInner, b::StateLabel, k::StateLabel) = inner_rule(p, b, k)
eval_inner_rule(p::AbstractInner, b, k) = inner_rule(p, StateLabel(b), StateLabel(k))

inner_rule{P<:AbstractInner}(::P, b, k) = error("inner_rule(::$P, ::StateLabel, ::StateLabel) must be defined to evaluate inner products with type $P")
inner_rule(::UndefinedInner, b, k) = InnerExpr(InnerProduct(b, k))
inner_rule(::KroneckerDelta, b, k) = b == k ? 1 : 0

inner_rettype{P<:AbstractInner}(::P) = first(Base.return_types(inner_rule, (P, StateLabel, StateLabel)))
inner_rettype(::UndefinedInner) = InnerExpr
inner_rettype(::KroneckerDelta) = Int64

# very common operation for state/operator inner products
inner_mul(v,c,prodtype,b,k) = v * c * eval_inner_rule(prodtype, b, k)

##############
# InnerExpr #
##############
# A InnerExpr is a type that wraps arthimetic expressions
# performed with InnerExprs. The allows storage and
# delayed evaluation of expressions. For example, this
# expression:
#
#   (< a | b >^2 + < c | d >^2 - 3.13+im) / 2
#
# is representable as a InnerExpr.

immutable InnerExpr <: Number
    ex::Expr
end

InnerExpr(iex::InnerExpr) = InnerExpr(iex.ex)
InnerExpr{N<:Number}(n::N) = convert(InnerExpr, n)

Base.convert(::Type{InnerExpr}, iex::InnerExpr) = iex
Base.convert{N<:Number}(::Type{InnerExpr}, n::N) = InnerExpr(Expr(:call, +, n))

Base.:(==)(a::InnerExpr, b::InnerExpr) = a.ex == b.ex
same_num(a::Number, b::Number) = a == b
same_num(a::InnerExpr, b::InnerExpr) = a == b
same_num(iex::InnerExpr, i::InnerProduct) = iex == InnerExpr(i)
same_num(i::InnerProduct, iex::InnerExpr) = iex == i
same_num(iex::InnerExpr, n::Number) = iex == InnerExpr(n)
same_num(n::Number, iex::InnerExpr) = iex == n

Base.hash(iex::InnerExpr) = hash(iex.ex)
Base.hash(iex::InnerExpr, h::UInt64) = hash(hash(iex), h)

Base.one(::InnerExpr) = InnerExpr(1)
Base.zero(::InnerExpr) = InnerExpr(0)

Base.promote_rule{N<:Number}(::Type{InnerExpr}, ::Type{N}) = InnerExpr

Base.length(iex::InnerExpr) = length(iex.ex.args)
Base.getindex(iex::InnerExpr, i::Integer) = iex.ex.args[i]

##############
# inner_eval #
##############
inner_eval(f, iex::InnerExpr) = eval(inner_reduce!(f, copy(iex.ex)))
inner_eval(f, c) = c

inner_reduce!(f, iex::InnerExpr) = inner_reduce!(f, copy(iex.ex))
inner_reduce!(f::Function, i::InnerProduct) = f(blabel(i), klabel(i))
inner_reduce!(p::AbstractInner, i::InnerProduct) = inner_rule(p, blabel(i), klabel(i))
inner_reduce!(f, n) = n

function inner_reduce!(f, ex::Expr)
    for i=1:length(ex.args)
        ex.args[i] = inner_reduce!(f, ex.args[i])
    end
    return ex
end

######################
# Printing Functions #
######################
Base.show(io::IO, iex::InnerExpr) = print(io, repr(iex.ex)[2:end])

##################
# Exponentiation #
##################
function iexpr_exp(a, b)
    if same_num(b, 1)
        return InnerExpr(a)
    elseif same_num(a, 0)
        return InnerExpr(0)
    elseif same_num(b, 0)
        return InnerExpr(1)
    else
        return InnerExpr(:($(s)^$(n)))
    end
end

Base.:^(a::InnerExpr, b::Integer) = iexpr_exp(a, b)
Base.:^(a::InnerExpr, b::Rational) = iexpr_exp(a, b)
Base.:^(a::InnerExpr, b::InnerExpr) = iexpr_exp(a, b)
Base.:^(a::InnerExpr, b::Number) = iexpr_exp(a, b)
#Base.:^(a::MathConst{:e}, b::InnerExpr) = iexpr_exp(a, b)
Base.:^(a::Number, b::InnerExpr) = iexpr_exp(a, b)

Base.exp(iex::InnerExpr) = InnerExpr(:(exp($(iex))))
Base.exp2(iex::InnerExpr) = InnerExpr(:(exp2($(iex))))

Base.sqrt(iex::InnerExpr) = InnerExpr(:(sqrt($(iex))))

Base.log(iex::InnerExpr) = length(iex)==2 && iex[1]==:exp ? iex[2] : InnerExpr(:(log($(iex))))
#Base.log(a::MathConst{:e}, b::InnerExpr) = InnerExpr(:(log($(a),$(b))))
Base.log(a::InnerExpr, b::InnerExpr) = InnerExpr(:(log($(a),$(b))))
Base.log(a::InnerExpr, b::Number) = InnerExpr(:(log($(a),$(b))))
Base.log(a::Number, b::InnerExpr) = InnerExpr(:(log($(a),$(b))))
Base.log2(iex::InnerExpr) = length(iex)==2 && iex[1]==:exp2 ? iex[2] : InnerExpr(:(log2($(iex))))

##################
# Multiplication #
##################
function iexpr_mul(a,b)
    if same_num(b, 1)
        return InnerExpr(a)
    elseif same_num(a, 1)
        return InnerExpr(b)
    elseif same_num(b, 0) || same_num(a, 0)
        return InnerExpr(0)
    else
        return InnerExpr(:($(a)*$(b)))
    end
end

Base.:*(a::InnerExpr, b::InnerExpr) =iexpr_mul(a,b)
Base.:*(a::Bool, b::InnerExpr) = iexpr_mul(a,b)
Base.:*(a::InnerExpr, b::Bool) = iexpr_mul(a,b)
Base.:*(a::InnerExpr, b::Number) = iexpr_mul(a,b)
Base.:*(a::Number, b::InnerExpr) = iexpr_mul(a,b)

############
# Division #
############
function iexpr_div(a, b)
    if same_num(a, b)
        return InnerExpr(1)
    elseif same_num(a, 0)
        return InnerExpr(0)
    elseif same_num(b, 0)
        return InnerExpr(Inf)
    elseif same_num(b, 1)
        return InnerExpr(a)
    else
        return InnerExpr(:($(a)/$(b)))
    end
end

Base.:/(a::InnerExpr, b::InnerExpr) = iexpr_div(a, b)
Base.:/(a::InnerExpr, b::Complex) = iexpr_div(a, b)
Base.:/(a::InnerExpr, b::Number) = iexpr_div(a, b)
Base.:/(a::Number, b::InnerExpr) = iexpr_div(a, b)

############
# Addition #
############
function iexpr_add(a, b)
    if same_num(a, 0)
        return InnerExpr(b)
    elseif same_num(b, 0)
        return InnerExpr(a)
    else
        return InnerExpr(:($(a)+$(b)))
    end
end

Base.:+(a::InnerExpr, b::InnerExpr) = iexpr_add(a, b)
Base.:+(a::InnerExpr, b::Number) = iexpr_add(a, b)
Base.:+(a::Number, b::InnerExpr) = iexpr_add(a, b)

###############
# Subtraction #
###############
Base.:-(iex::InnerExpr) = length(iex)==2 && iex[1]==:- ? InnerExpr(iex[2]) :  InnerExpr(:(-$(iex)))

function iexpr_subtract(a, b)
    if same_num(a, b)
        return InnerExpr(0)
    elseif same_num(a, 0)
        return InnerExpr(-b)
    elseif same_num(b, 0)
        return InnerExpr(a)
    else
        return InnerExpr(:($(a)-$(b)))
    end
end

Base.:-(a::InnerExpr, b::InnerExpr) = iexpr_subtract(a, b)
Base.:-(a::InnerExpr, b::Number) = iexpr_subtract(a, b)
Base.:-(a::Number, b::InnerExpr) = iexpr_subtract(a, b)

##################
# Absolute Value #
##################
Base.abs(iex::InnerExpr) = length(iex)==2 && iex[1]==:abs ? iex :  InnerExpr(:(abs($(iex))))
Base.abs2(iex::InnerExpr) = InnerExpr(:(abs2($iex)))

#####################
# Complex Conjugate #
#####################
Base.conj(iex::InnerExpr) = length(iex)==2 && iex[1]==:conj ? InnerExpr(iex[2]) :  InnerExpr(:(conj($(iex))))
Base.ctranspose(iex::InnerExpr) = conj(iex)

##########################
# Elementwise Operations #
##########################
for op=(:*,:-,:+,:/,:^)
    elop = Symbol(string(:.) * string(op))
    @eval begin
        ($elop)(a::InnerExpr, b::InnerExpr) = ($op)(a,b)
        ($elop)(a::InnerExpr, b::Number) = ($op)(a,b)
        ($elop)(a::Number, b::InnerExpr) = ($op)(a,b)
    end
end

######################
# Printing Functions #
######################
function Base.repr(iex::InnerExpr)
    if iex[1]==+
        if length(iex)==2
            return repr(iex[2])
        else
            return repr(Expr(:call, iex.ex[2:end]...))[3:end-1]
        end
    else
        return repr(iex.ex)[2:end]
    end
end

Base.show(io::IO, iex::InnerExpr) = print(io, repr(iex))

export InnerExpr,
    inner_eval,
    inner_rule,
    inner_rettype
