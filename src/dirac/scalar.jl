##############
# InnerLabel #
##############
immutable InnerProduct{N} <: DiracScalar
    b::StateLabel{N}
    k::StateLabel{N}
end

blabel(i::InnerProduct) = i.b
klabel(i::InnerProduct) = i.k

Base.repr(i::InnerProduct) = brstr(blabel(i))*ktstr(klabel(i))[2:end]
Base.show(io::IO, i::InnerProduct) = print(io, repr(i))

Base.(:(==))(a::InnerProduct, b::InnerProduct) = blabel(a) == blabel(b) && klabel(a) == klabel(b)                                             
Base.(:(==))(::InnerProduct, ::Number) = false
Base.(:(==))(::Number, ::InnerProduct) = false

Base.hash(i::InnerProduct) = hash(blabel(i), hash(klabel(i)))
Base.hash(i::InnerProduct, h::Uint64) = hash(hash(i), h)

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

##############
# InnerExpr #
##############
# A InnerExpr is a type that wraps arthimetic expressions
# performed with DiracScalars. The allows storage and 
# delayed evaluation of expressions. For example, this 
# expression:
#   
#   (< a | b >^2 + < c | d >^2 - 3.13+im) / 2
#
# is representable as a InnerExpr.
#
# One can then use the inner_eval(::Function, ::InnerExpr) 
# function to map an evaluation function onto all InnerProducts
# contained in the InnerExpr, and evaluate the expression
# arthimetically.

immutable InnerExpr <: DiracScalar
    ex::Expr
end

InnerExpr(s::InnerExpr) = InnerExpr(s.ex)
InnerExpr{N<:Number}(n::N) = convert(InnerExpr, n)

Base.convert(::Type{InnerExpr}, s::InnerExpr) = s
Base.convert{N<:Number}(::Type{InnerExpr}, n::N) = InnerExpr(Expr(:call, +, n))

Base.(:(==))(a::InnerExpr, b::InnerExpr) = a.ex == b.ex
Base.(:(==))(s::InnerExpr, i::InnerProduct) = s == InnerExpr(i)
Base.(:(==))(i::InnerProduct, s::InnerExpr) = s == i
Base.(:(==))(::InnerExpr, ::Number) = false
Base.(:(==))(::Number, ::InnerExpr) = false

Base.hash(s::InnerExpr) = hash(s.ex)
Base.hash(s::InnerExpr, h::Uint64) = hash(hash(s), h)

Base.one(::InnerExpr) = InnerExpr(1)
Base.zero(::InnerExpr) = InnerExpr(0)

Base.promote_rule{N<:Number}(::Type{InnerExpr}, ::Type{N}) = Number

Base.length(s::InnerExpr) = length(s.ex.args)
Base.getindex(s::InnerExpr, i) = s.ex.args[i]

##############
# inner_eval #
##############
inner_eval(f, s::InnerExpr) = eval(inner_reduce!(f, copy(s.ex)))

inner_reduce!(f, s::InnerExpr) = inner_reduce!(f, copy(s.ex))
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
Base.show(io::IO, s::InnerExpr) = print(io, repr(s.ex)[2:end])

##################
# Exponentiation #
##################
exponentiate(a::DiracScalar, b::DiracScalar) = InnerExpr(:($(a)^$(b)))

function exponentiate(s::DiracScalar, n::Number)
    if n==1
        return InnerExpr(s)
    elseif n==0
        return InnerExpr(1)
    else
        return InnerExpr(:($(s)^$(n)))
    end
end

Base.(:^)(s::DiracScalar, n::Integer) = exponentiate(s,n)
Base.(:^)(s::DiracScalar, n::Rational) = exponentiate(s,n)
Base.(:^)(s::DiracScalar, n::Number) = exponentiate(s,n)

Base.exp(s::DiracScalar) = InnerExpr(:(exp($(s))))
Base.exp2(s::DiracScalar) = InnerExpr(:(exp2($(s))))

Base.sqrt(s::DiracScalar) = InnerExpr(:(sqrt($(s))))

# The reason we don't actually implement the below comment
# out method for exp() is that we don't know for sure that 
# the log is actually natural (base e), and it's probably 
# not worth it in most cases to check. 
# exp(s::DiracScalar) = length(s)==2 && s[1]==:log ? s[2] : InnerExpr(:(exp($(s))))

Base.log(s::DiracScalar) = length(s)==2 && s[1]==:exp ? s[2] : InnerExpr(:(log($(s))))
Base.log2(s::DiracScalar) = length(s)==2 && s[1]==:exp2 ? s[2] : InnerExpr(:(log2($(s))))

Base.log(a::MathConst{:e}, b::DiracScalar) = InnerExpr(:(log($(a),$(b))))

Base.log(a::DiracScalar, b::DiracScalar) = InnerExpr(:(log($(a),$(b))))
Base.log(a::DiracScalar, b::Number) = InnerExpr(:(log($(a),$(b))))
Base.log(a::Number, b::DiracScalar) = InnerExpr(:(log($(a),$(b))))

##################
# Multiplication #
##################
Base.(:*)(a::DiracScalar, b::DiracScalar) = InnerExpr(:($(a)*$(b)))
Base.(:*)(a::Bool, b::DiracScalar) = a ? *(1,b) : *(0,b)
Base.(:*)(a::DiracScalar, b::Bool) = b ? *(a,1) : *(a,0)

function Base.(:*)(a::DiracScalar, b::Number)
    if b==1
        return InnerExpr(a)
    elseif b==0
        return InnerExpr(0)
    else
        return InnerExpr(:($(a)*$(b)))
    end
end

function Base.(:*)(a::Number, b::DiracScalar)
    if a==1
        return InnerExpr(b)
    elseif a==0
        return InnerExpr(0)
    else
        return InnerExpr(:($(a)*$(b)))
    end
end

##############
## Division ##
##############
Base.(:/)(a::DiracScalar, b::DiracScalar) = a==b ? InnerExpr(1) : InnerExpr(:($(a)/$(b)))

# the below is only implemented to prevent
# ambiguity warnings
function Base.(:/)(a::DiracScalar, b::Complex)
    if b==0
        return InnerExpr(Inf)
    elseif b==1
        return InnerExpr(a)
    else
        return InnerExpr(:($(a)/$(b)))
    end
end

function Base.(:/)(a::DiracScalar, b::Number)
    if b==0
        return InnerExpr(Inf)
    elseif b==1
        return InnerExpr(a)
    else
        return InnerExpr(:($(a)/$(b)))
    end
end

function Base.(:/)(a::Number, b::DiracScalar)
    if a==0
        return InnerExpr(0)
    else
        return InnerExpr(:($(a)/$(b)))
    end
end

##############
## Addition ##
##############
Base.(:+)(a::DiracScalar, b::DiracScalar) = InnerExpr(:($(a)+$(b)))

function Base.(:+)(a::DiracScalar, b::Number)
    if b==0
        return InnerExpr(a)
    else
        return InnerExpr(:($(a)+$(b)))
    end
end

function Base.(:+)(a::Number, b::DiracScalar)
    if a==0
        return InnerExpr(b)
    else
        return InnerExpr(:($(a)+$(b)))
    end
end

#################
## Subtraction ##
#################
Base.(:-)(s::InnerExpr) = length(s)==2 && s[1]==:- ? InnerExpr(s[2]) :  InnerExpr(:(-$(s)))
Base.(:-)(s::DiracScalar) = InnerExpr(:(-$(s)))

Base.(:-)(a::DiracScalar, b::DiracScalar) = a==b ? InnerExpr(0) : InnerExpr(:($(a)-$(b)))

function Base.(:-)(a::DiracScalar, b::Number)
    if b==0
        return InnerExpr(a)
    else
        return InnerExpr(:($(a)-$(b)))
    end
end

function Base.(:-)(a::Number, b::DiracScalar)
    if a==0
        return InnerExpr(-b)
    else
        return InnerExpr(:($(a)+$(b)))
    end
end

####################
## Absolute Value ##
####################
Base.abs(s::InnerExpr) = length(s)==2 && s[1]==:abs ? s :  InnerExpr(:(abs($(s))))
Base.abs(s::DiracScalar) = InnerExpr(:(abs($s)))

Base.abs2(s::DiracScalar) = InnerExpr(:(abs2($s)))

#######################
## Complex Conjugate ##
#######################
Base.conj(s::InnerExpr) = length(s)==2 && s[1]==:conj ? InnerExpr(s[2]) :  InnerExpr(:(conj($(s))))
Base.conj(s::DiracScalar) = InnerExpr(:(conj($(s))))
Base.ctranspose(s::DiracScalar) = conj(s)

############################
## Elementwise Operations ##
############################
for op=(:*,:-,:+,:/,:^)
    elop = symbol(string(:.) * string(op))
    @eval begin
        ($elop)(a::DiracScalar, b::DiracScalar) = ($op)(a,b)
        ($elop)(a::DiracScalar, b::Number) = ($op)(a,b)
        ($elop)(a::Number, b::DiracScalar) = ($op)(a,b)
    end
end

########################
## Printing Functions ##
########################
function Base.repr(s::InnerExpr)
    if s[1]==+
        if length(s)==2
            return repr(s[2])
        else
            return repr(Expr(:call, s.ex[2:end]...))[3:end-1]
        end
    else
        return repr(s.ex)[2:end]
    end
end

Base.show(io::IO, s::InnerExpr) = print(io, repr(s))

export inner_eval,
    inner_rule
