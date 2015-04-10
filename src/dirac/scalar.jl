##############
# InnerLabel #
##############
immutable InnerProduct{P<:AbstractInner,N} <: DiracScalar
    ptype::P
    b::StateLabel{N}
    k::StateLabel{N}
end

InnerProduct(ptype, b, k) = InnerProduct(ptype, StateLabel(k), StateLabel(b))

blabel(i::InnerProduct) = i.b
klabel(i::InnerProduct) = i.k

Base.repr(i::InnerProduct) = brstr(blabel(i))*ktstr(klabel(i))[2:end]
Base.show(io::IO, i::InnerProduct) = print(io, repr(i))

Base.(:(==))(a::InnerProduct, b::InnerProduct) = a.ptype == b.ptype && 
                                             blabel(a) == blabel(b) &&
                                             klabel(a) == klabel(b) 
Base.(:(==))(::InnerProduct, ::Number) = false
Base.(:(==))(::Number, ::InnerProduct) = false

Base.conj(i::InnerProduct) = InnerProduct(i.ptype, klabel(i), blabel(i))

##############
# inner_rule #
##############
inner_rule(p::AbstractInner, b, k) = ScalarExpr(InnerProduct(p, b, k))
inner_rule(::KroneckerDelta, b, k) = b == k ? 1 : 0

inner_type{P}(::P, N) = first(Base.return_types(inner_rule, (P, StateLabel{N}, StateLabel{N})))
inner_type(::UndefinedInner, N) = ScalarExpr
inner_type(::KroneckerDelta, N) = Int64

##############
# ScalarExpr #
##############
# A ScalarExpr is a type that wraps arthimetic expressions
# performed with DiracScalars. The allows storage and 
# delayed evaluation of expressions. For example, this 
# expression:
#   
#   (< a | b >^2 + < c | d >^2 - 3.13+im) / 2
#
# is representable as a ScalarExpr.
#
# One can then use the inner_eval(::Function, ::ScalarExpr) 
# function to map an evaluation function onto all InnerProducts
# contained in the ScalarExpr, and evaluate the expression
# arthimetically.

immutable ScalarExpr <: DiracScalar
    ex::Expr
end

ScalarExpr(s::ScalarExpr) = ScalarExpr(s.ex)
ScalarExpr{N<:Number}(n::N) = convert(ScalarExpr, n)

Base.(:(==))(a::ScalarExpr, b::ScalarExpr) = a.ex == b.ex
Base.(:(==))(a::InnerProduct, b::ScalarExpr) = ScalarExpr(a) == b
Base.(:(==))(a::ScalarExpr, b::InnerProduct) = a == ScalarExpr(b)
Base.(:(==))(::ScalarExpr, ::Number) = false
Base.(:(==))(::Number, ::ScalarExpr) = false

Base.convert(::Type{ScalarExpr}, s::ScalarExpr) = s
Base.convert{N<:Number}(::Type{ScalarExpr}, n::N) = ScalarExpr(Expr(:call, +, n))

Base.one(::ScalarExpr) = ScalarExpr(1)
Base.zero(::ScalarExpr) = ScalarExpr(0)

Base.promote_rule{N<:Number}(::Type{ScalarExpr}, ::Type{N}) = Number

Base.length(s::ScalarExpr) = length(s.ex.args)
Base.getindex(s::ScalarExpr, i) = s.ex.args[i]

##############
# inner_eval #
##############
inner_eval(i::InnerProduct) = inner_rule(i.ptype, blabel(i), klabel(i))
inner_eval(s::ScalarExpr) = eval(inner_reduce!(copy(s.ex)))
inner_eval(f::Function, s::ScalarExpr) = eval(f_reduce!(f, copy(s.ex)))
inner_eval(n::Number) = n

inner_reduce!(s::ScalarExpr) = inner_reduce!(copy(s.ex))
inner_reduce!(i::InnerProduct) = inner_eval(i)
inner_reduce!(n) = n

function inner_reduce!(ex::Expr)
    for i=1:length(ex.args)
        ex.args[i] = inner_reduce!(ex.args[i])
    end
    return ex
end

f_reduce!(f::Function, s::ScalarExpr) = f_reduce!(f, copy(s.ex))
f_reduce!(f::Function, i::InnerProduct) = f(blabel(i), klabel(i))
f_reduce!(f::Function, n) = n

function f_reduce!(f::Function, ex::Expr)
    for i=1:length(ex.args)
        ex.args[i] = f_reduce!(f, ex.args[i])
    end
    return ex
end

######################
# Printing Functions #
######################
Base.show(io::IO, s::ScalarExpr) = print(io, repr(s.ex)[2:end])

##################
# Exponentiation #
##################
exponentiate(a::DiracScalar, b::DiracScalar) = ScalarExpr(:($(a)^$(b)))

function exponentiate(s::DiracScalar, n::Number)
    if n==1
        return ScalarExpr(s)
    elseif n==0
        return ScalarExpr(1)
    else
        return ScalarExpr(:($(s)^$(n)))
    end
end

Base.(:^)(s::DiracScalar, n::Integer) = exponentiate(s,n)
Base.(:^)(s::DiracScalar, n::Rational) = exponentiate(s,n)
Base.(:^)(s::DiracScalar, n::Number) = exponentiate(s,n)

Base.exp(s::DiracScalar) = ScalarExpr(:(exp($(s))))
Base.exp2(s::DiracScalar) = ScalarExpr(:(exp2($(s))))

Base.sqrt(s::DiracScalar) = ScalarExpr(:(sqrt($(s))))

# The reason we don't actually implement the below comment
# out method for exp() is that we don't know for sure that 
# the log is actually natural (base e), and it's probably 
# not worth it in most cases to check. 
# exp(s::DiracScalar) = length(s)==2 && s[1]==:log ? s[2] : ScalarExpr(:(exp($(s))))

Base.log(s::DiracScalar) = length(s)==2 && s[1]==:exp ? s[2] : ScalarExpr(:(log($(s))))
Base.log2(s::DiracScalar) = length(s)==2 && s[1]==:exp2 ? s[2] : ScalarExpr(:(log2($(s))))

Base.log(a::MathConst{:e}, b::DiracScalar) = ScalarExpr(:(log($(a),$(b))))

Base.log(a::DiracScalar, b::DiracScalar) = ScalarExpr(:(log($(a),$(b))))
Base.log(a::DiracScalar, b::Number) = ScalarExpr(:(log($(a),$(b))))
Base.log(a::Number, b::DiracScalar) = ScalarExpr(:(log($(a),$(b))))

##################
# Multiplication #
##################
Base.(:*)(a::DiracScalar, b::DiracScalar) = ScalarExpr(:($(a)*$(b)))
Base.(:*)(a::Bool, b::DiracScalar) = a ? *(1,b) : *(0,b)
Base.(:*)(a::DiracScalar, b::Bool) = b ? *(a,1) : *(a,0)

function Base.(:*)(a::DiracScalar, b::Number)
    if b==1
        return ScalarExpr(a)
    elseif b==0
        return ScalarExpr(0)
    else
        return ScalarExpr(:($(a)*$(b)))
    end
end

function Base.(:*)(a::Number, b::DiracScalar)
    if a==1
        return ScalarExpr(b)
    elseif a==0
        return ScalarExpr(0)
    else
        return ScalarExpr(:($(a)*$(b)))
    end
end

##############
## Division ##
##############
Base.(:/)(a::DiracScalar, b::DiracScalar) = a==b ? ScalarExpr(1) : ScalarExpr(:($(a)/$(b)))

# the below is only implemented to prevent
# ambiguity warnings
function Base.(:/)(a::DiracScalar, b::Complex)
    if b==0
        return ScalarExpr(Inf)
    elseif b==1
        return ScalarExpr(a)
    else
        return ScalarExpr(:($(a)/$(b)))
    end
end

function Base.(:/)(a::DiracScalar, b::Number)
    if b==0
        return ScalarExpr(Inf)
    elseif b==1
        return ScalarExpr(a)
    else
        return ScalarExpr(:($(a)/$(b)))
    end
end

function Base.(:/)(a::Number, b::DiracScalar)
    if a==0
        return ScalarExpr(0)
    else
        return ScalarExpr(:($(a)/$(b)))
    end
end

##############
## Addition ##
##############
Base.(:+)(a::DiracScalar, b::DiracScalar) = ScalarExpr(:($(a)+$(b)))

function Base.(:+)(a::DiracScalar, b::Number)
    if b==0
        return ScalarExpr(a)
    else
        return ScalarExpr(:($(a)+$(b)))
    end
end

function Base.(:+)(a::Number, b::DiracScalar)
    if a==0
        return ScalarExpr(b)
    else
        return ScalarExpr(:($(a)+$(b)))
    end
end

#################
## Subtraction ##
#################
Base.(:-)(s::ScalarExpr) = length(s)==2 && s[1]==:- ? ScalarExpr(s[2]) :  ScalarExpr(:(-$(s)))
Base.(:-)(s::DiracScalar) = ScalarExpr(:(-$(s)))

Base.(:-)(a::DiracScalar, b::DiracScalar) = a==b ? ScalarExpr(0) : ScalarExpr(:($(a)-$(b)))

function Base.(:-)(a::DiracScalar, b::Number)
    if b==0
        return ScalarExpr(a)
    else
        return ScalarExpr(:($(a)-$(b)))
    end
end

function Base.(:-)(a::Number, b::DiracScalar)
    if a==0
        return ScalarExpr(-b)
    else
        return ScalarExpr(:($(a)+$(b)))
    end
end

####################
## Absolute Value ##
####################
Base.abs(s::ScalarExpr) = length(s)==2 && s[1]==:abs ? s :  ScalarExpr(:(abs($(s))))
Base.abs(s::DiracScalar) = ScalarExpr(:(abs($s)))

Base.abs2(s::DiracScalar) = ScalarExpr(:(abs2($s)))

#######################
## Complex Conjugate ##
#######################
Base.conj(s::ScalarExpr) = length(s)==2 && s[1]==:conj ? ScalarExpr(s[2]) :  ScalarExpr(:(conj($(s))))
Base.conj(s::DiracScalar) = ScalarExpr(:(conj($(s))))
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
function Base.repr(s::ScalarExpr)
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

Base.show(io::IO, s::ScalarExpr) = print(io, repr(s))

export inner_eval
