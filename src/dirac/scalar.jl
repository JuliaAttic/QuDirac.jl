################
# InnerProduct #
################
# An InnerProduct is a type that
# represents an abstract scalar formulated in
# Dirac notation - a br-kt product.
immutable InnerProduct{P<:AbstractInner} <: DiracScalar
    brlabel::Vector{Any}
    ktlabel::Vector{Any}
    InnerProduct(b::Array, k::Array) = new(b,k)
    InnerProduct(b, k) = InnerProduct{P}([b],[k])
end

inner_eval{A,B}(::Type{A}, ::Type{B}, k, b) = inner_rule(typejoin(A, B), k, b)

inner_rule{P<:AbstractInner}(::Type{P}, k, b) = ScalarExpr(InnerProduct{P}(k, b))
inner_rule{O<:Orthonormal}(::Type{O}, k, b) = k == b ? 1 : 0

######################
# Accessor Functions #
######################
brlabel(i::InnerProduct) = i.brlabel
ktlabel(i::InnerProduct) = i.ktlabel

######################
# Printing Functions #
######################
Base.repr(i::InnerProduct) = brstr(brlabel(i))*ktstr(ktlabel(i))[2:end]
Base.show(io::IO, i::InnerProduct) = print(io, repr(i))

###########################
# Mathematical Operations #
###########################
Base.(:(==))(a::InnerProduct, b::InnerProduct) = brlabel(a) == brlabel(b) && ktlabel(a) == ktlabel(b) 
Base.(:(==))(::InnerProduct, ::Number) = false
Base.(:(==))(::Number, ::InnerProduct) = false

Base.conj{P}(i::InnerProduct{P}) = InnerProduct{P}(getktlabel(i), getbrlabel(i))

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
# One can then use the queval(::Function, ::ScalarExpr) 
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

Base.convert(::Type{ScalarExpr}, s::ScalarExpr) = s
Base.convert{N<:Number}(::Type{ScalarExpr}, n::N) = ScalarExpr(Expr(:call, +, n))

Base.one(::ScalarExpr) = ScalarExpr(1)
Base.zero(::ScalarExpr) = ScalarExpr(0)

Base.promote_rule{N<:Number}(::Type{ScalarExpr}, ::Type{N}) = ScalarExpr

Base.length(s::ScalarExpr) = length(s.ex.args)
Base.getindex(s::ScalarExpr, i) = s.ex.args[i]

##########
# queval #
##########
queval(f::Function, s::ScalarExpr) = eval(qureduce!(f, copy(s.ex)))
queval(f::Function, i::InnerProduct) = f(i)
queval(f::Function, n::Number) = n

qureduce!(f::Function, s::ScalarExpr) = qureduce!(f, copy(s.ex))
qureduce!(f::Function, i::InnerProduct) = f(i)
qureduce!(f::Function, n) = n

function qureduce!(f::Function, ex::Expr)
    for i=1:length(ex.args)
        ex.args[i] = qureduce!(f, ex.args[i])
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
Base.sqrt(s::DiracScalar) = ScalarExpr(:(sqrt($(s))))

# The reason we don't actually implement the below comment
# out method for exp() is that we don't know for sure that 
# the log is actually natural (base e), and it's probably 
# not worth it in most cases to check. 
# exp(s::DiracScalar) = length(s)==2 && s[1]==:log ? s[2] : ScalarExpr(:(exp($(s))))

Base.log(s::DiracScalar) = length(s)==2 && s[1]==:exp ? s[2] : ScalarExpr(:(log($(s))))
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
Base.abs(s::DiracScalar) = ScalarExpr(:(abs($i)))

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

export ScalarExpr,
    ktlabel,
    brlabel,
    queval