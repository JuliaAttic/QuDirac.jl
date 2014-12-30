import Base: conj,
    repr,
    show,
    getindex,
    length,
    convert,
    one,
    zero,
    promote_rule,
    exp,
    log,
    abs,
    ctranspose,
    ^, .^,
    -, .-,
    /, ./,
    *, .*

################
# InnerProduct #
################
    # An InnerProduct is a type that
    # represents an abstract scalar formulated in
    # Dirac notation - a bra-ket product.
    #
    # The inner() function is defined to evaluate inner products
    # of AbstractStates. It defaults to constructing an InnerProduct
    # type from the factor states. To define an inner product of states 
    # with structure S<:AbstractStructure, one can overload the inner()
    # function as follows:
    #
    # inner(bra::DiracBra{S}, ket::DiracKet{S}) = # custom inner product for structure type S

    immutable InnerProduct{S<:AbstractStructure} <: DiracScalar
        bra::DiracBra{S}
        ket::DiracKet{S}
    end

    InnerProduct{A,B}(b::DiracBra{A}, k::DiracKet{B}) = InnerProduct{typejoin(A,B)}(b, k)

    @defstructure InnerProduct
    ######################
    # Accessor Functions #
    ######################
    getbra(i::InnerProduct) = i.bra
    getket(i::InnerProduct) = i.ket

    ######################
    # Printing Functions #
    ######################
    repr(i::InnerProduct) = repr(getbra(i))*repr(getket(i))[2:end]
    show(io::IO, i::InnerProduct) = print(io, repr(i))

    ###########################
    # Mathematical Operations #
    ###########################
    inner(b::DiracBra, k::DiracKet) = InnerProduct(b,k)
    inner(b::AbstractBra, k::AbstractKet) = (coeff(b)*coeff(k))*inner(state(b), state(k))

    conj{S}(i::InnerProduct{S}) = InnerProduct{S}(getket(i)', getbra(i)')

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

    convert(::Type{ScalarExpr}, s::ScalarExpr) = s
    convert{N<:Number}(::Type{ScalarExpr}, n::N) = ScalarExpr(:(1*$(n)))

    one(::ScalarExpr) = ScalarExpr(1)
    zero(::ScalarExpr) = ScalarExpr(0)

    promote_rule{N<:Number}(::Type{ScalarExpr}, ::Type{N}) = ScalarExpr

    length(s::ScalarExpr) = length(s.ex.args)
    getindex(s::ScalarExpr, i) = s.ex.args[i]

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
    show(io::IO, s::ScalarExpr) = print(io, repr(s.ex)[2:end])

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

    ^(s::DiracScalar, n::Integer) = exponentiate(s,n)
    ^(s::DiracScalar, n::Rational) = exponentiate(s,n)
    ^(s::DiracScalar, n::Number) = exponentiate(s,n)

    exp(s::DiracScalar) = ScalarExpr(:(exp($(s))))

    # The reason we don't actually implement the below comment
    # out method for exp() is that we don't know for sure that 
    # the log is actually natural (base e), and it's probably 
    # not worth it in most cases to check. 
    # exp(s::DiracScalar) = length(s)==2 && s[1]==:log ? s[2] : ScalarExpr(:(exp($(s))))

    log(s::DiracScalar) = length(s)==2 && s[1]==:exp ? s[2] : ScalarExpr(:(log($(s))))
    log(a::MathConst{:e}, b::DiracScalar) = ScalarExpr(:(log($(a),$(b))))

    log(a::DiracScalar, b::DiracScalar) = ScalarExpr(:(log($(a),$(b))))
    log(a::DiracScalar, b::Number) = ScalarExpr(:(log($(a),$(b))))
    log(a::Number, b::DiracScalar) = ScalarExpr(:(log($(a),$(b))))

    ##################
    # Multiplication #
    ##################
    *(a::DiracScalar, b::DiracScalar) = ScalarExpr(:($(a)*$(b)))
    *(a::Bool, b::DiracScalar) = a ? *(1,b) : *(0,b)
    *(a::DiracScalar, b::Bool) = b ? *(a,1) : *(a,0)

    function *(a::DiracScalar, b::Number)
        if b==1
            return ScalarExpr(a)
        elseif b==0
            return ScalarExpr(0)
        else
            return ScalarExpr(:($(a)*$(b)))
        end
    end

    function *(a::Number, b::DiracScalar)
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
    /(a::DiracScalar, b::DiracScalar) = a==b ? ScalarExpr(1) : ScalarExpr(:($(a)/$(b)))

    # the below is only implemented to prevent
    # ambiguity warnings
    function /(a::DiracScalar, b::Complex)
        if b==0
            return ScalarExpr(Inf)
        elseif b==1
            return ScalarExpr(a)
        else
            return ScalarExpr(:($(a)/$(b)))
        end
    end

    function /(a::DiracScalar, b::Number)
        if b==0
            return ScalarExpr(Inf)
        elseif b==1
            return ScalarExpr(a)
        else
            return ScalarExpr(:($(a)/$(b)))
        end
    end

    function /(a::Number, b::DiracScalar)
        if a==0
            return ScalarExpr(0)
        else
            return ScalarExpr(:($(a)/$(b)))
        end
    end

    ##############
    ## Addition ##
    ##############
    +(a::DiracScalar, b::DiracScalar) = ScalarExpr(:($(a)+$(b)))

    function +(a::DiracScalar, b::Number)
        if b==0
            return ScalarExpr(a)
        else
            return ScalarExpr(:($(a)+$(b)))
        end
    end

    function +(a::Number, b::DiracScalar)
        if a==0
            return ScalarExpr(b)
        else
            return ScalarExpr(:($(a)+$(b)))
        end
    end

    #################
    ## Subtraction ##
    #################
    -(s::ScalarExpr) = length(s)==2 && s[1]==:- ? ScalarExpr(s[2]) :  ScalarExpr(:(-$(s)))
    -(s::DiracScalar) = ScalarExpr(:(-$(s)))

    -(a::DiracScalar, b::DiracScalar) = a==b ? ScalarExpr(0) : ScalarExpr(:($(a)-$(b)))

    function -(a::DiracScalar, b::Number)
        if b==0
            return ScalarExpr(a)
        else
            return ScalarExpr(:($(a)-$(b)))
        end
    end

    function -(a::Number, b::DiracScalar)
        if a==0
            return ScalarExpr(-b)
        else
            return ScalarExpr(:($(a)+$(b)))
        end
    end

    ####################
    ## Absolute Value ##
    ####################
    abs(s::ScalarExpr) = length(s)==2 && s[1]==:abs ? s :  ScalarExpr(:(abs($(s))))
    abs(s::DiracScalar) = ScalarExpr(:(abs($i)))

    #######################
    ## Complex Conjugate ##
    #######################
    conj(s::ScalarExpr) = length(s)==2 && s[1]==:conj ? ScalarExpr(s[2]) :  ScalarExpr(:(conj($(s))))
    conj(s::DiracScalar) = ScalarExpr(:(conj($(s))))
    ctranspose(s::DiracScalar) = conj(s)

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

export InnerProduct,
    structure,
    getket,
    getbra,
    inner,
    ScalarExpr,
    queval