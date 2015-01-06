import Base: show,
    repr,
    convert,
    ctranspose,
    Ac_mul_B, #TODO: Define other members of the A_mul_B family as necessary 
    +,
    -,
    *,.*

################################
# Abstract types and Functions #
################################
    # An AbstractOperator is a type that
    # represents an operator formulated as a 
    # ket-bra projector, optionally multiplied 
    # by a scalar value. 
    #
    # Concrete subtypes of AbstractOperator
    # implement the following methods (T<:AbstractOperator):
    #
    # coeff(op::T) -> returns the scalar coefficient of the operator
    # operator(op::T) -> returns only the bra-ket projector, without a scalar coefficient, as a DiracOperator
    # getbra(op::T) -> returns the bra component of the operator as a DiracBra
    # getket(op::T) -> returns the ket component of the operator as a DiracKet
    # label(op::T) -> returns the operator's labels as a tuple of type (StateLabel, StateLabel)
    # repr(op::T) -> returns a string representation of the operator
    # coefftype(op::T), 
    # coefftype(::Type{T}) -> returns the coefficient type for a given operator/type

    show(io::IO, op::AbstractOperator) = print(io, repr(op))
    samelabels(a::AbstractOperator, b::AbstractOperator) = samelabels(getket(a), getket(b)) && samelabels(getbra(a), getbra(b))

#################
# DiracOperator #
#################
    # A DiracOperator is the type representation of
    # an unscaled abstract operator, formulated in 
    # Dirac notation as a ket-bra projector.
    immutable DiracOperator{S} <: AbstractOperator{S}
        ket::DiracKet{S}
        bra::DiracBra{S}
    end

    ######################
    # Accessor Functions #
    ######################
    structure{S}(::Type{DiracOperator{S}}) = S

    getket(op::DiracOperator) = op.ket
    getbra(op::DiracOperator) = op.bra

    coefftype(op::DiracOperator) = Complex128
    coefftype{O<:DiracOperator}(::Type{O}) = Complex128
    
    coeff(op::DiracOperator) = one(coefftype(op))
    operator(op::DiracOperator) = op
    label(op::DiracOperator) = (label(getket(op)), label(getbra(op)))

    ######################
    # Printing Functions #
    ######################
    repr(op::DiracOperator) = repr(getket(op))*repr(getbra(op))

    ###########################
    # Mathematical Operations #
    ###########################
    outer(k::DiracKet, b::DiracBra) = DiracOperator(k,b)
    ctranspose(op::DiracOperator) = outer(getbra(op)', getket(op)')

##################
# ScaledOperator #
##################
    # A ScaledOperator is the type representation of
    # a DiracOperator multiplied by a scalar.
    immutable ScaledOperator{S<:AbstractStructure, T} <: AbstractOperator{S}
        coeff::T
        operator::DiracOperator{S}
    end

    ######################
    # Accessor Functions #
    ######################
    coefftype{S,T}(::ScaledOperator{S,T}) = T
    coefftype(::Type{ScaledOperator}) = Any
    coefftype{S,T}(::Type{ScaledOperator{S,T}}) = T

    structure(::Type{ScaledOperator}) = AbstractStructure
    structure{S,T}(::Type{ScaledOperator{S,T}}) = S

    coeff(op::ScaledOperator) = op.coeff
    operator(op::ScaledOperator) = op.operator
    label(op::ScaledOperator) = label(operator(op))

    getket(op::ScaledOperator) = getket(operator(op))
    getbra(op::ScaledOperator) = getbra(operator(op))

    ######################
    # Printing Functions #
    ######################
    repr(op::ScaledOperator) = "$(coeff(op)) $(repr(operator(op)))"

    ###########################
    # Mathematical Operations #
    ###########################
    ctranspose(op::ScaledOperator) = ScaledOperator(coeff(op)', operator(op)')
    outer(k::AbstractKet, b::AbstractBra) = ScaledOperator(kron(coeff(k), coeff(b)), outer(state(k), state(b)))

#######################
# Both Operator Types #
#######################
    ###########################
    # Mathematical Operations #
    ###########################
    tensor(a::DiracOperator, b::DiracOperator) = outer(tensor(getket(a), getket(b)), tensor(getbra(a), getbra(b)))
    tensor(a::AbstractOperator, b::AbstractOperator) = ScaledOperator(kron(coeff(a),coeff(b)), tensor(operator(a), operator(b))) 

    tensor(op::DiracOperator, s::DiracBra) = outer(getket(op), tensor(getbra(op), s))
    tensor(s::DiracBra, op::DiracOperator) = outer(getket(op), tensor(s, getbra(op)))
    tensor(op::DiracOperator, s::DiracKet) = outer(tensor(getket(op), s), getbra(op))
    tensor(s::DiracKet, op::DiracOperator) = outer(tensor(s, getket(op)), getbra(op))
    tensor(op::AbstractOperator, s::AbstractState) = ScaledOperator(kron(coeff(op),coeff(s)), tensor(operator(op), s)) 
    tensor(s::AbstractState, op::AbstractOperator) = ScaledOperator(kron(coeff(s),coeff(op)), tensor(s, operator(op)))

    for op=(:*, :.*)
        @eval begin
            ($op)(s::AbstractKet, op::AbstractOperator) = tensor(s,op)
            ($op)(op::AbstractOperator, s::AbstractBra) = tensor(op,s)
            ($op)(a::AbstractOperator, b::AbstractOperator) = ScaledOperator(inner(getbra(a), getket(b)), tensor(getket(a), getbra(b))) 
            ($op)(s::AbstractBra, op::AbstractOperator) = ($op)(inner(s, getket(op)), getbra(op))
            ($op)(op::AbstractOperator, s::AbstractKet) = ($op)(inner(getbra(op), s), getket(op))
            ($op)(c::Number, op::AbstractOperator) = ScaledOperator(($op)(c,coeff(op)), operator(op))
            ($op)(op::AbstractOperator, c::Number) = ScaledOperator(($op)(coeff(op),c), operator(op))
        end
    end

    Ac_mul_B(a::AbstractOperator, b::AbstractOperator) = inner(a', b)
    Ac_mul_B(a::AbstractBra, b::AbstractOperator) = tensor(a', b)
    Ac_mul_B(a::AbstractOperator, b::AbstractBra) = tensor(a', b)
    Ac_mul_B(a::AbstractKet, b::AbstractOperator) = inner(a', getket(b)) * getbra(a)
    Ac_mul_B(a::AbstractOperator, b::AbstractKet) = inner(getket(a)', b) * getbra(a)'

    -(op::AbstractOperator) = ScaledOperator(-(coeff(op)), operator(op)) 
    +(op::AbstractOperator) = op

export DiracOperator,
    ScaledOperator,
    outer,
    operator,
    getket,
    getbra,
    coeff,
    coefftype,
    label,
    samelabels,
    structure,
    tensor