#############
# Projector #
#############
    type Projector{S} <: AbstractOperator{S}
        scalar::Number
        ket::Ket{S}
        bra::Bra{S}
    end
        
    Base.copy{S}(op::Projector{S}) = Projector{S}(copy(coeffs(op)))

    Base.convert{S}(::Type{DiracOp{S}}, op::Projector{S}) = DiracOp(op.ket, op.bra)
    Base.promote_rule{S}(::Type{DiracOp{S}}, ::Type{Projector{S}}) = DiracOp{S}

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){S}(a::Projector{S}, b::Projector{S}) = a.scalar == b.scalar && a.ket == b.ket && a.bra == b.bra

    Base.hash{S}(op::Projector{S}) = hash(S, hash(op.scalar, hash(op.ket, hash(op.bra)))

    Base.length(op::Projector) = length(op.ket)*length(op.bra)

    Base.getindex(op::Projector, k::Tuple, b::Tuple) = op.scalar * op.ket[k] * op.bra[b]
    Base.getindex(op::Projector, label::(Tuple,Tuple)) = op[label[1], label[2]]
    Base.getindex{S}(op::Projector{S}, ::Colon, ::Colon) = convert(DiracOp{S}, op)
    Base.getindex(op::Projector, k, ::Colon) = (op.scalar * op.ket[k]) * op.bra
    Base.getindex(op::Projector, ::Colon, b) = (op.scalar * op.bra[b]) * op.ket

    Base.haskey(op::Projector, label::(Tuple,Tuple)) = hasket(op, label[1]) && hasbra(op, label[2])
    hasket(op::Projector, label::Tuple) = haskey(op.ket, label)
    hasbra(op::Projector, label::Tuple) = haskey(op.bra, label)

    Base.get(op::Projector, label::(Tuple,Tuple), default) = haskey(op, label) ? op[label] : default

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!{S}(f::Function, op::Projector{S}) = filter!(f, convert(DiracOp{S}, op))
    Base.filter{S}(f::Function, op::Projector{S}) = filter(f, convert(DiracOp{S}, op))
    
    Base.map{S}(f::Function, op::Projector{S}) = map(f, convert(DiracOp{S}, op))

    mapcoeffs{S}(f::Function, op::Projector{S}) = mapcoeffs(f, convert(DiracOp{S}, op))
    maplabels{S}(f::Function, op::Projector{S}) = maplabels(f, convert(DiracOp{S}, op))

##########################
# Mathematical Functions #
##########################
    Base.ctranspose{S}(op::Projector{S}) = Projector{S}(op.scalar', op.bra', op.ket')

    QuBase.tensor(ket::Ket, bra::Bra) = Projector(1, ket, bra)
    QuBase.tensor(bra::Bra, ket::Ket) = tensor(ket, bra)

    Base.scale!(c::Number, op::Projector) = (op.scalar = c*op.scalar)
    Base.scale!(op::Projector, c::Number) = (op.scalar = op.scalar*c)

    Base.scale(c::Number, op::Projector) = scale!(c,copy(op))
    Base.scale(op::Projector, c::Number) = scale!(copy(op),c)

    Base.(:-)(op::Projector) = (op.scalar = -op.scalar)

export hasket,
    hasbra,
    mapcoeffs,
    maplabels
