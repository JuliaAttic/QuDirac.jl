#############
# Projector #
#############
    type Projector{S} <: AbstractOperator{S}
        scalar::Number
        ket::Ket{S}
        bra::Bra{S}
    end
        
    Base.copy{S}(op::Projector{S}) = Projector{S}(copy(op.scalar), copy(op.ket), copy(op.bra))

    Base.convert{S}(::Type{DiracOp{S}}, op::Projector{S}) = DiracOp(op.ket, op.bra)
    Base.promote_rule{S}(::Type{DiracOp{S}}, ::Type{Projector{S}}) = DiracOp{S}

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){S}(a::Projector{S}, b::Projector{S}) = a.scalar == b.scalar && a.ket == b.ket && a.bra == b.bra

    Base.hash{S}(op::Projector{S}) = hash(S, hash(op.scalar, hash(op.ket, hash(op.bra))))

    Base.length(op::Projector) = length(op.ket)*length(op.bra)

    Base.getindex(op::Projector, k::Tuple, b::Tuple) = op.scalar * op.ket[k] * op.bra[b]
    Base.getindex(op::Projector, label::(Tuple,Tuple)) = op[label[1], label[2]]
    Base.getindex{S}(op::Projector{S}, ::Colon, ::Colon) = convert(DiracOp{S}, op)
    Base.getindex(op::Projector, k, ::Colon) = (op.scalar * op.ket[k]) * op.bra
    Base.getindex(op::Projector, ::Colon, b) = (op.scalar * op.bra[b]) * op.ket
    Base.getindex(op::Projector, k, b) = op[tuple(k), tuple(b)]

    Base.haskey(op::Projector, label::(Tuple,Tuple)) = hasket(op, label[1]) && hasbra(op, label[2])
    hasket(op::Projector, label::Tuple) = haskey(op.ket, label)
    hasbra(op::Projector, label::Tuple) = haskey(op.bra, label)

    Base.get(op::Projector, label::(Tuple,Tuple), default) = haskey(op, label) ? op[label] : default

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter{S}(f::Function, op::Projector{S}) = filter(f, convert(DiracOp{S}, op))    
    Base.map{S}(f::Function, op::Projector{S}) = map(f, convert(DiracOp{S}, op))

    mapcoeffs{S}(f::Function, op::Projector{S}) = mapcoeffs(f, convert(DiracOp{S}, op))
    maplabels{S}(f::Function, op::Projector{S}) = maplabels(f, convert(DiracOp{S}, op))

##########################
# Mathematical Functions #
##########################
    Base.ctranspose{S}(op::Projector{S}) = Projector{S}(op.scalar', op.bra', op.ket')

    Base.scale!(c::Number, op::Projector) = (op.scalar = c*op.scalar; return op)
    Base.scale!(op::Projector, c::Number) = (op.scalar = op.scalar*c; return op)

    Base.scale(c::Number, op::Projector) = scale!(c,copy(op))
    Base.scale(op::Projector, c::Number) = scale!(copy(op),c)

    Base.(:-)(op::Projector) = (op.scalar = -op.scalar)
    Base.(:+){S}(a::Projector{S}, b::Projector{S}) = convert(DiracOp{S}, a) + convert(DiracOp{S}, b)

    function Base.norm(op::Projector)
        result = 0
        for k in keys(coeffs(op.ket))
            for b in keys(coeffs(op.bra))
                result += op[k,b]^2
            end
        end
        return sqrt(result)
    end

    function Base.trace(op::Projector)
        result = 0
        for k in keys(coeffs(op.ket))
            for b in keys(coeffs(op.bra))
                if b==k
                    result += op[k,b]
                end
            end
        end
        return result
    end

    inner(bra::Bra, op::Projector) = op.scalar * inner(bra, op.ket) * op.bra
    inner(op::Projector, ket::Ket) = op.scalar * op.ket * inner(op.bra, ket)
    inner(a::Projector, b::Projector) = Projector(a.scalar * b.scalar * inner(a.bra,b.ket), a.ket, b.bra)
    inner(a::Projector, b::GenericOperator) = a.scalar * a.ket * inner(a.bra, b)
    inner(a::GenericOperator, b::Projector) = inner(a, b.ket) * b.bra * b.scalar

    QuBase.tensor(ket::Ket, bra::Bra) = Projector(1, ket, bra)
    QuBase.tensor(bra::Bra, ket::Ket) = tensor(ket, bra)
    QuBase.tensor(a::Projector, b::Projector) = Projector(a.scalar * b.scalar, tensor(a.ket,b.ket), tensor(a.bra, b.bra))
    QuBase.tensor(op::Projector, bra::Bra) = Projector(op.scalar, op.ket, tensor(op.bra, bra))
    QuBase.tensor(bra::Bra, op::Projector) = Projector(op.scalar, op.ket, tensor(bra, op.bra))
    QuBase.tensor(op::Projector, ket::Ket) = Projector(op.scalar, tensor(op.ket, ket), op.bra)
    QuBase.tensor(ket::Ket, op::Projector) = Projector(op.scalar, tensor(ket, op.ket), op.bra)

    xsubspace(op::Projector,x) = xsubspace(convert(DiracOp{S}, op), x)
    filternz(op::Projector) = filternz(convert(DiracOp{S}, op))

######################
# Printing Functions #
######################
    function Base.show(io::IO, op::Projector)
        print(io, "$(summary(op)) with $(length(op)) operator(s):")
        pad = "  "
        maxlen = 30
        i = 1
        for k in keys(coeffs(op.ket))
            for b in keys(coeffs(op.bra))
                if i <= maxlen
                    println(io)
                    print(io, "$pad$(op[k,b]) $(statestr(k,Ket))$(statestr(b,Bra))")
                    i = i + 1
                else  
                    println(io)
                    print(io, "$pad$vdots")
                    break
                end
            end
        end
    end

export hasket,
    hasbra,
    mapcoeffs,
    maplabels
