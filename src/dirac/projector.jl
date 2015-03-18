#############
# Projector #
#############
    type Projector{S} <: AbstractOperator{S}
        scalar::Number
        ket::Ket{S}
        bra::Bra{S}
    end
        
    Base.copy{S}(op::Projector{S}) = Projector{S}(copy(op.scalar), copy(op.ket), copy(op.bra))

    Base.convert{S}(::Type{DiracOp{S}}, op::Projector{S}) = scale!(op.scalar, DiracOp(op.ket, op.bra))
    Base.promote_rule{S}(::Type{DiracOp{S}}, ::Type{Projector{S}}) = DiracOp{S}

    to_diracop{S}(op::Projector{S}) = convert(DiracOp{S}, op)

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){S}(a::Projector{S}, b::Projector{S}) = a.scalar == b.scalar && a.ket == b.ket && a.bra == b.bra

    Base.hash{S}(op::Projector{S}) = hash(S, hash(op.scalar, hash(op.ket, hash(op.bra))))

    Base.length(op::Projector) = length(op.ket)*length(op.bra)

    Base.getindex(op::Projector, k::Array, b::Array) = op.scalar * op.ket[k] * op.bra[b]
    Base.getindex(op::Projector, label::Array) = op[label[1], label[2]]
    Base.getindex(op::Projector, k, b) = op[[k],[b]]

    # would be great if the below worked with normal indexing
    # notation (e.g. op[k,:]) but slice notation is apparently
    # special and doesn't dispatch directly to getindex...
    # Base.getindex(op::Projector, k, ::Colon) = (op.scalar * op.ket[k]) * op.bra
    # Base.getindex(op::Projector, ::Colon, b) = (op.scalar * op.bra[b]) * op.ket
    # Base.getindex(op::Projector, ::Colon, ::Colon) = to_diracop(op)

    getbra(op::Projector, k::Array) = (op.scalar * op.ket[k]) * op.bra
    getket(op::Projector, b::Array) = (op.scalar * op.bra[b]) * op.ket

    Base.haskey(op::Projector, label::Array) = hasket(op, label[1]) && hasbra(op, label[2])
    hasket(op::Projector, label::Array) = haskey(op.ket, label)
    hasbra(op::Projector, label::Array) = haskey(op.bra, label)

    Base.get(op::Projector, label::Array, default) = haskey(op, label) ? op[label] : default

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter(f::Function, op::Projector) = filter(f, to_diracop(op))
    Base.map(f::Function, op::Projector) = map(f, to_diracop(op))

    mapcoeffs(f::Function, op::Projector) = mapcoeffs(f, to_diracop(op))
    maplabels(f::Function, op::Projector) = maplabels(f, to_diracop(op))

##########################
# Mathematical Functions #
##########################
    Base.ctranspose{S}(op::Projector{S}) = Projector{S}(op.scalar', op.bra', op.ket')

    Base.scale!(c::Number, op::Projector) = (op.scalar = c*op.scalar; return op)
    Base.scale!(op::Projector, c::Number) = (op.scalar = op.scalar*c; return op)

    Base.scale(c::Number, op::Projector) = scale!(c,copy(op))
    Base.scale(op::Projector, c::Number) = scale!(copy(op),c)

    Base.(:-)(op::Projector) = (op.scalar = -op.scalar)
    Base.(:+)(a::Projector, b::Projector) = to_diracop(a) + to_diracop(b)

    function Base.norm(op::Projector)
        result = 0
        for v in values(dict(op.ket))
            for c in values(dict(op.bra))
                result += (v*c)^2
            end
        end
        return sqrt(result)
    end

    function Base.trace(op::Projector)
        result = 0
        for k in keys(dict(op.ket))
            for b in keys(dict(op.bra))
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

    xsubspace(op::Projector,x) = xsubspace(to_diracop(op), x)
    filternz(op::Projector) = filternz(to_diracop(op))

    function ptrace{S<:Orthonormal}(op::Projector{S}, over::Integer)
        return DiracOp{S}(ptrace_proj!(OpDict(), op, over))
    end

    function ptrace_proj!(result::OpDict, op::Projector, over)
        for k in keys(dict(op.ket)), b in keys(dict(op.bra))
            if k[over] == b[over]
                add_to_dict!(result,
                             veclabel(except(k, over), except(b, over)),
                             op[k,b])
            end
        end
        return result
    end

######################
# Printing Functions #
######################
    function Base.show(io::IO, op::Projector)
        print(io, "$(summary(op)) with $(length(op)) operator(s):")
        pad = "  "
        maxlen = 30
        i = 1
        for k in keys(dict(op.ket)), b in keys(dict(op.bra))
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

export getbra,
    getket,
    hasket,
    hasbra,
    mapcoeffs,
    maplabels,
    ptrace
