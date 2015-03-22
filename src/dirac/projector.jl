#############
# Projector #
#############
    type Projector{P,N} <: AbstractOperator{P,N}
        scalar::Number
        kt::Ket{P,N}
        br::Bra{P,N}
    end
    
    Projector{P,N}(::Type{P}, scalar, kt::Ket{P,N}, br::Bra{P,N}) = Projector{P,N}(scalar, kt, br)

    Base.copy{P}(op::Projector{P}) = Projector(P, copy(op.scalar), copy(op.kt), copy(op.br))

    Base.convert(::Type{DiracOp}, op::Projector) = scale!(op.scalar, DiracOp(op.kt, op.br))
    Base.convert{P}(::Type{DiracOp{P}}, op::Projector{P}) = convert(DiracOp, op)
    Base.convert{P,N}(::Type{DiracOp{P,N}}, op::Projector{P,N}) = convert(DiracOp, op)

    Base.promote_rule(::Type{DiracOp}, ::Type{Projector}) = DiracOp
    Base.promote_rule{P}(::Type{DiracOp{P}}, ::Type{Projector{P}}) = DiracOp{P}
    Base.promote_rule{P,N}(::Type{DiracOp{P,N}}, ::Type{Projector{P,N}}) = DiracOp{P,N}

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){P}(a::Projector{P}, b::Projector{P}) = a.scalar == b.scalar && a.kt == b.kt && a.br == b.br

    Base.hash(op::Projector) = hash(op.scalar, hash(op.kt, hash(op.br)))
    Base.length(op::Projector) = length(op.kt)*length(op.br)

    Base.getindex(op::Projector, k::Array, b::Array) = op.scalar * op.kt[k] * op.br[b]
    Base.getindex(op::Projector, label::OpLabel) = op[ktlabel(label), brlabel(label)]
    Base.getindex(op::Projector, k, b) = op[[k],[b]]

    # would be great if the below worked with normal indexing
    # notation (e.g. op[k,:]) but slice notation is apparently
    # special and doesn't dispatch directly to getindex...
    # Base.getindex(op::Projector, k, ::Colon) = (op.scalar * op.kt[k]) * op.br
    # Base.getindex(op::Projector, ::Colon, b) = (op.scalar * op.br[b]) * op.kt
    # Base.getindex(op::Projector, ::Colon, ::Colon) = convert(DiracOp, op)

    getbra(op::Projector, k::Array) = (op.scalar * op.kt[k]) * op.br
    getket(op::Projector, b::Array) = (op.scalar * op.br[b]) * op.kt

    Base.haskey(op::Projector, k::Array, b::Array) = hasket(op,k) && hasbra(op, b)
    Base.haskey(op::Projector, label::OpLabel) = haskey(op, ktlabel(label), brlabel(label))
    hasket(op::Projector, label::Array) = haskey(op.kt, label)
    hasbra(op::Projector, label::Array) = haskey(op.br, label)

    Base.get(op::Projector, label::OpLabel, default) = get(op, ktlabel(label), brlabel(label), default)
    Base.get(op::Projector, k::Array, b::Array, default) = haskey(op, k, b) ? op[k,b] : default

    labels(op::Projector) = imap(pair->OpLabel(pair[1],pair[2]), product(labels(op.kt), labels(op.br)))
    coeffs(op::Projector) = imap(v->op.scalar*v[1]*v[2], product(coeffs(op.kt), coeffs(op.br)))

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter(f::Function, op::Projector) = filter(f, convert(DiracOp, op))
    Base.map(f::Function, op::Projector) = map(f, convert(DiracOp, op))

    mapcoeffs(f::Function, op::Projector) = mapcoeffs(f, convert(DiracOp, op))
    maplabels(f::Function, op::Projector) = maplabels(f, convert(DiracOp, op))

##########################
# Mathematical Functions #
##########################
    nfactors{P,N}(op::Projector{P,N}) = N

    Base.ctranspose{P}(op::Projector{P}) = Projector(P, op.scalar', op.br', op.kt')

    Base.scale!(c::Number, op::Projector) = (op.scalar = c*op.scalar; return op)
    Base.scale!(op::Projector, c::Number) = (op.scalar = op.scalar*c; return op)

    Base.scale(c::Number, op::Projector) = scale!(c,copy(op))
    Base.scale(op::Projector, c::Number) = scale!(copy(op),c)

    Base.(:-)(op::Projector) = (op.scalar = -op.scalar)
    Base.(:+)(a::Projector, b::Projector) = convert(DiracOp, a) + convert(DiracOp, b)

    function Base.norm(op::Projector)
        result = 0
        for v in values(dict(op.kt)), c in values(dict(op.br))
            result += (v*c)^2
        end
        return sqrt(result)
    end

    function Base.trace{O<:Orthonormal}(op::Projector{O})
        result = 0
        for k in keys(dict(op.kt)), b in keys(dict(op.br))
            if b==k
                result += op[k,b]
            end
        end
        return result
    end

    function Base.trace{P}(op::Projector{P})
        result = 0
        for k in keys(dict(op.kt)), b in keys(dict(op.br))
            result += op[k,b] * inner_rule(P, b, k)
        end
        return result
    end

    inner(br::Bra, op::Projector) = op.scalar * inner(br, op.kt) * op.br
    inner(op::Projector, kt::Ket) = op.scalar * op.kt * inner(op.br, kt)
    inner(a::Projector, b::Projector) = Projector(a.scalar * b.scalar * inner(a.br,b.kt), a.kt, b.br)
    inner(a::Projector, b::GenericOperator) = a.scalar * a.kt * inner(a.br, b)
    inner(a::GenericOperator, b::Projector) = inner(a, b.kt) * b.br * b.scalar

    QuBase.tensor(kt::Ket, br::Bra) = Projector(1, kt, br)
    QuBase.tensor(br::Bra, kt::Ket) = tensor(kt, br)
    QuBase.tensor(a::Projector, b::Projector) = Projector(a.scalar * b.scalar, tensor(a.kt,b.kt), tensor(a.br, b.br))

    xsubspace(op::Projector,x) = xsubspace(convert(DiracOp, op), x)
    filternz(op::Projector) = filternz(convert(DiracOp, op))
    purity(op::Projector) = trace(op^2)

    ptrace{P}(op::Projector{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())
    ptrace(op::Projector, over) = ptrace_proj!(op, over)

    function ptrace_proj!{O<:Orthonormal,N}(op::Projector{O,N}, over)
        result = OpDict()
        for k in keys(dict(op.kt)), b in keys(dict(op.br))
            if k[over] == b[over]
                add_to_dict!(result,
                             OpLabel(except(k, over), except(b, over)),
                             op[k,b])
            end
        end
        return DiracOp(O,result,Factors{N-1}())
    end

    function ptrace_proj!{P,N}(op::Projector{P,N}, over)
        result = OpDict()
        for k in keys(dict(op.kt)), b in keys(dict(op.br))
            add_to_dict!(result,
                         OpLabel(except(k, over), except(b, over)),
                         op[k,b]*inner_rule(P, k[over], b[over]))
        end
        return DiracOp(P,result,Factors{N-1}())
    end

######################
# Printing Functions #
######################
    labelrepr(op::Projector, k, b, pad) = "$pad$(op[k,b]) $(ktstr(k))$(brstr(b))"

    function Base.show(io::IO, op::Projector)
        print(io, summary(op)*":")
        pad = "  "
        maxlen = 4
        for k in take(keys(dict(op.kt)), maxlen),
            b in take(keys(dict(op.br)), maxlen)
            println(io)
            print(io, labelrepr(op, k, b, pad))
        end
        if length(op) > maxlen^2
            println(io)
            print(io, "$pad$vdots")
        end
    end

export getbra,
    getket,
    hasket,
    hasbra,
    mapcoeffs,
    maplabels,
    ptrace,
    purity,
    labels,
    coeffs
