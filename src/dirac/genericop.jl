###########
# OpLabel #
###########
    type OpLabel
        ktlabel::Vector{Any}
        brlabel::Vector{Any}
    end

    ktlabel(label::OpLabel) = label.ktlabel
    brlabel(label::OpLabel) = label.brlabel

    Base.reverse(label::OpLabel) = OpLabel(brlabel(label), ktlabel(label))
    Base.(:(==))(a::OpLabel, b::OpLabel) = ktlabel(a)==ktlabel(b) && brlabel(a)==brlabel(b)
    Base.hash(label::OpLabel) = hash(ktlabel(label), hash(brlabel(label)))

####################
# GenericOp/DualOp #
####################
    abstract GeneralOp{P,N} <: AbstractOperator{P,N}

    typealias OpDict Dict{OpLabel,Number}

    type GenericOp{P,N} <: GeneralOp{P,N}
        dict::OpDict
        fact::Factors{N}
        GenericOp(dict,fact) = new(dict,fact)
        GenericOp(dict,::Factors{0}) = error("Cannot construct a 0-factor operator; did you mean to construct a scalar?")
    end

    GenericOp{P,N}(::Type{P}, dict, fact::Factors{N}) = GenericOp{P,N}(dict, fact)

    type DualOp{P,N} <: GeneralOp{P,N}
        op::GenericOp{P,N}
        DualOp(op::GenericOp{P,N}) = new(op)
        DualOp(items...) = new(GenericOp{P,N}(items...))
    end

    DualOp{P,N}(op::GenericOp{P,N}) = DualOp{P,N}(op)
    DualOp(items...) = DualOp(GenericOp(items...))

    Base.convert(::Type{GenericOp}, opc::DualOp) = eager_ctran(opc.op)
    Base.convert{P}(::Type{GenericOp{P}}, opc::DualOp{P}) = convert(GenericOp, opc)
    Base.convert{P,N}(::Type{GenericOp{P,N}}, opc::DualOp{P,N}) = convert(GenericOp, opc)
    
    Base.promote_rule(::Type{GenericOp}, ::Type{DualOp}) = GenericOp
    Base.promote_rule{P}(::Type{GenericOp{P}}, ::Type{DualOp{P}}) = GenericOp{P}
    Base.promote_rule{P,N}(::Type{GenericOp{P,N}}, ::Type{DualOp{P,N}}) = GenericOp{P,N}

    dict(op::GenericOp) = op.dict
    dict(opc::DualOp) = dict(opc.op)

    fact(op::GenericOp) = op.fact
    fact(opc::DualOp) = fact(opc.op)

################
# Constructors #
################
    function GenericOp{P,N}(f::Function, kt::Ket{P,N})
        result = OpDict()
        for i in labels
            for j in labels 
                (c, new_j) = f(j)
                if length(new_j) == N 
                    result[OpLabel(i,j)] = c * inner_rule(S, i, new_j)
                else
                    throw(BoundsError())
                end
            end
        end
        return GenericOp(P,result,fact(kt))
    end

    function GenericOp{A,B,N}(kt::Ket{A,N}, br::Bra{B,N})
        result = OpDict()
        for (k,kc) in dict(kt)
            for (b,bc) in dict(br)
                result[OpLabel(k,b)] = kc * bc'
            end
        end
        return GenericOp(typejoin(A,B), result, fact(kt))
    end

    Base.copy(op::GeneralOp) = typeof(op)(copy(dict(op)), fact(op))
    Base.similar(op::GeneralOp, d=similar(dict(op))) = typeof(op)(d, fact(op))

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){P,N}(a::GenericOp{P,N}, b::GenericOp{P,N}) = dict(a) == dict(b)
    Base.(:(==)){P,N}(a::DualOp{P,N}, b::DualOp{P,N}) = a.op == b.op
    Base.(:(==))(a::AbstractOperator, b::AbstractOperator) = ==(promote(a,b)...)

    Base.hash(op::GeneralOp) = hash(dict(op), hash(typeof(op)))

    Base.length(op::GeneralOp) = length(dict(op))

    Base.getindex(op::GenericOp, label::OpLabel) = dict(op)[label]
    Base.getindex(op::GenericOp, k::Array, b::Array) = op[OpLabel(k,b)]
    Base.getindex(opc::DualOp, label::OpLabel) = opc.op[reverse(label)]'
    Base.getindex(opc::DualOp, k::Array, b::Array) = opc.op[OpLabel(b,k)]'
    Base.getindex(op::GeneralOp, k, b) = op[[k],[b]]

    function Base.setindex!{P,N}(s::GeneralOp{P,N}, c, label::OpLabel)
        if length(ktlabel(label)) == N && length(brlabel(label)) == N
            return _setindex!(s, c, label)
        else
            throw(BoundsError())
        end
    end

    function Base.setindex!{P,N}(s::GeneralOp{P,N}, c, k::Array, b::Array)
        if length(k) == N && length(b) == N
            return _setindex!(s, c, k, b)
        else
            throw(BoundsError())
        end
    end

    Base.setindex!(op::GeneralOp, c, k, b) = setindex!(op, c, [k], [b])

    Base.haskey(op::GenericOp, label::OpLabel) = haskey(dict(op), label)
    Base.haskey(opc::DualOp, label::OpLabel) = haskey(opc.op, reverse(label))
    Base.haskey(op::GeneralOp, k::Array, b::Array) = haskey(op, OpLabel(k,b))

    Base.get(op::GenericOp, label::OpLabel, default) = get(dict(op), label, default)
    Base.get(opc::DualOp, label::OpLabel, default) = haskey(opc, label) ? opc[label] : default
    Base.get(op::GeneralOp, k::Array, b::Array, default) = get(op, OpLabel(k,b))

    Base.delete!(op::GenericOp, label::OpLabel) = (delete!(dict(op), label); return op)
    Base.delete!(opc::DualOp, label::OpLabel) = delete!(opc.op, reverse(label))
    Base.delete!(op::GeneralOp, k::Array, b::Array) = delete!(op, OpLabel(k,b))

    labels(op::GenericOp) = keys(dict(op))
    labels(opc::DualOp) = imap(reverse, labels(opc.op))
    coeffs(op::GenericOp) = values(dict(op))
    coeffs(opc::DualOp) = imap(conj, coeffs(opc.op))

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!(f::Function, op::GenericOp) = (filter!(f, dict(op)); return op)
    Base.filter!(f::Function, opc::DualOp) = (filter!((k,v)->f(reverse(k),v'), dict(opc)); return opc)

    Base.filter(f::Function, op::GenericOp) = similar(op, filter(f, dict(op)))
    Base.filter(f::Function, opc::DualOp) = similar(opc, filter((k,v)->f(reverse(k),v'), dict(opc)))

    labelcheck(label::OpLabel, N) = length(ktlabel(label)) == length(brlabel(label)) == N ? label : throw(BoundsError())
    opprune(kv,N) = length(ktlabel(kv[1])) == length(brlabel(kv[1])) == N ? kv : throw(BoundsError())
    dualprune(kv,N) = (labelcheck(reverse(kv[1]),N),kv[2]')

    Base.map{P,N}(f::Function, op::GenericOp{P,N}) = similar(op, mapkv((k,v)->opprune(f(k,v),N), dict(op)))
    Base.map{P,N}(f::Function, opc::DualOp{P,N}) = similar(opc, mapkv((k,v)->dualprune(f(reverse(k),v'),N), dict(opc)))
    
    mapcoeffs!(f::Function, op::GenericOp) = (mapvals!(f, dict(op)); return op)
    mapcoeffs!(f::Function, opc::DualOp) = (mapvals!(v->f(v')', dict(opc)); return opc)
    mapcoeffs(f::Function, op::GenericOp) = similar(op, mapvals(f, dict(op)))
    mapcoeffs(f::Function, opc::DualOp) = similar(opc, mapvals(v->f(v')', dict(opc)))

    maplabels!{P,N}(f::Function, op::GenericOp{P,N}) = (mapkeys!(k->labelcheck(f(k), N), dict(op)); return op)
    maplabels!{P,N}(f::Function, opc::DualOp{P,N}) = (mapkeys!(k->reverse(labelcheck(f(reverse(k)),N)), dict(opc)); return opc)
    maplabels{P,N}(f::Function, op::GenericOp{P,N}) = similar(op, mapkeys(k->labelcheck(f(k),N), dict(op)))
    maplabels{P,N}(f::Function, opc::DualOp{P,N}) = similar(opc, mapkeys(k->reverse(labelcheck(f(reverse(k)),N)), dict(opc)))

##########################
# Mathematical Functions #
##########################
    nfactors{P,N}(::GeneralOp{P,N}) = N

    eager_ctran(op::GenericOp) = map((k,v)->(reverse(k),v'), op)
    
    Base.ctranspose(op::GenericOp) = DualOp(op)
    Base.ctranspose(opc::DualOp) = opc.op

    function inner{A,B,N}(br::Bra{A,N}, op::GenericOp{B,N})
        result = StateDict()
        for (label,v) in dict(op)
            coeff = 0
            for (b,c) in dict(br)
                coeff += c'*v*inner_eval(A,B,ktlabel(label),b) 
            end
            add_to_dict!(result, brlabel(label), coeff')
        end
        return Bra(typejoin(A,B), result, fact(op))
    end

    function inner{A,B,N}(op::GenericOp{A,N}, kt::Ket{B,N})
        result = StateDict()
        for (label,c) in dict(op)
            coeff = 0
            for (k,v) in dict(kt)
                coeff += c*v*inner_eval(A,B,brlabel(label),k) 
            end
            add_to_dict!(result, ktlabel(label), coeff)
        end
        return Ket(typejoin(A,B), result, fact(op))
    end

    function inner{A,B,N}(a::GenericOp{A,N}, b::GenericOp{B,N})
        result = OpDict()
        for (label1,v) in dict(a), (label2,c) in dict(b)
            add_to_dict!(result,
                         OpLabel(ktlabel(label1),brlabel(label2)),
                         v*c*inner_eval(A,B,brlabel(label1),ktlabel(label2)))
        end
        return GenericOp(typejoin(A,B), result, fact(a))
    end

    function inner{A,B,N}(a::GenericOp{A,N}, b::DualOp{B,N})
        result = OpDict()
        for (label1,v) in dict(a), (label2,c) in dict(b)
            add_to_dict!(result,
                         OpLabel(ktlabel(label1),ktlabel(label2)),
                         v*c'*inner_eval(A,B,brlabel(label1),brlabel(label2)))
        end
        return GenericOp(typejoin(A,B), result, fact(b))
    end

    function inner{A,B,N}(a::DualOp{A,N}, b::GenericOp{B,N})
        result = OpDict()
        for (label1,v) in dict(a), (label2,c) in dict(b)
            add_to_dict!(result,
                         OpLabel(brlabel(label1),brlabel(label2)),
                         v'*c*inner_eval(A,B,ktlabel(label1),ktlabel(label2)))
        end
        return GenericOp(typejoin(A,B), result, fact(b))
    end

    inner(br::Bra, opc::DualOp) = inner(opc.op, br')'
    inner(opc::DualOp, kt::Ket) = inner(kt', opc.op)'
    inner(a::DualOp, b::DualOp) = inner(a.op, b.op)'

    Base.(:*)(br::Bra, op::AbstractOperator) = inner(br,op)
    Base.(:*)(op::AbstractOperator, kt::Ket) = inner(op,kt)
    Base.(:*)(a::AbstractOperator, b::AbstractOperator) = inner(a,b)
    Base.(:*)(kt::Ket, br::Bra) = tensor(kt,br)

    Base.scale!(op::GenericOp, c::Number) = (dscale!(dict(op), c); return op)
    Base.scale!(c::Number, op::GenericOp) = scale!(op, c)
    Base.scale!(opc::DualOp, c::Number) = DualOp(scale!(opc.op, c'))
    Base.scale!(c::Number, opc::DualOp) = scale!(opc, c)

    Base.scale(op::GenericOp, c::Number) = similar(op, dscale(dict(op), c))
    Base.scale(c::Number, op::GenericOp) = scale(op, c)
    Base.scale(opc::DualOp, c::Number) = DualOp(scale(opc.op, c'))
    Base.scale(c::Number, opc::DualOp) = scale(opc, c)

    Base.(:*)(c::Number, op::AbstractOperator) = scale(c, op)
    Base.(:*)(op::AbstractOperator, c::Number) = scale(op, c)
    Base.(:/)(op::AbstractOperator, c::Number) = scale(op, 1/c)

    Base.(:+){P,N}(a::GenericOp{P,N}, b::GenericOp{P,N}) = similar(a,mergef(+, dict(a), dict(b)))
    Base.(:+){P,N}(a::DualOp{P,N}, b::DualOp{P,N}) = DualOp(a.op + b.op)
    Base.(:+)(a::AbstractOperator, b::AbstractOperator) = +(promote(a,b)...)

    Base.(:-)(a::AbstractOperator, b::AbstractOperator) = a + (-b)
    Base.(:-)(op::GenericOp) = mapcoeffs(-, op)
    Base.(:-)(opc::DualOp) = DualOp(-opc.op)

    Base.norm(op::GenericOp) = sqrt(sum(v->v^2, values(dict(op))))
    Base.norm(opc::DualOp) = norm(opc.op)
    
    QuBase.normalize(op::AbstractOperator) = scale(1/norm(op), op)
    QuBase.normalize!(op::AbstractOperator) = scale!(1/norm(op), op)

    function Base.trace{O<:Orthonormal}(op::GenericOp{O})
        result = 0
        for label in labels(op)
            if ktlabel(label)==brlabel(label)
                result += op[label]
            end
        end
        return result
    end

    function Base.trace{P}(op::GenericOp{P})
        result = 0
        for i in distinct(map(ktlabel, labels(op))), (label,v) in dict(op)
            result += v * inner_rule(P, i, ktlabel(label)) * inner_rule(P, brlabel(label), i)
        end
        return result
    end

    Base.trace(opc::DualOp) = trace(opc.op)'

    QuBase.tensor{P}(ops::GenericOp{P}...) = GenericOp(P,mergecart!(tensor_op, OpDict(), map(dict, ops)),mapreduce(fact, +, ops))
    QuBase.tensor(opcs::DualOp...) = tensor(map(opc->opc.op, opcs)...)'
    QuBase.tensor(ops::AbstractOperator...) = tensor(promote(ops...)...)

    xsubspace(op::GeneralOp, x) = filter((k,v)->sum(ktlabel(k))==x && sum(brlabel(k))==x, op)

    filternz!(op::GeneralOp) = filter!((k, v) -> v != 0, op)
    filternz(op::GeneralOp) = filter((k, v) -> v != 0, op)
    purity(op::GeneralOp) = trace(op^2)

    queval(f, op::AbstractOperator) = mapcoeffs(x->queval(f,x),op)

#################
# Partial trace #
#################
    ptrace{P}(op::GeneralOp{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())
    ptrace(op::GeneralOp, over) = ptrace_op(op, over)

    traced_label(::GenericOp, label, over) = OpLabel(except(ktlabel(label), over), except(brlabel(label), over))
    traced_label(::DualOp, label, over) = OpLabel(except(brlabel(label), over), except(ktlabel(label), over))

    traced_coeff(::GenericOp, i, label, v, over) = v * inner_rule(P, i[over], ktlabel(label)[over]) * inner_rule(P, brlabel(label)[over], i[over])
    traced_coeff(::DualOp, i, label, v, over) = v' * inner_rule(P, i[over], brlabel(label)[over]) * inner_rule(P, ktlabel(label)[over], i[over])

    function ptrace_op{O<:Orthonormal,N}(op::GeneralOp{O,N}, over)
        result = OpDict()
        for label in labels(op)
            if ktlabel(label)[over] == brlabel(label)[over]
                add_to_dict!(result, traced_label(op,label,over), op[label])
            end
        end
        return GenericOp(O,result,Factors{N-1}())
    end

    function ptrace_op{P,N}(op::GeneralOp{P,N}, over)
        result = OpDict()
        for i in distinct(map(ktlabel, labels(op))), (label,v) in dict(op)
            add_to_dict!(result, traced_label(op,label,over), traced_coeff(op,i,label,v,over))
        end
        return GenericOp(P,result,Factors{N-1}())
    end

######################
# Printing Functions #
######################
    labelrepr(op::GenericOp, label, pad) = "$pad$(op[label]) $(ktstr(ktlabel(label)))$(brstr(brlabel(label)))"
    labelrepr(opc::DualOp, label, pad) = "$pad$(opc[reverse(label)]) $(ktstr(brlabel(label)))$(brstr(ktlabel(label)))"

    Base.summary(op::AbstractOperator) = "$(typeof(op)) with $(length(op)) operator(s)"
    Base.show(io::IO, op::GeneralOp) = dirac_show(io, op)
    Base.showcompact(io::IO, op::GeneralOp) = dirac_showcompact(io, op)
    Base.repr(op::GeneralOp) = dirac_repr(op)

####################
# Helper Functions #
####################
    function tensor_op(pairs)
        #pairs structure is: (((op1ktlabel, op1brlabel), op1value), ((op2ktlabel, op2brlabel), op2value)...,)
        labels = map(first, pairs)
        return (OpLabel(vcat(map(ktlabel, labels)...), vcat(map(brlabel, labels)...)), prod(second, pairs))
    end

export ptrace,
    xsubspace,
    nfactors,
    maplabels!,
    mapcoeffs!,
    mapcoeffs,
    maplabels,
    filternz,
    filternz!,
    purity,
    labels,
    coeffs,
    OpLabel,
    ktlabel,
    brlabel
