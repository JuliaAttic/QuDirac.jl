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

##################
# DiracOp/DualOp #
##################
    abstract GenericOperator{P,N} <: AbstractOperator{P,N}

    typealias OpDict Dict{OpLabel,Number}

    type DiracOp{P,N} <: GenericOperator{P}
        dict::OpDict
        fact::Factors{N}
        DiracOp(dict,fact) = new(dict,fact)
        DiracOp(dict,::Factors{0}) = error("Cannot construct a 0-factor operator; did you mean to construct a scalar?")
    end

    DiracOp{P,N}(::Type{P}, dict, fact::Factors{N}) = DiracOp{P,N}(dict, fact)

    type DualOp{P,N} <: GenericOperator{P}
        op::DiracOp{P,N}
    end

    DualOp(items...) = DualOp(DiracOp(items...))

    Base.convert(::Type{DiracOp}, opc::DualOp) = eager_ctran(opc.op)
    Base.convert{P}(::Type{DiracOp{P}}, opc::DualOp{P}) = convert(DiracOp, opc)
    Base.convert{P,N}(::Type{DiracOp{P,N}}, opc::DualOp{P,N}) = convert(DiracOp, opc)
    
    Base.promote_rule(::Type{DiracOp}, ::Type{DualOp}) = DiracOp
    Base.promote_rule{P}(::Type{DiracOp{P}}, ::Type{DualOp{P}}) = DiracOp{P}
    Base.promote_rule{P,N}(::Type{DiracOp{P,N}}, ::Type{DualOp{P,N}}) = DiracOp{P,N}

    dict(op::DiracOp) = op.dict
    dict(opc::DualOp) = dict(opc.op)

    fact(op::DiracOp) = op.fact
    fact(opc::DualOp) = fact(opc.op)

################
# Constructors #
################
    function DiracOp{P,N}(f::Function, kt::Ket{P,N})
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
        return DiracOp(P,result,fact(kt))
    end

    function DiracOp{A,B,N}(kt::Ket{A,N}, br::Bra{B,N})
        result = OpDict()
        for (k,kc) in dict(kt)
            for (b,bc) in dict(br)
                result[OpLabel(k,b)] = kc * bc'
            end
        end
        return DiracOp(typejoin(A,B), result, fact(kt))
    end

    Base.copy(op::GenericOperator) = typeof(op)(copy(dict(op)), fact(op))
    Base.similar(op::GenericOperator, d=similar(dict(op))) = typeof(op)(d, fact(op))

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){P,N}(a::DiracOp{P,N}, b::DiracOp{P,N}) = dict(a) == dict(b)
    Base.(:(==)){P,N}(a::DualOp{P,N}, b::DualOp{P,N}) = a.op == b.op
    Base.(:(==))(a::AbstractOperator, b::AbstractOperator) = ==(promote(a,b)...)

    Base.hash(op::GenericOperator) = hash(dict(op), hash(typeof(op)))

    Base.length(op::GenericOperator) = length(dict(op))

    Base.getindex(op::DiracOp, label::OpLabel) = dict(op)[label]
    Base.getindex(op::DiracOp, k::Array, b::Array) = op[OpLabel(k,b)]
    Base.getindex(opc::DualOp, label::OpLabel) = opc.op[reverse(label)]'
    Base.getindex(opc::DualOp, k::Array, b::Array) = opc.op[OpLabel(b,k)]'
    Base.getindex(op::GenericOperator, k, b) = op[[k],[b]]

    _setindex!(op::DiracOp, c, label::OpLabel) = setindex!(dict(op), c, label)
    _setindex!(op::DiracOp, c, k::Array, b::Array) = setindex!(op, c, OpLabel(k,b))
    _setindex!(opc::DualOp, c, label::OpLabel) = setindex!(opc.op, c', reverse(label))
    _setindex!(opc::DualOp, c, k::Array, b::Array) = setindex!(opc.op, c', OpLabel(b,k))
    _setindex!(op::GenericOperator, c, k, b) = setindex!(op, c, [k], [b])

    Base.haskey(op::DiracOp, label::OpLabel) = haskey(dict(op), label)
    Base.haskey(opc::DualOp, label::OpLabel) = haskey(opc.op, reverse(label))
    Base.haskey(op::GenericOperator, k::Array, b::Array) = haskey(op, OpLabel(k,b))

    Base.get(op::DiracOp, label::OpLabel, default) = get(dict(op), label, default)
    Base.get(opc::DualOp, label::OpLabel, default) = haskey(opc, label) ? opc[label] : default
    Base.get(op::GenericOperator, k::Array, b::Array, default) = get(op, OpLabel(k,b))

    Base.delete!(op::DiracOp, label::OpLabel) = (delete!(dict(op), label); return op)
    Base.delete!(opc::DualOp, label::OpLabel) = delete!(opc.op, reverse(label))
    Base.delete!(op::GenericOperator, k::Array, b::Array) = delete!(op, OpLabel(k,b))

    labels(op::DiracOp) = keys(dict(op))
    labels(opc::DualOp) = imap(reverse, labels(opc.op))
    coeffs(op::DiracOp) = values(dict(op))
    coeffs(opc::DualOp) = imap(conj, coeffs(opc.op))

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!(f::Function, op::DiracOp) = (filter!(f, dict(op)); return op)
    Base.filter!(f::Function, opc::DualOp) = (filter!((k,v)->f(reverse(k),v'), opc.op); return opc)

    Base.filter(f::Function, op::DiracOp) = similar(op, filter(f, dict(op)))
    Base.filter(f::Function, opc::DualOp) = DualOp(filter((k,v)->f(reverse(k),v'), opc.op))

    Base.map(f::Function, op::DiracOp) = similar(op, mapkv(f, dict(op)))
    Base.map(f::Function, opc::DualOp) = mapkv!((k,v)->f(reverse(k),v'), similar(opc), opc.op)
    
    mapcoeffs(f::Function, op::DiracOp) = similar(op, mapvals(f, dict(op)))
    mapcoeffs(f::Function, opc::DualOp) = mapvals!(v->f(v'), similar(opc), opc.op)

    maplabels(f::Function, op::DiracOp) = similar(op, mapkeys(f, dict(op)))
    maplabels(f::Function, opc::DualOp) = mapkeys!(k->f(reverse(k)), similar(opc), opc.op)

##########################
# Mathematical Functions #
##########################
    nfactors{P,N}(::GenericOperator{P,N}) = N

    eager_ctran(op::DiracOp) = map((k,v)->(reverse(k),v'), op)
    
    Base.ctranspose(op::DiracOp) = DualOp(op)
    Base.ctranspose(opc::DualOp) = opc.op

    function inner{A,B,N}(br::Bra{A,N}, op::DiracOp{B,N})
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

    function inner{A,B,N}(op::DiracOp{A,N}, kt::Ket{B,N})
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

    function inner{A,B,N}(a::DiracOp{A,N}, b::DiracOp{B,N})
        result = OpDict()
        for (label1,v) in dict(a), (label2,c) in dict(b)
            add_to_dict!(result,
                         OpLabel(ktlabel(label1),brlabel(label2)),
                         v*c*inner_eval(A,B,brlabel(label1),ktlabel(label2)))
        end
        return DiracOp(typejoin(A,B), result, fact(a))
    end

    function inner{A,B,N}(a::DiracOp{A,N}, b::DualOp{B,N})
        result = OpDict()
        for (label1,v) in dict(a), (label2,c) in dict(b)
            add_to_dict!(result,
                         OpLabel(ktlabel(label1),ktlabel(label2)),
                         v*c'*inner_eval(A,B,brlabel(label1),brlabel(label2)))
        end
        return DiracOp(typejoin(A,B), result, fact(b))
    end

    function inner{A,B,N}(a::DualOp{A,N}, b::DiracOp{B,N})
        result = OpDict()
        for (label1,v) in dict(a), (label2,c) in dict(b)
            add_to_dict!(result,
                         OpLabel(brlabel(label1),brlabel(label2)),
                         v'*c*inner_eval(A,B,ktlabel(label1),ktlabel(label2)))
        end
        return DiracOp(typejoin(A,B), result, fact(b))
    end

    inner(br::Bra, opc::DualOp) = inner(opc.op, br')'
    inner(opc::DualOp, kt::Ket) = inner(kt', opc.op)'
    inner(a::DualOp, b::DualOp) = inner(a.op, b.op)'

    Base.(:*)(br::Bra, op::AbstractOperator) = inner(br,op)
    Base.(:*)(op::AbstractOperator, kt::Ket) = inner(op,kt)
    Base.(:*)(a::AbstractOperator, b::AbstractOperator) = inner(a,b)
    Base.(:*)(kt::Ket, br::Bra) = tensor(kt,br)

    Base.scale!(c::Number, op::DiracOp) = (castvals!(*, c, dict(op)); return op)
    Base.scale!(op::DiracOp, c::Number) = (castvals!(*, dict(op), c); return op)
    Base.scale!(c::Number, opc::DualOp) = DualOp(scale!(c',opc.op))
    Base.scale!(opc::DualOp, c::Number) = DualOp(scale!(opc.op,c'))

    Base.scale(c::Number, op::DiracOp) = similar(op, castvals(*, c, dict(op)))
    Base.scale(op::DiracOp, c::Number) = similar(op, castvals(*, dict(op), c))
    Base.scale(c::Number, opc::DualOp) = DualOp(scale(c',opc.op))
    Base.scale(opc::DualOp, c::Number) = DualOp(scale(opc.op,c'))

    Base.(:*)(c::Number, op::AbstractOperator) = scale(c, op)
    Base.(:*)(op::AbstractOperator, c::Number) = scale(op, c)
    Base.(:/)(op::AbstractOperator, c::Number) = scale(op, 1/c)

    Base.(:+){P,N}(a::DiracOp{P,N}, b::DiracOp{P,N}) = similar(a,mergef(+, dict(a), dict(b)))
    Base.(:+){P,N}(a::DualOp{P,N}, b::DualOp{P,N}) = DualOp(a.op + b.op)
    Base.(:+)(a::AbstractOperator, b::AbstractOperator) = +(promote(a,b)...)

    Base.(:-)(a::AbstractOperator, b::AbstractOperator) = a + (-b)
    Base.(:-)(op::DiracOp) = mapcoeffs(-, op)
    Base.(:-)(opc::DualOp) = DualOp(-opc.op)

    Base.norm(op::DiracOp) = sqrt(sum(v->v^2, values(dict(op))))
    Base.norm(opc::DualOp) = norm(opc.op)
    
    QuBase.normalize(op::AbstractOperator) = scale(1/norm(op), op)
    QuBase.normalize!(op::AbstractOperator) = scale!(1/norm(op), op)

    function Base.trace{O<:Orthonormal}(op::DiracOp{O})
        result = 0
        for label in keys(dict(op))
            if ktlabel(label)==brlabel(label)
                result += op[label]
            end
        end
        return result
    end

    function Base.trace{P}(op::DiracOp{P})
        result = 0
        for (label,v) in dict(op)
            result += v * inner_rule(P, brlabel(label), ktlabel(label))
        end
        return result
    end

    Base.trace(opc::DualOp) = trace(opc.op)'

    QuBase.tensor{P}(ops::DiracOp{P}...) = DiracOp(P,mergecart!(tensor_op, OpDict(), map(dict, ops)),mapreduce(fact, +, ops))
    QuBase.tensor(opcs::DualOp...) = tensor(map(opc->opc.op, opcs)...)'
    QuBase.tensor(ops::AbstractOperator...) = tensor(promote(ops...)...)

    xsubspace(op::GenericOperator, x) = filter((k,v)->sum(ktlabel(k))==x && sum(brlabel(k))==x, op)

    filternz!(op::GenericOperator) = filter!((k, v) -> v != 0, op)
    filternz(op::GenericOperator) = filter((k, v) -> v != 0, op)
    purity(op::GenericOperator) = trace(op^2)

#################
# Partial trace #
#################
    ptrace(opc::DualOp, over::Integer) = ptrace(opc.op, over)'
    ptrace{P}(op::DiracOp{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())
    ptrace(op::DiracOp, over) = ptrace_op!(op, over)

    function ptrace_op!{O<:Orthonormal,N}(op::DiracOp{O,N}, over)
        result = OpDict()
        for label in keys(dict(op))
            if ktlabel(label)[over] == brlabel(label)[over]
                add_to_dict!(result,
                             OpLabel(except(ktlabel(label), over), except(brlabel(label), over)),
                             op[label])
            end
        end
        return DiracOp(O,result,Factors{N-1}())
    end

    function ptrace_op!{P,N}(op::DiracOp{P,N}, over)
        result = OpDict()
        for (label,v) in dict(op)
            add_to_dict!(result,
                         OpLabel(except(ktlabel(label), over), except(brlabel(label), over)),
                         v*inner_rule(P, ktlabel(label)[over], brlabel(label)[over]))
        end
        return DiracOp(P,result,Factors{N-1}())
    end

######################
# Printing Functions #
######################
    labelrepr(op::DiracOp, label, pad) = "$pad$(op[label]) $(ktstr(ktlabel(label)))$(brstr(brlabel(label)))"
    labelrepr(opc::DualOp, label, pad) = "$pad$(opc[label]) $(ktstr(brlabel(label)))$(brstr(ktlabel(label)))"

    Base.summary(op::AbstractOperator) = "$(typeof(op)) with $(length(op)) operator(s)"
    Base.show(io::IO, op::GenericOperator) = dirac_show(io, op)
    Base.showcompact(io::IO, op::GenericOperator) = dirac_showcompact(io, op)
    Base.repr(op::GenericOperator) = dirac_repr(op)

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
    mapcoeffs,
    maplabels,
    filternz,
    filternz!,
    purity,
    labels,
    coeffs
