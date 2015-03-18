###########
# OpLabel #
###########
    type OpLabel
        ketlabel::Vector{Any}
        bralabel::Vector{Any}
    end

    ketlabel(label::OpLabel) = label.ketlabel
    bralabel(label::OpLabel) = label.bralabel

    Base.reverse(label::OpLabel) = OpLabel(bralabel(label), ketlabel(label))
    Base.(:(==))(a::OpLabel, b::OpLabel) = ketlabel(a)==ketlabel(b) && bralabel(a)==bralabel(b)
    Base.hash(label::OpLabel) = hash(ketlabel(label), hash(bralabel(label)))

##################
# DiracOp/DualOp #
##################
    abstract GenericOperator{P} <: AbstractOperator{P}

    typealias OpDict Dict{OpLabel,Number}

    type DiracOp{P<:AbstractInner} <: GenericOperator{P}
        dict::OpDict
        DiracOp() = new(OpDict())
        DiracOp(dict) = new(dict)
    end

    type DualOp{P} <: GenericOperator{P}
        op::DiracOp{P}
        DualOp(items...) = new(DiracOp{P}(items...))
        DualOp(op::DiracOp{P}) = new(op)
    end

    DualOp{P}(op::DiracOp{P}) = DualOp{P}(op)
    DualOp(items...) = DualOp(DiracOp(items...))

    Base.convert{P}(::Type{DiracOp{P}}, opc::DualOp{P}) = eager_ctran(opc.op)
    Base.promote_rule{P}(::Type{DiracOp{P}}, ::Type{DualOp{P}}) = DiracOp{P}

################
# Constructors #
################
    function DiracOp{P}(f::Function, ket::Ket{P})
        return DiracOp(f, S, keys(dict((ket))))
    end

    # f(label) -> (newval, newlabel)
    function DiracOp{P}(f::Function, ::Type{P}, labels)
        result = OpDict()
        for i in labels
            for j in labels 
                (c, new_j) = f(j)
                result[OpLabel(i,j)] = c * inner_rule(S, i, new_j)
            end
        end
        return DiracOp{P}(result)
    end

    function DiracOp{A,B}(ket::Ket{A}, bra::Bra{B})
        result = OpDict()
        for (k,kc) in dict(ket)
            for (b,bc) in dict(bra)
                result[OpLabel(k,b)] = kc * bc'
            end
        end
        return DiracOp{typejoin(A,B)}(result)
    end

    dict(op::DiracOp) = op.dict
    dict(opc::DualOp) = dict(opc.op)

    Base.copy(op::GenericOperator) = typeof(op)(copy(dict(op)))
    Base.similar(op::GenericOperator) = typeof(op)(similar(dict(op)))

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){P}(a::DiracOp{P}, b::DiracOp{P}) = dict(a) == dict(b)
    Base.(:(==)){P}(a::DualOp{P}, b::DualOp{P}) = a.op == b.op
    Base.(:(==)){P}(a::AbstractOperator{P}, b::AbstractOperator{P}) = ==(promote(a,b)...)

    Base.hash(op::GenericOperator) = hash(dict(op))

    Base.length(op::GenericOperator) = length(dict(op))

    Base.getindex(op::DiracOp, label::OpLabel) = dict(op)[label]
    Base.getindex(op::DiracOp, k::Array, b::Array) = op[OpLabel(k,b)]
    Base.getindex(opc::DualOp, label::OpLabel) = opc.op[reverse(label)]'
    Base.getindex(opc::DualOp, k::Array, b::Array) = opc.op[OpLabel(b,k)]'
    Base.getindex(op::GenericOperator, k, b) = op[[k],[b]]

    Base.setindex!(op::DiracOp, c, label::OpLabel) = setindex!(dict(op), c, label)
    Base.setindex!(op::DiracOp, c, k::Array, b::Array) = setindex!(op, c, OpLabel(k,b))
    Base.setindex!(opc::DualOp, c, label::OpLabel) = setindex!(opc.op, c', reverse(label))
    Base.setindex!(opc::DualOp, c, k::Array, b::Array) = setindex!(opc.op, c', OpLabel(b,k))
    Base.setindex!(op::GenericOperator, c, k, b) = setindex!(op, c, [k], [b])

    Base.haskey(op::DiracOp, label::OpLabel) = haskey(dict(op), label)
    Base.haskey(opc::DualOp, label::OpLabel) = haskey(opc.op, reverse(label))
    Base.haskey(op::GenericOperator, k::Array, b::Array) = haskey(op, OpLabel(k,b))

    Base.get(op::DiracOp, label::OpLabel, default) = get(dict(op), label, default)
    Base.get(opc::DualOp, label::OpLabel, default) = haskey(opc, label) ? opc[label] : default
    Base.get(op::GenericOperator, k::Array, b::Array, default) = get(op, OpLabel(k,b))

    Base.delete!(op::DiracOp, label::OpLabel) = (delete!(dict(op), label); return op)
    Base.delete!(opc::DualOp, label::OpLabel) = delete!(opc.op, reverse(label))
    Base.delete!(op::GenericOperator, k::Array, b::Array) = delete!(op, OpLabel(k,b))

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!(f::Function, op::DiracOp) = (filter!(f, dict(op)); return op)
    Base.filter!(f::Function, opc::DualOp) = (filter!((k,v)->f(reverse(k),v'), opc.op); return opc)

    Base.filter{P}(f::Function, op::DiracOp{P}) = DiracOp{P}(filter(f, dict(op)))
    Base.filter(f::Function, opc::DualOp) = DualOp(filter((k,v)->f(reverse(k),v'), opc.op))

    Base.map{P}(f::Function, op::DiracOp{P}) = DiracOp{P}(mapkv(f, dict(op)))
    Base.map(f::Function, opc::DualOp) = mapkv!((k,v)->f(reverse(k),v'), similar(opc), opc.op)
    
    mapcoeffs{P}(f::Function, op::DiracOp{P}) = DiracOp{P}(mapvals(f, dict(op)))
    mapcoeffs(f::Function, opc::DualOp) = mapvals!(v->f(v'), similar(opc), opc.op)

    maplabels{P}(f::Function, op::DiracOp{P}) = DiracOp{P}(mapkeys(f, dict(op)))
    maplabels(f::Function, opc::DualOp) = mapkeys!(k->f(reverse(k)), similar(opc), opc.op)

##########################
# Mathematical Functions #
##########################
    eager_ctran(op::DiracOp) = map((k,v)->(reverse(k),v'), op)
    
    Base.ctranspose(op::DiracOp) = DualOp(op)
    Base.ctranspose(opc::DualOp) = opc.op

    function inner{A,B}(bra::Bra{A}, op::DiracOp{B})
        result = Bra{typejoin(A,B)}()
        for (label,v) in dict(op)
            if !haskey(result, bralabel(label))
                result[bralabel(label)] = 0
            end
            coeff = 0
            for (b,c) in dict(bra)
                coeff += c'*v*inner_eval(A,B,ketlabel(label),b) 
            end
            result[bralabel(label)] += coeff
        end
        return result
    end

    function inner{A,B}(op::DiracOp{A}, ket::Ket{B})
        result = Ket{typejoin(A,B)}()
        for (label,c) in dict(op)
            if !haskey(result, ketlabel(label))
                result[ketlabel(label)] = 0
            end
            coeff = 0
            for (k,v) in dict(ket)
                coeff += c*v*inner_eval(A,B,bralabel(label),k) 
            end
            result[ketlabel(label)] += coeff
        end
        return result
    end

    function inner{A,B}(a::DiracOp{A}, b::DiracOp{B})
        result = DiracOp{typejoin(A,B)}()
        for (label1,v) in dict(a)
            for (label2,c) in dict(b)
                if !haskey(result, OpLabel(ketlabel(label1),bralabel(label2)))
                    result[ketlabel(label1),bralabel(label2)] = 0
                end
                result[ketlabel(label1),bralabel(label2)] += v*c*inner_eval(A,B,bralabel(label1),ketlabel(label2)) 
            end
        end
        return result
    end

    function inner{A,B}(a::DiracOp{A}, b::DualOp{B})
        result = DiracOp{typejoin(A,B)}()
        for (label1,v) in dict(a)
            for (label2,c) in dict(b)
                if !haskey(result, OpLabel(ketlabel(label1),ketlabel(label2)))
                    result[ketlabel(label1),ketlabel(label2)] = 0
                end
                result[ketlabel(label1),ketlabel(label2)] += v*c'*inner_eval(A,B,bralabel(label1),bralabel(label2)) 
            end
        end
        return result
    end

    function inner{A,B}(a::DualOp{A}, b::DiracOp{B})
        result = DiracOp{typejoin(A,B)}()
        for (label1,v) in dict(a)
            for (label2,c) in dict(b)
                if !haskey(result, OpLabel(bralabel(label1),bralabel(label2)))
                    result[bralabel(label1),bralabel(label2)] = 0
                end
                result[bralabel(label1),bralabel(label2)] += v'*c*inner_eval(A,B,ketlabel(label1),ketlabel(label2)) 
            end
        end
        return result
    end

    inner(bra::Bra, opc::DualOp) = inner(opc.op, bra')'
    inner(opc::DualOp, ket::Ket) = inner(ket', opc.op)'
    inner(a::DualOp, b::DualOp) = inner(a.op, b.op)'

    Base.(:*)(bra::Bra, op::AbstractOperator) = inner(bra,op)
    Base.(:*)(op::AbstractOperator, ket::Ket) = inner(op,ket)
    Base.(:*)(ket::Ket, op::AbstractOperator) = tensor(ket,op)
    Base.(:*)(op::AbstractOperator, bra::Bra) = tensor(op,bra)
    Base.(:*)(a::AbstractOperator, b::AbstractOperator) = inner(a,b)
    Base.(:*)(ket::Ket, bra::Bra) = tensor(ket,bra)

    Base.scale!(c::Number, op::GenericOperator) = (castvals!(*, c, dict(op)); return op)
    Base.scale!(op::GenericOperator, c::Number) = (castvals!(*, dict(op), c); return op)

    Base.scale(c::Number, op::GenericOperator) = typeof(op)(castvals(*, c, dict(op)))
    Base.scale(op::GenericOperator, c::Number) = typeof(op)(castvals(*, dict(op), c))

    Base.(:*)(c::Number, op::AbstractOperator) = scale(c, op)
    Base.(:*)(op::AbstractOperator, c::Number) = scale(op, c)
    Base.(:/)(op::AbstractOperator, c::Number) = scale(op, 1/c)

    Base.(:+){P}(a::DiracOp{P}, b::DiracOp{P}) = DiracOp{P}(mergef(+, dict(a), dict(b)))
    Base.(:+){P}(a::DualOp{P}, b::DualOp{P}) = DualOp(a.op + b.op)
    Base.(:+){P}(a::AbstractOperator{P}, b::AbstractOperator{P}) = +(promote(a,b)...)

    Base.(:-){P}(a::AbstractOperator{P}, b::AbstractOperator{P}) = a + (-b)
    Base.(:-)(op::DiracOp) = mapcoeffs(-, op)
    Base.(:-)(opc::DualOp) = DualOp(-opc.op)

    Base.norm(op::DiracOp) = sqrt(sum(v->v^2, values(dict(op))))
    Base.norm(opc::DualOp) = norm(opc.op)
    
    QuBase.normalize(op::AbstractOperator) = scale(1/norm(op), op)
    QuBase.normalize!(op::AbstractOperator) = scale!(1/norm(op), op)

    function Base.trace{O<:Orthonormal}(op::DiracOp{O})
        result = 0
        for label in keys(dict(op))
            if ketlabel(label)==bralabel(label)
                result += op[label]
            end
        end
        return result
    end
    Base.trace(opc::DualOp) = trace(opc.op)'

    QuBase.tensor{P}(ops::DiracOp{P}...) = DiracOp{P}(mergecart!(tensor_op, OpDict(), map(dict, ops)))
    QuBase.tensor{P}(ket::Ket{P}, op::DiracOp{P}) = DiracOp{P}(mergecart!(tensor_ket_to_op, OpDict(), dict(ket), dict(op)))
    QuBase.tensor{P}(op::DiracOp{P}, ket::Ket{P}) = DiracOp{P}(mergecart!(tensor_op_to_ket, OpDict(), dict(op), dict(ket)))
    QuBase.tensor{P}(op::DiracOp{P}, bra::Bra{P}) = DiracOp{P}(mergecart!(tensor_bra_to_op, OpDict(), dict(op), mapvals(ctranspose, dict(bra))))
    QuBase.tensor{P}(bra::Bra{P}, op::DiracOp{P}) = DiracOp{P}(mergecart!(tensor_op_to_bra, OpDict(), mapvals(ctranspose, dict(bra)), dict(op)))
    QuBase.tensor(ops::AbstractOperator...) = tensor(promote(ops...)...)
    QuBase.tensor(ket::Ket, opc::DualOp) = tensor(ket', opc.op)'
    QuBase.tensor(opc::DualOp, ket::Ket) = tensor(opc.op, ket')'
    QuBase.tensor(opc::DualOp, bra::Bra) = tensor(opc.op, bra')'
    QuBase.tensor(bra::Bra, opc::DualOp) = tensor(bra', opc.op)'
    QuBase.tensor(a::DualOp, b::DualOp) = tensor(a.op, b.op)'

    xsubspace(op::GenericOperator, x) = filter((k,v)->sum(ketlabel(k))==x && sum(bralabel(k))==x, op)

    filternz!(op::GenericOperator) = filter!((k, v) -> v != 0, op)
    filternz(op::GenericOperator) = filter((k, v) -> v != 0, op)
    purity(op::GenericOperator) = trace(op^2)

#################
# Partial trace #
#################
    ptrace(opc::DualOp, over::Integer) = DualOp(ptrace(opc.op, over))

    function ptrace{O<:Orthonormal}(op::DiracOp{O}, over::Integer)
        return DiracOp{O}(ptrace_op!(OpDict(), op, over))
    end

    function ptrace_op!(result::OpDict, op::DiracOp, over)
        for label in keys(dict(op))
            if ketlabel(label)[over] == bralabel(label)[over]
                add_to_dict!(result,
                             OpLabel(except(ketlabel(label), over), except(bralabel(label), over)),
                             op[label])
            end
        end
        return result
    end


######################
# Printing Functions #
######################
    labelrepr(op::DiracOp, label, pad) = "$pad$(op[label]) $(ketstr(ketlabel(label)))$(brastr(bralabel(label)))"
    labelrepr(opc::DualOp, label, pad) = "$pad$(opc[label]) $(ketstr(bralabel(label)))$(brastr(ketlabel(label)))"

    Base.summary(op::AbstractOperator) = "$(typeof(op)) with $(length(op)) operator(s)"
    Base.show(io::IO, op::GenericOperator) = dirac_show(io, op)
    Base.showcompact(io::IO, op::GenericOperator) = dirac_showcompact(io, op)
    Base.repr(op::GenericOperator) = dirac_repr(op)

####################
# Helper Functions #
####################
    function tensor_op(pairs)
        #pairs structure is: (((op1ketlabel, op1bralabel), op1value), ((op2ketlabel, op2bralabel), op2value)...,)
        labels = map(first, pairs)
        return (OpLabel(vcat(map(ketlabel, labels)...), vcat(map(bralabel, labels)...)), prod(second, pairs))
    end

    function tensor_ket_to_op(pairs)
        #pairs structure is: ((ketlabel, ketvalue), ((opketlabel, opbralabel), opvalue))
        return (OpLabel(vcat(pairs[1][1], ketlabel(pairs[2][1])), bralabel(pairs[2][1])), prod(second, pairs))
    end

    function tensor_op_to_ket(pairs)
        #pairs structure is: (((opketlabel, opbralabel), (ketlabel, ketvalue), opvalue))
        return (OpLabel(vcat(ketlabel(pairs[1][1]), pairs[2][1]), bralabel(pairs[1][1])), prod(second, pairs))
    end

    function tensor_bra_to_op(pairs)
        #pairs structure is: (((opketlabel, opbralabel), opvalue), (bralabel, bravalue))
        return (OpLabel(ketlabel(pairs[1][1]), vcat(bralabel(pairs[1][1]), pairs[2][1])), prod(second, pairs))
    end

    function tensor_op_to_bra(pairs)
        #pairs structure is: ((bralabel, bravalue), ((opketlabel, opbralabel), opvalue))
        return (OpLabel(ketlabel(pairs[2][1]), vcat(pairs[1][1], bralabel(pairs[2][1]))), prod(second, pairs))
    end

export ptrace,
    xsubspace,
    mapcoeffs,
    maplabels,
    filternz,
    filternz!,
    purity
