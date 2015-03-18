##################
# DiracOp/DualOp #
##################
    abstract GenericOperator{S} <: AbstractOperator{S}

    typealias OpDict Dict{(Tuple,Tuple),Number}

    type DiracOp{S} <: GenericOperator{S}
        dict::OpDict
        DiracOp() = new(OpDict())
        DiracOp(dict) = new(dict)
    end

    type DualOp{S} <: GenericOperator{S}
        op::DiracOp{S}
        DualOp(items...) = new(DiracOp{S}(items...))
        DualOp(op::DiracOp{S}) = new(op)
    end

    DualOp{S}(op::DiracOp{S}) = DualOp{S}(op)
    DualOp(items...) = DualOp(DiracOp(items...))

    Base.convert{S}(::Type{DiracOp{S}}, opc::DualOp{S}) = eager_ctran(opc.op)
    Base.promote_rule{S}(::Type{DiracOp{S}}, ::Type{DualOp{S}}) = DiracOp{S}

################
# Constructors #
################
    function DiracOp{S}(f::Function, ket::Ket{S})
        return DiracOp(f, S, keys(dict((ket))))
    end

    # f(label) -> (newval, newlabel)
    function DiracOp{S}(f::Function, ::Type{S}, labels)
        result = OpDict()
        for i in labels
            for j in labels 
                (c, new_j) = f(j)
                result[i,j] = c * inner_rule(S, i, new_j)
            end
        end
        return DiracOp{S}(result)
    end

    function DiracOp{A,B}(ket::Ket{A}, bra::Bra{B})
        result = OpDict()
        for (k,kc) in dict(ket)
            for (b,bc) in dict(bra)
                result[k,b] = kc * bc'
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
    Base.(:(==)){S}(a::DiracOp{S}, b::DiracOp{S}) = dict(a) == dict(b)
    Base.(:(==)){S}(a::DualOp{S}, b::DualOp{S}) = a.op == b.op
    Base.(:(==)){S}(a::AbstractOperator{S}, b::AbstractOperator{S}) = ==(promote(a,b)...)

    Base.hash(op::GenericOperator) = hash(dict(op), hash(typeof(op)))

    Base.length(op::GenericOperator) = length(dict(op))

    Base.getindex(op::DiracOp, label::(Tuple,Tuple)) = dict(op)[label]
    Base.getindex(op::DiracOp, k::Tuple, b::Tuple) = op[(k,b)]
    Base.getindex(opc::DualOp, label::(Tuple,Tuple)) = opc.op[reverse(label)]'
    Base.getindex(opc::DualOp, k::Tuple, b::Tuple) = opc.op[(b,k)]'
    Base.getindex(op::GenericOperator, k, b) = op[tuple(k),tuple(b)]

    Base.setindex!(op::DiracOp, c, label::(Tuple,Tuple)) = setindex!(dict(op), c, label)
    Base.setindex!(op::DiracOp, c, k::Tuple, b::Tuple) = setindex!(op, c, (k,b))
    Base.setindex!(opc::DualOp, c, label::(Tuple,Tuple)) = setindex!(opc.op, c', reverse(label))
    Base.setindex!(opc::DualOp, c, k::Tuple, b::Tuple) = setindex!(opc.op, c', (b,k))
    Base.setindex!(op::GenericOperator, c, k, b) = setindex!(op, c, tuple(k), tuple(b))

    Base.haskey(op::DiracOp, label::(Tuple,Tuple)) = haskey(dict(op), label)
    Base.haskey(opc::DualOp, label::(Tuple,Tuple)) = haskey(opc.op, reverse(label))

    Base.get(op::DiracOp, label::(Tuple,Tuple), default) = get(dict(op), label, default)
    Base.get(opc::DualOp, label::(Tuple,Tuple), default) = haskey(opc, label) ? opc[label] : default

    Base.delete!(op::DiracOp, label::(Tuple,Tuple)) = (delete!(dict(op), label); return op)
    Base.delete!(opc::DualOp, label::(Tuple,Tuple)) = delete!(opc.op, reverse(label))

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!(f::Function, op::DiracOp) = (filter!(f, dict(op)); return op)
    Base.filter!(f::Function, opc::DualOp) = (filter!((k,v)->f(reverse(k),v'), opc.op); return opc)

    Base.filter{S}(f::Function, op::DiracOp{S}) = DiracOp{S}(filter(f, dict(op)))
    Base.filter(f::Function, opc::DualOp) = DualOp(filter((k,v)->f(reverse(k),v'), opc.op))

    Base.map{S}(f::Function, op::DiracOp{S}) = DiracOp{S}(mapkv(f, dict(op)))
    Base.map(f::Function, opc::DualOp) = mapkv!((k,v)->f(reverse(k),v'), similar(opc), opc.op)
    
    mapcoeffs{S}(f::Function, op::DiracOp{S}) = DiracOp{S}(mapvals(f, dict(op)))
    mapcoeffs(f::Function, opc::DualOp) = mapvals!(v->f(v'), similar(opc), opc.op)

    maplabels{S}(f::Function, op::DiracOp{S}) = DiracOp{S}(mapkeys(f, dict(op)))
    maplabels(f::Function, opc::DualOp) = mapkeys!(k->f(reverse(k)), similar(opc), opc.op)

##########################
# Mathematical Functions #
##########################
    eager_ctran(op::DiracOp) = map((k,v)->(reverse(k),v'), op)
    
    Base.ctranspose(op::DiracOp) = DualOp(op)
    Base.ctranspose(opc::DualOp) = opc.op

    function inner{A,B}(bra::Bra{A}, op::DiracOp{B})
        result = Bra{typejoin(A,B)}()
        for ((ok,ob),oc) in dict(op)
            if !haskey(result, ob)
                result[ob] = 0
            end
            coeff = 0
            for (label,c) in dict(bra)
                coeff += c'*oc*inner_eval(A,B,label,ok) 
            end
            result[ob] += coeff
        end
        return result
    end

    function inner{A,B}(op::DiracOp{A}, ket::Ket{B})
        result = Ket{typejoin(A,B)}()
        for ((ok,ob),oc) in dict(op)
            if !haskey(result, ok)
                result[ok] = 0
            end
            coeff = 0
            for (label,v) in dict(ket)
                coeff += oc*v*inner_eval(A,B,ob,label) 
            end
            result[ok] += coeff
        end
        return result
    end

    function inner{A,B}(a::DiracOp{A}, b::DiracOp{B})
        result = DiracOp{typejoin(A,B)}()
        for ((ak,ab),ac) in dict(a)
            for ((bk,bb),bc) in dict(b)
                if !haskey(result, (ak,bb))
                    result[ak,bb] = 0
                end
                result[ak,bb] += ac*bc*inner_eval(A,B,ab,bk) 
            end
        end
        return result
    end

    function inner{A,B}(a::DiracOp{A}, b::DualOp{B})
        result = DiracOp{typejoin(A,B)}()
        for ((ak,ab),ac) in dict(a)
            for ((bb,bk),bc) in dict(b)
                if !haskey(result, (ak,bb))
                    result[ak,bb] = 0
                end
                result[ak,bb] += ac*bc'*inner_eval(A,B,ab,bk) 
            end
        end
        return result
    end

    function inner{A,B}(a::DualOp{A}, b::DiracOp{B})
        result = DiracOp{typejoin(A,B)}()
        for ((ab,ak),ac) in dict(a)
            for ((bk,bb),bc) in dict(b)
                if !haskey(result, (ak,bb))
                    result[ak,bb] = 0
                end
                result[ak,bb] += ac'*bc*inner_eval(A,B,ab,bk) 
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

    Base.(:+){S}(a::DiracOp{S}, b::DiracOp{S}) = DiracOp{S}(mergef(+, dict(a), dict(b)))
    Base.(:+){S}(a::DualOp{S}, b::DualOp{S}) = DualOp(a.op + b.op)
    Base.(:+){S}(a::AbstractOperator{S}, b::AbstractOperator{S}) = +(promote(a,b)...)

    Base.(:-){S}(a::AbstractOperator{S}, b::AbstractOperator{S}) = a + (-b)
    Base.(:-)(op::DiracOp) = mapcoeffs(-, op)
    Base.(:-)(opc::DualOp) = DualOp(-opc.op)

    Base.norm(op::DiracOp) = sqrt(sum(v->v^2, values(dict(op))))
    Base.norm(opc::DualOp) = norm(opc.op)
    
    QuBase.normalize(op::AbstractOperator) = scale(1/norm(op), op)
    QuBase.normalize!(op::AbstractOperator) = scale!(1/norm(op), op)

    Base.trace(op::DiracOp) = sum(k->op[k], filter(k->k[1]==k[2], keys(dict(op))))
    Base.trace(opc::DualOp) = trace(opc.op)'

    QuBase.tensor{S}(ops::DiracOp{S}...) = DiracOp{S}(mergecart!(tensor_op, OpDict(), map(dict, ops)))
    QuBase.tensor{S}(ket::Ket{S}, op::DiracOp{S}) = DiracOp{S}(mergecart!(tensor_ket_to_op, OpDict(), dict(ket), dict(op)))
    QuBase.tensor{S}(op::DiracOp{S}, ket::Ket{S}) = DiracOp{S}(mergecart!(tensor_op_to_ket, OpDict(), dict(op), dict(ket)))
    QuBase.tensor{S}(op::DiracOp{S}, bra::Bra{S}) = DiracOp{S}(mergecart!(tensor_bra_to_op, OpDict(), dict(op), mapvals(ctranspose, dict(bra))))
    QuBase.tensor{S}(bra::Bra{S}, op::DiracOp{S}) = DiracOp{S}(mergecart!(tensor_op_to_bra, OpDict(), mapvals(ctranspose, dict(bra)), dict(op)))
    QuBase.tensor(ops::AbstractOperator...) = tensor(promote(ops...)...)
    QuBase.tensor(ket::Ket, opc::DualOp) = tensor(ket', opc.op)'
    QuBase.tensor(opc::DualOp, ket::Ket) = tensor(opc.op, ket')'
    QuBase.tensor(opc::DualOp, bra::Bra) = tensor(opc.op, bra')'
    QuBase.tensor(bra::Bra, opc::DualOp) = tensor(bra', opc.op)'
    QuBase.tensor(a::DualOp, b::DualOp) = tensor(a.op, b.op)'

    xsubspace(op::GenericOperator, x) = filter((k,v)->sum(k[1])==x && sum(k[2])==x, op)

    filternz!(op::GenericOperator) = filter!((k, v) -> v != 0, op)
    filternz(op::GenericOperator) = filter((k, v) -> v != 0, op)

#################
# Partial trace #
#################
    ptrace(opc::DualOp, over::Integer...) = DualOp(ptrace(opc.op, over...))

    function ptrace{S<:Orthonormal}(op::DiracOp{S}, over::Integer...)
        result = dict(op)
        for i in over
            result = ptrace_dict(result, i)
        end
        return DiracOp{S}(result)
    end

    function ptrace_dict(d, over::Integer)
        result = OpDict()
        for (k, v) in d
            if k[1][over]==k[2][over]
                new_label = (except(k[1], over), except(k[2], over))
                if haskey(result, new_label)
                    result[new_label] += v
                else
                    result[new_label] = v
                end
            end
        end
        return result
    end

######################
# Printing Functions #
######################
    opstr(::DiracOp, k, v, pad) = "$pad$v $(statestr(k[1],Ket))$(statestr(k[2],Bra))"
    opstr(::DualOp, k, v, pad) = "$pad$(v') $(statestr(k[2],Ket))$(statestr(k[1],Bra))"

    function Base.show(io::IO, op::GenericOperator)
        print(io, "$(summary(op)) with $(length(op)) operator(s):")
        pad = "  "
        maxlen = 30
        i = 1
        for (k,v) in dict(op)
            if i <= maxlen
                println(io)
                print(io, opstr(op, k, v, pad))
                i = i + 1
            else  
                println(io)
                print(io, "$pad$vdots")
                break
            end
        end
    end

####################
# Helper Functions #
####################
    function tensor_op(pairs)
        #pairs structure is: (((op1ketlabel, op1bralabel), op1value), ((op2ketlabel, op2bralabel), op2value)...,)
        k = map(first, pairs)
        return ((join_tup(map(first, k)), join_tup(map(second, k))), prod(second, pairs))
    end

    function tensor_ket_to_op(pairs)
        #pairs structure is: ((ketlabel, ketvalue), ((opketlabel, opbralabel), opvalue))
        return ((join_tup(pairs[1][1], pairs[2][1][1]), pairs[2][1][2]), prod(second, pairs))
    end

    function tensor_op_to_ket(pairs)
        #pairs structure is: (((opketlabel, opbralabel), (ketlabel, ketvalue), opvalue))
        return ((join_tup(pairs[1][1][1], pairs[2][1]), pairs[1][1][2]), prod(second, pairs))
    end

    function tensor_bra_to_op(pairs)
        #pairs structure is: (((opketlabel, opbralabel), opvalue), (bralabel, bravalue))
        return ((pairs[1][1][1], join_tup(pairs[1][1][2], pairs[2][1])), prod(second, pairs))
    end

    function tensor_op_to_bra(pairs)
        #pairs structure is: ((bralabel, bravalue), ((opketlabel, opbralabel), opvalue))
        return ((pairs[2][1][1], join_tup(pairs[1][1], pairs[2][1][2])), prod(second, pairs))
    end

export ptrace,
    xsubspace,
    mapcoeffs,
    maplabels,
    filternz,
    filternz!
