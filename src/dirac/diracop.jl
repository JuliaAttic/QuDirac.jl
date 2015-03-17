##################
# DiracOp/DualOp #
##################
    abstract GenericOperator{S} <: AbstractOperator{S}

    typealias OpCoeffs Dict{(Tuple,Tuple),Number}

    type DiracOp{S} <: GenericOperator{S}
        coeffs::OpCoeffs
        DiracOp() = new(OpCoeffs())
        DiracOp(coeffs) = new(coeffs)
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
        return DiracOp(f, S, keys(ket))
    end

    # f(label) -> (newval, newlabel)
    function DiracOp{S}(f::Function, ::Type{S}, labels)
        coeffs = OpCoeffs()
        for i in labels
            for j in labels 
                (c, new_j) = f(j)
                coeffs[i,j] = c * inner_rule(S, i, new_j)
            end
        end
        return DiracOp{S}(coeffs)
    end

    function DiracOp{A,B}(ket::Ket{A}, bra::Bra{B})
        coeffs = OpCoeffs()
        for (k,kc) in ket
            for (b,bc) in bra.ket
                coeffs[k,b] = kc * bc'
            end
        end
        return DiracOp{typejoin(A,B)}(coeffs)
    end

    coeffs(op::DiracOp) = op.coeffs
    coeffs(opc::DualOp) = coeffs(opc.op)

    Base.copy(op::GenericOperator) = typeof(op)(copy(coeffs(op)))
    Base.similar(op::GenericOperator) = typeof(op)(similar(coeffs(op)))

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){S}(a::DiracOp{S}, b::DiracOp{S}) = coeffs(a) == coeffs(b)
    Base.(:(==)){S}(a::DualOp{S}, b::DualOp{S}) = a.op == b.op
    Base.(:(==)){S}(a::AbstractOperator{S}, b::AbstractOperator{S}) = ==(promote(a,b)...)

    Base.hash(op::GenericOperator) = hash(coeffs(op), hash(typeof(op)))

    Base.length(op::GenericOperator) = length(coeffs(op))

    Base.getindex(op::DiracOp, label::(Tuple,Tuple)) = coeffs(op)[label]
    Base.getindex(op::DiracOp, k::Tuple, b::Tuple) = op[(k,b)]
    Base.getindex(opc::DualOp, label::(Tuple,Tuple)) = opc.op[reverse(label)]'
    Base.getindex(opc::DualOp, k::Tuple, b::Tuple) = opc.op[(b,k)]'
    Base.getindex(op::GenericOperator, k, b) = op[tuple(k),tuple(b)]

    Base.setindex!(op::DiracOp, c, label::(Tuple,Tuple)) = setindex!(coeffs(op), c, label)
    Base.setindex!(op::DiracOp, c, k::Tuple, b::Tuple) = setindex!(op, c, (k,b))
    Base.setindex!(opc::DualOp, c, label::(Tuple,Tuple)) = setindex!(opc.op, c', reverse(label))
    Base.setindex!(opc::DualOp, c, k::Tuple, b::Tuple) = setindex!(opc.op, c', (b,k))
    Base.setindex!(op::GenericOperator, c, k, b) = setindex!(op, c, tuple(k), tuple(b))

    Base.keys(op::DiracOp) = keys(coeffs(op))
    Base.values(op::DiracOp) = values(coeffs(op))

    Base.start(op::DiracOp) = start(coeffs(op))
    Base.done(op::DiracOp, state) = done(coeffs(op), state)
    Base.next(op::DiracOp, state) = next(coeffs(op), state)

    Base.haskey(op::DiracOp, label::(Tuple,Tuple)) = haskey(coeffs(op), label)
    Base.haskey(opc::DualOp, label::(Tuple,Tuple)) = haskey(opc.op, reverse(label))

    Base.get(op::DiracOp, label::(Tuple,Tuple), default) = get(coeffs(op), label, default)
    Base.get(opc::DualOp, label::(Tuple,Tuple), default) = haskey(opc, label) ? opc[label] : default

    Base.delete!(op::DiracOp, label::(Tuple,Tuple)) = (delete!(coeffs(op), label); return op)
    Base.delete!(opc::DualOp, label::(Tuple,Tuple)) = delete!(opc.op, reverse(label))

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!(f::Function, op::DiracOp) = (filter!(f, coeffs(op)); return op)
    Base.filter!(f::Function, opc::DualOp) = (filter!((k,v)->f(reverse(k),v'), opc.op); return opc)

    Base.filter{S}(f::Function, op::DiracOp{S}) = DiracOp{S}(filter(f, coeffs(op)))
    Base.filter(f::Function, opc::DualOp) = DualOp(filter((k,v)->f(reverse(k),v'), opc.op))

    Base.map{S}(f::Function, op::DiracOp{S}) = DiracOp{S}(mapkv(f,coeffs(op)))
    Base.map(f::Function, opc::DualOp) = mapkv!((k,v)->f(reverse(k),v'), similar(opc), opc.op)
    
    mapcoeffs{S}(f::Function, op::DiracOp{S}) = DiracOp{S}(mapvals(f, coeffs(op)))
    mapcoeffs(f::Function, opc::DualOp) = mapvals!(v->f(v'), similar(opc), opc.op)

    maplabels{S}(f::Function, op::DiracOp{S}) = DiracOp{S}(mapkeys(f, coeffs(op)))
    maplabels(f::Function, opc::DualOp) = mapkeys!(k->f(reverse(k)), similar(opc), opc.op)

##########################
# Mathematical Functions #
##########################
    eager_ctran(op::DiracOp) = map((k,v)->(reverse(k),v'), op)
    
    Base.ctranspose(op::DiracOp) = DualOp(op)
    Base.ctranspose(opc::DualOp) = opc.op

    function inner{A,B}(bra::Bra{A}, op::DiracOp{B})
        result = Bra{typejoin(A,B)}()
        for ((ok,ob),oc) in op            
            if !haskey(result, ob)
                result[ob] = 0
            end
            coeff = 0
            for (label,c) in bra.ket
                coeff += c'*oc*inner_eval(A,B,label,ok) 
            end
            result[ob] += coeff
        end
        return result
    end

    function inner{A,B}(op::DiracOp{A}, ket::Ket{B})
        result = Ket{typejoin(A,B)}()
        for ((ok,ob),oc) in op            
            if !haskey(result, ok)
                result[ok] = 0
            end
            coeff = 0
            for (label,v) in ket
                coeff += oc*v*inner_eval(A,B,ob,label) 
            end
            result[ok] += coeff
        end
        return result
    end

    function inner{A,B}(a::DiracOp{A}, b::DiracOp{B})
        result = DiracOp{typejoin(A,B)}()
        for ((ak,ab),ac) in a
            for ((bk,bb),bc) in b
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
        for ((ak,ab),ac) in a
            for ((bb,bk),bc) in b.op
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
        for ((ab,ak),ac) in a.op
            for ((bk,bb),bc) in b
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

    Base.scale!(c::Number, op::GenericOperator) = (castvals!(*, c, coeffs(op)); return op)
    Base.scale!(op::GenericOperator, c::Number) = (castvals!(*, coeffs(op), c); return op)

    Base.scale(c::Number, op::GenericOperator) = typeof(op)(castvals(*, c, coeffs(op)))
    Base.scale(op::GenericOperator, c::Number) = typeof(op)(castvals(*, coeffs(op), c))

    Base.(:*)(c::Number, op::AbstractOperator) = scale(c, op)
    Base.(:*)(op::AbstractOperator, c::Number) = scale(op, c)
    Base.(:/)(op::AbstractOperator, c::Number) = scale(op, 1/c)

    Base.(:+){S}(a::DiracOp{S}, b::DiracOp{S}) = mergelabels(+, a, b)
    Base.(:+){S}(a::DualOp{S}, b::DualOp{S}) = DualOp(a.op + b.op)
    Base.(:+){S}(a::AbstractOperator{S}, b::AbstractOperator{S}) = +(promote(a,b)...)

    Base.(:-){S}(a::AbstractOperator{S}, b::AbstractOperator{S}) = a + (-b)
    Base.(:-)(op::DiracOp) = mapcoeffs(-, op)
    Base.(:-)(opc::DualOp) = DualOp(-opc.op)

    Base.norm(op::DiracOp) = sqrt(sum(v->v^2, values(op)))
    Base.norm(opc::DualOp) = norm(opc.op)
    
    QuBase.normalize(op::AbstractOperator) = scale(1/norm(op), op)
    QuBase.normalize!(op::AbstractOperator) = scale!(1/norm(op), op)

    Base.trace(op::DiracOp) = sum(k->op[k], filter(k->k[1]==k[2], keys(op)))
    Base.trace(opc::DualOp) = trace(opc.op)'

    QuBase.tensor{S}(ops::DiracOp{S}...) = DiracOp{S}(mergecart!(tensor_op, OpCoeffs(), ops))
    QuBase.tensor{S}(ket::Ket{S}, op::DiracOp{S}) = DiracOp{S}(mergecart!(tensor_ket_to_op, OpCoeffs(), ket, op))
    QuBase.tensor{S}(op::DiracOp{S}, ket::Ket{S}) = DiracOp{S}(mergecart!(tensor_op_to_ket, OpCoeffs(), op, ket))
    QuBase.tensor{S}(op::DiracOp{S}, bra::Bra{S}) = DiracOp{S}(mergecart!(tensor_bra_to_op, OpCoeffs(), op, mapvals(ctranspose, coeffs(bra))))
    QuBase.tensor{S}(bra::Bra{S}, op::DiracOp{S}) = DiracOp{S}(mergecart!(tensor_op_to_bra, OpCoeffs(), mapvals(ctranspose, coeffs(bra)), op))
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
        result = coeffs(op)
        for i in over
            result = ptrace_coeffs(result, i)
        end
        return DiracOp{S}(result)
    end

    function ptrace_coeffs(coeffs, over::Integer)
        result = OpCoeffs()
        for tr_label in factor_labels(coeffs, over) # labels to be traced over
            for key in filter_at_labels(coeffs, tr_label, over)
                new_label = (except(key[1], over), except(key[2], over))
                if haskey(result, new_label)
                    result[new_label] += coeffs[key]
                else
                    result[new_label] = coeffs[key]
                end
            end
        end  
        return result
    end

    filter_at_labels(coeffs, tr_label, i) = filter(k -> tr_label==k[1][i] && tr_label==k[2][i], keys(coeffs))
    factor_labels(coeffs, factor) = distinct(imap(i->i[1][factor], keys(coeffs)))

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
        for (k,v) in coeffs(op)
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

    function mergelabels{S}(f::Function, a::DiracOp{S}, b::DiracOp{S})
        return DiracOp{S}(mergef(f, a.coeffs, b.coeffs))
    end

export ptrace,
    xsubspace,
    mapcoeffs,
    maplabels,
    filternz,
    filternz!
