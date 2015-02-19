###########
# DiracOp #
###########
    type DiracOp{S} <: AbstractOperator{S}
        coeffs::ObjectIdDict
        DiracOp() = new(ObjectIdDict())
        DiracOp(coeffs) = new(coeffs)
    end

################
# Constructors #
################
    function DiracOp{S}(f::Function, ket::DiracKet{S})
        return DiracOp(f, S, keys(ket))
    end

    # f(label) -> (eigval, eiglabel)
    function DiracOp{S}(f::Function, ::Type{S}, labels)
        coeffs = ObjectIdDict()
        for i in labels
            for j in labels 
                (c, new_j) = f(j)
                coeffs[i,j] = c * inner_rule(S, i, new_j)
            end
        end
        return DiracOp{S}(coeffs)
    end

    function DiracOp{A,B}(ket::DiracKet{A}, bra::DiracBra{B})
        coeffs = ObjectIdDict()
        for (k,kc) in ket
            for (b,bc) in bra
                coeffs[k,b] = kc * bc
            end
        end
        return DiracOp{typejoin(A,B)}(coeffs)
    end

    copy_type{S}(::DiracOp{S}, coeffs) = DiracOp{S}(coeffs)

    Base.copy(o::DiracOp) = copy_type(o, copy(o.coeffs))
    Base.similar(o::DiracOp) = copy_type(o, ObjectIdDict())

#######################
# Dict-Like Functions #
#######################
    Base.length(o::DiracOp) = length(o.coeffs)

    Base.getindex(o::DiracOp, labels::(Tuple,Tuple)) = o.coeffs[labels]
    Base.getindex(o::DiracOp, labels...) = o[labels]

    Base.setindex!(o::DiracOp, c, labels::(Tuple,Tuple)) = setindex!(o.coeffs, c, labels)
    Base.setindex!(o::DiracOp, c, k::Tuple, b::Tuple) = setindex!(o, c, (k,b))

    Base.start(o::DiracOp) = start(o.coeffs)
    Base.done(o::DiracOp, state) = done(o.coeffs, state)
    Base.next(o::DiracOp, state) = next(o.coeffs, state)
    Base.endof(o::DiracOp) = endof(o.coeffs)
    Base.last(o::DiracOp) = last(o.coeffs)
    Base.first(o::DiracOp) = first(o.coeffs)
    Base.collect(o::DiracOp) = collect(o.coeffs)

    Base.keys(o::DiracOp) = keys(o.coeffs)
    Base.values(o::DiracOp) = values(o.coeffs)

    Base.get(o::DiracOp, label, default) = get(o.coeffs, label, default)
    Base.haskey(o::DiracOp, label) = haskey(o.coeffs, label)
    Base.filter!(f::Function, o::DiracOp) = (filter!(f, o.coeffs); return o)
    Base.filter(f::Function, o::DiracOp) = copy_type(o, filter(f, o.coeffs))
    Base.map(f::Function, o::DiracOp) = mapkv(f, o)
    Base.delete!(o::DiracOp, label) = copy_type(o, delete!(o.coeffs, label))

    filternz!(o::DiracOp) = filter!((k, v) -> v != 0, o)
    filternz(o::DiracOp) = filter((k, v) -> v != 0, o)


##########################
# Mathematical Functions #
##########################
    function inner{A,B}(s::DiracBra{A}, o::DiracOp{B})
        result = DiracBra{typejoin(A,B)}()
        for ((ok,ob),oc) in o            
            if !haskey(result, ob)
                result[ob] = 0
            end
            coeff = 0
            for (sb,sc) in s
                coeff += sc*oc*inner_eval(A,B,sb,ok) 
            end
            result[ob] += coeff
        end
        return result
    end

    function inner{A,B}(o::DiracOp{A}, s::DiracKet{B})
        result = DiracKet{typejoin(A,B)}()
        for ((ok,ob),oc) in o            
            if !haskey(result, ok)
                result[ok] = 0
            end
            coeff = 0
            for (sk,sc) in s
                coeff += oc*sc*inner_eval(A,B,ob,sk) 
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

    mult_state_op(b::DiracBra, o::DiracOp) = inner(b,o)
    mult_state_op(o::DiracOp, k::DiracKet) = inner(o,k)
    mult_state_op{S}(k::DiracKet{S}, o::DiracOp{S}) = tensor(k,o)
    mult_state_op{S}(o::DiracOp{S}, b::DiracBra{S}) = tensor(o,b)

    *(s::DiracState, o::DiracOp) = mult_state_op(s, o)
    *(o::DiracOp, s::DiracState) = mult_state_op(o, s)
    *(a::DiracOp, b::DiracOp) = inner(a,b)
    *(k::DiracKet, b::DiracBra) = DiracOp(k,b)
    *(c, o::DiracOp) = copy_type(o, castvals(*, c, o.coeffs))
    *(o::DiracOp, c) = copy_type(o, castvals(*, o.coeffs, c))

    +{S}(a::DiracOp{S}, b::DiracOp{S}) = mergelabels(+, a, b)
    -{S}(a::DiracOp{S}, b::DiracOp{S}) = a + (-b)
    -(o::DiracOp) = copy_type(o, mapvals(-, o.coeffs))

    /(o::DiracOp, c) = copy_type(o, castvals(/, o.coeffs, c))

    Base.conj(o::DiracOp) = mapvals(conj, o)
    Base.ctranspose(o::DiracOp) = copy_type(o, mapkv((k,v)->(reverse(k), conj(v)), o.coeffs))

    Base.norm(o::DiracOp) = sqrt(sum(v->v^2, values(o)))
    Base.trace(o::DiracOp) = sum(k->o[k], filter(k->k[1]==k[2], keys(o)))

    QuBase.tensor{S}(ops::DiracOp{S}...) = DiracOp{S}(mergecart!(tensor_op, ObjectIdDict(), ops))
    QuBase.tensor{S}(k::DiracKet{S}, o::DiracOp{S}) = DiracOp{S}(mergecart!(tensor_ket_to_op, ObjectIdDict(), k, o))
    QuBase.tensor{S}(o::DiracOp{S}, k::DiracKet{S}) = DiracOp{S}(mergecart!(tensor_op_to_ket, ObjectIdDict(), o, k))
    QuBase.tensor{S}(o::DiracOp{S}, b::DiracBra{S}) = DiracOp{S}(mergecart!(tensor_bra_to_op, ObjectIdDict(), o, b))
    QuBase.tensor{S}(b::DiracBra{S}, o::DiracOp{S}) = DiracOp{S}(mergecart!(tensor_op_to_bra, ObjectIdDict(), b, o))

    normalize(o::DiracOp) = (1/norm(o))*o

    function ptrace{S<:Orthonormal}(o::DiracOp{S}, over::Integer...)
        result = o.coeffs
        for i in over
            result = ptrace_coeffs(result, i)
        end
        return DiracOp{S}(result)
    end

    xsubspace(o::DiracOp, x, i=1) = filter((k,v)->sum(k[i])==x, o)
    matchlabel_at(o::DiracOp, x, y, i) = filter((k,v)-> x==k[i][y], o)
    matchlabel_at(o::DiracOp, x, y) = filter((k,v)-> x==k[1][y] && x==k[2][y], o)
    matchlabel_in(o::DiracOp, x, i) = filter((k,v)-> x in k[i], o)
    matchlabel_in(o::DiracOp, x) = filter((k,v)-> x in k[1] && x in k[2], o)

######################
# Printing Functions #
######################
    function Base.show(io::IO, o::DiracOp)
        print(io, "$(summary(o)) with $(length(o)) operator(s):")
        pad = "  "
        maxlen = 30
        i = 1
        for (k,v) in o
            if i <= maxlen
                println(io)
                print(io, "$pad$v $(statestr(first(k),Ket))$(statestr(second(k),Bra))")
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
    function ptrace_coeffs(coeffs, over::Integer)
        result = ObjectIdDict()
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

    mergekets(a::DiracOp, b::DiracOp) = merge(a.ketlabels, b.ketlabels)
    mergebras(a::DiracOp, b::DiracOp) = merge(a.bralabels, b.bralabels)

    mapkv(f::Function, o::DiracOp) = copy_type(o, mapkv(f, o.coeffs))
    mapvals(f::Function, o::DiracOp) = copy_type(o, mapvals(f, o.coeffs))
    mapkeys(f::Function, o::DiracOp) = copy_type(o, mapkeys(f, o.coeffs))

export DiracOp,
    ptrace,
    normalize,
    xsubspace,
    matchlabel_at,
    matchlabel_in
