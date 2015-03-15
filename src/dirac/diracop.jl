###########
# DiracOp #
###########
    typealias OpCoeffs Dict{(Tuple,Tuple),Number}

    type DiracOp{S} <: AbstractOperator{S}
        coeffs::OpCoeffs
        DiracOp() = new(OpCoeffs())
        DiracOp(coeffs) = new(coeffs)
    end

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
            for (b,bc) in bra
                coeffs[k,b] = kc * bc
            end
        end
        return DiracOp{typejoin(A,B)}(coeffs)
    end

    coeffs(o::DiracOp) = coeffs(o)

    Base.copy{S}(o::DiracOp{S}) = DiracOp{S}(copy(coeffs(o)))
    Base.similar{S}(o::DiracOp{S}) = DiracOp{S}(similar(coeffs(o)))

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){S}(a::DiracOp{S}, b::DiracOp{S}) = coeffs(a) == coeffs(b)
    Base.hash(o::DiracOp) = hash(coeffs(o), hash(typeof(o)))

    Base.length(o::DiracOp) = length(coeffs(o))

    Base.getindex(o::DiracOp, labels::(Tuple,Tuple)) = coeffs(o)[labels]
    Base.getindex(o::DiracOp, k::Tuple, b::Tuple) = o[(k,b)]

    Base.setindex!(o::DiracOp, c, labels::(Tuple,Tuple)) = setindex!(coeffs(o), c, labels)
    Base.setindex!(o::DiracOp, c, k::Tuple, b::Tuple) = setindex!(o, c, (k,b))

    Base.keys(o::DiracOp) = keys(coeffs(o))
    Base.values(o::DiracOp) = values(coeffs(o))

    Base.start(o::DiracOp) = start(coeffs(o))
    Base.done(o::DiracOp, state) = done(coeffs(o), state)
    Base.next(o::DiracOp, state) = next(coeffs(o), state)

    Base.get(o::DiracOp, label, default) = get(coeffs(o), label, default)
    Base.haskey(o::DiracOp, label) = haskey(coeffs(o), label)

    Base.delete!(o::DiracOp, label) = (delete!(coeffs(o), label); return o)

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!(f::Function, o::DiracOp) = (filter!(f, coeffs(o)); return o)
    Base.filter{S}(f::Function, o::DiracOp{S}) = DiracOp{S}(filter(f, coeffs(o)))
    Base.map(f::Function, o::DiracOp{S}) = DiracOp{S}(mapkv(f,o))
    
    mapcoeffs{S}(f::Function, o::DiracOp{S}) = DiracOp{S}(mapvals(f, coeffs(o)))
    maplabels{S}(f::Function, o::DiracOp{S}) = DiracOp{S}(mapkeys(f, coeffs(o)))

##########################
# Mathematical Functions #
##########################
    function inner{A,B}(bra::Bra{A}, o::DiracOp{B})
        result = Bra{typejoin(A,B)}()
        for ((ok,ob),oc) in o            
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

    function inner{A,B}(o::DiracOp{A}, ket::Ket{B})
        result = Ket{typejoin(A,B)}()
        for ((ok,ob),oc) in o            
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

    Base.(:*)(b::Bra, o::DiracOp) = inner(b,o)
    Base.(:*)(o::DiracOp, k::Ket) = inner(o,k)
    Base.(:*)(k::Ket, o::DiracOp) = tensor(k,o)
    Base.(:*)(o::DiracOp, b::Bra) = tensor(o,b)
    Base.(:*)(a::DiracOp, b::DiracOp) = inner(a,b)
    Base.(:*)(k::Ket, b::Bra) = DiracOp(k,b)

    Base.scale!(c::Number, o::DiracOp) = (castvals!(*, c, coeffs(o)); return o)
    Base.scale!(o::DiracOp, c::Number) = (castvals!(*, coeffs(o), c); return o)

    Base.scale{S}(c::Number, o::DiracOp{S}) = DiracOp{S}(castvals(*, c, coeffs(o)))
    Base.scale{S}(o::DiracOp{S}, c::Number) = DiracOp{S}(castvals(*, coeffs(o), c))

    Base.(:*)(c::Number, o::DiracOp) = scale(c, o)
    Base.(:*)(o::DiracOp, c::Number) = scale(o, c)
    Base.(:/)(o::DiracOp, c::Number) = scale(o, 1/c)

    Base.(:+){S}(a::DiracOp{S}, b::DiracOp{S}) = mergelabels(+, a, b)
    Base.(:-){S}(a::DiracOp{S}, b::DiracOp{S}) = a + (-b)
    Base.(:-){S}(o::DiracOp{S}) = mapcoeffs(-, o)

    Base.ctranspose{S}(o::DiracOp{S}) = map((k,v)->(reverse(k), v'), o)

    Base.norm(o::DiracOp) = sqrt(sum(v->v^2, values(o)))
    QuBase.normalize(o::DiracOp) = (1/norm(o))*o
    QuBase.normalize!(o::DiracOp) = scale!(1/norm(o), o)

    Base.trace(o::DiracOp) = sum(k->o[k], filter(k->k[1]==k[2], keys(o)))

    QuBase.tensor{S}(ops::DiracOp{S}...) = DiracOp{S}(mergecart!(tensor_op, OpCoeffs(), ops))
    QuBase.tensor{S}(k::Ket{S}, o::DiracOp{S}) = DiracOp{S}(mergecart!(tensor_ket_to_op, OpCoeffs(), k, o))
    QuBase.tensor{S}(o::DiracOp{S}, k::Ket{S}) = DiracOp{S}(mergecart!(tensor_op_to_ket, OpCoeffs(), o, k))
    QuBase.tensor{S}(o::DiracOp{S}, b::Bra{S}) = DiracOp{S}(mergecart!(tensor_bra_to_op, OpCoeffs(), o, mapvals(ctranspose, coeffs(b.ket)))
    QuBase.tensor{S}(b::Bra{S}, o::DiracOp{S}) = DiracOp{S}(mergecart!(tensor_op_to_bra, OpCoeffs(), mapvals(ctranspose, coeffs(b.ket)), o))

    xsubspace(o::DiracOp, x, i=1) = filter((k,v)->sum(k[i])==x, o)

    filternz!(o::DiracOp) = filter!((k, v) -> v != 0, o)
    filternz(o::DiracOp) = filter((k, v) -> v != 0, o)

#################
# Partial trace #
#################
    function ptrace{S<:Orthonormal}(o::DiracOp{S}, over::Integer...)
        result = coeffs(o)
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

export DiracOp,
    ptrace,
    xsubspace,
    mapcoeffs,
    maplabels,
    filternz,
    filternz!
