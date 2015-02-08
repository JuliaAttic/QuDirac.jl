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
    function DiracOp{S}(ket::DiracKet{S}, bra::DiracBra{S})
        coeffs = ObjectIdDict()
        for (k,kc) in ket
            for (b,bc) in bra
                coeffs[k,b] = kc * bc
            end
        end
        return DiracOp{S}(coeffs)
    end

    copy_type{S}(::DiracOp{S}, coeffs) = DiracState{S}(coeffs)

    Base.copy(o::DiracOp) = copy_type(o, copy(o.coeffs))
    Base.similar(o::DiracOp) = copy_type(o, ObjectIdDict())

#######################
# Dict-Like Functions #
#######################
    Base.length(o::DiracOp) = length(o.coeffs)

    Base.getindex(o::DiracOp, labels::(Tuple,Tuple)) = o.coeffs[labels]
    Base.getindex(o::DiracOp, labels...) = o[labels]

    Base.setindex!(o::DiracOp, c, labels::(Tuple,Tuple)) = setindex!(o.coeffs, c, labels)
    Base.setindex!(o::DiracOp, c, labels...) = setindex!(o, c, labels)

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

##########################
# Mathematical Functions #
##########################
    function inner{S}(s::DiracBra{S}, o::DiracOp{S})
        result = DiracBra{S}()
        for ((ok,ob),oc) in o            
            if !haskey(result, ob)
                result[ob] = 0
            end
            coeff = 0
            for (sb,sc) in s
                coeff += sc*oc*inner(S,sb,ok) 
            end
            result[ob] += coeff
        end
        return result
    end

    function inner{S}(o::DiracOp{S}, s::DiracKet{S})
        result = DiracKet{S}()
        for ((ok,ob),oc) in o            
            if !haskey(result, ok)
                result[ok] = 0
            end
            coeff = 0
            for (sk,sc) in s
                coeff += oc*sc*inner(S,ob,sk) 
            end
            result[ok] += coeff
        end
        return result
    end

    function inner{S}(a::DiracOp{S}, b::DiracOp{S})
        result = DiracOp{S}()
        for ((ak,ab),ac) in a
            for ((bk,bb),bc) in b
                if !haskey(result, (ak,bb))
                    result[ak,bb] = 0
                end
                result[ak,bb] += ac*bc*inner(S,ab,bk) 
            end
        end
        return op
    end

    *{S}(a::DiracOp{S}, b::DiracOp{S}) = inner(a,b)

    *{S}(b::DiracBra{S}, o::DiracOp{S}) = inner(b,o)
    *{S}(o::DiracOp{S}, k::DiracKet{S}) = inner(o,k)

    *{S}(k::DiracKet{S}, o::DiracOp{S}) = tensor(k,o)
    *{S}(o::DiracOp{S}, b::DiracBra{S}) = tensor(o,b)

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

    normalize(o::DiracOp) = (1/norm(o))*o

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
        k = map(first, pairs)
        return ((tensor_tup(map(first, k)), tensor_tup(map(second, k))), prod(second, pairs))
    end

    function mergelabels{S}(f::Function, a::DiracOp{S}, b::DiracOp{S})
        return DiracOp{S}(mergef(f, a.coeffs, b.coeffs))
    end

    mergekets(a::DiracOp, b::DiracOp) = merge(a.ketlabels, b.ketlabels)
    mergebras(a::DiracOp, b::DiracOp) = merge(a.bralabels, b.bralabels)

    mapkv(f::Function, o::DiracOp) = copy_type(o, mapkv(f, o.coeffs))
    mapvals(f::Function, o::DiracOp) = copy_type(o, mapvals(f, o.coeffs))
    mapkeys(f::Function, o::DiracOp) = copy_type(o, mapkeys(f, o.coeffs))

export DiracOp
