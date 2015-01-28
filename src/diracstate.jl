import Base: 
    eltype,
    copy,
    +,
    *,
    -,
    show

abstract BasisType
abstract Orthonormal <: BasisType

immutable DiracState{D,B,N,T}
    labels::OrderedDict{StateLabel{N},T}   
end

function DiracState{D,N,T,B}(labels::OrderedDict{StateLabel{N},T}, ::Type{D}=Ket, ::Type{B}=Orthonormal) 
    return DiracState{D,B,N,T}(labels)
end

ket(label...) = DiracState(singlet_dict(StateLabel(label...), 1), Ket)
bra(label...) = DiracState(singlet_dict(StateLabel(label...), 1), Bra)

labels(ds::DiracState) = ds.labels
nfactors{D,B,N}(ds::DiracState{D,B,N}) = N
basistype{D,B}(ds::DiracState{D,B}) = B
dualtype{D}(ds::DiracState{D}) = D

eltype{D,B,N,T}(::DiracState{D,B,N,T}) = T

copy{D,B}(ds::DiracState{D,B}) = DiracState(copy(labels(ds)), D, B)

+{D,B,N}(a::DiracState{D,B,N}, b::DiracState{D,B,N}) = DiracState(mergef(+, labels(a), labels(b)), D, B)
-{D,B,N}(a::DiracState{D,B,N}, b::DiracState{D,B,N}) = a + (-b)
-(ds::DiracState) = DiracState(mapvals((k,v)->-v, labels(ds)), dualtype(ds), basistype(ds))

# *(a::DiracState, b::DiracState) = error("mismatched mutliplication")
*{D,B}(a::DiracState{D,B}, b::DiracState{D,B}) = DiracState(mergecart(tensor_reduce, labels(a), labels(b)), D, B)
*(c, ds::DiracState) = DiracState(castvals(*, labels(ds), c), dualtype(ds), basistype(ds))
*(ds::DiracState, c) = DiracState(castvals(*, c, labels(ds)), dualtype(ds), basistype(ds))

######################
# Printing Functions #
######################
statestr(c, label, ::Type{Ket}) = "$c | $(labelstr(label)) $rang"
statestr(c, label, ::Type{Bra}) = "$c $lang $(labelstr(label)) |"

function show{D}(io::IO, ds::DiracState{D})
    print(io, "$(summary(ds)):")
    pad = "  "
    maxlen = 30
    i = 1
    for (k,v) in labels(ds)
        if i <= maxlen
            println(io)
            print(io, "$pad$(statestr(v,k,D))")
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
assoc_nfactors{N}(::Associative{StateLabel{N}}) = N

function singlet_dict{N,T}(label::StateLabel{N}, c::T)
    od = OrderedDict(StateLabel{N}, T)
    od[label] = c
    return od
end

function value_type(d::Associative, others::Associative...)
    V = eltype(d)[2]
    for other in others
        V = promote_type(V, eltype(other)[2])
    end
    return V
end

function tensor_reduce(pairs)
    return reduce((kv1,kv2) -> (tensor(kv1[1],kv2[1]), kron(kv1[2],kv2[2])), pairs)
end

function mergecart(f::Function, d...)
    result = OrderedDict(StateLabel{sum(map(assoc_nfactors, d))}, value_type(d...))
    for pairs in product(d...)
        merged_pair = f(pairs)
        result[merged_pair[1]] = merged_pair[2]
    end
    return result
end

function mergef!(f::Function, d, others...)
    for other in others
        for (k,v) in other
            if haskey(d, k)
                d[k] = f(d[k], v)
            else   
                d[k] = v
            end
        end
    end
    return d
end

function mergef{K}(f::Function, d::Associative{K}, others::Associative{K}...)
    return mergef!(f, OrderedDict(K,value_type(d, others...)), d, others...)
end


function castvals{K,V}(f::Function, a::Associative{K,V}, b::Associative{K,V})
    return mergef(f, a, b)
end

function castvals{K,V,T}(f::Function, d::Associative{K,V}, c::T)
    return mapvals!((k,v)->f(v,c), d, OrderedDict(K, promote_type(V,T)))
end

function castvals{T,K,V}(f::Function, c::T, d::Associative{K,V})
    return mapvals!((k,v)->f(c,v), d, OrderedDict(K, promote_type(V,T)))
end

function mapvals!(f::Function, d)
    for (k,v) in d
        d[k] = f(k,v)
    end
    return d
end

function mapvals!(f::Function, d, result)
    for (k,v) in d
        result[k] = f(k,v)
    end
    return result
end

mapvals(f::Function, d) = mapvals!(f, copy(d))

export DiracState,
    ket,
    bra