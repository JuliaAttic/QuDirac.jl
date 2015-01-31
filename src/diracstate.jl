import Base: 
    eltype,
    copy,
    norm,
    ctranspose,
    conj,
    +,
    *,
    -,
    getindex,
    show,
    start,
    done,
    next,
    endof,
    last,
    first,
    collect,
    values,
    filter,
    filter!,
    haskey

abstract Orthogonal <: AbstractStructure
abstract Orthonormal <: Orthogonal

type DiracState{D,B,N,T}
    labels::Dict{StateLabel{N},T}   
end

typealias DiracKet{B,N,T} DiracState{Ket,B,N,T}
typealias DiracBra{B,N,T} DiracState{Bra,B,N,T}

function DiracState{D,N,T,B}(labels::Dict{StateLabel{N},T}, ::Type{D}, ::Type{B}) 
    return DiracState{D,B,N,T}(labels)
end

ket{B}(::Type{B}, label...) = DiracState(singlet_dict(StateLabel(label...), 1), Ket, B)
ket(label...) = ket(Orthonormal, label...)
bra{B}(::Type{B}, label...) = DiracState(singlet_dict(StateLabel(label...), 1), Bra, B)
bra(label...) = bra(Orthonormal, label...)

labels(ds::DiracState) = ds.labels
nfactors{D,B,N}(ds::DiracState{D,B,N}) = N
basistype{D,B}(ds::DiracState{D,B}) = B
dualtype{D}(ds::DiracState{D}) = D

eltype{D,B,N,T}(::DiracState{D,B,N,T}) = T

copy{D,B}(ds::DiracState{D,B}) = DiracState(copy(labels(ds)), D, B)

#######################
# Dict-Like Functions #
#######################

length(ds::DiracState) = length(labels(ds))

getindex(ds::DiracState, label::StateLabel) = ds.labels[label]
getindex(ds::DiracState, i...) = ds.labels[StateLabel(i...)]
getstate(ds::DiracState, i...) = ds[i...] * ket(i...)

start(ds::DiracState) = start(labels(ds))
done(ds::DiracState, state) = done(labels(ds), state)
next(ds::DiracState, state) = next(labels(ds), state)
endof(ds::DiracState) = endof(labels(ds))
last(ds::DiracState) = last(labels(ds))
first(ds::DiracState) = first(labels(ds))
collect{D,B}(ds::DiracState{D,B}) = DiracState(collect(labels(ds)), D, B)

haskey(ds::DiracState, label) = haskey(labels(ds), label)
values(ds::DiracState) = values(labels(ds))
filter!(f::Function, ds::DiracState) = filter!(f, labels(ds))
filter(f::Function, ds::DiracState) = filter!(f, copy(ds))

##########################
# Mathematical Functions #
##########################
function inner{B}(db::DiracBra{B}, dk::DiracKet{B})
    result = 0
    for (b,c) in labels(db)
        for (k,v) in labels(dk)
            result += c*v*inner(B,b,k) 
        end
    end
    return result
end

function inner{B<:Orthogonal}(db::DiracBra{B}, dk::DiracKet{B})
    if length(db) > length(dk)
        return ortho_inner(labels(dk), labels(db))
    else
        return ortho_inner(labels(db), labels(dk))
    end
end

+{D,B,N}(a::DiracState{D,B,N}, b::DiracState{D,B,N}) = DiracState(mergef(+, labels(a), labels(b)), D, B)
-{D,B,N}(a::DiracState{D,B,N}, b::DiracState{D,B,N}) = a + (-b)
-(ds::DiracState) = DiracState(mapvals(-, labels(ds)), dualtype(ds), basistype(ds))

*(a::DiracBra, b::DiracKet) = inner(a,b)
*{D,B}(a::DiracState{D,B}, b::DiracState{D,B}) = DiracState(mergecart(tensor_reduce, labels(a), labels(b)), D, B)
*(c, ds::DiracState) = DiracState(castvals(*, labels(ds), c), dualtype(ds), basistype(ds))
*(ds::DiracState, c) = DiracState(castvals(*, c, labels(ds)), dualtype(ds), basistype(ds))

conj(ds::DiracState) = mapvals(conj, ds)
ctranspose{D,B}(ds::DiracState{D,B}) = DiracState(labels(conj(ds)), D', B)

norm(ds::DiracState) = sqrt(sum(v->v^2, values(ds)))
normalize(ds::DiracState) = (1/norm(ds))*ds

######################
# Printing Functions #
######################
statestr(label, ::Type{Ket}) = "| $(labelstr(label)) $rang"
statestr(label, ::Type{Bra}) = "$lang $(labelstr(label)) |"

function show{D}(io::IO, ds::DiracState{D})
    print(io, "$(summary(ds)) with $(length(ds)) state(s):")
    pad = "  "
    maxlen = 30
    i = 1
    for (k,v) in labels(ds)
        if i <= maxlen
            println(io)
            print(io, "$pad$v $(statestr(k,D))")
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
function ortho_inner(a, b)
    result = 0
    for (label,c) in a
        if haskey(b, label)
            result += b[label]
        end
    end
    return result
end

assoc_nfactors{N}(::Associative{StateLabel{N}}) = N

function singlet_dict{N,T}(label::StateLabel{N}, c::T)
    od = Dict{StateLabel{N}, T}()
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
    result = Dict{StateLabel{sum(map(assoc_nfactors, d))}, value_type(d...)}()
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
    return mergef!(f, Dict{K,value_type(d, others...)}(), d, others...)
end


function castvals{K,V}(f::Function, a::Associative{K,V}, b::Associative{K,V})
    return mergef(f, a, b)
end

function castvals{K,V,T}(f::Function, d::Associative{K,V}, c::T)
    return mapvals!(v->f(v,c), d, Dict{K, promote_type(V,T)}())
end

function castvals{T,K,V}(f::Function, c::T, d::Associative{K,V})
    return mapvals!(v->f(c,v), d, Dict{K, promote_type(V,T)}())
end

function mapkv!(f::Function, d, result)
    for (k,v) in d
        (k0,v0) = f(k,v)
        delete!(d,k)
        result[k0] = v0
    end
    return result
end

function mapvals!(f::Function, d, result)
    for (k,v) in d
        result[k] = f(v)
    end
    return result
end

mapvals(f::Function, d) = mapvals!(f, d, copy(d))
mapvals{D,B}(f::Function, ds::DiracState{D,B}) = DiracState(mapvals(f, labels(ds)), D, B)

mapkv(f::Function, d) = mapkv!(f, d, copy(d))
mapkv{D,B}(f::Function, ds::DiracState{D,B}) = DiracState(mapkv(f, labels(ds)), D, B)

export DiracState,
    ket,
    bra,
    getstate,
    normalize