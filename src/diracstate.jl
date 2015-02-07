import Base: 
    eltype,
    copy,
    norm,
    ctranspose,
    conj,
    +,
    *,
    /,
    -,
    getindex,
    setindex!,
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

type DiracState{D,S,N}
    labels::ObjectIdDict
end

typealias DiracKet{S,N} DiracState{Ket,S,N}
typealias DiracBra{S,N} DiracState{Bra,S,N}

make_ket{S,N}(::Type{S}, label::NTuple{N}) = DiracKet{S,N}(singlet_dict(label, 1))
ket{S}(::Type{S}, label...) = make_ket(S, label)
ket(label...) = make_ket(Orthonormal, label)

make_bra{S,N}(::Type{S}, label::NTuple{N}) = DiracBra{S,N}(singlet_dict(label, 1))
bra{S}(::Type{S}, label...) = make_bra(S, label)
bra(label...) = make_bra(Orthonormal, label)

nfactors{D,B,N}(ds::DiracState{D,B,N}) = N
basistype{D,B}(ds::DiracState{D,B}) = B
dualtype{D}(ds::DiracState{D}) = D

copy_type{D,B,N}(::DiracState{D,B,N}, new_labels) = DiracState{D,B,N}(new_labels)
copy{D,B,N}(ds::DiracState{D,B,N}) = copy_type(ds, copy(ds.labels))

#######################
# Dict-Like Functions #
#######################

length(ds::DiracState) = length(ds.labels)

getindex{D,S,N}(ds::DiracState{D,S,N}, label::NTuple{N}) = ds.labels[label]
getindex{D,S,A,B}(ds::DiracState{D,S,A}, label::NTuple{B}) = error("length of input label must be $A; got label of length $B")
getindex(ds::DiracState, i...) = ds[i]

setindex!{D,S,N}(ds::DiracState{D,S,N}, c, label::NTuple{N}) = setindex!(ds.labels, c, label)
setindex!{D,S,A,B}(ds::DiracState{D,S,A}, c, label::NTuple{B}) = error("length of input label must be $A; got label of length $B")
setindex!(ds::DiracState, x, i...) = setindex!(ds, x, i)

getstate{D,S,N}(ds::DiracState{D,S,N}, label::Tuple) = ds[label] * DiracState{D,S,N}(singlet_dict(label, 1))
getstate(ds::DiracState, i...) = getstate(ds, i)

start(ds::DiracState) = start(ds.labels)
done(ds::DiracState, state) = done(ds.labels, state)
next(ds::DiracState, state) = next(ds.labels, state)
endof(ds::DiracState) = endof(ds.labels)
last(ds::DiracState) = last(ds.labels)
first(ds::DiracState) = first(ds.labels)
collect(ds::DiracState) = collect(ds.labels)

haskey(ds::DiracState, label) = haskey(ds.labels, label)
keys(ds::DiracState) = keys(ds.labels)
values(ds::DiracState) = values(ds.labels)
filter!(f::Function, ds::DiracState) = (filter!(f, ds.labels); return ds)
filter(f::Function, ds::DiracState) = copy_type(ds, filter(f, ds.labels))

##########################
# Mathematical Functions #
##########################
function inner{B}(db::DiracBra{B}, dk::DiracKet{B})
    result = 0
    for (b,c) in db.labels
        for (k,v) in dk.labels
            result += c*v*inner(B,b,k) 
        end
    end
    return result
end

function inner{B<:Orthogonal}(db::DiracBra{B}, dk::DiracKet{B})
    if length(db) > length(dk)
        return ortho_inner(dk.labels, db.labels)
    else
        return ortho_inner(db.labels, dk.labels)
    end
end

+{D,B,N}(a::DiracState{D,B,N}, b::DiracState{D,B,N}) = copy_type(a, mergef(+, a.labels, b.labels))
-{D,B,N}(a::DiracState{D,B,N}, b::DiracState{D,B,N}) = a + (-b)
-(ds::DiracState) = copy_type(ds, mapvals(-, ds.labels))

*(a::DiracBra, b::DiracKet) = inner(a,b)
*{D,B}(a::DiracState{D,B}, b::DiracState{D,B}) = tensor(a,b)
*(c, ds::DiracState) = castvals(*, ds, c)
*(ds::DiracState, c) = castvals(*, c, ds)

/(ds::DiracState, c) = castvals(/, ds, c)

conj(ds::DiracState) = mapvals(conj, ds)
ctranspose{D,B,N}(ds::DiracState{D,B,N}) = DiracState{D',B,N}(conj(ds.labels))
norm(ds::DiracState) = sqrt(sum(v->v^2, values(ds)))
normalize(ds::DiracState) = (1/norm(ds))*ds

xsubspace(ds::DiracState, x) = filter((k,v)->sum(k)==x, ds)
matchlabels_at(ds::DiracState, x, y) = filter((k,v)-> x==k[y], ds)
matchlabels_in(ds::DiracState, x) = filter((k,v)-> x in k, ds)

switch(ds::DiracState, i, j) = mapkeys(k->switch(k,i,j), ds)
permute(ds::DiracState, p) = mapkeys(k->permute(k,p), ds)

######################
# Printing Functions #
######################
labelstr(label) = strip(repr(label)[2:end-1], ',')
statestr(label, ::Type{Ket}) = "| $(labelstr(label)) $rang"
statestr(label, ::Type{Bra}) = "$lang $(labelstr(label)) |"

function show{D}(io::IO, ds::DiracState{D})
    print(io, "$(summary(ds)) with $(length(ds)) state(s):")
    pad = "  "
    maxlen = 30
    i = 1
    for (k,v) in ds
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

tensor_labels(labels) = apply(tuple, labels...)

function tensor_reduce(pairs)
    return (tensor_labels(map(first, pairs)), prod(map(p->p[2], pairs)))
end

function tensor{D,B}(states::DiracState{D,B}...)
    return DiracState{D,B,sum(nfactors,states)}(mergecart!(tensor_reduce, ObjectIdDict(), states))
end

castvals(f::Function, a::DiracState, b::DiracState) = copy_type(ds, dict_castvals(f, a.labels, b.labels))
castvals(f::Function, c, ds::DiracState) = copy_type(ds, dict_castvals(f, c, ds.labels))
castvals(f::Function, ds::DiracState, c) = copy_type(ds, dict_castvals(f, ds.labels, c))
mapkv(f::Function, ds::DiracState) = copy_type(ds, mapkv(f, ds.labels))
mapvals(f::Function, ds::DiracState) = copy_type(ds, mapvals(f, ds.labels))
mapkeys(f::Function, ds::DiracState) = copy_type(ds, mapkeys(f, ds.labels))

export DiracState,
    ket,
    bra,
    getstate,
    normalize,
    xsubspace,
    matchlabels_at,
    matchlabels_in,
    switch,
    permute