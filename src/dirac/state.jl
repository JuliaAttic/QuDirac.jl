###########
# Ket/Bra #
###########
typealias StateDict{N,T} Dict{StateLabel{N},T}

type Ket{P,N,T} <: DiracState{P,N}
    ptype::P
    dict::StateDict{N,T}
    Ket(ptype, dict) = new(ptype, dict)
    Ket(ptype, dict::StateDict{0}) = error("Cannot construct a 0-factor state; did you mean to construct a scalar?")
end

Ket{P,N,T}(ptype::P, dict::StateDict{N,T}) = Ket{P,N,T}(ptype, dict)

ket(ptype::AbstractInner, label::StateLabel) = Ket(ptype, [label => 1])
ket(ptype::AbstractInner, items...) = ket(ptype, StateLabel(items))

function default_inner(ptype::AbstractInner)
    QuDirac.ket(items...) = ket(ptype, StateLabel(items))
end

default_inner(KroneckerDelta())

type Bra{P,N,T} <: DiracState{P,N}
    kt::Ket{P,N,T}
end

Bra{P,N,T}(kt::Ket{P,N,T}) = Bra{P,N,T}(kt)
Bra(items...) = Bra(Ket(items...))
bra(items...) = Bra(ket(items...))

######################
# Accessor functions #
######################
dict(k::Ket) = k.dict
dict(b::Bra) = dict(b.kt)

ptype(k::Ket) = k.ptype
ptype(b::Bra) = ptype(b.kt)

#######################
# Dict-Like Functions #
#######################
Base.eltype{P,N,T}(::Ket{P,N,T}) = T
Base.eltype{P,N,T}(::Bra{P,N,T}) = T

Base.copy(kt::Ket) = Ket(ptype(kt), copy(dict(kt)))
Base.copy(br::Bra) = Bra(copy(br.kt))

Base.similar(kt::Ket, d=similar(dict(kt)); P=ptype(kt)) = Ket(P, d)
Base.similar(br::Bra, d=similar(dict(br)); P=ptype(br)) = Bra(P, d)

Base.(:(==)){P,N}(a::Ket{P,N}, b::Ket{P,N}) = dict(filternz(a)) == dict(filternz(b))
Base.(:(==)){P,N}(a::Bra{P,N}, b::Bra{P,N}) = a.kt == b.kt
Base.hash(s::DiracState) = hash(dict(filternz(s)), hash(ptype(s)))
Base.hash(s::DiracState, h::Uint64) = hash(hash(s), h)

Base.length(s::DiracState) = length(dict(s))

Base.getindex(k::Ket, label::StateLabel) = getindex(dict(k), label)
Base.getindex(b::Bra, label::StateLabel) = b.kt[label]'
Base.getindex(s::DiracState, tup::Tuple) = s[StateLabel(tup)]
Base.getindex(s::DiracState, i...) = s[StateLabel(i)]

Base.setindex!(k::Ket, c, label::StateLabel) = setindex!(dict(k), c, label)
Base.setindex!(b::Bra, c, label::StateLabel) = setindex!(dict(b), c', label)
Base.setindex!(s::DiracState, c, tup::Tuple) = setindex!(s, c, StateLabel(tup))
Base.setindex!(s::DiracState, c, i...) = setindex!(s, c, StateLabel(i))

Base.haskey(s::DiracState, label::StateLabel) = haskey(dict(s), label)
Base.haskey(s::DiracState, label) = haskey(s, StateLabel(label))

Base.get(k::Ket, label::StateLabel, default=0) = get(dict(k), label, default)
Base.get(b::Bra, label::StateLabel, default=0) = get(dict(b), label, default')'
Base.get(s::DiracState, label, default=0) = get(s, StateLabel(label), default)

Base.delete!(s::DiracState, label::StateLabel) = (delete!(dict(s), label); return s)
Base.delete!(s::DiracState, label) = delete!(s, StateLabel(label))

##############
# ctranspose #
##############
Base.ctranspose(k::Ket) = Bra(k)
Base.ctranspose(b::Bra) = b.kt

#########
# inner #
#########
inner(br::Bra, kt::Ket) = error("inner(b::Bra,k::Ket) is only defined when nfactors(b) == nfactors(k)")

function inner{P,N}(br::Bra{P,N}, kt::Ket{P,N})
    result = 0
    prodtype = ptype(kt)
    for (b,c) in dict(br), (k,v) in dict(kt)
        result += c'*v*eval_inner_rule(prodtype,b,k)
    end
    return result  
end

function ortho_inner(a::DiracState{KroneckerDelta}, b::DiracState{KroneckerDelta})
    result = 0
    for label in keys(dict(b))
        if haskey(a, label)
            result += a[label]*b[label]
        end
    end
    return result
end

function inner{N}(br::Bra{KroneckerDelta,N}, kt::Ket{KroneckerDelta,N})
    if length(br) < length(kt)
        return ortho_inner(kt, br)
    else
        return ortho_inner(br, kt)
    end
end

Base.(:*)(br::Bra, kt::Ket) = inner(br,kt)

##########
# act_on #
##########

act_on(br::Bra, kt::Ket, i) = error("inner(b::Bra,k::Ket,i) is only defined when nfactors(b) == 1")

function act_on{P}(br::Bra{P,1}, kt::Ket{P,1}, i)
    if i==1
        return inner(br, kt)
    else
        throw(BoundsError())
    end
end

function act_on{P,N,A,B}(br::Bra{P,1,A}, kt::Ket{P,N,B}, i)
    prodtype = ptype(br)
    result = StateDict{N-1, promote_type(A, B, inner_rettype(prodtype))}()
    return Ket(prodtype, act_on_dict!(result, br, kt, i, prodtype))
end

function act_on_dict!(result, br::Bra, kt::Ket, i, prodtype)
    for (b,c) in dict(br), (k,v) in dict(kt)
        add_to_dict!(result, except(k,i), c'*v*eval_inner_rule(prodtype, b[1], k[i]))
    end
    return result
end

###########
# Scaling #
###########
Base.scale!(k::Ket, c::Number) = (dscale!(dict(k), c); return k)
Base.scale!(c::Number, k::Ket) = scale!(k,c)
Base.scale!(b::Bra, c::Number) = Bra(scale!(b.kt, c'))
Base.scale!(c::Number, b::Bra) = scale!(b,c)

Base.scale(k::Ket, c::Number) = similar(k, dscale(dict(k), c))
Base.scale(c::Number, k::Ket) = scale(k,c)
Base.scale(b::Bra, c::Number) = Bra(scale(b.kt, c'))
Base.scale(c::Number, b::Bra) = scale(b,c)

Base.(:*)(c::Number, s::DiracState) = scale(c, s)
Base.(:*)(s::DiracState, c::Number) = scale(s, c)
Base.(:/)(s::DiracState, c::Number) = scale(s, 1/c)

###########
# + and - #
###########
Base.(:-){P,N}(kt::Ket{P,N}) = -1 * kt
Base.(:-)(br::Bra) = Bra(-br.kt)

Base.(:+){P,N}(a::Ket{P,N}, b::Ket{P,N}) = similar(b, add_merge(dict(a), dict(b)))
Base.(:-){P,N}(a::Ket{P,N}, b::Ket{P,N}) = similar(b, sub_merge(dict(a), dict(b)))

Base.(:+)(a::Bra, b::Bra) = Bra(a.kt + b.kt)
Base.(:-)(a::Bra, b::Bra) = Bra(a.kt - b.kt)

##########
# tensor #
##########
QuBase.tensor{P}(a::Ket{P}, b::Ket{P}) = Ket(ptype(b), tensordict(dict(a), dict(b)))
QuBase.tensor(a::Bra, b::Bra) = tensor(a.kt, b.kt)'

Base.(:*)(a::Ket, b::Ket) = tensor(a,b)
Base.(:*)(a::Bra, b::Bra) = tensor(a,b)

#################
# Normalization #
#################
Base.norm(s::DiracState) = sqrt(sum(abs2, values(dict(s))))
QuBase.normalize(s::DiracState) = (1/norm(s))*s
QuBase.normalize!(s::DiracState) = scale!(1/norm(s), s)

########################
# Misc. Math Functions #
########################
nfactors{P,N}(::DiracState{P,N}) = N
xsubspace(s::DiracState, x) = similar(s, filter((k,v)->is_sum_x(k,x), dict(s)))
switch(s::DiracState, i, j) = similar(s, mapkeys(label->switch(label,i,j), dict(s)))
permute(s::DiracState, perm::Vector) = similar(s, mapkeys(label->permute(label,perm), dict(s)))

filternz!(s::DiracState) = (filter!(nzcoeff, dict(s)); return s)
filternz(s::DiracState) = similar(s, filter(nzcoeff, dict(s)))

function wavefunc(f::Function, kt::Ket)
    return (args...) -> sum(pair->pair[2]*f(pair[1])(args...), dict(kt))
end

# should always be pure, of course,
# but makes a good sanity check function
purity(kt::Ket) = purity(kt*kt')
purity(br::Bra) = purity(br.kt)

inner_eval(f, s::DiracState) = mapcoeffs(x->inner_eval(f,x),s)

######################
# Printing Functions #
######################
ktstr(label) = "| $(labelstr(label)) $rang"
brstr(label) = "$lang $(labelstr(label)) |"

labelrepr(kt::Ket, label, pad) = "$pad$(kt[label]) $(ktstr(label))"
labelrepr(br::Bra, label, pad) = "$pad$(br[label]) $(brstr(label))"

Base.summary(s::DiracState) = "$(typeof(s)) with $(length(s)) state(s)"
Base.show(io::IO, s::DiracState) = dirac_show(io, s)
Base.showcompact(io::IO, s::DiracState) = dirac_showcompact(io, s)
Base.repr(s::DiracState) = dirac_repr(s)

export ket,
    bra,
    nfactors,
    xsubspace,
    permute,
    switch,
    switch!,
    permute!,
    filternz!,
    filternz,
    purity,
    wavefunc,
    act_on,
    default_inner,
    inner_eval