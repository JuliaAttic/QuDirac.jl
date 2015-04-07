###########
# Ket/Bra #
###########
typealias StateDict{N} Dict{NTuple{N},Number}

type Ket{P,N} <: DiracState{P,N}
    ptype::P
    dict::StateDict{N}
    Ket(ptype, dict) = new(ptype, dict)
    Ket(ptype, dict::StateDict{0}) = error("Cannot construct a 0-factor state; did you mean to construct a scalar?")
end

Ket{P,N}(ptype::P, dict::StateDict{N}) = Ket{P,N}(ptype,dict)

ket{N}(ptype::AbstractInner, label::NTuple{N}) = Ket(ptype, single_dict(StateDict{N}(), label, 1))
ket(ptype::AbstractInner, items...) = ket(ptype, items)

macro default_inner(ptype)
    @eval begin
        ket(items...) = ket(($ptype)(), items)
    end
end

@default_inner Orthonormal

type Bra{P,N} <: DiracState{P,N}
    kt::Ket{P,N}
end

Bra{P,N}(kt::Ket{P,N}) = Bra{P,N}(kt)
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
Base.copy(kt::Ket) = Ket(ptype(s), copy(dict(s)))
Base.copy(br::Bra) = Bra(copy(br.ket))

Base.similar(kt::Ket, d::StateDict=similar(dict(kt)); P=ptype(kt)) = Ket(P, d)
Base.similar(br::Bra, d::StateDict=similar(dict(br)); P=ptype(br)) = Bra(P, d)

Base.(:(==)){P,N}(a::Ket{P,N}, b::Ket{P,N}) = dict(a) == dict(b)
Base.(:(==)){P,N}(a::Bra{P,N}, b::Bra{P,N}) = a.kt == b.kt
Base.hash(s::DiracState) = hash(dict(s), hash(typeof(s)))

Base.length(s::DiracState) = length(dict(s))

Base.getindex(k::Ket, label::Tuple) = getindex(dict(k), label)
Base.getindex(b::Bra, label::Tuple) = b.kt[label]'
Base.getindex(s::DiracState, i...) = s[i]

Base.setindex!(k::Ket, c, label::Tuple) = setindex!(dict(k), c, label)
Base.setindex!(b::Bra, c, label::Tuple) = setindex!(dict(b), c', label)
Base.setindex!(s::DiracState, c, i...) = setindex!(s,c,i)

Base.haskey(s::DiracState, label) = haskey(dict(s), label)

Base.get(k::Ket, label::Tuple, default=0) = get(dict(k), label, default)
Base.get(b::Bra, label::Tuple, default=0) = haskey(b, label) ? b[label] : default
Base.get(k::Ket, label, default=0) = get(k, tuple(label), default)
Base.get(b::Bra, label, default=0) = get(b, tuple(label), default)

Base.delete!(s::DiracState, label::Tuple) = (delete!(dict(s), label); return s)
Base.delete!(s::DiracState, label) = delete!(s, tuple(label))

labels(s::DiracState) = keys(dict(s))
QuBase.coeffs(kt::Ket) = values(dict(kt))
QuBase.coeffs(br::Bra) = imap(ctranspose, coeffs(br.kt))

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
        result += c'*v*inner_rule(prodtype,b,k)
    end
    return result  
end

function ortho_inner(a::DiracState{Orthonormal}, b::DiracState{Orthonormal})
    result = 0
    prodtype = ptype(a)
    for label in keys(dict(b))
        if haskey(a, label)
            result += a[label]*b[label]*inner_rule(prodtype,label,label)
        end
    end
    return result
end

function inner{N}(br::Bra{Orthonormal,N}, kt::Ket{Orthonormal,N})
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

function act_on{P,N}(br::Bra{P,1}, kt::Ket{P,N}, i)
    result = StateDict{N-1}()
    prodtype = ptype(kt)
    for (b,c) in dict(br), (k,v) in dict(kt)
        add_to_dict!(result, 
                     except(k,i),
                     c'*v*inner_rule(prodtype,b[1],k[i]))
    end
    return Ket(prodtype, result)
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
QuBase.tensor{P,A,B}(a::Ket{P,A}, b::Ket{P,B}) = Ket(ptype(b), tensorstate!(StateDict{A+B}(), dict(a), dict(b)))
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
xsubspace(s::DiracState, x) = similar(s, filter((k,v)->isx(k,x), dict(s)))
switch(s::DiracState, i, j) = similar(s, mapkeys(label->switch(label,i,j), dict(s)))
permute(s::DiracState, perm::Vector) = similar(s, mapkeys(label->permute(label,perm), dict(s)))
switch!(s::DiracState, i, j) = (mapkeys!(label->switch(label,i,j), dict(s)); return s)
Base.permute!(s::DiracState, perm::Vector) = (mapkeys!(label->permute(label,perm), dict(s)); return s)

filternz!(s::DiracState) = (filter!(nzcoeff, dict(s)); return s)
filternz(s::DiracState) = similar(s, filter(nzcoeff, dict(s)))

function wavefunc(f::Function, kt::Ket)
    return (args...) -> sum(pair->pair[2]*f(pair[1])(args...), dict(kt))
end

# should always be pure, of course,
# but makes a good sanity check function
purity(kt::Ket) = purity(kt*kt')
purity(br::Bra) = purity(br.kt)

inner_eval(f::Function, s::DiracState) = mapcoeffs(x->inner_eval(f,x),s)
inner_eval(s::DiracState) = mapcoeffs(inner_eval,s)

######################
# Printing Functions #
######################
labelstr(label) = join(map(repr, label), ',')
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
    labels,
    act_on,
    @default_inner,
    inner_eval