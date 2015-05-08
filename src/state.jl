###########
# Ket/Bra #
###########
abstract Ket{P,N,T} <: DiracState{P,N}

type SingleKet{P,N,T} <: Ket{P,N,T}
    coeff::T
    label::StateLabel{N}
end

SingleKet{P,N,T}(::Type{P}, coeff::T, label::StateLabel{N}) = SingleKet{P,N,T}(coeff, label)

typealias StateDict{N,T} Dict{StateLabel{N},T}

type MultiKet{P,N,T} <: Ket{P,N,T}
    dict::StateDict{N,T}
    MultiKet(dict) = new(dict)
    MultiKet(dict::StateDict{0}) = error("Cannot construct a 0-factor state.")
end

MultiKet{P,N,T}(::Type{P}, dict::StateDict{N,T}) = MultiKet{P,N,T}(dict)

ket{P<:AbstractInner}(::Type{P}, label::StateLabel) = SingleKet(P,1,label)
ket{P<:AbstractInner}(::Type{P}, items...) = ket(P, StateLabel(items))

type Bra{P,N,T,K<:Ket} <: DiracState{P,N}
    kt::Ket{P,N,T}
    Bra(kt::Ket{P,N,T}) = new(kt)
end

Bra{P,N,T}(kt::Ket{P,N,T}) = Bra{P,N,T,typeof(kt)}(kt)
bra(items...) = Bra(ket(items...))

typealias MultiBra{P,N,T,K<:MultiKet} Bra{P,N,T,K}
typealias SingleBra{P,N,T,K<:SingleKet} Bra{P,N,T,K}
typealias SingleState Union(SingleKet, SingleBra)
typealias MultiState Union(MultiKet, MultiBra)

########################
# Conversion/Promotion #
########################
Base.convert{P}(::Type{MultiKet}, k::SingleKet{P}) = MultiKet(P, dict(k)) 

######################
# Accessor functions #
######################
coeff(k::SingleKet) = k.coeff
label(k::SingleKet) = k.label
coeff(b::SingleBra) = coeff(b.kt)'
label(b::SingleBra) = label(b.kt)

dict(k::SingleKet) = @compat(Dict(label(k) => coeff(k)))
dict(k::MultiKet) = k.dict
dict(b::Bra) = dict(b.kt)

#######################
# Dict-Like Functions #
#######################
Base.eltype{P,N,T}(::Ket{P,N,T}) = T
Base.eltype{P,N,T}(::Bra{P,N,T}) = T

Base.copy{P,N,T}(kt::SingleKet{P,N,T}) = SingleKet{P,N,T}(copy(coeff(kt)), copy(label(kt)))
Base.copy{P,N,T}(kt::MultiKet{P,N,T}) = MultiKet{P,N,T}(copy(dict(kt)))
Base.copy(br::Bra) = Bra(copy(br.kt))

Base.similar{P,N,T}(kt::SingleKet{P,N,T}, d=StateDict{N,T}()) = MultiKet(P, d)
Base.similar{P}(kt::MultiKet{P}, d=similar(dict(kt))) = MultiKet(P, d)
Base.similar(br::Bra, args...) = Bra(similar(br.kt, args...))

Base.(:(==)){P,N}(a::Ket{P,N}, b::Ket{P,N}) = dict(filternz(a)) == dict(filternz(b))
Base.(:(==)){P,N}(a::Bra{P,N}, b::Bra{P,N}) = a.kt == b.kt
Base.hash{P}(s::DiracState{P}) = hash(dict(filternz(s)), hash(P))
Base.hash(s::DiracState, h::Uint64) = hash(hash(s), h)

Base.length(s::MultiKet) = length(dict(s))
Base.length(::SingleKet) = 1
Base.length(b::Bra) = length(b.kt)

Base.getindex(k::SingleKet, sl::StateLabel) = sl == label(k) ? coeff(k) : throw(KeyError(sl))
Base.getindex(k::MultiKet, sl::StateLabel) = getindex(dict(k), sl)
Base.getindex(b::Bra, sl::StateLabel) = b.kt[sl]'
Base.getindex(s::DiracState, tup::Tuple) = s[StateLabel(tup)]
Base.getindex(s::DiracState, i...) = s[StateLabel(i)]

function Base.setindex!(k::SingleKet, c, sl::StateLabel)
    if sl == label(k)
        coeff(k) = c
        return c
    else
        return setindex!(convert(MultiKet, k), c, sl)
    end
end

Base.setindex!(k::MultiKet, c, sl::StateLabel) = setindex!(dict(k), c, sl)
Base.setindex!(b::Bra, c, sl::StateLabel) = setindex!(dict(b), c', sl)
Base.setindex!(s::DiracState, c, tup::Tuple) = setindex!(s, c, StateLabel(tup))
Base.setindex!(s::DiracState, c, i...) = setindex!(s, c, StateLabel(i))

Base.haskey(k::MultiKet, sl::StateLabel) = haskey(dict(k), sl)
Base.haskey(k::SingleKet, sl::StateLabel) = label(k) == sl
Base.haskey(b::Bra, sl::StateLabel) = haskey(b.kt, sl)
Base.haskey(s::DiracState, sl) = haskey(s, StateLabel(sl))

function Base.get(k::SingleKet, sl::StateLabel, default=predict_zero(eltype(k)))
    if sl == label(k)
        return coeff(k)
    else
        return default
    end
end

Base.get(k::MultiKet, sl::StateLabel, default=predict_zero(eltype(k))) = get(dict(k), sl, default)
Base.get(b::Bra, sl::StateLabel, default=predict_zero(eltype(b))) = get(b.kt, sl, default')'
Base.get(s::DiracState, sl, default=predict_zero(eltype(s))) = get(s, StateLabel(sl), default)

function Base.delete!(k::SingleKet, sl::StateLabel)
    if sl == label(k)
        coeff(k) = predict_zero(eltype(k))
    end
    return k
end

Base.delete!(k::MultiKet, sl::StateLabel) = (delete!(dict(k), sl); return k)
Base.delete!(b::Bra, sl::StateLabel) = delete!(b.kt, sl)
Base.delete!(s::DiracState, sl) = delete!(s, StateLabel(sl))

########################
# Iteration/Collection #
########################
iter(k::SingleKet) = tuple(tuple(label(k), coeff(k)))
iter(k::MultiKet) = dict(k)
iter(b::Bra) = iter(b.kt)

coeffs(k::SingleKet) = tuple(coeff(k))
labels(k::SingleKet) = tuple(label(k))
coeffs(k::MultiKet) = values(dict(k))
labels(k::MultiKet) = keys(dict(k))
coeffs(b::Bra) = coeffs(b.kt)
labels(b::Bra) = labels(b.kt)

Base.collect(s::SingleState) = [(label(s), coeff(s))]
Base.collect(kt::MultiKet) = collect(dict(kt))
Base.collect{P,N,T}(br::MultiBra{P,N,T}) = collect_pairs!(Array((StateLabel{N}, T), length(br)), br)

function collect_pairs!(result, br::MultiBra)
    i = 1
    for (k,v) in iter(br)
        result[i] = (k, v')
        i += 1
    end
    return result
end

abstract IterFlag{b}

Base.start(state::DiracState) = start(iter(state))
Base.start(::SingleState) = IterFlag{true}

Base.next(kt::Ket, i) = next(iter(kt), i)
Base.next(s::SingleState, ::Type{IterFlag{true}}) = ((label(s), coeff(s)), IterFlag{false})

function Base.next(br::MultiBra, i)
    (k,v), n = next(iter(br), i)
    return ((k,v'), n)
end

Base.done(kt::MultiKet, i) = done(iter(kt), i)
Base.done(kt::SingleKet, ::Type{IterFlag{false}}) = true
Base.done(kt::SingleKet, ::Type{IterFlag{true}}) = false
Base.done(br::Bra, i) = done(br.kt, i)

Base.first(state::DiracState) = next(state, start(state))

##############
# ctranspose #
##############
Base.ctranspose(k::Ket) = Bra(k)
Base.ctranspose(b::Bra) = b.kt

#########
# inner #
#########
inner{P,N}(br::SingleBra{P,N}, kt::SingleKet{P,N}) = coeff(br) * coeff(kt) * P(label(br), label(kt))

function inner{P,N}(br::SingleBra{P,N}, kt::MultiKet{P,N})
    result = predict_zero(inner_coefftype(br, kt))
    c = coeff(br)
    b = label(br)
    for (k,v) in iter(kt)
        result += c * v * P(b, k)
    end
    return result
end

function inner{P,N}(br::MultiBra{P,N}, kt::SingleKet{P,N})
    result = predict_zero(inner_coefftype(br, kt))
    v = coeff(kt)
    k = label(kt)
    for (b,c) in iter(br)
        result += c' * v * P(b, k)
    end
    return result
end

function inner{P,N}(br::MultiBra{P,N}, kt::MultiKet{P,N})
    result = predict_zero(inner_coefftype(br, kt))
    for (b,c) in iter(br), (k,v) in iter(kt)
        result += c' * v * P(b, k)
    end
    return result  
end

inner{N}(br::MultiBra{KroneckerDelta,N}, kt::SingleKet{KroneckerDelta,N}) = get(br, label(kt)) * coeff(kt)
inner{N}(br::SingleBra{KroneckerDelta,N}, kt::MultiKet{KroneckerDelta,N}) = get(kt, label(br)) * coeff(br)

function inner{N}(br::MultiBra{KroneckerDelta,N}, kt::MultiKet{KroneckerDelta,N})
    if length(br) < length(kt)
        return ortho_inner(kt, br)
    else
        return ortho_inner(br, kt)
    end
end

function ortho_inner(a::DiracState{KroneckerDelta}, b::DiracState{KroneckerDelta})
    result = predict_zero(inner_coefftype(a, b))
    for sl in labels(b)
        if haskey(a, sl)
            result += a[sl]*b[sl]
        end
    end
    return result
end

Base.(:*)(br::Bra, kt::Ket) = inner(br,kt)

inner_eval(f, s::DiracState) = mapcoeffs(x->inner_eval(f,x),s)

##########
# act_on #
##########
act_on(kt::Ket, br::Bra, i) = act_on(kt', br', i)'
act_on{P}(br::Bra{P,1}, kt::Ket{P,1}, i) = i==1 ? inner(br, kt) : throw(BoundsError())

function act_on{P,N}(br::SingleBra{P,1}, kt::SingleKet{P,N}, i)
    return (br * ket(P, label(kt)[i])) * ket(P, except(label(kt), i))
end

function act_on{P,N}(br::Bra{P,1}, kt::Ket{P,N}, i)
    result = StateDict{N-1, inner_coefftype(br,kt)}()
    return MultiKet(P, act_on_dict!(result, br, kt, i, P))
end

function act_on_dict!{P}(result, br::MultiBra, kt::MultiKet, i, ::Type{P})
    for (b,c) in iter(br), (k,v) in iter(kt)
        add_to_dict!(result, except(k,i), c'*v*P(b[1], k[i]))
    end
    return result
end

function act_on_dict!{P}(result, br::SingleBra, kt::MultiKet, i, ::Type{P})
    b = label(br)[1]
    c = coeff(br)
    for (k,v) in iter(kt)
        add_to_dict!(result, except(k,i), c'*v*P(b, k[i]))
    end
    return result
end

function act_on_dict!{P}(result, br::MultiBra, kt::SingleKet, i, ::Type{P})
    k = label(kt)[i]
    new_k = except(label(kt),i)
    v = coeff(kt)
    for (b,c) in iter(br)
        add_to_dict!(result, new_k, c'*v*P(b[1], k))
    end
    return result
end

###########
# Scaling #
###########
Base.scale!(k::MultiKet, c::Number) = (dscale!(dict(k), c); return k)
Base.scale!(k::SingleKet, c::Number) = (coeff(k) = c; return k)
Base.scale!(c::Number, k::Ket) = scale!(k,c)
Base.scale!(b::Bra, c::Number) = scale!(b', c')'
Base.scale!(c::Number, b::Bra) = scale!(b,c)

Base.scale(k::MultiKet, c::Number) = similar(k, dscale(dict(k), c))
Base.scale{P}(k::SingleKet{P}, c::Number) = SingleKet(P, coeff(k) * c, label(k))
Base.scale(c::Number, k::Ket) = scale(k,c)
Base.scale(b::Bra, c::Number) = scale(b', c')'
Base.scale(c::Number, b::Bra) = scale(b,c)

Base.(:*)(c::Number, s::DiracState) = scale(c, s)
Base.(:*)(s::DiracState, c::Number) = scale(s, c)
Base.(:/)(s::DiracState, c::Number) = scale(s, 1/c)

###########
# + and - #
###########
Base.(:-){P,N}(kt::Ket{P,N}) = -1 * kt
Base.(:-)(br::Bra) = ctranspose(-(br'))

function Base.(:+){P,N,T,V}(a::SingleKet{P,N,T}, b::SingleKet{P,N,V})
    result = StateDict{N, promote_type(T,V)}()
    add_to_dict!(result, label(a), coeff(a))
    add_to_dict!(result, label(b), coeff(b))
    return MultiKet(P, result)
end

Base.(:+){P,N}(a::MultiKet{P,N}, b::SingleKet{P,N}) = similar(a, add_merge(dict(a), b))
Base.(:+){P,N}(a::SingleKet{P,N}, b::MultiKet{P,N}) = +(b, a)
Base.(:+){P,N}(a::MultiKet{P,N}, b::MultiKet{P,N}) = similar(b, add_merge(dict(a), dict(b)))

function Base.(:-){P,N,T,V}(a::SingleKet{P,N,T}, b::SingleKet{P,N,V})
    result = StateDict{N, promote_type(T,V)}()
    add_to_dict!(result, label(a), coeff(-a))
    add_to_dict!(result, label(b), coeff(-b))
    return MultiKet(P, result)
end

Base.(:-){P,N}(a::MultiKet{P,N}, b::SingleKet{P,N}) = similar(a, sub_merge(dict(a), b))
Base.(:-){P,N}(a::SingleKet{P,N}, b::MultiKet{P,N}) = -(b, a)
Base.(:-){P,N}(a::MultiKet{P,N}, b::MultiKet{P,N}) = similar(b, sub_merge(dict(a), dict(b)))

Base.(:+)(a::Bra, b::Bra) = ctranspose(a' + b')
Base.(:-)(a::Bra, b::Bra) = ctranspose(a' - b')

##########
# tensor #
##########
tensor{P}(a::SingleKet{P}, b::SingleKet{P}) = SingleKet(P, coeff(a)*coeff(b), tensor(label(a), label(b)))
tensor{P}(a::SingleKet{P}, b::MultiKet{P}) = MultiKet(P, tensor_merge(a, dict(b)))
tensor{P}(a::MultiKet{P}, b::SingleKet{P}) = MultiKet(P, tensor_merge(dict(b), b))
tensor{P}(a::MultiKet{P}, b::MultiKet{P}) = MultiKet(P, tensor_merge(dict(a), dict(b)))
tensor(a::Bra, b::Bra) = tensor(a', b')'

Base.(:*)(a::Ket, b::Ket) = tensor(a,b)
Base.(:*)(a::Bra, b::Bra) = tensor(a,b)

#################
# Normalization #
#################
Base.norm(k::SingleKet) = abs(coeff(k))
Base.norm(k::MultiKet) = sqrt(sum(abs2, coeffs(k)))
Base.norm(b::Bra) = norm(b')
normalize(s::DiracState) = (1/norm(s))*s
normalize!(s::DiracState) = scale!(1/norm(s), s)

####################
# Raising/Lowering #
####################
ladderdict{P,N}(state::DiracState{P,N}) = StateDict{N, promote_type(Float64, eltype(state))}()

lower{P}(state::DiracState{P,1}) = lower(state, 1)
lower(k::MultiKet, i) = similar(k, lowerdict!(ladderdict(k), dict(k), i))

function lower(k::SingleKet, i)
    c = sqrt(label(k)[i])*coeff(k) 
    sl = setindex(label(k), label(k)[i] - 1, i)
    return SingleKet(P, c, sl)
end

lower(b::Bra, i) = lower(b', i)'

raise{P}(state::DiracState{P,1}) = raise(state, 1)
raise(k::MultiKet, i) = similar(k, raisedict!(ladderdict(k), dict(k), i))

function raise(k::SingleKet, i)
    c = sqrt(label(k)[i]+1)*coeff(k) 
    sl = setindex(label(k), label(k)[i] + 1, i)
    return SingleKet(P, c, sl)
end

raise(b::Bra, i) = raise(b', i)'

function lowerdict!(result, d, i)
    for (k,v) in d
        add_to_dict!(result, setindex(k, k[i] - 1, i), sqrt(k[i])*v)
    end
    return result
end

function raisedict!(result, d, i)
    for (k,v) in d
        add_to_dict!(result, setindex(k, k[i] + 1, i), sqrt(k[i]+1)*v)
    end
    return result
end

########################
# Misc. Math Functions #
########################
nfactors{P,N}(::DiracState{P,N}) = N
xsubspace(s::DiracState, x) = similar(s, filter((k,v)->is_sum_x(k,x), dict(s)))

switch(s::MultiKet, i, j) = maplabels(sl->switch(sl, i, j), s)
switch(s::SingleKet, i, j) = SingleKet(coeff(s), switch(label(s), i, j))
switch(b::Bra, i, j) = switch(b', i, j)'

permute(s::MultiKet, perm::Vector) = maplabels(sl->permute(sl,perm), s)
permute{P}(s::SingleKet{P}, perm::Vector) = SingleKet(P, coeff(s), permute(label(s), perm))
permute(b::Bra, perm::Vector) = permute(b', perm)'

filternz!(k::MultiKet) = (filter!(nzcoeff, dict(k)); return k)
filternz!(k::SingleKet) = k # no-op
filternz!(b::Bra) = (filternz!(b.kt); return b)

filternz(s::DiracState) = similar(s, filter(nzcoeff, dict(s)))

function represent{P}(kt::Ket{P}, basis)
    T = promote_type(return_type(P), eltype(kt))
    return T[bra(i) * kt for i in basis]
end

function represent{P}(kt::Ket{P}, basis...)
    prodbasis = product(basis...)
    T = promote_type(return_type(P), eltype(kt))
    return T[bra(i...) * kt for i in prodbasis]
end

represent(br::Bra, basis...) = represent(br', basis...)'

# should always be pure, of course,
# but makes a good sanity check function
purity(kt::Ket) = purity(kt*kt')
purity(br::Bra) = purity(br')

######################
# Printing Functions #
######################
ktstr(sl) = "| $(labelstr(sl)) $rang"
brstr(sl) = "$lang $(labelstr(sl)) |"

labelrepr(kt::Ket, sl, pad) = "$pad$(kt[sl]) $(ktstr(sl))"
labelrepr(br::Bra, sl, pad) = "$pad$(br[sl]) $(brstr(sl))"

Base.summary(s::DiracState) = "$(typeof(s)) with $(length(s)) state(s)"
Base.show(io::IO, s::DiracState) = dirac_show(io, s)
Base.showcompact(io::IO, s::DiracState) = dirac_showcompact(io, s)
Base.repr(s::DiracState) = dirac_repr(s)

export Ket,
    Bra,
    ket,
    bra,
    represent,
    nfactors,
    xsubspace,
    permute,
    switch,
    switch!,
    permute!,
    filternz!,
    filternz,
    purity,
    lower,
    raise,
    act_on,
    inner_eval