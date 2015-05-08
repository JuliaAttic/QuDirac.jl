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

typealias MultiBra{P,N,T} Bra{P,N,T,MultiKet{P,N,T}}
typealias SingleBra{P,N,T} Bra{P,N,T,SingleKet{P,N,T}}

########################
# Conversion/Promotion #
########################
Base.convert{P}(::Type{MultiKet}, k::SingleKet{P}) = MultiKet(P, dict(k)) 

######################
# Accessor functions #
######################
dict(k::SingleKet) = @compat(Dict(k.label => k.coeff))
dict(k::MultiKet) = k.dict
dict(b::Bra) = dict(b.kt)

#######################
# Dict-Like Functions #
#######################
Base.eltype{P,N,T}(::Ket{P,N,T}) = T
Base.eltype{P,N,T}(::Bra{P,N,T}) = T

Base.copy{P,N,T}(kt::SingleKet{P,N,T}) = SingleKet{P,N,T}(copy(kt.coeff), copy(kt.label))
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

Base.getindex(k::SingleKet, label::StateLabel) = label == k.label ? k.coeff : throw(KeyError(label))
Base.getindex(k::MultiKet, label::StateLabel) = getindex(dict(k), label)
Base.getindex(b::Bra, label::StateLabel) = b.kt[label]'
Base.getindex(s::DiracState, tup::Tuple) = s[StateLabel(tup)]
Base.getindex(s::DiracState, i...) = s[StateLabel(i)]

function Base.setindex!(k::SingleKet, c, label::StateLabel)
    if label == k.label
        k.coeff = c
        return c
    else
        return setindex!(convert(MultiKet, k), c, label)
    end
end

Base.setindex!(k::MultiKet, c, label::StateLabel) = setindex!(dict(k), c, label)
Base.setindex!(b::Bra, c, label::StateLabel) = setindex!(dict(b), c', label)
Base.setindex!(s::DiracState, c, tup::Tuple) = setindex!(s, c, StateLabel(tup))
Base.setindex!(s::DiracState, c, i...) = setindex!(s, c, StateLabel(i))

Base.haskey(k::MultiKet, label::StateLabel) = haskey(dict(k), label)
Base.haskey(k::SingleKet, label::StateLabel) = k.label == label
Base.haskey(b::Bra, label::StateLabel) = haskey(b.kt, label)
Base.haskey(s::DiracState, label) = haskey(s, StateLabel(label))

function Base.get(k::SingleKet, label::StateLabel, default=predict_zero(eltype(k)))
    if label == k.label
        return k.coeff
    else
        return default
    end
end

Base.get(k::MultiKet, label::StateLabel, default=predict_zero(eltype(k))) = get(dict(k), label, default)
Base.get(b::Bra, label::StateLabel, default=predict_zero(eltype(b))) = get(b.kt, label, default')'
Base.get(s::DiracState, label, default=predict_zero(eltype(s))) = get(s, StateLabel(label), default)

function Base.delete!(k::SingleKet, label::StateLabel)
    if label == k.label
        k.coeff = predict_zero(eltype(k))
    end
    return k
end

Base.delete!(k::MultiKet, label::StateLabel) = (delete!(dict(k), label); return k)
Base.delete!(b::Bra, label::StateLabel) = delete!(b.kt, label)
Base.delete!(s::DiracState, label) = delete!(s, StateLabel(label))

Base.collect(kt::SingleKet) = [(kt.label, kt.coeff)]
Base.collect(kt::MultiKet) = collect(dict(kt))
Base.collect{P,N,T}(br::SingleBra{P,N,T}) = [(br.kt.label, br.kt.coeff')]
Base.collect{P,N,T}(br::MultiBra{P,N,T}) = collect_pairs!(Array((StateLabel{N}, T), length(br)), br)

function collect_pairs!{P,N,T}(result, br::MultiBra{P,N,T})
    i = 1
    for (k,v) in dict(br)
        result[i] = (k, v')
        i += 1
    end
    return result
end

Base.start(state::DiracState) = start(dict(state))
Base.next(kt::Ket, i) = next(dict(kt), i)

function Base.next(br::Bra, i)
    (k,v), n = next(dict(br), i)
    return ((k,v'), n)
end

Base.done(state::DiracState, i) = done(dict(state), i)
Base.first(state::DiracState) = next(state, start(state))

##############
# ctranspose #
##############
Base.ctranspose(k::Ket) = Bra(k)
Base.ctranspose(b::Bra) = b.kt

#########
# inner #
#########
inner{P,N,T}(br::SingleBra{P,N,T}, kt::SingleKet{P,N}) = br.kt.coeff' * kt.coeff * P(br.kt.label, kt.label)

function inner{P,N,T}(br::SingleBra{P,N,T}, kt::MultiKet{P,N})
    result = predict_zero(inner_coefftype(br, kt))
    c = br.kt.coeff'
    b = br.kt.label
    for (k,v) in dict(kt)
        result += c * v * P(b, k)
    end
    return result
end

function inner{P,N,T}(br::MultiBra{P,N,T}, kt::SingleKet{P,N})
    result = predict_zero(inner_coefftype(br, kt))
    v = kt.coeff
    k = kt.label
    for (b,c) in dict(br)
        result += c' * v * P(b, k)
    end
    return result
end

inner{N,T}(br::MultiBra{KroneckerDelta,N,T}, kt::SingleKet{KroneckerDelta,N}) = get(br, kt.label) * kt.coeff
inner{N,T}(br::SingleBra{KroneckerDelta,N,T}, kt::MultiKet{KroneckerDelta,N}) = get(kt, br.kt.label) * br.kt.coeff'

function inner{P,N,T}(br::MultiBra{P,N,T}, kt::MultiKet{P,N})
    result = predict_zero(inner_coefftype(br, kt))
    for (b,c) in dict(br), (k,v) in dict(kt)
        result += c' * v * P(b, k)
    end
    return result  
end

function inner{N,T}(br::MultiBra{KroneckerDelta,N,T}, kt::MultiKet{KroneckerDelta,N})
    if length(br) < length(kt)
        return ortho_inner(kt, br)
    else
        return ortho_inner(br, kt)
    end
end

function ortho_inner(a::DiracState{KroneckerDelta}, b::DiracState{KroneckerDelta})
    result = predict_zero(inner_coefftype(a, b))
    for label in keys(dict(b))
        if haskey(a, label)
            result += a[label]*b[label]
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

function act_on{P,N,T}(br::SingleBra{P,1,T}, kt::SingleKet{P,N}, i)
    return (br * ket(P, kt.label[i])) * ket(P, except(kt.label, i))
end

function act_on{P,N}(br::Bra{P,1}, kt::Ket{P,N}, i)
    result = StateDict{N-1, inner_coefftype(br,kt)}()
    return MultiKet(P, act_on_dict!(result, br, kt, i, P))
end

function act_on_dict!{P,N,T}(result, br::MultiBra{P,N,T}, kt::MultiKet, i, ::Type{P})
    for (b,c) in dict(br), (k,v) in dict(kt)
        add_to_dict!(result, except(k,i), c'*v*P(b[1], k[i]))
    end
    return result
end

function act_on_dict!{P,N,T}(result, br::SingleBra{P,N,T}, kt::MultiKet, i, ::Type{P})
    b = br.kt.label[1]
    c = br.kt.coeff'
    for (k,v) in dict(kt)
        add_to_dict!(result, except(k,i), c'*v*P(b, k[i]))
    end
    return result
end

function act_on_dict!{P,N,T}(result, br::MultiBra{P,N,T}, kt::SingleKet, i, ::Type{P})
    k = kt.label[i]
    new_k = except(kt.label,i)
    v = kt.coeff
    for (b,c) in dict(br)
        add_to_dict!(result, new_k, c'*v*P(b[1], k))
    end
    return result
end

###########
# Scaling #
###########
Base.scale!(k::MultiKet, c::Number) = (dscale!(dict(k), c); return k)
Base.scale!(k::SingleKet, c::Number) = (k.coeff = c; return k)
Base.scale!(c::Number, k::Ket) = scale!(k,c)
Base.scale!(b::Bra, c::Number) = scale!(b', c')'
Base.scale!(c::Number, b::Bra) = scale!(b,c)

Base.scale(k::MultiKet, c::Number) = similar(k, dscale(dict(k), c))
Base.scale{P}(k::SingleKet{P}, c::Number) = SingleKet(P, k.coeff * c, k.label)
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
    add_to_dict!(result, a.label, a.coeff)
    add_to_dict!(result, b.label, b.coeff)
    return MultiKet(P, result)
end

Base.(:+){P,N}(a::MultiKet{P,N}, b::SingleKet{P,N}) = similar(a, add_merge(dict(a), b))
Base.(:+){P,N}(a::SingleKet{P,N}, b::MultiKet{P,N}) = +(b, a)
Base.(:+){P,N}(a::MultiKet{P,N}, b::MultiKet{P,N}) = similar(b, add_merge(dict(a), dict(b)))

function Base.(:-){P,N,T,V}(a::SingleKet{P,N,T}, b::SingleKet{P,N,V})
    result = StateDict{N, promote_type(T,V)}()
    add_to_dict!(result, a.label, -a.coeff)
    add_to_dict!(result, b.label, -b.coeff)
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
tensor{P}(a::SingleKet{P}, b::SingleKet{P}) = SingleKet(P, a.coeff*b.coeff, tensor(a.label, b.label))
tensor{P}(a::SingleKet{P}, b::MultiKet{P}) = MultiKet(P, tensor_merge(a, dict(b)))
tensor{P}(a::MultiKet{P}, b::SingleKet{P}) = MultiKet(P, tensor_merge(dict(b), b))
tensor{P}(a::MultiKet{P}, b::MultiKet{P}) = MultiKet(P, tensor_merge(dict(a), dict(b)))
tensor(a::Bra, b::Bra) = tensor(a', b')'

Base.(:*)(a::Ket, b::Ket) = tensor(a,b)
Base.(:*)(a::Bra, b::Bra) = tensor(a,b)

#################
# Normalization #
#################
Base.norm(k::SingleKet) = abs(k.coeff)
Base.norm(k::MultiKet) = sqrt(sum(abs2, values(dict(k))))
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
    coeff = sqrt(k.label[i])*k.coeff 
    label = setindex(k.label, k.label[i] - 1, i)
    return SingleKet(P, coeff, label)
end

lower(b::Bra, i) = lower(b', i)'

raise{P}(state::DiracState{P,1}) = raise(state, 1)
raise(k::MultiKet, i) = similar(k, raisedict!(ladderdict(k), dict(k), i))

function raise(k::SingleKet, i)
    coeff = sqrt(k.label[i]+1)*k.coeff 
    label = setindex(k.label, k.label[i] + 1, i)
    return SingleKet(P, coeff, label)
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

switch(s::MultiKet, i, j) = maplabels(label->switch(label, i, j), s)
switch(s::SingleKet, i, j) = SingleKet(s.coeff, switch(s.label, i, j))
switch(b::Bra, i, j) = switch(b', i, j)'

permute(s::MultiKet, perm::Vector) = maplabels(label->permute(label,perm), s)
permute{P}(s::SingleKet{P}, perm::Vector) = SingleKet(P, s.coeff, permute(s.label, perm))
permute(b::Bra, perm::Vector) = permute(b', perm)'

filternz!(k::MultiKet) = (filter!(nzcoeff, dict(k)); return k)
filternz!(k::SingleKet) = k # no-op
filternz!(b::Bra) = (filternz!(b.kt); return b)

filternz(s::DiracState) = similar(s, filter(nzcoeff, dict(s)))

function represent{P}(kt::Ket{P}, labels)
    T = promote_type(return_type(P), eltype(kt))
    return T[bra(i) * kt for i in labels]
end

function represent{P}(kt::Ket{P}, labels...)
    iter = product(labels...)
    T = promote_type(return_type(P), eltype(kt))
    return T[bra(i...) * kt for i in iter]
end

represent(br::Bra, labels...) = represent(br', labels...)'

# should always be pure, of course,
# but makes a good sanity check function
purity(kt::Ket) = purity(kt*kt')
purity(br::Bra) = purity(br')

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