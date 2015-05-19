Base.eltype(s::DiracState) = eltype(typeof(s))
labeltype(s::DiracState) = labeltype(typeof(s))

#######
# Ket #
#######
type Ket{P,N,S<:SumAssoc} <: DiracState{P,N}
    data::S
    Ket{L,T}(data::SumAssoc{StateLabel{N,L},T}) = new(data)
end

Ket{P,N,L,T}(::Type{P}, data::SumAssoc{StateLabel{N,L},T}) = Ket{P,N,typeof(data)}(data)

ket{P<:AbstractInner}(::Type{P}, label::StateLabel) = Ket(P, SumTerm(label, 1))
ket{P<:AbstractInner}(::Type{P}, args...) = ket(P, StateLabel(args))

data(k::Ket) = k.data

Base.convert{P,N,S<:SumAssoc}(::Type{Ket{P,N,S}}, k::Ket{P}) = Ket(P, convert(S, data(k)))
Base.convert{P,N,S<:SumAssoc}(::Type{Ket{P,N,S}}, k::Ket{P,N,S}) = k

Base.promote_rule{P,N,A,B}(::Type{Ket{P,N,A}}, ::Type{Ket{P,N,B}}) = Ket{P,N,promote_type(A,B)}

Base.eltype{P,N,S}(::Type{Ket{P,N,S}}) = eltype(S)
labeltype{P,N,S}(::Type{Ket{P,N,S}}) = labeltype(S)

Base.hash{P}(k::Ket{P}) = hash(data(k), hash(P))
Base.hash(k::Ket, h::Uint64) = hash(hash(k), h)
Base.copy{P}(k::Ket{P}) = Ket(P, copy(data(k)))
Base.(:(==)){P,N}(a::Ket{P,N}, b::Ket{P,N}) = data(a) == data(b)
Base.length(k::Ket) = length(data(k))

Base.getindex(k::Ket, x::StateLabel) = k.data[x]
Base.setindex!(k::Ket, x, y::StateLabel) = setindex!(data(k), x, y)

Base.haskey(k::Ket, x::StateLabel) = haskey(data(k), x)
Base.get(k::Ket, x::StateLabel, default=predict_zero(eltype(k))) = get(data(k), x, default)

Base.start(k::Ket) = start(data(k))
Base.next(k::Ket, i) = next(data(k), i)
Base.done(k::Ket, i) = done(data(k), i)
Base.collect(k::Ket) = collect(data(k))

#######
# Bra #
#######
type Bra{P,N,S} <: DiracState{P,N}
    kt::Ket{P,N,S}
end

const bra_hash = hash(Bra)

Bra{P,N,S}(kt::Ket{P,N,S}) = Bra{P,N,S}(kt)
bra(args...) = Bra(ket(args...))

data(b::Bra) = data(b.kt)

Base.convert{P,N,S}(::Type{Bra{P,N,S}}, b::Bra{P}) = Bra(convert(Ket{P,N,S}, b.kt))
Base.convert{P,N,S}(::Type{Bra{P,N,S}}, b::Bra{P,N,S}) = b

Base.promote_rule{P,N,A,B}(::Type{Bra{P,N,A}}, ::Type{Bra{P,N,B}}) = Bra{P,N,promote_type(A,B)}

Base.eltype{P,N,S}(::Type{Bra{P,N,S}}) = eltype(S)
labeltype{P,N,S}(::Type{Bra{P,N,S}}) = labeltype(S)

Base.hash(b::Bra) = hash(b.kt, bra_hash)
Base.hash(b::Bra, h::Uint64) = hash(hash(b), h)
Base.copy(b::Bra) = Bra(copy(b.kt))
Base.(:(==)){P,N}(a::Bra{P,N}, b::Bra{P,N}) = a.kt == b.kt
Base.length(b::Bra) = length(b.kt)

Base.getindex(b::Bra, x::StateLabel) = b.kt[x]'
Base.setindex!(b::Bra, x, y::StateLabel) = setindex!(b.kt, x', y)

Base.haskey(b::Bra, x::StateLabel) = haskey(b.kt, x)
Base.get(b::Bra, x::StateLabel, default=predict_zero(eltype(b))) = get(b.kt, x, default)

Base.start(b::Bra) = start(b.kt)

function Base.next(b::Bra, i)
    (k,v), n = next(b.kt, i)
    return ((k,v'), n)
end

Base.done(b::Bra, i) = done(b.kt, i)

Base.collect{P,N,S}(b::Bra{P,N,S}) = @compat collect_pairs!(Array(Tuple{labeltype(b), eltype(b)}, length(b)), b)

function collect_pairs!(result, b::Bra)
    i = 1
    for (k,v) in data(br)
        result[i] = (k, v')
        i += 1
    end
    return result
end

##############
# DiracState #
##############
Base.getindex(s::DiracState, i) = s[StateLabel(i)]
Base.getindex(s::DiracState, i, j...) = s[StateLabel(i, j...)]
Base.setindex!(s::DiracState, x, i) = setindex!(s, x, StateLabel(i))
Base.setindex!(s::DiracState, x, i, j...) = setindex!(s, x, StateLabel(i, j...))

################
# Type Aliases #
################
typealias KetSum{P,N,S<:SumDict} Ket{P,N,S}
typealias SingleKet{P,N,S<:SumTerm} Ket{P,N,S}
typealias BraSum{P,N,S<:SumDict} Bra{P,N,S}
typealias SingleBra{P,N,S<:SumTerm} Bra{P,N,S}

coeff(k::SingleKet) = val(data(k))
coeff(b::SingleBra) = val(data(b))'
label(k::SingleKet) = key(data(k))
label(b::SingleBra) = key(data(b))

##############
# ctranspose #
##############
Base.ctranspose(k::Ket) = Bra(k)
Base.ctranspose(b::Bra) = b.kt

#########
# inner #
#########
function inner{P,N}(br::SingleBra{P,N}, kt::SingleKet{P,N})
    return coeff(br) * coeff(kt) * inner(P, label(br), label(kt))
end

function inner{P,N}(br::SingleBra{P,N}, kt::KetSum{P,N})
    result = predict_zero(inner_rettype(br, kt))
    c = coeff(br)
    b = label(br)
    for (k,v) in data(kt)
        result += c * v * inner(P, b, k)
    end
    return result
end

function inner{P,N}(br::BraSum{P,N}, kt::SingleKet{P,N})
    result = predict_zero(inner_rettype(br, kt))
    v = coeff(kt)
    k = label(kt)
    for (b,c) in data(br)
        result += c' * v * inner(P, b, k)
    end
    return result
end

function inner{P,N}(br::BraSum{P,N}, kt::KetSum{P,N})
    result = predict_zero(inner_rettype(br, kt))
    for (b,c) in data(br), (k,v) in data(kt)
        result += c' * v * inner(P, b, k)
    end
    return result  
end

inner{N}(br::BraSum{KronDelta,N}, kt::SingleKet{KronDelta,N}) = get(br, label(kt)) * coeff(kt)
inner{N}(br::SingleBra{KronDelta,N}, kt::KetSum{KronDelta,N}) = get(kt, label(br)) * coeff(br)

function inner{N}(br::BraSum{KronDelta,N}, kt::KetSum{KronDelta,N})
    if length(br) < length(kt)
        return ortho_inner(kt, br)
    else
        return ortho_inner(br, kt)
    end
end

function ortho_inner(a::DiracState{KronDelta}, b::DiracState{KronDelta})
    result = predict_zero(inner_rettype(a, b))
    for l in keys(data(b))
        if haskey(a, l)
            result += a[l]*b[l]
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
# redundant definitions to resolve ambiguity warnings
act_on{P,L,T}(br::Bra{P,1}, kt::Ket{P,1,SumDict{StateLabel{1,L},T}}, i) = i==1 ? inner(br, kt) : throw(BoundsError())
act_on{P}(br::SingleBra{P,1}, kt::SingleKet{P,1}, i) = i==1 ? inner(br, kt) : throw(BoundsError())
act_on{P}(br::BraSum{P,1}, kt::SingleKet{P,1}, i) = i==1 ? inner(br, kt) : throw(BoundsError())

function act_on{P,N}(br::SingleBra{P,1}, kt::SingleKet{P,N}, i)
    return (br * ket(P, label(kt)[i])) * ket(P, except(label(kt), i))
end

function act_on{P,N,L,T}(br::Bra{P,1}, kt::Ket{P,N,SumDict{StateLabel{N,L},T}}, i)
    result = SumDict{StateLabel{N-1, L}, inner_rettype(br,kt)}()
    return Ket(P, act_on_dict!(result, br, kt, i))
end

function act_on_dict!{P}(result::SumDict, br::BraSum{P}, kt::KetSum{P}, i)
    for (b,c) in data(br), (k,v) in data(kt)
        add_to_sum!(result, except(k,i), c'*v*P(b[1], k[i]))
    end
    return result
end

function act_on_dict!{P}(result, br::SingleBra{P}, kt::KetSum{P}, i,)
    b = label(br)[1]
    c = coeff(br)
    for (k,v) in data(kt)
        add_to_sum!(result, except(k,i), c'*v*P(b, k[i]))
    end
    return result
end

function act_on_dict!{P}(result, br::BraSum{P}, kt::SingleKet{P}, i,)
    k = label(kt)[i]
    new_k = except(label(kt),i)
    v = coeff(kt)
    for (b,c) in data(br)
        add_to_sum!(result, new_k, c'*v*P(b[1], k))
    end
    return result
end

###########
# Scaling #
###########
Base.scale!(k::Ket, c::Number) = (scale!(data(k), c); return k)
Base.scale!(c::Number, k::Ket) = scale!(k, c)
Base.scale!(b::Bra, c::Number) = (scale!(b', c'); return b)
Base.scale!(c::Number, b::Bra) = scale!(b, c)

Base.scale{P}(k::Ket{P}, c::Number) = Ket(P, scale(data(k), c))
Base.scale(c::Number, k::Ket) = scale(k, c)
Base.scale(b::Bra, c::Number) = scale(b', c')'
Base.scale(c::Number, b::Bra) = scale(b, c)

Base.(:*)(c::Number, s::DiracState) = scale(c, s)
Base.(:*)(s::DiracState, c::Number) = scale(s, c)
Base.(:/)(s::DiracState, c::Number) = scale(s, 1/c)

###########
# + and - #
###########
Base.(:-){P}(k::Ket{P}) = Ket(P, -data(k))
Base.(:-)(b::Bra) = ctranspose(-(b'))

Base.(:+){P,N}(a::Ket{P,N}, b::Ket{P,N}) = Ket(P, data(a) + data(b))
Base.(:-){P,N}(a::Ket{P,N}, b::Ket{P,N}) = Ket(P, data(a) - data(b))

Base.(:+)(a::Bra, b::Bra) = ctranspose(a' + b')
Base.(:-)(a::Bra, b::Bra) = ctranspose(a' - b')

##########
# tensor #
##########
tensor{P}(a::Ket{P}, b::Ket{P}) = Ket(P, tensor(data(a), data(b)))
tensor(a::Bra, b::Bra) = tensor(a', b')'

Base.(:*)(a::Ket, b::Ket) = tensor(a,b)
Base.(:*)(a::Bra, b::Bra) = tensor(a,b)

#################
# Normalization #
#################
Base.norm(k::SingleKet) = abs(coeff(k))
Base.norm(k::KetSum) = sqrt(sum(abs2, values(data(k))))
Base.norm(b::Bra) = norm(b')
normalize(s::DiracState) = (1/norm(s))*s
normalize!(s::DiracState) = scale!(1/norm(s), s)

####################
# Raising/Lowering #
####################
lower{P}(s::DiracState{P,1}) = lower(s, 1)
raise{P}(s::DiracState{P,1}) = raise(s, 1)

function lower{P}(k::SingleKet{P}, i)
    l = setindex(label(k), label(k)[i] - 1, i)
    c = sqrt(label(k)[i])*coeff(k) 
    return Ket(P, SumTerm(l, c))
end

function raise{P}(k::SingleKet{P}, i)
    l = setindex(label(k), label(k)[i] + 1, i)
    c = sqrt(label(k)[i]+1)*coeff(k) 
    return Ket(P, SumTerm(l, c))
end

function ladder_result(k::Ket)
    T = promote_type(Float64, eltype(labeltype(k)), eltype(k))
    @compat sizehint!(SumDict{labeltype(k), T}(), length(k))
end

lower{P}(k::KetSum{P}, i) = Ket(P, lowerdict!(ladder_result(k), data(k), i))
raise{P}(k::KetSum{P}, i) = Ket(P, raisedict!(ladder_result(k), data(k), i))

function lowerdict!(result, d, i)
    for (k,v) in d
        add_to_sum!(result, setindex(k, k[i] - 1, i), sqrt(k[i])*v)
    end
    return result
end

function raisedict!(result, d, i)
    for (k,v) in d
        add_to_sum!(result, setindex(k, k[i] + 1, i), sqrt(k[i]+1)*v)
    end
    return result
end

lower(b::Bra, i) = lower(b', i)'
raise(b::Bra, i) = raise(b', i)'

########################
# Misc. Math Functions #
########################
nfactors{P,N}(::DiracState{P,N}) = N

xsubspace(s::DiracState, x) = filter((k,v)->is_sum_x(k,x), s)
switch(s::DiracState, i, j) = maplabels(l->switch(l, i, j), s)
permute(s::DiracState, perm::Vector) = maplabels(l->permute(l, perm), k)

filternz{P}(k::Ket{P}) = Ket(P, filternz(data(k)))
filternz(b::Bra) = filternz(b')'
filternz!(s::DiracState) = (filternz!(data(s)); return s)

represent{P}(kt::Ket{P}, basis) = [bra(P, i) * kt for i in basis]

function represent{P}(kt::Ket{P}, basis...)
    prodbasis = product(basis...)
    return [bra(P, StateLabel(i)) * kt for i in prodbasis]
end

represent(br::Bra, basis...) = represent(br', basis...)'

purity(s::DiracState) = 1
purity(kt::Ket, i) = purity(kt*kt', i)
purity(br::Bra, i) = purity(br', i)

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
    switch,
    switch!,
    permute,
    permute!,
    filternz,
    filternz!,
    purity,
    lower,
    raise,
    act_on,
    inner_eval