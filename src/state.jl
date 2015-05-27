#################
# Ket/Bra Types #
#################
abstract Ket{P,N,T,L} <: DiracState{P,N}
abstract Bra{P,N,T,L} <: DiracState{P,N}

immutable SingleKet{P,N,T,L} <: Ket{P,N,T,L}
    data::SumTerm{StateLabel{N,L},T}
end

SingleKet{P,N,T,L}(::Type{P}, data::SumTerm{StateLabel{N,L},T}) = SingleKet{P,N,T,L}(data)

type KetSum{P,N,T,L} <: Ket{P,N,T,L}
    data::SumDict{StateLabel{N,L},T}
end

KetSum{P,N,T,L}(::Type{P}, data::SumDict{StateLabel{N,L},T}) = KetSum{P,N,T,L}(data)

immutable SingleBra{P,N,T,L} <: Bra{P,N,T,L}
    kt::SingleKet{P,N,T,L}
end

type BraSum{P,N,T,L} <: Bra{P,N,T,L}
    kt::KetSum{P,N,T,L}
end

typealias SingleState{P,N,T,L} Union(SingleKet{P,N,T,L}, SingleBra{P,N,T,L})
typealias StateSum{P,N,T,L} Union(KetSum{P,N,T,L}, BraSum{P,N,T,L})

################
# Constructors #
################
make_kt{P}(::Type{P}, data::SumTerm) = SingleKet(P,data)
make_kt{P}(::Type{P}, data::SumDict) = KetSum(P,data)

ket{P}(::Type{P}, lbl::StateLabel) = SingleKet(P, SumTerm(lbl, 1))
ket{P}(::Type{P}, i...) = ket(P, StateLabel(i))

bra(i...) = SingleBra(ket(i...))

##############
# ctranspose #
##############
ctranspose(br::SingleBra) = br.kt
ctranspose(br::BraSum) = br.kt
ctranspose(kt::SingleKet) = SingleBra(kt)
ctranspose(kt::KetSum) = BraSum(kt)

######################
# Accessor Functions #
######################
data(kt::SingleKet) = kt.data
data(kt::KetSum) = kt.data
data(br::Bra) = data(br')

coeff(kt::SingleKet) = val(data(kt))
coeff(br::SingleBra) = val(data(br))'
label(state::SingleState) = key(data(state))

Base.eltype{S<:DiracState}(::S) = eltype(S)
labeltype{S<:DiracState}(::S) = labeltype(S)
nfactors{S<:DiracState}(::S) = nfactors(S)

for t in (:Ket, :SingleKet, :KetSum, :Bra, :SingleBra, :BraSum)
    @eval begin
        Base.eltype{P,N,T,L}(::Type{($t){P,N,T,L}}) = T
        labeltype{P,N,T,L}(::Type{($t){P,N,T,L}}) = L
        nfactors{P,N,T,L}(::Type{($t){P,N,T,L}}) = N
    end
end

########################
# Conversion/Promotion #
########################
issimilar{P,N,T,L}(::Ket{P,N,T,L}, ::Ket{P,N,T,L}) = true
issimilar{P,N,T,L}(::Bra{P,N,T,L}, ::Bra{P,N,T,L}) = true

Base.convert{P,N,T,L}(::Type{SingleKet{P,N,T,L}}, k::SingleKet{P,N,T,L}) = k 
Base.convert{P,N,T,L}(::Type{KetSum{P,N,T,L}}, k::KetSum{P,N,T,L}) = k
Base.convert{P,N,T,L}(::Type{SingleBra{P,N,T,L}}, b::SingleBra{P,N,T,L}) = b
Base.convert{P,N,T,L}(::Type{BraSum{P,N,T,L}}, b::BraSum{P,N,T,L}) = b

Base.convert{P,N,T,L}(::Type{SingleKet{P,N,T,L}}, kt::SingleKet) = SingleKet{P,N,T,L}(data(kt))
Base.convert{P,N,T,L}(::Type{KetSum{P,N,T,L}}, kt::SingleKet) = KetSum{P,N,T,L}(data(kt))
Base.convert{P,N,T,L}(::Type{KetSum{P,N,T,L}}, kt::KetSum) = KetSum{P,N,T,L}(data(kt))
Base.convert{P,N,T,L}(::Type{SingleBra{P,N,T,L}}, br::Bra) = SingleBra(convert(SingleKet{P,N,T,L}, br'))
Base.convert{P,N,T,L}(::Type{BraSum{P,N,T,L}}, br::Bra) = BraSum(convert(KetSum{P,N,T,L}, br'))

function Base.promote_rule{P,N,T1,T2,L1,L2}(::Type{SingleKet{P,N,T1,L1}}, ::Type{SingleKet{P,N,T2,L2}})
    return SingleKet{P,N,promote_type(T1,T2),label_promote(L1,L2)}
end

function Base.promote_rule{P,N,T1,T2,L1,L2}(::Type{KetSum{P,N,T1,L1}}, ::Type{KetSum{P,N,T2,L2}})
    return KetSum{P,N,promote_type(T1,T2),label_promote(L1,L2)}
end

function Base.promote_rule{P,N,T1,T2,L1,L2}(::Type{SingleKet{P,N,T1,L1}}, ::Type{KetSum{P,N,T2,L2}})
    return KetSum{P,N,promote_type(T1,T2),label_promote(L1,L2)}
end

function Base.promote_rule{P,N,T1,T2,L1,L2}(::Type{SingleBra{P,N,T1,L1}}, ::Type{SingleBra{P,N,T2,L2}})
    return SingleBra{P,N,promote_type(T1,T2),label_promote(L1,L2)}
end

function Base.promote_rule{P,N,T1,T2,L1,L2}(::Type{BraSum{P,N,T1,L1}}, ::Type{BraSum{P,N,T2,L2}})
    return BraSum{P,N,promote_type(T1,T2),label_promote(L1,L2)}
end

function Base.promote_rule{P,N,T1,T2,L1,L2}(::Type{SingleBra{P,N,T1,L1}}, ::Type{BraSum{P,N,T2,L2}})
    return BraSum{P,N,promote_type(T1,T2),label_promote(L1,L2)}
end

############################
# Hashing/Equality/Copying #
############################
const kt_hash = hash(Ket)
const br_hash = hash(Bra)

Base.hash{P}(kt::Ket{P}) = hash(P, hash(data(kt), kt_hash))
Base.hash{P}(br::Bra{P}) = hash(P, hash(data(br), br_hash))
Base.hash(state::DiracState, h::Uint64) = hash(hash(state), h)

Base.(:(==)){P,N}(a::Ket{P,N}, b::Ket{P,N}) = data(a) == data(b)
Base.(:(==)){P,N}(a::Bra{P,N}, b::Bra{P,N}) = data(a) == data(b)

Base.copy{P}(kt::SingleKet{P}) = SingleKet(P,copy(data(kt)))
Base.copy{P}(kt::KetSum{P}) = KetSum(P,copy(data(kt)))
Base.copy(br::Bra) = ctranspose(copy(br'))

Base.similar{P}(kt::SingleKet{P}) = SingleKet(P,similar(data(kt)))
Base.similar{P}(kt::KetSum{P}) = KetSum(P,similar(data(kt)))
Base.similar(br::Bra) = ctranspose(similar(br'))

#######################
# Dict-like Functions #
#######################
Base.length(state::DiracState) = length(data(state))

Base.getindex(s::DiracState, i) = s[StateLabel(i)]
Base.getindex(s::DiracState, i, j...) = s[StateLabel(i, j...)]
Base.getindex(kt::Ket, x::StateLabel) = getindex(data(kt), x)
Base.getindex(br::Bra, x::StateLabel) = getindex(br', x)'

Base.setindex!(s::StateSum, x, i) = setindex!(s, x, StateLabel(i))
Base.setindex!(s::StateSum, x, i, j...) = setindex!(s, x, StateLabel(i, j...))
Base.setindex!(kt::KetSum, x, y::StateLabel) = setindex!(data(kt), x, y)
Base.setindex!(br::BraSum, x, y::StateLabel) = setindex!(br', x', y)

Base.haskey(state::DiracState, x::StateLabel) = haskey(data(state), x)

Base.get(kt::Ket, x::StateLabel, default=predict_zero(eltype(kt))) = get(data(kt), x, default)
Base.get(br::Bra, x::StateLabel, default=predict_zero(eltype(br))) = haskey(br', x) ? br[x] : default

#############
# Iteration #
#############
Base.start(state::DiracState) = start(data(state))

Base.next(kt::Ket, i) = next(data(kt), i)

function Base.next(br::Bra, i)
    (k,v), n = next(br', i)
    return ((k,v'), n)
end

Base.done(state::DiracState, i) = done(data(state), i)

###########
# collect #
###########
Base.collect(kt::Ket) = collect(data(kt))

function Base.collect(br::Bra)
    arr = @compat Array(Tuple{StateLabel{nfactors(br), labeltype(br)}, eltype(br)}, length(br))
    return collect_pairs!(arr, br)
end

function collect_pairs!(result, b::Bra)
    i = 1
    for (k,v) in data(br)
        result[i] = (k, v')
        i += 1
    end
    return result
end

#########
# inner #
#########
function inner{P,N}(br::SingleBra{P,N}, kt::SingleKet{P,N})
    return coeff(br) * coeff(kt) * inner(P, label(br), label(kt))
end

function inner{P,N}(br::SingleBra{P,N}, kt::KetSum{P,N})
    c = coeff(br)
    b = label(br)
    (k0,v0) = first(data(kt))
    result = c * v0 * inner(P, b, k0)

    for (k,v) in drop(data(kt), 1)
        result += c * v * inner(P, b, k)
    end

    return result
end

function inner{P,N}(br::BraSum{P,N}, kt::SingleKet{P,N})
    v = coeff(kt)
    k = label(kt)
    (b0,c0) = first(data(br))
    result = c0' * v * inner(P, b0, k)

    for (b,c) in drop(data(br), 1)
        result += c' * v * inner(P, b, k)
    end

    return result
end

function inner{P,N}(br::BraSum{P,N}, kt::KetSum{P,N})
    (b0,c0) = first(data(br))
    (k0,v0) = first(data(kt))
    result = c0' * v0 * inner(P, b0, k0)
    result -= result 

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
    result = predict_zero(inner_rettype(a,b))
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
act_on{P<:ProvidedInner,T,L}(br::Bra{P,1}, kt::Ket{P,1,T,L}, i) = i==1 ? inner(br, kt) : throw(BoundsError())
act_on{P,T,L}(br::Bra{P,1}, kt::Ket{P,1,T,L}, i) = i==1 ? inner(br, kt) : throw(BoundsError())
act_on{P}(br::SingleBra{P,1}, kt::SingleKet{P,1}, i) = i==1 ? inner(br, kt) : throw(BoundsError())
act_on{P}(br::BraSum{P,1}, kt::SingleKet{P,1}, i) = i==1 ? inner(br, kt) : throw(BoundsError())

function act_on{P,N}(br::SingleBra{P,1}, kt::SingleKet{P,N}, i)
    return (br * ket(P, label(kt)[i])) * ket(P, except(label(kt), i))
end

function act_on{P,N}(br::Bra{P,1}, kt::Ket{P,N}, i)
    (k0,v0) = first(data(kt))
    (b0,c0) = first(data(br))
    result = act_on(c0*bra(P,b0), v0*ket(P,k0), i)
    result -= result

    for (k,v) in data(kt), (b,c) in data(br)
        tmp = act_on(c*bra(P,b), v*ket(P,k), i)
        if issimilar(result, tmp)
            add!(result, tmp)
        else
            result += tmp
        end
    end

    return result
end

function act_on{P<:ProvidedInner,N,T,L}(br::Bra{P,1}, kt::Ket{P,N,T,L}, i)
    result = SumDict{StateLabel{N-1, L}, inner_rettype(br,kt)}()
    return KetSum(P, act_on_dict!(result, br, kt, i))
end

function act_on_dict!{P<:ProvidedInner}(result::SumDict, br::BraSum{P}, kt::KetSum{P}, i)
    for (b,c) in data(br), (k,v) in data(kt)
        add_to_sum!(result, except(k,i), c'*v*P(b[1], k[i]))
    end
    return result
end

function act_on_dict!{P<:ProvidedInner}(result, br::SingleBra{P}, kt::KetSum{P}, i,)
    b = label(br)[1]
    c = coeff(br)
    for (k,v) in data(kt)
        add_to_sum!(result, except(k,i), c'*v*P(b, k[i]))
    end
    return result
end

function act_on_dict!{P<:ProvidedInner}(result, br::BraSum{P}, kt::SingleKet{P}, i,)
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
Base.scale!(kt::KetSum, c::Number) = (scale!(data(kt), c); return kt)
Base.scale!(c::Number, kt::KetSum) = scale!(kt, c)
Base.scale!(br::BraSum, c::Number) = (scale!(br', c'); return br)
Base.scale!(c::Number, br::BraSum) = scale!(br, c)

Base.scale{P}(kt::Ket{P}, c::Number) = make_kt(P, scale(data(kt), c))
Base.scale(c::Number, kt::Ket) = scale(kt, c)
Base.scale(br::Bra, c::Number) = scale(br', c')'
Base.scale(c::Number, br::Bra) = scale(br, c)

Base.(:*)(c::Number, state::DiracState) = scale(c, state)
Base.(:*)(state::DiracState, c::Number) = scale(state, c)
Base.(:/)(state::DiracState, c::Number) = scale(state, 1/c)

###########
# + and - #
###########
add!{P,N}(a::KetSum{P,N}, b::Ket{P,N}) = (add!(data(a), data(b)); return a)
add!{P,N}(a::BraSum{P,N}, b::Bra{P,N}) = (add!(data(a), data(b)); return a)

sub!{P,N}(a::KetSum{P,N}, b::Ket{P,N}) = (sub!(data(a), data(b)); return a)
sub!{P,N}(a::BraSum{P,N}, b::Bra{P,N}) = (sub!(data(a), data(b)); return a)

Base.(:-){P}(kt::Ket{P}) = make_kt(P, -data(kt))
Base.(:-)(br::Bra) = ctranspose(-(br'))

Base.(:+){P,N}(a::Ket{P,N}, b::Ket{P,N}) = make_kt(P, data(a) + data(b))
Base.(:-){P,N}(a::Ket{P,N}, b::Ket{P,N}) = make_kt(P, data(a) - data(b))

Base.(:+)(a::Bra, b::Bra) = ctranspose(a' + b')
Base.(:-)(a::Bra, b::Bra) = ctranspose(a' - b')

##########
# tensor #
##########
tensor{P}(a::Ket{P}, b::Ket{P}) = make_kt(P, tensor(data(a), data(b)))
tensor(a::Bra, b::Bra) = tensor(a', b')'

Base.(:*)(a::Ket, b::Ket) = tensor(a,b)
Base.(:*)(a::Bra, b::Bra) = tensor(a,b)

#################
# Normalization #
#################
Base.norm(kt::SingleKet) = abs(coeff(kt))
Base.norm(kt::KetSum) = sqrt(sum(abs2, values(data(kt))))
Base.norm(br::Bra) = norm(br')
normalize(state::DiracState) = (1/norm(state))*state
normalize!(state::DiracState) = scale!(1/norm(state), state)

####################
# Raising/Lowering #
####################
lower{P}(state::DiracState{P,1}) = lower(state, 1)
raise{P}(state::DiracState{P,1}) = raise(state, 1)

function lower{P}(kt::SingleKet{P}, i)
    lbl = setindex(label(kt), label(kt)[i] - 1, i)
    c = sqrt(label(kt)[i])*coeff(kt) 
    return SingleKet(P, SumTerm(lbl, c))
end

function raise{P}(kt::SingleKet{P}, i)
    lbl = setindex(label(kt), label(kt)[i] + 1, i)
    c = sqrt(label(kt)[i]+1)*coeff(kt) 
    return SingleKet(P, SumTerm(lbl, c))
end

function ladder_result(kt::Ket)
    T = promote_type(Float64, labeltype(kt), eltype(kt))
    @compat sizehint!(SumDict{StateLabel{nfactors(kt), labeltype(kt)}, T}(), length(kt))
end

lower{P}(kt::KetSum{P}, i) = KetSum(P, lowerdict!(ladder_result(kt), data(kt), i))
raise{P}(kt::KetSum{P}, i) = KetSum(P, raisedict!(ladder_result(kt), data(kt), i))

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

lower(br::Bra, i) = lower(br', i)'
raise(br::Bra, i) = raise(br', i)'

########################
# Misc. Math Functions #
########################
xsubspace(state::DiracState, x) = filter((k,v)->is_sum_x(k,x), state)
switch(state::DiracState, i, j) = maplabels(l->switch(l, i, j), state)
permute(state::DiracState, perm::Vector) = maplabels(l->permute(l, perm), state)

represent{P}(kt::Ket{P}, basis) = [bra(P, i) * kt for i in basis]

function represent{P}(kt::Ket{P}, basis...)
    prodbasis = product(basis...)
    return [bra(P, StateLabel(i)) * kt for i in prodbasis]
end

represent(br::Bra, basis...) = represent(br', basis...)'

purity(s::DiracState) = 1
purity(kt::Ket, i) = purity(kt*kt', i)
purity(br::Bra, i) = purity(br', i)

export Ket,
    SingleKet,
    KetSum,
    Bra,
    SingleBra,
    BraSum,
    ket,
    bra,
    represent,
    nfactors,
    xsubspace,
    switch,
    permute,
    purity,
    lower,
    raise,
    act_on,
    inner_eval,
    add!,
    sub!