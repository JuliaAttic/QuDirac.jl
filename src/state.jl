#################
# Ket/Bra Types #
#################
abstract Ket{P,T,L} <: DiracState{P}
abstract Bra{P,T,L} <: DiracState{P}

immutable SingleKet{P,T,L} <: Ket{P,T,L}
    data::SumTerm{StateLabel{L},T}
end

SingleKet{P,T,L}(::Type{P}, data::SumTerm{StateLabel{L},T}) = SingleKet{P,T,L}(data)
SingleKet{P}(::Type{P}, label::StateLabel, coeff) = SingleKet(P, SumTerm(label, coeff))

type KetSum{P,T,L} <: Ket{P,T,L}
    data::SumDict{StateLabel{L},T}
    nfactors::Int
    function KetSum(data::SumDict{StateLabel{L},T})        
        n = isempty(data) ? 0 : nfactors(first(keys(data)))
        return new(data, n)
    end
end

KetSum{P,T,L}(::Type{P}, data::SumDict{StateLabel{L},T}) = KetSum{P,T,L}(data)

immutable SingleBra{P,T,L} <: Bra{P,T,L}
    kt::SingleKet{P,T,L}
end

type BraSum{P,T,L} <: Bra{P,T,L}
    kt::KetSum{P,T,L}
end

typealias SingleState{P,T,L} Union(SingleKet{P,T,L}, SingleBra{P,T,L})
typealias StateSum{P,T,L} Union(KetSum{P,T,L}, BraSum{P,T,L})

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
coeff(br::SingleBra) = coeff(br')'
label(state::SingleState) = key(data(state))

Base.eltype{S<:DiracState}(::S) = eltype(S)
labeltype{S<:DiracState}(::S) = labeltype(S)

for t in (:Ket, :SingleKet, :KetSum, :Bra, :SingleBra, :BraSum)
    @eval begin
        Base.eltype{P,T,L}(::Type{($t){P,T,L}}) = T
        labeltype{P,T,L}(::Type{($t){P,T,L}}) = L
    end
end

nfactors(s::SingleState) = nfactors(label(s))
nfactors(kt::KetSum) = kt.nfactors
nfactors(br::BraSum) = nfactors(br')

########################
# Conversion/Promotion #
########################
Base.convert{P,T,L}(::Type{SingleKet{P,T,L}}, k::SingleKet{P,T,L}) = k 
Base.convert{P,T,L}(::Type{KetSum{P,T,L}}, k::KetSum{P,T,L}) = k
Base.convert{P,T,L}(::Type{SingleBra{P,T,L}}, b::SingleBra{P,T,L}) = b
Base.convert{P,T,L}(::Type{BraSum{P,T,L}}, b::BraSum{P,T,L}) = b

Base.convert{P,T,L}(::Type{SingleKet{P,T,L}}, kt::SingleKet) = SingleKet{P,N,T,L}(data(kt))
Base.convert{P,T,L}(::Type{KetSum{P,T,L}}, kt::SingleKet) = KetSum{P,N,T,L}(data(kt))
Base.convert{P,T,L}(::Type{KetSum{P,T,L}}, kt::KetSum) = KetSum{P,N,T,L}(data(kt))
Base.convert{P,T,L}(::Type{SingleBra{P,T,L}}, br::Bra) = SingleBra(convert(SingleKet{P,T,L}, br'))
Base.convert{P,T,L}(::Type{BraSum{P,T,L}}, br::Bra) = BraSum(convert(KetSum{P,T,L}, br'))

function Base.promote_rule{P,T1,T2,L1,L2}(::Type{SingleKet{P,T1,L1}}, ::Type{SingleKet{P,T2,L2}})
    return SingleKet{P,promote_type(T1,T2),promote_type(L1,L2)}
end

function Base.promote_rule{P,T1,T2,L1,L2}(::Type{KetSum{P,T1,L1}}, ::Type{KetSum{P,T2,L2}})
    return KetSum{P,promote_type(T1,T2),promote_type(L1,L2)}
end

function Base.promote_rule{P,T1,T2,L1,L2}(::Type{SingleKet{P,T1,L1}}, ::Type{KetSum{P,T2,L2}})
    return KetSum{P,promote_type(T1,T2),promote_type(L1,L2)}
end

function Base.promote_rule{P,T1,T2,L1,L2}(::Type{SingleBra{P,T1,L1}}, ::Type{SingleBra{P,T2,L2}})
    return SingleBra{P,promote_type(T1,T2),promote_type(L1,L2)}
end

function Base.promote_rule{P,T1,T2,L1,L2}(::Type{BraSum{P,T1,L1}}, ::Type{BraSum{P,T2,L2}})
    return BraSum{P,promote_type(T1,T2),promote_type(L1,L2)}
end

function Base.promote_rule{P,T1,T2,L1,L2}(::Type{SingleBra{P,T1,L1}}, ::Type{BraSum{P,T2,L2}})
    return BraSum{P,promote_type(T1,T2),promote_type(L1,L2)}
end

############################
# Hashing/Equality/Copying #
############################
const kt_hash = hash(Ket)
const br_hash = hash(Bra)

Base.hash{P}(kt::Ket{P}) = hash(P, hash(data(kt), kt_hash))
Base.hash{P}(br::Bra{P}) = hash(P, hash(data(br), br_hash))
Base.hash(state::DiracState, h::Uint64) = hash(hash(state), h)

Base.(:(==)){P}(a::Ket{P}, b::Ket{P}) = data(a) == data(b)
Base.(:(==)){P}(a::Bra{P}, b::Bra{P}) = data(a) == data(b)

Base.copy{P}(kt::SingleKet{P}) = SingleKet(P,copy(data(kt)))
Base.copy{P}(kt::KetSum{P}) = KetSum(P,copy(data(kt)))
Base.copy(br::Bra) = ctranspose(copy(br'))

Base.similar{P}(kt::Ket{P}) = KetSum(P,similar(data(kt)))
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

Base.isempty(state::DiracState) = isempty(data(state))
Base.empty!(state::DiracState) = (empty!(data(state)); return s)

#############
# Iteration #
#############
Base.start(state::SingleState) = false
Base.next(state::SingleState, i) = tuple(state, true)
Base.done(state::SingleState, i) = i

Base.start(state::StateSum) = start(data(state))

function Base.next{P,T,L}(kt::KetSum{P,T,L}, i)
    (lab,c), n = next(data(kt), i)
    return tuple(SingleKet{P,T,L}(SumTerm(lab, c)), n)
end

function Base.next(br::BraSum, i)
    kt, n = next(br', i)
    return tuple(kt', n)
end

Base.done(state::StateSum, i) = done(data(state), i)

Base.collect(state::DiracState) = [i for i in state]

#########
# inner #
#########
function execute_inner{P}(br::SingleBra{P}, kt::SingleKet{P})
    return coeff(br) * coeff(kt) * inner(P, label(br), label(kt))
end

function execute_inner{P}(br::SingleBra{P}, kt::KetSum{P})
    c = coeff(br)
    b = label(br)
    k0,v0 = first(data(kt))
    result = c * v0 * inner(P, b, k0)

    for (k,v) in drop(data(kt), 1)
        result += c * v * inner(P, b, k)
    end

    return result
end

function execute_inner{P}(br::BraSum{P}, kt::SingleKet{P})
    v = coeff(kt)
    k = label(kt)
    b0,c0 = first(data(br))
    result = c0' * v * inner(P, b0, k)

    for (b,c) in drop(data(br), 1)
        result += c' * v * inner(P, b, k)
    end

    return result
end

function execute_inner{P}(br::BraSum{P}, kt::KetSum{P})
    b0,c0 = first(data(br))
    k0,v0 = first(data(kt))
    result = c0' * v0 * inner(P, b0, k0)
    result -= result 

    for (b,c) in data(br), (k,v) in data(kt)
        result += c' * v * inner(P, b, k)
    end

    return result  
end

execute_inner(br::BraSum{KronDelta}, kt::SingleKet{KronDelta}) = get(br, label(kt)) * coeff(kt)
execute_inner(br::SingleBra{KronDelta}, kt::KetSum{KronDelta}) = get(kt, label(br)) * coeff(br)

function execute_inner(br::BraSum{KronDelta}, kt::KetSum{KronDelta})
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

##################
# execute_act_on #
##################
execute_act_on(kt::Ket, br::Bra, i) = act_on(kt', br', i)'

function execute_act_on{P}(br::SingleBra{P}, kt::SingleKet{P}, i)
    return (br * ket(P, label(kt)[i])) * SingleKet(P, except(label(kt), i), coeff(kt))
end

function execute_act_on{P}(br::Bra{P}, kt::Ket{P}, i)
    k0,v0 = first(data(kt))
    b0,c0 = first(data(br))
    result = act_on(SingleKet(P,b0,c0)', SingleKet(P,k0,v0), i)
    result -= result

    for (k,v) in data(kt), (b,c) in data(br)
        result = unsafe_add!(result, act_on(SingleKet(P,b,c)', SingleKet(P,k,v), i))
    end

    return result
end

function execute_act_on{P<:ProvidedInner,T,L}(br::Bra{P}, kt::Ket{P,T,L}, i)
    result = SumDict{StateLabel{L}, inner_rettype(br,kt)}()
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
Base.(:-){P}(kt::Ket{P}) = make_kt(P, -data(kt))
Base.(:-)(br::Bra) = ctranspose(-(br'))

Base.(:+){P}(a::Ket{P}, b::Ket{P}) = (@assert matching_nfactors(a, b); make_kt(P, data(a) + data(b)))
Base.(:-){P}(a::Ket{P}, b::Ket{P}) = (@assert matching_nfactors(a, b); make_kt(P, data(a) - data(b)))

Base.(:+)(a::Bra, b::Bra) = ctranspose(a' + b')
Base.(:-)(a::Bra, b::Bra) = ctranspose(a' - b')

# The below methods are unsafe because:
# 1. The first argument is potentially not mutated
# 2. No nfactors check is performed
unsafe_add!(a, b) = a + b
unsafe_add!{P,T,L}(a::KetSum{P,T,L}, b::Ket{P,T,L}) = (add!(data(a), data(b)); a)
unsafe_add!{P,T,L}(a::BraSum{P,T,L}, b::Bra{P,T,L}) = (add!(data(a), data(b)); a)

unsafe_sub!(a, b) = a - b
unsafe_sub!{P,T,L}(a::KetSum{P,T,L}, b::Ket{P,T,L}) = (sub!(data(a), data(b)); a)
unsafe_sub!{P,T,L}(a::BraSum{P,T,L}, b::Bra{P,T,L}) = (sub!(data(a), data(b)); a)

add!{P}(a::KetSum{P}, b::Ket{P}) = (@assert matching_nfactors(a, b); unsafe_add!(a, b))
add!{P}(a::BraSum{P}, b::Bra{P}) = (@assert matching_nfactors(a, b); unsafe_add!(a, b))

sub!{P}(a::KetSum{P}, b::Ket{P}) = (@assert matching_nfactors(a, b); unsafe_sub!(a, b))
sub!{P}(a::BraSum{P}, b::Bra{P}) = (@assert matching_nfactors(a, b); unsafe_sub!(a, b))

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
lower{P}(state::DiracState{P}) = (@assert nfactors(state) == 1; lower(state, 1))
raise{P}(state::DiracState{P}) = (@assert nfactors(state) == 1; raise(state, 1))

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
    @compat sizehint!(SumDict{StateLabel{labeltype(kt)}, T}(), length(kt))
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

represent(kt::Ket, bras) = [b * kt for b in bras]
represent(br::Bra, kets) = [br * k for k in kets]

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