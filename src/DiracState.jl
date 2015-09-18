#################
# Ket/Bra Types #
#################

# Abstract Types #
#----------------#
abstract DiracState{P,L,T} <: AbstractDirac{P}
abstract AbstractKet{P,L,T} <: DiracState{P,L,T}
abstract AbstractBra{P,L,T} <: DiracState{P,L,T}

# Ket Types #
#-----------#
immutable BasisKet{P,L,T} <: AbstractKet{P,L,T}
    data::LabelTerm{StateLabel{L},T}
end

BasisKet{P,L,T}(::Type{P}, data::LabelTerm{StateLabel{L},T}) = BasisKet{P,L,T}(data)
BasisKet{P}(::Type{P}, label::StateLabel, coeff) = BasisKet(P, LabelTerm(label, coeff))

type KetSum{P,L,T} <: AbstractKet{P,L,T}
    data::LabelSum{StateLabel{L},T}
end

KetSum{P,L,T}(::Type{P}, data::LabelSum{StateLabel{L},T}) = KetSum{P,L,T}(data)

# Bra Types #
#-----------#
immutable BasisBra{P,L,T} <: AbstractBra{P,L,T}
    kt::BasisKet{P,L,T}
end

BasisBra(args...) = BasisBra(BasisKet(args...))

type BraSum{P,L,T} <: AbstractBra{P,L,T}
    kts::KetSum{P,L,T}
end

BraSum(args...) = BraSum(KetSum(args...))

# Type Aliases #
#--------------#
typealias BasisState{P,T,L} Union(BasisKet{P,T,L}, BasisBra{P,T,L})
typealias StateSum{P,T,L} Union(KetSum{P,T,L}, BraSum{P,T,L})

################
# Constructors #
################
ket{P}(::Type{P}, data::LabelTerm) = BasisKet(P, data)
ket{P}(::Type{P}, data::LabelSum) = KetSum(P, data)
ket{P}(::Type{P}, x::StateLabel) = BasisKet(P, term(x))
ket{P}(::Type{P}, i...) = ket(P, StateLabel(i))

bra(i...) = BasisBra(ket(i...))

##############
# ctranspose #
##############
Base.ctranspose(br::BasisBra) = br.kt
Base.ctranspose(brs::BraSum) = brs.kts
Base.ctranspose(kt::BasisKet) = BasisBra(kt)
Base.ctranspose(kts::KetSum) = BraSum(kts)

######################
# Property Functions #
######################
data(kt::BasisKet) = kt.data
data(kts::KetSum) = kts.data
data(br::AbstractBra) = data(br')

coeff(kt::BasisKet) = coeff(data(kt))
coeff(br::BasisBra) = coeff(br')'
label(s::BasisState) = label(data(s))

Base.eltype(s::StateSum) = eltype(typeof(s))
Base.eltype{P,L,T}(::Type{KetSum{P,L,T}}) = BasisKet{P,L,T}
Base.eltype{P,L,T}(::Type{BraSum{P,L,T}}) = BasisBra{P,L,T}

coefftype(s::DiracState) = coefftype(typeof(s))
labeltype(s::DiracState) = labeltype(typeof(s))
nfactors(s::DiracState) = nfactors(typeof(s))

for S in (:AbstractKet, :BasisKet, :KetSum, 
          :AbstractBra, :BasisBra, :BraSum)
    @eval begin
        innertype{P,L,T}(::Type{($S){P,L,T}}) = P
        labeltype{P,L,T}(::Type{($S){P,L,T}}) = L
        coefftype{P,L,T}(::Type{($S){P,L,T}}) = T
        nfactors{P,L,T}(::Type{($S){P,L,T}}) = nfactors(L)
    end
end

########################
# Conversion/Promotion #
########################
for S in (:BasisKet, :KetSum, :BasisBra, :BraSum)
    @eval begin
        if $S == :KetSum || $S == :BraSum
            Base.convert{P,L,T}(::Type{($S){P,L,T}}, s::($S)) = ($S)(P, convert(LabelSum{StateLabel{L},T}, data(s)))
        else
            Base.convert{P,L,T}(::Type{($S){P,L,T}}, s::($S)) = ($S)(P, convert(LabelTerm{StateLabel{L},T}, data(s)))
        end
        Base.convert{P,L,T}(::Type{($S){P,L,T}}, s::($S){P,L,T}) = s
    end
end

Base.promote_rule{P,L,A,B}(::Type{BasisKet{P,L,A}}, ::Type{BasisKet{P,L,B}}) = BasisKet{P,L,promote_type(A, B)}
Base.promote_rule{P,L,A,B}(::Type{BasisKet{P,L,A}}, ::Type{KetSum{P,L,B}}) = KetSum{P,L,promote_type(A, B)}
Base.promote_rule{P,L,A,B}(::Type{KetSum{P,L,A}}, ::Type{KetSum{P,L,B}}) = KetSum{P,L,promote_type(A, B)}
Base.promote_rule{P,L,A,B}(::Type{BasisBra{P,L,A}}, ::Type{BasisBra{P,L,B}}) = BasisBra{P,L,promote_type(A, B)}
Base.promote_rule{P,L,A,B}(::Type{BasisBra{P,L,A}}, ::Type{BraSum{P,L,B}}) = BraSum{P,L,promote_type(A, B)}
Base.promote_rule{P,L,A,B}(::Type{BraSum{P,L,A}}, ::Type{BraSum{P,L,B}}) = BraSum{P,L,promote_type(A, B)}

############################
# Hashing/Equality/Copying #
############################
const kt_hash = hash(AbstractKet)
const br_hash = hash(AbstractBra)

Base.hash{P,L}(kt::AbstractKet{P,L}) = hash(L, hash(P, hash(data(kt), kt_hash)))
Base.hash{P,L}(br::AbstractBra{P,L}) = hash(L, hash(P, hash(data(br), br_hash)))
Base.hash(s::DiracState, h::Uint64) = hash(hash(s), h)

Base.(:(==)){P,L}(a::AbstractKet{P,L}, b::AbstractKet{P,L}) = data(a) == data(b)
Base.(:(==)){P,L}(a::AbstractBra{P,L}, b::AbstractBra{P,L}) = data(a) == data(b)

Base.copy(s::BasisState) = s
Base.copy{P}(kts::KetSum{P}) = KetSum(P, data(kts))
Base.copy{P}(brs::BraSum{P}) = BraSum(P, data(brs))

#######################
# Dict-like Functions #
#######################
Base.length(s::StateSum) = length(data(s))

Base.getindex(s::StateSum, i...) = s[StateLabel(i...)]
Base.getindex(kts::KetSum, x::StateLabel) = getindex(data(kts), x)
Base.getindex(brs::BraSum, x::StateLabel) = getindex(brs', x)'

Base.setindex!(s::StateSum, x, i...) = setindex!(s, x, StateLabel(i...))
Base.setindex!(kts::KetSum, x, y::StateLabel) = setindex!(data(kts), x, y)
Base.setindex!(brs::BraSum, x, y::StateLabel) = setindex!(brs', x', y)

Base.get(kts::KetSum, x::StateLabel, default=any_zero(coefftype(kts))) = get(data(kts), x, default)
Base.get(brs::BraSum, x::StateLabel, default=any_zero(coefftype(brs))) = ifelse(haslabel(brs, x), brs[x], default)
Base.get(kts::KetSum, x::Tuple, default=any_zero(coefftype(kts))) = get(kts, StateLabel(x), default)
Base.get(brs::BraSum, x::Tuple, default=any_zero(coefftype(brs))) = get(kts, StateLabel(x), default)

Base.isempty(s::StateSum) = isempty(data(s))
Base.empty!(s::StateSum) = (empty!(data(s)); return s)
Base.delete!(s::StateSum, x::StateLabel) = (delete!(data(s), x); return s)
Base.delete!(s::StateSum, x::Tuple) = delete!(s, StateLabel(x))

LabelSums.haslabel(s::StateSum, x::StateLabel) = haslabel(data(s), x)

#############
# Iteration #
#############
Base.start(s::StateSum) = start(data(s))

function Base.next{P,L,T}(kt::KetSum{P,L,T}, i)
    item, n = next(data(kt), i)
    return tuple(BasisKet{P,L,T}(item), n)
end

function Base.next{P,L,T}(br::BraSum{P,L,T}, i)
    item, n = next(data(br), i)
    return tuple(BasisBra{P,L,T}(item), n)
end

Base.done(s::StateSum, i) = done(data(s), i)

Base.collect(s::StateSum) = [i for i in s]

#########
# inner #
#########
Base.(:*)(br::AbstractBra, kt::AbstractKet) = inner(br, kt)

# Generalized inner #
#-------------------#
function inner{P}(br::BasisBra{P}, kt::BasisKet{P})
    return coeff(br) * coeff(kt) * inner(P, label(br), label(kt))
end

function inner{P}(brs::BraSum{P}, kt::BasisKet{P})
    result = inner(first(brs), kt)
    result -= result

    for br in brs
        result += inner(br, kt)
    end

    return result
end

function inner{P}(br::BasisBra{P}, kts::KetSum{P})
    result = inner(br, first(kt))
    result -= result

    for kt in kts
        result += inner(br, kt)
    end

    return result
end

function inner{P}(brs::BraSum{P}, kts::KetSum{P})
    result = inner(first(brs), first(kts))
    result -= result

    for br in brs, kt in kts
        result += inner(br, kt)
    end

    return result 
end

# Optimized inner for KronDelta #
#-------------------------------#
inner(brs::BraSum{KronDelta}, kt::BasisKet{KronDelta}) = get(brs, label(kt)) * coeff(kt)
inner(br::BasisBra{KronDelta}, kts::KetSum{KronDelta}) = coeff(br) * get(kts, label(br))

function inner(brs::BraSum{KronDelta}, kts::KetSum{KronDelta})
    if length(brs) < length(kts)
        return ortho_inner(kts, brs)
    else
        return ortho_inner(brs, kts)
    end
end

function ortho_inner(long_state::DiracState{KronDelta}, short_state::DiracState{KronDelta})
    T = promote_type(coefftype(long_state), coefftype(short_state), rettype(KronDelta))
    result = any_zero(T)
    for basis_state in short_state
        result += get(long_state, label(basis_state)) * coeff(basis_state)
    end
    return result
end

#######
# act #
#######
act{i}(kt::AbstractKet, br::AbstractBra, idx::Type{Val{i}}) = act(kt', br', idx)'

function act{P,i}(br::BasisBra{P}, kt::BasisKet{P}, idx::Type{Val{i}})
    return inner(br, ket(P, label(kt)[idx])) * BasisKet(P, except(label(kt), idx), coeff(kt))
end

function act{P,i}(br::BasisBra{P}, kts::KetSum{P},  idx::Type{Val{i}})
    result = act(br, first(kts), idx)
    result -= result

    for kt in kts
        add!(result, act(br, kt, idx))
    end

    return result 
end

function act{P,i}(brs::BraSum{P}, kt::BasisKet{P},  idx::Type{Val{i}})
    result = act(first(brs), kt, idx)
    result -= result

    for br in brs
        add!(result, act(br, kt, idx))
    end

    return result 
end

function act{P,i}(brs::BraSum{P}, kts::KetSum{P},  idx::Type{Val{i}})
    result = act(first(brs), first(kts), idx)
    result -= result

    for br in brs, kt in kts
        add!(result, act(br, kt, idx))
    end

    return result 
end

##############
# Arithmetic #
##############

# Multiplication #
#----------------#
Base.scale!(kts::KetSum, c::Number) = (scale!(data(kts), c); return kts)
Base.scale!(c::Number, kts::KetSum) = scale!(kts, c)
Base.scale!(brs::BraSum, c::Number) = (scale!(brs', c'); return brs)
Base.scale!(c::Number, brs::BraSum) = scale!(brs, c)

Base.scale(kt::AbstractKet, c::Number) = ket(innertype(kt), scale(data(kt), c))
Base.scale(c::Number, kt::AbstractKet) = scale(kt, c)
Base.scale(br::AbstractBra, c::Number) = scale(br', c')'
Base.scale(c::Number, br::AbstractBra) = scale(br, c)

Base.(:*)(c::Number, state::DiracState) = scale(c, state)
Base.(:*)(state::DiracState, c::Number) = scale(state, c)
Base.(:/)(state::DiracState, c::Number) = scale(state, inv(c))

tensor{P}(a::AbstractKet{P}, b::AbstractKet{P}) = ket(P, data(a)*data(b))
tensor(a::AbstractBra, b::AbstractBra) = tensor(a', b')'

Base.(:*)(a::AbstractKet, b::AbstractKet) = tensor(a, b)
Base.(:*)(a::AbstractBra, b::AbstractBra) = tensor(a, b)

# Addition/Subtraction #
#----------------------#
Base.(:-){P}(kt::AbstractKet{P}) = ket(P, -data(kt))
Base.(:-)(br::AbstractBra) = ctranspose(-(br'))

Base.(:+){P,L}(a::AbstractKet{P,L}, b::AbstractKet{P,L}) = ket(P, data(a) + data(b))
Base.(:-){P,L}(a::AbstractKet{P,L}, b::AbstractKet{P,L}) = ket(P, data(a) - data(b))

Base.(:+)(a::AbstractBra, b::AbstractBra) = ctranspose(a' + b')
Base.(:-)(a::AbstractBra, b::AbstractBra) = ctranspose(a' - b')

add!{P,L}(a::KetSum{P,L}, b::AbstractKet{P,L}) = (add!(data(a), data(b)); a)
add!{P,L}(a::BraSum{P,L}, b::AbstractBra{P,L}) = (add!(data(a), data(b)); a)

sub!{P,L}(a::KetSum{P,L}, b::AbstractKet{P,L}) = (sub!(data(a), data(b)); a)
sub!{P,L}(a::BraSum{P,L}, b::AbstractBra{P,L}) = (sub!(data(a), data(b)); a)

#################
# Normalization #
#################
norm_term(kt::BasisKet) = abs2(coeff(kt))

Base.norm(kt::BasisKet) = abs(coeff(kt))
Base.norm(kts::KetSum) = sqrt(sum(norm_term, kts))
Base.norm(br::AbstractBra) = norm(br')

normalize(state::DiracState) = inv(norm(state))*state
normalize!(state::DiracState) = scale!(inv(norm(state)), state)

####################
# Raising/Lowering #
####################
lower(state::DiracState) = lower(state, Val{1})
raise(state::DiracState) = raise(state, Val{1})

lower{i}(br::AbstractBra, idx::Type{Val{i}}) = lower(br', idx)'
raise{i}(br::AbstractBra, idx::Type{Val{i}}) = raise(br', idx)'

function lower{P,i}(kt::BasisKet{P}, idx::Type{Val{i}})
    idx_lbl = label(kt)[idx]
    lbl = setindex(label(kt),  idx_lbl - 1, idx)
    c = sqrt(idx_lbl)*coeff(kt) 
    return BasisKet(P, lbl, c)
end

function raise{P,i}(kt::BasisKet{P}, idx::Type{Val{i}})
    idx_lbl = label(kt)[idx] + 1
    lbl = setindex(label(kt),  idx_lbl, idx)
    c = sqrt(idx_lbl)*coeff(kt)
    return BasisKet(P, lbl, c)
end

for f in (:raise, :lower)
    @eval function ($f){i}(kts::KetSum, idx::Type{Val{i}})
        result = ($f)(first(kts), idx)
        result -= result

        for kt in kts
            add!(result, ($f)(kt, idx))
        end

        return result
    end
end

########################
# Misc. Math Functions #
########################
label_sums_to_x(state::BasisState, x) = sum(label(state)) == x

xsubspace(state::StateSum, x) = filter(s->label_sums_to_x(s, x), state)
xsubspace!(state::StateSum, x) = filter!(s->label_sums_to_x(s, x), state)

function switch{i,j}(state::StateSum, idxi::Type{Val{i}}, idxj::Type{Val{j}})
    P = innertype(state)
    return map(s->ket(P, switch(label(s), idxi, idxj), coeff(s)), state)
end

function switch{i,j}(state::BasisState, idxi::Type{Val{i}}, idxj::Type{Val{j}})
    t = LabelTerm(switch(label(state), idxi, idxj), coeff(state))
    return ket(innertype(state), t)
end

function permute{T<:Tuple}(state::StateSum, perm::Type{T})
    return map(s->permute(label(s), perm), state)
end

function permute{T<:Tuple}(state::BasisState, perm::Type{T})
    t = LabelTerm(permute(label(state), perm), coeff(state))
    return ket(innertype(state), t)
end

represent(kt::AbstractKet, bras) = [b * kt for b in bras]
represent(br::AbstractBra, kets) = [br * k for k in kets]

purity(kt::AbstractKet, i) = purity(kt*kt', i)
purity(br::AbstractBra, i) = purity(br', i)
