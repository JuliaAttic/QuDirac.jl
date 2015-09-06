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
ket{P}(::Type{P}, x::StateLabel) = BasisKet(P, term(x))
ket{P}(::Type{P}, i...) = ket(P, StateLabel(i))

bra(i...) = BasisBra(ket(i...))

##############
# ctranspose #
##############
Base.ctranspose(br::BasisBra) = br.kt
Base.ctranspose(brs::BraSum) = br.kts
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
        coefftype{P,L,T}(::Type{($S){P,L,T}}) = T
        labeltype{P,L,T}(::Type{($S){P,L,T}}) = L
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

for S in (:BasisKet, :KetSum, :BasisBra, :BraSum)
    @eval begin
        Base.copy{P}(s::($S){P}) = ($S)(P, copy(data(s)))
    end
end

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

Base.get(kts::KetSum, x::StateLabel, default=any_zero(eltype(kts))) = get(data(kts), x, default)
Base.get(brs::BraSum, x::StateLabel, default=any_zero(eltype(brs))) = ifelse(haslabel(brs, x), br[x], default)
Base.get(kts::KetSum, x::Tuple, default=any_zero(eltype(kts))) = get(kts, StateLabel(x), default)
Base.get(brs::BraSum, x::Tuple, default=any_zero(eltype(brs))) = get(kts, StateLabel(x), default)

Base.isempty(s::StateSum) = isempty(data(s))
Base.empty!(s::StateSum) = (empty!(data(s)); return s)
Base.delete!(s::StateSum, x::StateLabel) = (delete!(data(s), x); return s)
Base.delete!(s::StateSum, x::Tuple) = delete!(s, StateLabel(x))

haslabel(s::StateSum, x::StateLabel) = haslabel(data(s), x)

#############
# Iteration #
#############
labels(s::StateSum) = labels(data(ts))
coeffs(s::StateSum) = coeffs(data(ts))

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

# Optimizations for KronDelta #
#-----------------------------#
inner(brs::BraSum{KronDelta}, kt::BasisKet{KronDelta}) = get(brs, label(kt)) * coeff(kt)
inner(br::BasisBra{KronDelta}, kts::KetSum{KronDelta}) = coeff(br) * get(kts, label(br))

function inner(br::BraSum{KronDelta}, kt::KetSum{KronDelta})
    if length(br) < length(kt)
        return ortho_inner(kt, br)
    else
        return ortho_inner(br, kt)
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
