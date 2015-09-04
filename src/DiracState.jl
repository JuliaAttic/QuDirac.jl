#################
# Ket/Bra Types #
#################

# Abstract Types #
#----------------#
abstract AbstractKet{P,L,T} <: DiracState{P,L,T}
abstract AbstractBra{P,L,T} <: DiracState{P,L,T}

# Ket Types #
#-----------#
immutable BasisKet{P,L,T} <: AbstractKet{P,L,T}
    data::Term{StateLabel{L},T}
end

BasisKet{P,L,T}(::Type{P}, data::Term{StateLabel{L},T}) = BasisKet{P,L,T}(data)

type KetSum{P,L,T} <: AbstractKet{P,L,T}
    data::TermSum{StateLabel{L},T}
end

KetSum{P,L,T}(::Type{P}, data::TermSum{StateLabel{L},T}) = KetSum{P,L,T}(data)

# Bra Types #
#-----------#
immutable BasisBra{P,L,T} <: AbstractBra{P,L,T}
    kt::BasisKet{P,L,T}
end

BasisBra(args...) = BasisBra(BasisKet(args...))

type BraSum{P,L,T} <: AbstractBra{P,L,T}
    kt::KetSum{P,L,T}
end

BraSum(args...) = BraSum(KetSum(args...))

# Type Aliases #
#--------------#
typealias BasisState{P,T,L} Union(BasisKet{P,T,L}, BasisBra{P,T,L})
typealias StateSum{P,T,L} Union(KetSum{P,T,L}, BraSum{P,T,L})

################
# Constructors #
################
ket{P}(::Type{P}, x::StateLabel) = BasisKet(P, Term(x))
ket{P}(::Type{P}, i...) = ket(P, StateLabel(i))

bra(i...) = BasisBra(ket(i...))

##############
# ctranspose #
##############
ctranspose(br::BasisBra) = br.kt
ctranspose(br::BraSum) = br.kt
ctranspose(kt::BasisKet) = BasisBra(kt)
ctranspose(kt::KetSum) = BraSum(kt)

######################
# Property Functions #
######################
data(kt::BasisKet) = kt.data
data(kt::KetSum) = kt.data
data(br::AbstractBra) = data(br')

Base.eltype(s::DiracState) = eltype(typeof(s))
labeltype(s::DiracState) = labeltype(typeof(s))
nfactors(s::DiracState) = nfactors(typeof(s))

for S in (:AbstractKet, :BasisKet, :KetSum, 
          :AbstractBra, :BasisBra, :BraSum)
    @eval begin
        Base.eltype{P,L,T}(::Type{($S){P,L,T}}) = T
        labeltype{P,L,T}(::Type{($S){P,L,T}}) = L
        nfactors{P,L,T}(::Type{($S){P,L,T}}) = nfactors(L)
    end
end

########################
# Conversion/Promotion #
########################
for S in (:BasisKet, :KetSum, :BasisBra, :BraSum)
    @eval begin
        if S == :KetSum || S == :BraSum
            Base.convert{P,L,T}(::Type{($S){P,L,T}}, s::($S)) = ($S)(P, convert(TermSum{StateLabel{L},T}, data(s)))
        else
            Base.convert{P,L,T}(::Type{($S){P,L,T}}, s::($S)) = ($S)(P, convert(Term{StateLabel{L},T}, data(s)))
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
const kt_hash = hash(Ket)
const br_hash = hash(Bra)

Base.hash{P,L}(kt::AbstractKet{P,L}) = hash(L, hash(P, hash(data(kt), kt_hash)))
Base.hash{P,L}(br::AbstractBra{P,L}) = hash(L, hash(P, hash(data(br), br_hash)))
Base.hash(s::DiracState, h::Uint64) = hash(hash(s), h)

Base.(:(==)){P,L}(a::Ket{P,L}, b::Ket{P,L}) = data(a) == data(b)
Base.(:(==)){P,L}(a::Bra{P,L}, b::Bra{P,L}) = data(a) == data(b)

for S in (:BasisKet, :KetSum, :BasisBra, :BraSum)
    @eval begin
        Base.copy{P}(s::($S){P}) = ($S)(P, copy(data(s)))
    end
end

#######################
# Dict-like Functions #
#######################
Base.length(s::DiracState) = length(data(s))

Base.getindex(s::DiracState, i...) = s[StateLabel(i...)]
Base.getindex(kt::AbstractKet, x::StateLabel) = getindex(data(kt), x)
Base.getindex(br::AbstractBra, x::StateLabel) = getindex(br', x)'

Base.setindex!(s::StateSum, x, i...) = setindex!(s, x, StateLabel(i...))
Base.setindex!(kt::KetSum, x, y::StateLabel) = setindex!(data(kt), x, y)
Base.setindex!(br::BraSum, x, y::StateLabel) = setindex!(br', x', y)

Base.haskey(s::DiracState, x::StateLabel) = haskey(data(s), x)

Base.get(kt::AbstractKet, x::StateLabel, default=any_zero(eltype(kt))) = get(data(kt), x, default)
Base.get(br::AbstractBra, x::StateLabel, default=any_zero(eltype(br))) = haskey(br', x) ? br[x] : default
Base.get(kt::AbstractKet, x::Tuple, default=any_zero(eltype(kt))) = get(kt, StateLabel(x), default)
Base.get(br::AbstractBra, x::Tuple, default=any_zero(eltype(br))) = get(kt, StateLabel(x), default)

Base.isempty(s::StateSum) = isempty(data(s))
Base.empty!(s::StateSum) = (empty!(data(s)); return s)
Base.delete!(s::StateSum, x::StateLabel) = (delete!(data(s), x); return s)
Base.delete!(s::StateSum, x::Tuple) = delete!(s, StateLabel(x))

#############
# Iteration #
#############
Base.start(s::BasisState) = false
Base.next(s::BasisState, i) = tuple(s, true)
Base.done(s::BasisState, i) = i

Base.start(s::StateSum) = start(data(s))

function Base.next{P,L,T}(kt::KetSum{P,L,T}, i)
    (x, c), n = next(data(kt), i)
    return tuple(BasisKet{P,L,T}(Term(x, c)), n)
end

function Base.next{P,L,T}(br::BraSum{P,L,T}, i)
    (x, c), n = next(data(br), i)
    return tuple(BasisBra{P,L,T}(Term(x, c)), n)
end

Base.done(s::StateSum, i) = done(data(s), i)

Base.collect(s::DiracState) = [i for i in s]

