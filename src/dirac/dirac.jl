import Base: ctranspose

####################
# Type Definitions #
####################
    abstract DualType
    abstract Ket <: DualType
    abstract Bra <: DualType

    abstract AbstractDirac{S<:AbstractStructure} <: AbstractQuantum{S}
    abstract AbstractOperator{S<:AbstractStructure} <: AbstractDirac{S}
    abstract AbstractState{D<:DualType, S<:AbstractStructure} <: AbstractDirac{S}

    abstract DiracScalar <: Number

    typealias AbstractKet{S<:AbstractStructure} AbstractState{Ket, S}
    typealias AbstractBra{S<:AbstractStructure} AbstractState{Bra, S}

#############
# Functions #
#############
    structure{D,S}(::Type{AbstractState{D,S}}) = S
    structure(::Type{AbstractState}) = AbstractStructure
    structure{S}(::Type{AbstractDirac{S}}) = S
    structure(::Type{AbstractDirac}) = AbstractStructure
    structure{S}(::Type{AbstractOperator{S}}) = S
    structure(::Type{AbstractOperator}) = AbstractStructure

    dualtype{D,S}(::AbstractState{D,S}) = D
    dualtype{D,S}(::Type{AbstractState{D,S}}) = D
    dualtype(::Type{AbstractState}) = DualType

    ctranspose(::Type{Ket}) = Bra
    ctranspose(::Type{Bra}) = Ket

######################
# Include Statements #
######################
    include("statelabel.jl")
    include("diracstates.jl")
    include("diracoperators.jl")
    include("scalar.jl")

export Ket,
    Bra,
    AbstractDirac,
    AbstractOperator,
    AbstractState,
    AbstractKet,
    AbstractBra,
    dualtype,
    structure
