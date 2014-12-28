import Base: filter,
    map

####################
# Type Definitions #
####################
    abstract AbstractLabelBasis{S<:AbstractStructure} <: AbstractFiniteBasis{S}

######################
# Include Statements #
######################
    include("fockbasis.jl")
    include("labelbasis.jl")

#############
# Functions #
#############
    structure{S}(::Type{AbstractLabelBasis{S}}) = S
    structure(::Type{AbstractLabelBasis}) = AbstractStructure

    function getstate{D<:DualType, S}(basis::AbstractLabelBasis{S}, 
                                      i, 
                                      ::Type{D}=Ket)
        return DiracState{D,S}(basis[i])
    end

    function getstate{D<:DualType}(basis::AbstractLabelBasis, 
                                   arr::AbstractArray, 
                                   ::Type{D}=Ket)
        return [getstate(basis, i, D) for i in arr]
    end

    filter{S}(f::Function, basis::AbstractLabelBasis{S}) = LabelBasis{S}(filter(f, labelvec(basis)), BypassFlag)
    map{S}(f::Function, basis::AbstractLabelBasis{S}) = LabelBasis{S}(map(f, labelvec(basis)))
    
    xsubspace(basis::AbstractLabelBasis, x::Int) = filter(s->sum(s)==x, basis)

export structure,
    getstate,
    xsubspace


