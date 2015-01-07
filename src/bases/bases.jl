import Base: filter,
    map,
    find

####################
# Type Definitions #
####################
    # All B<:AbstractLabelBasis types should implement the following:
    # 
    #   getindex(basis::B, i) -> the StateLabel at index `i`
    #   getindex(basis::B, label::StateLabel) -> the index at which `label` resides in `basis`
    #   in(label::StateLabel, basis::B) -> checks if the given `label` is an element of `basis`
    #   labelvec(basis::B) -> returns a Vector{StateLabel} of the labels in this basis.
    #   samelabels(a::B, b::B) -> returns true if the labels of `a` are the same as 
    #                             those in `b`, and in the same order. This
    #                             should be implmented to run in constant time/low-input 
    #                             linear time if at all possible, as this function is used 
    #                             by DiracArrays to check whether or not an array 
    #                             operation can be performed solely with coefficients or 
    #                             requires the use of bases.

    abstract AbstractLabelBasis{S<:AbstractStructure,N} <: AbstractFiniteBasis{S}

######################
# Include Statements #
######################
    include("fockbasis.jl")
    include("labelbasis.jl")

#############
# Functions #
#############
    structure{S}(::Type{AbstractLabelBasis{S}}) = S
    structure{S,N}(::Type{AbstractLabelBasis{S,N}}) = S

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

    filter{S}(f::Function, basis::AbstractLabelBasis{S}) = LabelBasis(S, filter(f, labelvec(basis)))
    map{S}(f::Function, basis::AbstractLabelBasis{S}) = LabelBasis(S, map(f, labelvec(basis)))
    find{S}(f::Function, basis::AbstractLabelBasis{S}) = find(f, labelvec(basis))

    xsubspace(basis::AbstractLabelBasis, x::Int) = filter(s->sum(s)==x, basis)

export AbstractLabelBasis,
    structure,
    getstate,
    xsubspace


