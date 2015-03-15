module QuDirac
    
    using QuBase
    using Iterators

    import Base: +, .+,
                 -, .-,
                 /, ./,
                 *, .*,
                 ^, .^

    ####################
    # String Constants #
    ####################
        const lang = "\u27E8"
        const rang = "\u27E9"
        const otimes = "\u2297"
        const vdots ="\u205E"

    ##################
    # Abstract Types #
    ##################
        # Various constructor methods in this repo allow an argument 
        # of type Type{BypassFlag} to be passed in in order to 
        # circumvent value precalculation/checking. This is useful for
        # conversion methods and the like, where you know the input 
        # has already been vetted elsewhere. Don't use this unless
        # you're sure of what you're doing, and don't export this.
        abstract BypassFlag

        # abstract DualType
        # abstract Ket <: DualType
        # abstract Bra <: DualType

        abstract AbstractOperator{S<:AbstractStructure} #<: Associative{(Tuple,Tuple),Number}
        abstract AbstractState{S<:AbstractStructure} #<: Associative{Tuple,Number}

        abstract DiracScalar <: Number

        # typealias AbstractKet{S<:AbstractStructure} AbstractState{Ket, S}
        # typealias AbstractBra{S<:AbstractStructure} AbstractState{Bra, S}

    #############
    # Functions #
    #############
        # dualtype{D,S}(::AbstractState{D,S}) = D

        QuBase.tensor() = error("Cannot call tensor function without arguments")
        # Base.ctranspose(::Type{Ket}) = Bra
        # Base.ctranspose(::Type{Bra}) = Ket
    
    ######################
    # Include Statements #
    ######################
        include("helperfuncs/tuplefuncs.jl")
        include("helperfuncs/dictfuncs.jl")
        
        include("dirac/scalar.jl")
        include("dirac/state.jl")
        # include("dirac/diracstate.jl")
        # include("dirac/diracop.jl")

    export Ket,
        Bra,
        dualtype
end 

using QuBase
