module QuDirac
    
    using QuBase
    using Iterators

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
        abstract AbstractInner
        abstract Orthonormal <: AbstractInner
        
        abstract AbstractOperator{P<:AbstractInner}
        abstract AbstractState{P<:AbstractInner}

        abstract DiracScalar <: Number

    #############
    # Functions #
    #############
        QuBase.tensor() = error("Cannot call tensor function without arguments")

    ######################
    # Include Statements #
    ######################
        include("helperfuncs/arrfuncs.jl")
        include("helperfuncs/dictfuncs.jl")
        
        include("dirac/scalar.jl")
        include("dirac/state.jl")
        include("dirac/diracop.jl")
        include("dirac/projector.jl")

    export AbstractInner,
        Orthonormal
end 

using QuBase
