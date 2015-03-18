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
        abstract AbstractOperator{S<:AbstractStructure}
        abstract AbstractState{S<:AbstractStructure}

        abstract DiracScalar <: Number

    #############
    # Functions #
    #############
        QuBase.tensor() = error("Cannot call tensor function without arguments")

    ######################
    # Include Statements #
    ######################
        include("helperfuncs/tuplefuncs.jl")
        include("helperfuncs/dictfuncs.jl")
        
        include("dirac/scalar.jl")
        include("dirac/state.jl")
        include("dirac/diracop.jl")
        include("dirac/projector.jl")

    export Ket,
        Bra,
        dualtype
end 

using QuBase
