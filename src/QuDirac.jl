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
        
        abstract AbstractOperator{P<:AbstractInner,N}
        abstract AbstractState{P<:AbstractInner,N}

        abstract DiracScalar <: Number

    #############
    # Functions #
    #############
        QuBase.tensor() = error("Cannot call tensor function without arguments")

    ###########
    # Factors #
    ###########
        # could use Val{N} in julia v0.4,
        # but we're targeting v0.3...
        immutable Factors{N} end

        Base.copy{N}(::Factors{N}) = Factors{N}()
        Base.(:+){A,B}(::Factors{A}, ::Factors{B}) = Factors{A+B}()

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
