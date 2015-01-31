module QuDirac
    
    using QuBase
    using Iterators
    using DataStructures

    import QuBase: tensor,
        structure,
        checkcoeffs,
        nfactors,
        bases,
        coeffs,
        QuCoeffs,
        QuBra,
        QuKet

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

    ######################
    # Include Statements #
    ######################
        include("dirac.jl")
        include("statelabel.jl")
        include("diracstate.jl")
        include("scalar.jl")

end 

using QuBase
