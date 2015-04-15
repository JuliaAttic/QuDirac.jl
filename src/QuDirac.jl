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
    
    immutable UndefinedInner <: AbstractInner end 
    immutable KroneckerDelta <: AbstractInner end
    
    abstract AbstractDirac{P<:AbstractInner,N}
    abstract DiracOp{P,N} <: AbstractDirac{P,N}
    abstract DiracState{P,N} <: AbstractDirac{P,N}

    abstract DiracScalar <: Number

    #############
    # Functions #
    #############
    QuBase.tensor() = error("Cannot call tensor function without arguments")

    ######################
    # Include Statements #
    ######################
    include("labels.jl")
    
    include("dirac/scalar.jl")
    include("dirac/state.jl")
    include("dirac/genericop.jl")
    include("dirac/outerproduct.jl")

    include("printfuncs.jl")
    include("dictfuncs.jl")
    include("mapfuncs.jl")
    
    ########
    # @drc #
    ########
    const ktpat = r"\|.*?\>"
    const brpat = r"\<.*?\|"

    ktrep(str) = "ket("*str[2:end-1]*")"
    brrep(str) = "bra("*str[2:end-1]*")"

    macro d_str(str)
        result = replace(str, brpat, brrep)
        result = replace(result, ktpat, ktrep)
        return parse(result)
    end

    export AbstractInner,
        UndefinedInner,
        KroneckerDelta,
        @d_str,
        AbstractDirac,
        DiracState,
        DiracOp,
        DiracScalar
end 

using QuBase
