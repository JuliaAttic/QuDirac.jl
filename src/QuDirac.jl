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
    immutable Orthonormal <: AbstractInner end
    
    abstract AbstractDirac{P,N}
    abstract DiracOp{P<:AbstractInner,N} <: AbstractDirac{P,N}
    abstract DiracState{P<:AbstractInner,N} <: AbstractDirac{P,N}

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
    Base.(:-){A,B}(::Factors{A}, ::Factors{B}) = Factors{A-B}()
    decr{N}(::Factors{N}) = Factors{N-1}()

    ######################
    # Include Statements #
    ######################
    include("dirac/scalar.jl")
    include("dirac/state.jl")
    include("dirac/genericop.jl")
    include("dirac/outerproduct.jl")

    include("miscfuncs.jl")
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
        Orthonormal,
        @d_str,
        AbstractDirac,
        DiracState,
        DiracOp,
        DiracScalar
end 

using QuBase
