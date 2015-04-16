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
    
    ##########
    # @d_str #
    ##########
    const ktpat = r"\|.*?\>"
    const brpat = r"\<.*?\|"

    ktrep(str) = "ket("*str[2:end-1]*")"
    brrep(str) = "bra("*str[2:end-1]*")"

    function prune_dirac(str)
        return replace(replace(str, brpat, brrep), ktpat, ktrep)
    end

    macro d_str(str)
        return esc(parse(prune_dirac(str)))
    end

    macro d_mstr(str)
        return esc(quote
            local s;
            for s in filter(x -> ! isempty(x), split($str, '\n'))
                eval(parse(QuDirac.prune_dirac(s)))
            end
        end)
    end

    export AbstractInner,
        UndefinedInner,
        KroneckerDelta,
        @d_str,
        @d_mstr,
        AbstractDirac,
        DiracState,
        DiracOp,
        DiracScalar
end 

using QuBase
