module QuDirac
    
    importall LabelSums
    using Iterators

    if VERSION < v"0.4-"
        warn("QuDirac v0.2 only officially supports the v0.4 release of Julia. Your version of Julia is $VERSION.")
    end
    
    ####################
    # String Constants #
    ####################
    const lang = "⟨"
    const rang = "⟩"
    const otimes = "⊗"
    const vdots ="⋮"

    ##################
    # Abstract Types #
    ##################
    abstract AbstractInner
    abstract AbstractDirac{P<:AbstractInner}

    #############
    # Functions #
    #############
    tensor() = error("Cannot call tensor function without arguments")
    tensor(s...) = reduce(tensor, s)

    innertype{P}(::AbstractDirac{P}) = P

    any_zero{T}(::Type{T}) = zero(T)
    any_zero(::Type{Any}) = false

    ######################
    # Include Statements #
    ######################
    include("DiracLabels.jl")
    include("inner.jl")
    include("DiracState.jl")

    #################
    # default_inner #
    #################
    # Julia doesn't recompile functions within other 
    # previously compiled functions, so we can't have a 
    # get_default_inner() or something like that.
    
    # Also, global optimization is poor, so we don't want
    # to use that either for most things. Thus, we go for
    # a function that straight-up redefines the default
    # constructors for the relevant objects. This is hacky,
    # but works for now, seeing as how only a few functions
    # actually "use" the default ptype.
    function default_inner{P<:AbstractInner}(::Type{P})
        QuDirac.ket(label::StateLabel) = ket(P, label)
        QuDirac.ket(items...) = ket(P, StateLabel(items))
        QuDirac.bra(items...) = BasisBra(ket(items...))
        info("QuDirac's default inner product type is currently $P.")
    end

    default_inner(KronDelta);

    export default_inner,
           innertype,
           labeltype,
           coefftype,
           label,
           coeff,
           labels,
           coeffs,
           haslabel,
           tensor,
           ket,
           bra

end # module QuDirac
