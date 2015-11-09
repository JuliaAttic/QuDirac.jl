module DiracNotation

if VERSION < v"0.4"
    warn("DiracNotation only officially supports the v0.4 release of Julia. Your version of Julia is $VERSION.")
end

####################
# String Constants #
####################
const lang = "⟨"
const rang = "⟩"
const otimes = "⊗"
const vdots = "⋮"

##################
# Abstract Types #
##################
abstract InnerProduct
abstract AbstractDirac{P<:InnerProduct}

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
include("labels/label_types.jl")
include("labels/label_sums.jl")
include("labels/arithmetic.jl")
include("labels/functors.jl")
include("inner_products.jl")
include("states.jl")
include("outer_products.jl")
include("operators.jl")
include("functors.jl")

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
function default_inner_product{P<:InnerProduct}(::Type{P})
    DiracNotation.ket(label::StateLabel) = ket(P, label)
    DiracNotation.ket(items...) = ket(P, StateLabel(items))
    DiracNotation.bra(items...) = BasisBra(ket(items...))
    info("DiracNotation's default inner product type is currently $P.")
end

default_inner_product(KronDelta);

export innertype,
       labeltype,
       coefftype,
       nfactors,
       label,
       coeff,
       labels,
       coeffs,
       haslabel,
       tensor,
       ket,
       bra,
       act,
       raise,
       lower,
       switch,
       permute,
       normalize,
       normalize!

end # module DiracNotation
