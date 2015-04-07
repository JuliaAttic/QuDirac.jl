using QuDirac
using Base.Test

@default_inner UndefinedInner

include("constest.jl")
include("abstractinnertests.jl")
include("generaltests.jl")

@default_inner Orthonormal

include("constest.jl")
include("orthotests.jl")
include("generaltests.jl")
include("lowertests.jl")
