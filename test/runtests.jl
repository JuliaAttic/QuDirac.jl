using QuDirac
using Base.Test

QuDirac.set_default_inner(UndefinedInner())

include("constest.jl")
include("abstractinnertests.jl")
include("generaltests.jl")

QuDirac.set_default_inner(Orthonormal())

include("constest.jl")
include("orthotests.jl")
include("generaltests.jl")
include("lowertests.jl")
