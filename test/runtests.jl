using QuDirac
using Base.Test

default_inner(UndefinedInner)

include("abstractinnertests.jl")
include("generaltests.jl")

default_inner(KroneckerDelta)

include("generaltests.jl")
include("orthotests.jl")
include("laddertests.jl")
