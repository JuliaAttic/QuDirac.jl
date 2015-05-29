using QuDirac
using Base.Test

default_inner(UndefInner)

include("abstractinnertests.jl")
include("generaltests.jl")

default_inner(KronDelta)

include("generaltests.jl")
include("orthotests.jl")
include("laddertests.jl")

@definner TestOrtho
TestOrtho(a,b) = KronDelta(a,b)

default_inner(TestOrtho)

include("generaltests.jl")
include("orthotests.jl")
include("laddertests.jl")
