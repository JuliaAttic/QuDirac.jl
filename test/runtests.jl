using QuDirac
using Base.Test

filt3s(k,c) = ! (in(3, klabel(k)) || in(3, blabel(k)))

default_inner(UndefInner)

include("generaltests.jl")
include("abstractinnertests.jl")

###########################################################

default_inner(KronDelta)

@def_op " H | n > = 1/√2 * ( | 0 > + (-1)^n * | 1 > ) "
@assert represent(H, map(ket,0:1)) == 1/√2 * [1 1 ; 1 -1]

include("generaltests.jl")
include("orthotests.jl")

@assert raise(d" -im * < 2,1 | ", 2) == d" (sqrt(2) * -im) * < 2,2 |"
@assert lower(d" -im * < 3,5 | ", 2) == d" (sqrt(5) * -im) * < 3,4 |"

@def_op "a | n > = sqrt(n) * | n-1 >"
@def_op "< n | a = sqrt(n + 1) * < n+1 |"
@assert d" act_on(a, | 1,2,3 > + 2| 3,3,1 >, 2) == sqrt(2)| 1,1,3 > + 2*sqrt(3)| 3,2,1 > " 
@assert d" act_on(a', | 1,2,3 > + 2| 3,2,1 >, 2) == sqrt(3)| 1,3,3 > + 2*sqrt(3)| 3,3,1 > "

k3 = sum(i->d" (i + (i * im)) * | i > ", 1:2)^3

@def_op "a_2 | x,y,z > = sqrt(y) * | x, y-1, z >"
@def_op "< x,y,z | a_2 = sqrt(y + 1) * < x, y+1, z |"

@assert a_2 * k3 == act_on(a, k3, 2) == lower(k3, 2)
@assert a_2' * k3 == act_on(a', k3, 2) == raise(k3, 2)
@assert k3' * a_2  == raise(k3', 2)
@assert k3' * a_2' == lower(k3', 2)

###########################################################

@definner TestOrtho
TestOrtho(a,b) = KronDelta(a,b)
default_inner(TestOrtho)

@def_op " Htest | n > = 1/√2 * ( | 0 > + (-1)^n * | 1 > ) "
@assert represent(Htest, map(ket,0:1)) == 1/√2 * [1 1 ; 1 -1]

include("generaltests.jl")
include("orthotests.jl")

@def_op "atest | n > = sqrt(n) * | n-1 >"
@def_op "< n | atest = sqrt(n + 1) * < n+1 |"
@assert d" act_on(atest, | 1,2,3 > + 2| 3,3,1 >, 2) == sqrt(2)| 1,1,3 > + 2*sqrt(3)| 3,2,1 > " 
@assert d" act_on(atest', | 1,2,3 > + 2| 3,2,1 >, 2) == sqrt(3)| 1,3,3 > + 2*sqrt(3)| 3,3,1 > "

k3 = sum(i->d" (i + (i * im)) * | i > ", 1:2)^3

@def_op "a_2test | x,y,z > = sqrt(y) * | x, y-1, z >"
@def_op "< x,y,z | a_2test = sqrt(y + 1) * < x, y+1, z |"

@assert a_2test * k3 == act_on(atest, k3, 2) == lower(k3, 2)
@assert a_2test' * k3 == act_on(atest', k3, 2) == raise(k3, 2)
@assert k3' * a_2test  == raise(k3', 2)
@assert k3' * a_2test' == lower(k3', 2)