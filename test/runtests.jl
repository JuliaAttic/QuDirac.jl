using QuDirac
using Base.Test

filt3s(k,c) = ! (in(3, klabel(k)) || in(3, blabel(k)))

default_inner(UndefInner)

include("generaltests.jl")
include("abstractinnertests.jl")

###########################################################

default_inner(KronDelta)

@defop " H | n > = 1/√2 * ( | 0 > + (-1)^n * | 1 > ) "
@test represent(H, map(ket,0:1)) == 1/√2 * [1 1 ; 1 -1]

include("generaltests.jl")
include("orthotests.jl")

@test raise(d" -im * < 2,1 | ", 2) == d" (sqrt(2) * -im) * < 2,2 |"
@test lower(d" -im * < 3,5 | ", 2) == d" (sqrt(5) * -im) * < 3,4 |"

@defop "a | n > = sqrt(n) * | n-1 >"
@defop "< n | a = sqrt(n + 1) * < n+1 |"
@test d" act_on(a, | 1,2,3 > + 2| 3,3,1 >, 2) == sqrt(2)| 1,1,3 > + 2*sqrt(3)| 3,2,1 > "
@test d" act_on(a', | 1,2,3 > + 2| 3,2,1 >, 2) == sqrt(3)| 1,3,3 > + 2*sqrt(3)| 3,3,1 > "

k3 = sum(i->d" (i + (i * im)) * | i > ", 1:2)^3

@defop "a_2 | x,y,z > = sqrt(y) * | x, y-1, z >"
@defop "< x,y,z | a_2 = sqrt(y + 1) * < x, y+1, z |"

@test a_2 * k3 == act_on(a, k3, 2) == lower(k3, 2)
@test a_2' * k3 == act_on(a', k3, 2) == raise(k3, 2)
@test k3' * a_2  == raise(k3', 2)
@test k3' * a_2' == lower(k3', 2)

###########################################################

@definner TestOrtho
TestOrtho(a,b) = KronDelta(a,b)
default_inner(TestOrtho)

@defop " Htest | n > = 1/√2 * ( | 0 > + (-1)^n * | 1 > ) "
@test represent(Htest, map(ket,0:1)) == 1/√2 * [1 1 ; 1 -1]

include("generaltests.jl")
include("orthotests.jl")

@defop "atest | n > = sqrt(n) * | n-1 >"
@defop "< n | atest = sqrt(n + 1) * < n+1 |"
@test d" act_on(atest, | 1,2,3 > + 2| 3,3,1 >, 2) == sqrt(2)| 1,1,3 > + 2*sqrt(3)| 3,2,1 > "
@test d" act_on(atest', | 1,2,3 > + 2| 3,2,1 >, 2) == sqrt(3)| 1,3,3 > + 2*sqrt(3)| 3,3,1 > "

k3 = sum(i->d" (i + (i * im)) * | i > ", 1:2)^3

@defop "a_2test | x,y,z > = sqrt(y) * | x, y-1, z >"
@defop "< x,y,z | a_2test = sqrt(y + 1) * < x, y+1, z |"

@test a_2test * k3 == act_on(atest, k3, 2) == lower(k3, 2)
@test a_2test' * k3 == act_on(atest', k3, 2) == raise(k3, 2)
@test k3' * a_2test  == raise(k3', 2)
@test k3' * a_2test' == lower(k3', 2)