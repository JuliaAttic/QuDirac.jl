@rep_op "a | n > = sqrt(n)| n-1 >" 0:4
@assert d" act_on(a, | 1,2,3 > + 2| 3,3,1 >, 2) == sqrt(2)| 1,1,3 > + 2*sqrt(3)| 3,2,1 > "
@assert d" act_on(a', | 1,2,3 > + 2| 3,2,1 >, 2) == sqrt(3)| 1,3,3 > + 2*sqrt(3)| 3,3,1 > "

k3 = sum(i->d" (i + (i * im)) * | i > ", 1:2)^3

@rep_op "a_2 | x,y,z > = sqrt(y) * | x, y-1, z >" 0:3 0:3 0:3
@rep_op "< x,y,z | a_2_c = sqrt(y+1) * < x, y+1, z |" 0:2 0:2 0:2

@def_op "lowf_2 | x,y,z > = sqrt(y) * | x, y-1, z >"
@def_op "< x,y,z | lowf_2 = sqrt(y + 1) * < x, y+1, z |"

filt3s(k,c) = ! (in(3, klabel(k)) || in(3, blabel(k)))

@assert filter(filt3s, a_2) == filter(filt3s, a_2_c)
@assert a_2 * k3 == act_on(a, k3, 2) == lowf_2 * k3  == lower(k3, 2)
@assert a_2' * k3 == act_on(a', k3, 2) == lowf_2' * k3 == raise(k3, 2)
@assert k3' * a_2  ==  k3' * lowf_2  == raise(k3', 2)
@assert k3' * a_2'  ==  k3' * lowf_2'  == lower(k3', 2)