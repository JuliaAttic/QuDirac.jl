@repr_op "a | n > = sqrt(n)| n-1 >" 1:3
@assert d" act_on(a, | 1,2,3 > + 2| 3,3,1 >, 2) == sqrt(2)| 1,1,3 > + 2*sqrt(3)| 3,2,1 > "
@assert d" act_on(a', | 1,2,3 > + 2| 3,2,1 >, 2) == sqrt(3)| 1,3,3 > + 2*sqrt(3)| 3,3,1 > "

k3 = sum(i->d" i * | i > ", 1:3)^3

@repr_op "a_2 | x,y,z > = sqrt(y)| x, y-1, z >" [(i,j,k) for i=0:3, j=0:3, k=0:3]
@def_op "lowf_2 | x,y,z > = sqrt(y)| x, y-1, z >"

@assert a_2 * k3 == act_on(a, k3, 2)
@assert a_2 * k3 == lowf_2 * k3
@assert lowf_2 * k3 == lower(k3, 2)


