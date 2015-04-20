k = sum(i->i * ket(i), 1:3)

@repr_op "a | n > = sqrt(n) * | n-1 >" 1:3

@assert a == sum(i->sqrt(i) * ket(i-1) * bra(i), 1:3)

lowk = map((label,v) -> (StateLabel(label[1]-1), v*sqrt(label[1])), k)

hik = map((label,v) -> (StateLabel(label[1]+1), v*sqrt(label[1]+1)), k)
filternz!(delete!(hik, StateLabel(4)))

@assert lowk == a*k
@assert hik == a'*k
@assert act_on(a, ket(1,2,3) + 2*ket(3,3,1), 2) == sqrt(2)*ket(1,1,3) + 2*sqrt(3)*ket(3,2,1)
@assert act_on(a', ket(1,2,3) + 2*ket(3,2,1), 2) == sqrt(3)*ket(1,3,3) + 2*sqrt(3)*ket(3,3,1)

k3 = k^3

@repr_op "a_2 | x,y,z > = sqrt(y) * | x, y-1, z >" [(i,j,k) for i=0:3, j=0:3, k=0:3]
@def_op "a_2f | x,y,z > = sqrt(y) * | x, y-1, z >"

@assert a_2 * k3 == act_on(a, k3, 2)
@assert a_2 * k3 == a_2f * k3