k = sum(i->i * ket(i), 1:5)

lower(label) = (sqrt(label[1]), StateLabel(label[1]-1))
a = func_permop(lower, k)

@assert a == sum(i->sqrt(i) * ket(i-1) * bra(i), 1:5)

lowk = map((label,v) -> (StateLabel(label[1]-1), v*sqrt(label[1])), k)

hik = map((label,v) -> (StateLabel(label[1]+1), v*sqrt(label[1]+1)), k)
filternz!(delete!(hik, StateLabel(6)))

@assert lowk == a*k
@assert hik == a'*k
@assert act_on(a, ket(1,2,3) + 2*ket(3,5,1), 2) == sqrt(2)*ket(1,1,3) + 2*sqrt(5)*ket(3,4,1)
@assert act_on(a', ket(1,2,3) + 2*ket(3,3,1), 2) == sqrt(3)*ket(1,3,3) + 2*sqrt(4)*ket(3,4,1)

k3 = k^3

lower2(label) = (sqrt(label[2]), StateLabel(label[1], label[2]-1, label[3]))

a_2 = func_permop(lower2, k3)

@assert a_2 * k3 == act_on(a, k3, 2)