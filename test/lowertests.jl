k = sum(i->i * ket(i), 1:5)

lower(label) = (sqrt(label), label-1)
a = GenericOp(lower, k, 1)

@assert a == sum(i->sqrt(i) * ket(i-1) * bra(i), 1:5)

lowk = map((label,v) -> ({label[1]-1}, v*sqrt(label[1])), k)

hik = map((label,v) -> ({label[1]+1}, v*sqrt(label[1]+1)), k)
filternz!(delete!(hik, {6}))

@assert lowk == a*k
@assert hik == a'*k
@assert act_on(a, ket(1,2,3) + 2*ket(3,5,1), 2) == sqrt(2)*ket(1,1,3) + 2*sqrt(5)*ket(3,4,1)
@assert act_on(a', ket(1,2,3) + 2*ket(3,3,1), 2) == sqrt(3)*ket(1,3,3) + 2*sqrt(4)*ket(3,4,1)

k3 = k^3

a_2 = GenericOp(lower, k3, 2)

@assert a_2 * k3 == act_on(a, k3, 2)