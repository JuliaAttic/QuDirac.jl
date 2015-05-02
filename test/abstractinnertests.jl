testf(b,k) = k[1] - b[1]

@assert inner_eval(testf, b*(k*b)*k) == inner_eval(testf, (b*k)*(b*k))

@def_inner SomeInner Any
SomeInner(b, k) = b * k
s = (e^(bra(SomeInner,1) * ket(SomeInner,2)) + (bra(SomeInner,3)*ket(SomeInner,4))im)^4
@assert s == (e^2 + 12im)^4