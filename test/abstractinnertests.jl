k = normalize!(sum(i -> rand(Complex{Float64}) * ket(i), 1:5))
b = normalize!(sum(i -> rand(Complex{Float64}) * bra(i), 1:5))

@test_approx_eq inner_eval(KronDelta, b*(k*b)*k) inner_eval(KronDelta, (b*k)*(b*k))

testf(b,k) = k[1] - b[1]
@test_approx_eq inner_eval(testf, b*(k*b)*k) inner_eval(testf, (b*k)*(b*k))

@definner SomeInner
SomeInner(b, k) = b * k
s = (e^(bra(SomeInner,1) * ket(SomeInner,2)) + (bra(SomeInner,3)*ket(SomeInner,4))im)^4
@assert s == (e^2 + 12im)^4