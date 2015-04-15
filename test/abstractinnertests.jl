testf(b,k) = k[1] - b[1]

@assert inner_eval(testf, b*(k*b)*k) == inner_eval(testf, (b*k)*(b*k))
@assert inner_eval(testf, ptrace(ptrace(bitdens, 2), 1)) == inner_eval(testf, ptrace(ptrace(bitdens, 1), 2))

immutable SomeInner <: AbstractInner end
QuDirac.inner_rule(::SomeInner, b, k) = b[1]*k[1]
s = (e^(bra(SomeInner(),1) * ket(SomeInner(),2)) + (bra(SomeInner(),3)*ket(SomeInner(),4))im)^4
@assert s == (e^2 + 12im)^4