testf(b,k) = k[1] - b[1]

@assert inner_eval(testf, b*(k*b)*k) == inner_eval(testf, (b*k)*(b*k))
@assert inner_eval(testf, ptrace(ptrace(bitdens, 2), 1)) == inner_eval(testf, ptrace(ptrace(bitdens, 1), 2))
