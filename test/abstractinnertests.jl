f(b,k) = k[1] - b[1]

@assert inner_eval(f, b*(k*b)*k) == inner_eval(f, (b*k)*(b*k))
@assert inner_eval(f, ptrace(ptrace(bitdens, 2), 1)) == inner_eval(f, ptrace(ptrace(bitdens, 1), 2))
