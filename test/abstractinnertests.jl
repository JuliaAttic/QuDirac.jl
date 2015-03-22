f(i) = ktlabel(i)[1] - brlabel(i)[1]

@assert queval(f, b*(k*b)*k) == queval(f, (b*k)*(b*k))
@assert queval(f, ptrace(ptrace(bitdens, 2), 1)) == queval(f, ptrace(ptrace(bitdens, 1), 2))
