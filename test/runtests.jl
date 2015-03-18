using QuDirac
using Base.Test

k = Ket([[i]=>i+(i*im) for i=0:3])
b = Ket([[i]=>i+(i*3*im) for i=0:3])'
op = k*b + Ket(3)*Bra(0)

@assert b*k == 56 - 28im
@assert b'*b'*b' == (b^3)'
@assert (op'*op') == (op*op)'
@assert maplabels(reverse, mapcoeffs(ctranspose, op)) == op'
@assert op[2,1] == op[Vector[[2],[1]]] == 8-4im
@assert k*k*b*b == tensor(k*b,k*b)
@assert b*op*k == 2352 - 3136im
@assert op+op == 2 * op
@assert op-op == 0 * op
@assert op-op+op == op
@assert filternz!(xsubspace(tensor(op,op,op),3)) == (16 - 88im) * Ket(1,1,1) * Bra(1,1,1) 

qubits = normalize!(sum(Ket, 0:1))^3
bitdens = qubits*qubits'
@test_approx_eq qubits'*qubits 1
@test_approx_eq norm(qubits) 1
@test_approx_eq norm(ptrace(bitdens, 2)) 1
@test_approx_eq trace(ptrace(bitdens, 2)) 1
@assert ptrace(ptrace(bitdens, 2), 1) == ptrace(ptrace(bitdens, 1), 2)

bell = 1/sqrt(2) * (Ket(1,1) + Ket(0,0))
belldens = bell * bell'

@assert ptrace(belldens, 1) == ptrace(belldens, 2)
@test_approx_eq trace(ptrace(belldens, 1)^2) .5
