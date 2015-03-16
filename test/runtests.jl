using QuDirac
using Base.Test

k = Ket([(i,)=>i+(i*im) for i=0:3])
b = Ket([(i,)=>i+(i*3*im) for i=0:3])'
op = k*b
@assert b*k == 56 - 28im
@assert b'*b'*b' == (b^3)'
@assert (op'*op') == (op*op)'
@assert (op'*op) == (op*op')
@assert mapcoeffs(ctranspose, op') == op
@assert op[2,1] == op[(2,),(1,)] == 8-4im
@assert k*op*b == tensor(op,op)
@assert b*op*k == 2352 - 3136im
@assert op+op == 2 * op
@assert op-op == 0 * op
@assert op-op+op == op
@assert filternz!(xsubspace(tensor(op,op,op),3)) == (16 - 88im) * Ket(1,1,1) * Bra(1,1,1) 

qubits = normalize!(sum(Ket, 0:1))^3
@test_approx_eq qubits'*qubits 1
@test_approx_eq norm(qubits) 1
@test_approx_eq norm(ptrace(qubits*qubits', 2)) 1
@test_approx_eq trace(ptrace(qubits*qubits', 2)) 1

bell1 = 1/sqrt(2) * (Ket(1,1) + Ket(0,0))
bell2 = 1/sqrt(2) * (Ket(1,0) + Ket(0,1))
dens1 = bell1 * bell1'
dens2 = bell2 * bell2'

@test_approx_eq trace(ptrace(dens1, 1)^2) .5
@test_approx_eq trace(ptrace(dens1, 2)^2) .5
@test_approx_eq trace(ptrace(dens2, 1)^2) .5
@test_approx_eq trace(ptrace(dens2, 2)^2) .5