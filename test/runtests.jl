using QuDirac
using Base.Test

k = (1+3im) * ket(1)
proj = k*k'
op = convert(QuDirac.DiracOp, proj)

@assert proj[1,1] == 10
@assert (k'*proj)[1] == 10 - 30im
@assert (proj*k)[1] == 10 + 30im
@assert k'*proj*k == (k'*k)*(k'*k) == 100
@assert (3-im) * k' ==  (3-im) * ((1-3im) * bra(1))
@assert (3-im) * k' == ((3+im)*(1+3im)*ket(1))'
@assert (-im * op')[1,1] == 0-10im

k = sum([(i+(i*im))*ket(i) for i=0:3])
b = sum([(i+(i*3*im))*ket(i) for i=0:3])'
op = (k*b) + (ket(3)*bra(0))

@assert b*k == 56 - 28im
@assert b'*b'*b' == (b^3)'
@assert (op'*op') == (op*op)'
@assert maplabels(reverse, mapcoeffs(ctranspose, op)) == op'
@assert op[2,1] == 8-4im
@assert (k*k)*(b*b) == tensor(k*b,k*b)
@assert b*(k*b)*k == (b*k)*(b*k)
@assert b*op*k == 2352 - 3136im
@assert op+op == 2 * op
@assert op-op == 0 * op
@assert op-op+op == op
@assert filternz!(xsubspace(tensor(op,op,op),3)) == (16 - 88im) * ket(1,1,1) * bra(1,1,1) 

qubits = normalize!(sum(ket, 0:1))^3
bitdens = qubits*qubits'
@test_approx_eq qubits'*qubits 1
@test_approx_eq norm(qubits) 1
@test_approx_eq norm(ptrace(bitdens, 2)) 1
@test_approx_eq trace(ptrace(bitdens, 2)) 1
@assert ptrace(ptrace(bitdens, 2), 1) == ptrace(ptrace(bitdens, 1), 2)

bell = 1/sqrt(2) * (ket(1,1) + ket(0,0))
belldens = bell * bell'

@assert ptrace(belldens, 1) == ptrace(belldens, 2)
@test_approx_eq trace(ptrace(belldens, 1)^2) .5
@test_approx_eq purity(bell) 1
