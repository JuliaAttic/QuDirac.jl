simpk = d" (1+3im)| 1 > "
proj = simpk*simpk'
projop = convert(OpSum, proj)

k = sum([d" (i+(i*im)) * | i >" for i=0:3])
b = sum([d" (i+(i*3*im)) * | i >" for i=0:3])'
op = (k*b) + d"| 3 >< 0 |"

qubits = normalize(sum(ket, 0:1))^3
bitdens = qubits*qubits'

bell = d" 1/√2 * (| 1,1 > + | 0,0 >) "
bell_unbal = d" normalize!( | 0,1 > + 2.0| 1,0 > ) "
belldens = bell * bell'

######################

@assert proj[1,1] == 10
@assert (simpk'*proj)[1] == 10 - 30im
@assert (proj*simpk)[1] == 10 + 30im
@assert simpk'*proj*simpk == (simpk'*simpk)*(simpk'*simpk) == 100

@assert b*(k*b)*k == (b*k)*(b*k)
@assert b*k == 56 - 28im
@assert b*op*k == 2352 - 3136im

@test_approx_eq qubits'*qubits 1
@test_approx_eq norm(ptrace(bitdens, 2)) 1
@test_approx_eq trace(ptrace(bitdens, 2)) 1
@test_approx_eq trace(ptrace(belldens, 1)^2) .5
@test_approx_eq purity(bell) 1

@rep_op " H | n > = 1/√2 * ( | 0 > + (-1)^n * | 1 > )" 0:1
@assert represent(H, 0:1) == 1/√2 * [1 1 ; 1 -1]
m = 0.5 * [1 1;1 -1]
@test_approx_eq convert(Matrix{Float64}, represent(tensor(H, H), 0:1, 0:1)) [m m ; m -m]

@assert ptrace(ptrace(bitdens, 2), 1) == ptrace(ptrace(bitdens, 1), 2)

@assert act_on(bra(0), bell_unbal, 2)[1] == bell_unbal[1,0]
@assert act_on(bra(0), bell_unbal, 1)[1] == bell_unbal[0,1]

@assert anticommute(op,op) == 2 * op^2
@assert commute(op,op) == op^2 - op^2

b = d" 1/√2 * (< 0 | + < 1 |) "
itest = act_on(b, bell_unbal, 2)

@assert itest[0] == b[1]*bell_unbal[0,1]
@assert itest[1] == b[0]*bell_unbal[1,0]