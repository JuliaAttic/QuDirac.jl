@assert proj[1,1] == 10
@assert (simpk'*proj)[1] == 10 - 30im
@assert (proj*simpk)[1] == 10 + 30im
@assert simpk'*proj*simpk == (simpk'*simpk)*(simpk'*simpk) == 100

@assert b*(k*b)*k == (b*k)*(b*k)
@assert b*k == 56 - 28im
@assert b*op*k == 2352 - 3136im

@test qubits'*qubits ≈ 1
@test norm(ptrace(bitdens, 2)) ≈ 1
@test trace(ptrace(bitdens, 2)) ≈ 1
@test trace(ptrace(belldens, 1)^2) ≈ .5
@test purity(bell) ≈ 1

@rep_op " H | n > = 1/√2 * ( | 0 > + (-1)^n * | 1 > )" 0:1
@assert matrep(H, 0:1) == 1/√2 * [1 1 ; 1 -1]
m = 0.5 * [1 1;1 -1]
@test matrep(tensor(H, H), 0:1, 0:1) ≈ [m m ; m -m]

@assert ptrace(ptrace(bitdens, 2), 1) == ptrace(ptrace(bitdens, 1), 2)

@assert act_on(bra(0), bell_unbal, 2)[1] == bell_unbal[1,0]
@assert act_on(bra(0), bell_unbal, 1)[1] == bell_unbal[0,1]

@assert anticommute(op,op) == 2 * op^2
@assert commute(op,op) == op^2 - op^2

b = d" 1/√2 * (< 0 | + < 1 |) "
itest = act_on(b, bell_unbal, 2)

@assert itest[0] == b[1]*bell_unbal[0,1]
@assert itest[1] == b[0]*bell_unbal[1,0]