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

@assert ptrace(ptrace(bitdens, 2), 1) == ptrace(ptrace(bitdens, 1), 2)

@assert innerat(bra(0), bell_unbal, 2)[1] == bell_unbal[1,0]
@assert innerat(bra(0), bell_unbal, 1)[1] == bell_unbal[0,1]

b = 1/√2 * (bra(0) + bra(1))
itest = innerat(b, bell_unbal, 2)

@assert itest[0] == b[1]*bell_unbal[0,1]
@assert itest[1] == b[0]*bell_unbal[1,0]