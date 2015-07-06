simpk = d" (1+3im)| 1 > "
proj = simpk*simpk'
projop = convert(OuterSum, proj)

k = sum([d" (i+(i*im)) * | i >" for i=0:3])
b = sum([d" (i+(i*3*im)) * | i >" for i=0:3])'
op = (k*b) + d"| 3 >< 0 |"

qubits = normalize(sum(ket, 0:1))^3
bitdens = qubits*qubits'

bell = d" 1/√2 * (| 1,1 > + | 0,0 >) "
bell_unbal = d" normalize!( | 0,1 > + 2.0| 1,0 > ) "
belldens = bell * bell'

######################

@test proj[1,1] == 10
@test (simpk'*proj)[1] == 10 - 30im
@test (proj*simpk)[1] == 10 + 30im
@test simpk'*proj*simpk == (simpk'*simpk)*(simpk'*simpk) == 100

@test b*(k*b)*k == (b*k)*(b*k)
@test b*k == 56 - 28im
@test b*op*k == 2352 - 3136im

@test_approx_eq qubits'*qubits 1
@test_approx_eq norm(ptrace(bitdens, 2)) 1
@test_approx_eq trace(ptrace(bitdens, 2)) 1
@test_approx_eq trace(ptrace(belldens, 1)^2) .5

@test ptrace(ptrace(bitdens, 2), 1) == ptrace(ptrace(bitdens, 1), 2)

@test act_on(bra(0), bell_unbal, 2)[1] == bell_unbal[1,0]
@test act_on(bra(0), bell_unbal, 1)[1] == bell_unbal[0,1]

@test anticommute(op,op) == 2 * op^2
@test commute(op,op) == op^2 - op^2

b = d" 1/√2 * (< 0 | + < 1 |) "
itest = act_on(b, bell_unbal, 2)

@test itest[0] == b[1]*bell_unbal[0,1]
@test itest[1] == b[0]*bell_unbal[1,0]

@test raise(d" -im * < 2,1 | ", 2) == d" (sqrt(2) * -im) * < 2,2 |"
@test lower(d" -im * < 3,5 | ", 2) == d" (sqrt(5) * -im) * < 3,4 |"
