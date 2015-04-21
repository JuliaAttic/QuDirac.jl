simpk = d" (1+3im)| 1 > "
proj = simpk*simpk'
projop = convert(QuDirac.OpSum, proj)

k = sum([d" (i+(i*im)) * | i >" for i=0:3])
b = sum([d" (i+(i*3*im)) * | i >" for i=0:3])'
op = (k*b) + d"| 3 >< 0 |"

qubits = normalize(sum(ket, 0:1))^3
bitdens = qubits*qubits'

bell = d" 1/√2 * (| 1,1 > + | 0,0 >) "
bell_unbal = d" normalize!( | 0,1 > + 2.0| 1,0 > ) "
belldens = bell * bell'