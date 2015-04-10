simpk = (1+3im) * ket(1)
proj = simpk*simpk'
projop = convert(QuDirac.GenericOp, proj)

k = sum([(i+(i*im))*ket(i) for i=0:3])
b = sum([(i+(i*3*im))*ket(i) for i=0:3])'
op = (k*b) + (ket(3)*bra(0))

qubits = normalize(sum(ket, 0:1))^3
bitdens = qubits*qubits'

bell = 1/sqrt(2) * (ket(1,1) + ket(0,0))
bell_unbal = normalize!(ket(0,1) + 2.0*ket(1,0))
belldens = bell * bell'