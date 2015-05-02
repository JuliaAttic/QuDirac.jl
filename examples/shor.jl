using QuDirac

# Shor's algorithm for N = 15

# n = 4 qubits
# y = 13 such that gcd(y, N) = 1

regs = normalize(sum(ket, 0:15)) * d" | 0 > "