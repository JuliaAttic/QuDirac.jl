@assert sum(i -> d" i * | 1,i > ", 1:3) == ket(1) * sum(i -> d" i * | i > ", 1:3)
@assert sum(i -> d" i * | i,1 > ", 1:3) == sum(i -> d" i * | i > ", 1:3) * ket(1)

a1,a2,a3,a4 = rand(Complex{Int64},4)
b1,b2,b3,b4 = conj(rand(Complex{Int64},4))
c1,c2 = rand(Complex{Int64},2)

d"""
@assert (a1 * | 1 >)' == a1' * < 1 | 
@assert (a1 * | 1 > + a2 * | :b >)' == a1' * < 1 | + a2' * < :b | 
@assert a1 * (a2 * | 1 > + a3 * | :b >)' ==  (a1*a2') * < 1 | + (a1*a3') * < :b |

kt = a1*| 1 > + a2*| 'b' > + a3*| :c > + a4*| \"d\" >
br = b1*< 1 | + b2*< 'b' | + b3*< :c | + b4*< \"d\" |
"""
@assert kt+kt == 2kt
@assert kt-kt == 0kt
@assert 3kt - 5kt == -2kt
@assert kt-kt+kt == kt
@assert (im*kt*br)' == (im'*br'*kt')
@assert (c1 * kt') * (c2' * br) == ((c1' * kt) * (c2 * br'))'
@assert (kt*kt)*(br*br) == tensor(kt*br,kt*br)

A = sum([d" c1*(i-j)*| i >< j |" for i in 1:3, j in 1:3])
B = sum([d" c2*(i-j)*| i >< j |" for i in 1:3, j in 1:3])

@assert A'*B' == (B*A)'
@assert A' == map((k,v)->(k', v'), A)
@assert A[2,3] == c1 * (2-3)
@assert A+A == 2A
@assert A-A == 0A
@assert A-B+B == A
@assert length(filter((o,v) -> klabel(o)==StateLabel(3), A)) == 2
@assert convert(OpSum, A') == maplabels(ctranspose, mapcoeffs(ctranspose, A))
@assert xsubspace(tensor(A,A,A),3) == d" A[1,1]^3 * | 1,1,1 >< 1,1,1 | "

C = tensor(A, B')
@assert C[(3,1),(1,2)] == (c1*(3-1)) * (c2 * (2-1))'
C[(3,1),(1,2)] = 1
@assert C[(3,1),(1,2)] == 1

bell = d" 1/√2 * (| 1,1 > + | 0,0 >) "
belldens = bell * bell'
@test_approx_eq norm(bell) 1
@assert ptrace(belldens, 1) == ptrace(belldens, 2)

@assert d" act_on(< 1 |, | 'a','b','c' >, 2) == < 1 | 'b' >| 'a','c' > "
@assert d" act_on(| 1 >, < 'a','b','c' |, 2) == conj(< 1 | 'b' >)< 'a','c' | "

tranop = d" c1 * | 'a','b','c' >< 'd','e','f' | + c2 * | 'i','j','k' >< 'l','m','n' | "
@assert ptranspose(tranop, 2) == d" c1 * | 'a','e','c' >< 'd','b','f' | + c2 * | 'i','m','k' >< 'l','j','n' | "
