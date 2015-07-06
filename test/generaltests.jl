@test sum(i -> d" i * | 1,i > ", 1:3) == ket(1) * sum(i -> d" i * | i > ", 1:3)
@test sum(i -> d" i * | i,1 > ", 1:3) == sum(i -> d" i * | i > ", 1:3) * ket(1)

a1,a2,a3,a4 = rand(-10:10,4) + rand(-10:10,4)im
b1,b2,b3,b4 = rand(-10:10,4) + rand(-10:10,4)im
c1,c2 = rand(-10:10,2) + rand(-10:10,2)im

d"""
@test (a1 * | 1 >)' == a1' * < 1 |
@test (a1 * | 1 > + a2 * | :b >)' == a1' * < 1 | + a2' * < :b |
@test a1 * (a2 * | 1 > + a3 * | :b >)' ==  (a1*a2') * < 1 | + (a1*a3') * < :b |

kt = a1*| 1 > + a2*| 'b' > + a3*| :c > + a4*| \"d\" >
br = b1*< 1 | + b2*< 'b' | + b3*< :c | + b4*< \"d\" |
"""

kt = d" a1*| 1 > + a2*| 'b' > + a3*| :c > + a4*| \"d\" > "
br = d" b1*< 1 | + b2*< 'b' | + b3*< :c | + b4*< \"d\" | "

@test kt+kt == 2kt
@test kt-kt == 0kt
@test 3kt - 5kt == -2kt
@test kt-kt+kt == kt
@test (im*kt*br)' == (im'*br'*kt')
@test (c1 * kt') * (c2' * br) == ((c1' * kt) * (c2 * br'))'
@test (kt*kt)*(br*br) == tensor(kt*br,kt*br)

A = sum([d" c1*(i-j)*| i >< j |" for i in 1:3, j in 1:3])
B = sum([d" c2*(i-j)*| i >< j |" for i in 1:3, j in 1:3])

@test A'*B' == (B*A)'
@test A' == map((k,v)->(k', v'), A)
@test A[2,3] == c1 * (2-3)
@test A+A == 2A
@test A-A == 0A
@test A-B+B == A
@test length(filter((o,v) -> klabel(o)==StateLabel(3), A)) == 2
@test convert(OuterSum, A') == maplabels(ctranspose, mapcoeffs(ctranspose, A))
@test xsubspace(tensor(A,A,A),3) == d" A[1,1]^3 * | 1,1,1 >< 1,1,1 | "

C = tensor(A, B')
@test C[(3,1),(1,2)] == (c1*(3-1)) * (c2 * (2-1))'
C[(3,1),(1,2)] = 1
@test C[(3,1),(1,2)] == 1

bell = d" 1/√2 * (| 1,1 > + | 0,0 >) "
belldens = bell * bell'
@test_approx_eq norm(bell) 1
@test ptrace(belldens, 1) == ptrace(belldens, 2)

@test d" act_on(< 1 |, | 'a','b','c' >, 2) == < 1 | 'b' >| 'a','c' > "
@test d" act_on(| 1 >, < 'a','b','c' |, 2) == conj(< 1 | 'b' >)< 'a','c' | "

tranop = d" c1 * | 'a','b','c' >< 'd','e','f' | + c2 * | 'i','j','k' >< 'l','m','n' | "
@test ptranspose(tranop, 2) == d" c1 * | 'a','e','c' >< 'd','b','f' | + c2 * | 'i','m','k' >< 'l','j','n' | "
