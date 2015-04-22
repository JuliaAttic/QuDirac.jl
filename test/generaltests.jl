@assert (3-im) * simpk' ==  d"(3-im) * (1-3im)< 1 |"
@assert (3-im) * simpk' == (d"(3+im)*(1+3im)| 1 >")'
@assert (-im * projop')[1,1] == 0-10im

@assert length(filter((o,v) -> klabel(o)==StateLabel(3), op')) == 3
@assert b'*b'*b' == (b^3)'
@assert (op'*op') == (op*op)'
@assert op' == maplabels(ctranspose, mapcoeffs(ctranspose, op)) 
@assert op' == map((k,v)->(k', v'), op)
@assert op[2,1] == 8-4im
@assert (k*k)*(b*b) == tensor(k*b,k*b)
@assert op+op == 2 * op
@assert op-op == 0 * op
@assert op-op+op == op
@assert filternz!(xsubspace(tensor(op,op,op),3)) == d"(16 - 88im)| 1,1,1 >< 1,1,1 |" 

op_copy = copy(op)
op_copy[3,0] = 32+im
@assert op_copy[3,0] == 32+im
@assert tensor(op_copy, op_copy')[(3,1),(1,2)] == 120

@test_approx_eq norm(qubits) 1
@assert ptrace(belldens, 1) == ptrace(belldens, 2)

@assert d" act_on(< 1 |, | 'a','b','c' >, 2) == < 1 | 'b' >| 'a','c' > "
@assert d" act_on(| 1 >, < 'a','b','c' |, 2) == conj(< 1 | 'b' >)< 'a','c' | "

tranop = d" | 'a','b','c' >< 'd','e','f' | + 2 * | 'i','j','k' >< 'l','m','n' | "
@assert ptranspose(tranop, 2) == d" | 'a','e','c' >< 'd','b','f' | + 2 * | 'i','m','k' >< 'l','j','n' | "