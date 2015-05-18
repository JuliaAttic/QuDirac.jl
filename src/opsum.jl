abstract AbsOpSum{P,N,K,B,T} <: DiracOp{P,N}
Base.eltype(op::AbsOpSum) = eltype(typeof(op))
ketlabeltype{P,N,K,B}(op::AbsOpSum{P,N,K,B}) = StateLabel{N,K}
bralabeltype{P,N,K,B}(op::AbsOpSum{P,N,K,B}) = StateLabel{N,B}

#########
# OpSum #
#########
type OpSum{P,N,K,B,T} <: AbsOpSum{P,N,K,B,T}
    data::SumDict{OpLabel{N,K,B},T}
end

OpSum{P,N,K,B,T}(::Type{P}, data::SumDict{OpLabel{N,K,B},T}) = OpSum{P,N,K,B,T}(data)

function OpSum{P,N}(kt::Ket{P,N}, br::Bra{P,N})
    T = promote_type(eltype(kt), eltype(br))
    K = eltype(labeltype(kt))
    B = eltype(labeltype(br))
    result = @compat sizehint!(SumDict{OpLabel{N,K,B},T}(), length(kt) * length(br))
    return OpSum(P, cons_outer!(result, kt, br))
end

function cons_outer!(result::SumDict, kt::KetSum, br::BraSum)
    for (k,kc) in data(kt)
        for (b,bc) in data(br)
            add_to_dict!(result, OpLabel(k, b), kc * bc')
        end
    end
    return result
end

function cons_outer!(result::SumDict, kt::SingleKet, br::BraSum)
    k = label(kt)
    kc = coeff(kt)
    for (b,bc) in data(br)
        add_to_dict!(result, OpLabel(k, b), kc * bc')
    end
    return result
end

function cons_outer!(result::SumDict, kt::KetSum, br::SingleBra)
    b = label(br)
    bc = coeff(br)
    for (k,kc) in data(kt)
        add_to_dict!(result, OpLabel(k, b), kc * bc)
    end
    return result
end

function cons_outer!(result::SumDict, kt::SingleKet, br::SingleBra)
    add_to_dict!(result, OpLabel(label(kt), label(br)), coeff(kt) * coeff(br))
    return result
end


data(op::OpSum) = op.data

Base.convert{P,N,K,B,T}(::Type{OpSum{P,N,K,B,T}}, op::OpSum{P}) = OpSum(P, convert(SumDict{OpLabel{N,K,B},T}, data(op)))

function Base.promote_rule{P,N,K1,B1,T1,K2,B2,T2}(::Type{OpSum{P,N,K1,B1,T1}}, ::Type{OpSum{P,N,K2,B2,T2}})
    return OpSum{P,N,promote_type(K1,K2),promote_type(B1,B2),promote_type(T1,T2)}
end

Base.eltype{P,N,K,B,T}(::Type{OpSum{P,N,K,B,T}}) = T
Base.hash{P}(op::OpSum{P}) = hash(data(op), hash(P))
Base.hash(op::OpSum, h::Uint64) = hash(hash(op), h)
Base.copy{P}(op::OpSum{P}) = OpSum(P, copy(data(op)))
Base.(:(==)){P,N}(a::OpSum{P,N}, b::OpSum{P,N}) = data(a) == data(b)
Base.length(op::OpSum) = length(data(op))

Base.getindex(op::OpSum, x::OpLabel) = op.data[x]
Base.setindex!(op::OpSum, x, y::OpLabel) = setindex!(data(op), x, y)
Base.haskey(op::OpSum, x::OpLabel) = haskey(data(op), x)
Base.get(op::OpSum, x::OpLabel, default=predict_zero(eltype(op))) = get(data(op), x, default)

Base.start(op::OpSum) = start(data(op))
Base.next(op::OpSum, i) = next(data(op), i)
Base.done(op::OpSum, i) = done(data(op), i)

Base.collect(op::OpSum) = collect(data(op))

#############
# DualOpSum #
#############
type DualOpSum{P,N,K,B,T} <: AbsOpSum{P,N,K,B,T}
    op::OpSum{P,N,B,K,T}
end

Base.convert{P,N,K,B,T}(::Type{OpSum{P,N,B,K,T}}, opc::DualOpSum{P,N,K,B,T}) = eager_ctranspose(opc.op)
function Base.convert{P,N,K,B,T}(::Type{DualOpSum{P,N,K,B,T}}, opc::DualOpSum)
    return DualOpSum(convert(OpSum{P,N,B,K,T}, opc.op))
end

Base.promote_rule{P,N,K,B,T}(::Type{OpSum{P,N,B,K,T}}, ::Type{DualOpSum{P,N,K,B,T}}) = OpSum{P,N,B,K,T}
function Base.promote_rule{P,N,K1,B1,T1,K2,B2,T2}(::Type{DualOpSum{P,N,K1,B1,T1}}, ::Type{DualOpSum{P,N,K2,B2,T2}})
    return DualOpSum{P,N,promote_type(K1,K2),promote_type(B1,B2),promote_type(T1,T2)}
end

data(opc::DualOpSum) = data(opc.op)

Base.eltype{P,N,K,B,T}(::Type{DualOpSum{P,N,K,B,T}}) = T
Base.hash(opc::DualOpSum) = hash(eager_ctranspose(opc.op))
Base.hash(opc::DualOpSum, h::Uint64) = hash(hash(opc), h)
Base.copy(opc::DualOpSum) = DualOpSum(copy(opc.op))

Base.(:(==)){P,N}(a::DualOpSum{P,N}, b::DualOpSum{P,N}) = a.op == b.op

Base.length(opc::DualOpSum) = length(opc.op)
Base.getindex(opc::DualOpSum, x::OpLabel) = opc.op[x']'
Base.setindex!(opc::DualOpSum, x, y::OpLabel) = setindex!(opc.op, x', y')
Base.haskey(opc::DualOpSum, x::OpLabel) = haskey(opc.op, x')
Base.get(opc::DualOpSum, x::OpLabel, default=predict_zero(eltype(opc))) = get(opc.op, x', default)

Base.collect{P,N,K,B,T}(opc::DualOpSum{P,N,K,B,T}) = collect_pairs!(Array(@compat(Tuple{OpLabel{N,K,B}, T}), length(opc)), opc)

ctpair(k,v) = (k', v')

function collect_pairs!(result, opc::DualOpSum)
    i = 1
    for (k,v) in data(opc)
        result[i] = ctpair(k,v)
        i += 1
    end
    return result
end

Base.start(opc::DualOpSum) = start(opc.op)

function Base.next(opc::DualOpSum, i)
    (k,v), n = next(dict(opc), i)
    return ((k',v'), n)
end

Base.done(opc::DualOpSum, i) = done(opc.op, i)

#############
# AbsOpSum #
#############
Base.getindex(op::AbsOpSum, k::StateLabel, b::StateLabel) = op[OpLabel(k,b)]
Base.getindex(op::AbsOpSum, k, b) = op[StateLabel(k),StateLabel(b)]

Base.setindex!(op::AbsOpSum, x, k::StateLabel, b::StateLabel) = setindex!(op, x, OpLabel(k,b))
Base.setindex!(op::AbsOpSum, x, k, b) = setindex!(op, x, StateLabel(k), StateLabel(b))

Base.haskey(op::AbsOpSum, k::StateLabel, b::StateLabel) = haskey(op, OpLabel(k,b))
Base.haskey(op::AbsOpSum, k, b) = haskey(op, StateLabel(k), StateLabel(b))

Base.get(op::AbsOpSum, k::StateLabel, b::StateLabel, default=predict_zero(eltype(op))) = get(op, OpLabel(k,b), default)
Base.get(op::AbsOpSum, k, b, default=predict_zero(eltype(op))) = get(op, StateLabel(k), StateLabel(b), default)

#############
# ctranpose #
#############
eager_ctranspose{P}(op::OpSum{P}) = OpSum(P, SumDict(collect(op')))

Base.ctranspose(op::OpSum) = DualOpSum(op)
Base.ctranspose(opc::DualOpSum) = opc.op

#########
# inner #
#########
function inner{P,N}(br::Bra{P,N}, op::OpSum{P,N})
    result = SumDict{bralabeltype(op), inner_rettype(br, op)}()
    return Bra(P, inner_load!(result, br, op))'
end

function inner_load!{P}(result::SumDict, br::SingleBra{P}, op::OpSum{P})
    c = coeff(br)
    b = label(br)
    for (o,v) in data(op)
        add_to_dict!(result, blabel(o), ctranspose(v*c*inner_rule(P, b, klabel(o))))
    end
    return result
end

function inner_load!{P}(result::SumDict, br::BraSum{P}, op::OpSum{P})
    for (o,v) in data(op)
        k = klabel(o)
        coeff = predict_zero(eltype(result))
        for (b,c) in data(br)
            coeff += c' * v * inner(P, b, k) 
        end
        return 
        add_to_dict!(result, blabel(o), ctranspose(coeff))
    end
    return result
end

function inner{P,N,A,B}(op::OpSum{P,N,A}, kt::Ket{P,N,B})
    result = SumDict{ketlabeltype(op), inner_rettype(op, kt)}()
    return Ket(P, inner_load!(result, op, kt))
end

function inner_load!{P}(result::SumDict, op::OpSum{P}, kt::SingleKet{P})
    c = coeff(kt)
    k = label(kt)
    for (o,v) in data(op)
        add_to_dict!(result, klabel(o), v*c*inner_rule(P, blabel(o), k))
    end
    return result
end

function inner_load!{P}(result::SumDict, op::OpSum{P}, kt::KetSum{P})
    for (o,v) in data(op)
        b = blabel(o)
        coeff = predict_zero(eltype(result))
        for (k,c) in data(kt)
            coeff += c * v * inner(P, b, k) 
        end
        return 
        add_to_dict!(result, klabel(o), coeff)
    end
    return result
end

function inner{P,N,K1,B1,K2,B2}(a::OpSum{P,N,K1,B1}, b::OpSum{P,N,K2,B2})
    result = SumDict{OpLabel{N,K1,B2}, inner_rettype(a, b)}()
    return OpSum(P, inner_load!(result, a, b))
end

function inner_load!{P}(result::SumDict, a::OpSum{P}, b::OpSum{P})
    for (o1,v) in data(a), (o2,c) in data(b)
        add_to_dict!(result, 
                     OpLabel(klabel(o1), blabel(o2)),
                     v*c*inner(P, blabel(o1), klabel(o2)))
    end
    return result
end

function inner{P,N,K1,B1,K2,B2}(op::OpSum{P,N,K1,B1}, opc::DualOpSum{P,N,K2,B2})
    result = SumDict{OpLabel{N,K1,B2}, inner_rettype(op, opc)}()
    return OpSum(P, inner_load!(result, op, opc))
end

function inner_load!{P}(result::SumDict, op::OpSum{P}, opc::DualOpSum{P})
    for (o,v) in data(op), (oc,c) in data(opc)
        add_to_dict!(result, 
                     OpLabel(klabel(o), klabel(oc)),
                     v*c'*inner(P, blabel(o), blabel(oc)))
    end
    return result
end

function inner{P,N,K1,B1,K2,B2}(opc::DualOpSum{P,N,K1,B1}, op::OpSum{P,N,K2,B2})
    result = SumDict{OpLabel{N,K1,B2}, inner_rettype(opc, op)}()
    return OpSum(P, inner_load!(result, opc, op))
end

function inner_load!{P}(result::SumDict, opc::DualOpSum{P}, op::OpSum{P})
    for (oc,v) in data(opc), (o,c) in data(op)
        add_to_dict!(result,
                     OpLabel(blabel(oc), blabel(o)),
                     v'*c*inner(P, klabel(oc), klabel(o)))
    end
    return result
end

inner(br::Bra, opc::DualOpSum) = inner(opc', br')'
inner(opc::DualOpSum, kt::Ket) = inner(kt', opc')'
inner(a::DualOpSum, b::DualOpSum) = inner(b', a')'

Base.(:*)(br::Bra, op::DiracOp) = inner(br,op)
Base.(:*)(op::DiracOp, kt::Ket) = inner(op,kt)
Base.(:*)(a::DiracOp, b::DiracOp) = inner(a,b)

inner_eval(f, op::DiracOp) = mapcoeffs(x->inner_eval(f,x),op)

##########
# act_on #
##########
act_on(op::AbsOpSum, br::Bra, i) = act_on(op', br', i)'

act_on{P}(op::AbsOpSum{P,1}, kt::Ket{P,1}, i) = i==1 ? inner(op, kt) : throw(BoundsError())

function act_result{P,N}(op::AbsOpSum{P,1}, kt::Ket{P,N})
    return SumDict{promote_type(ketlabeltype(op), labeltype(kt)), inner_rettype(op, kt)}()
end

act_on{P,N}(op::AbsOpSum{P,1}, kt::Ket{P,N}, i) = Ket(P, act_on_dict!(act_result(op,kt), op, kt, i))

function act_on_dict!{P}(result, op::OpSum{P}, kt::SingleKet{P}, i)
    k = label(kt)
    v = coeff(kt)
    k_i = k[i]
    for (o,c) in data(op)
        add_to_dict!(result, 
                     setindex(k, klabel(o)[1], i),
                     c*v*P(blabel(o)[1], k_i))
    end
    return result
end

function act_on_dict!{P}(result, op::OpSum{P}, kt::KetSum{P}, i)
    for (o,c) in data(op), (k,v) in data(kt)
        add_to_dict!(result, 
                     setindex(k, klabel(o)[1], i),
                     c*v*P(blabel(o)[1], k[i]))
    end
    return result
end

function act_on_dict!{P}(result, opc::DualOpSum{P}, kt::SingleKet{P}, i)
    k = label(kt)
    v = coeff(kt)
    k_i = k[i]
    for (o,c) in data(opc)
        add_to_dict!(result,
                     setindex(k, blabel(o)[1], i),
                     c'*v*P(klabel(o)[1], k_i))
    end
    return result
end

function act_on_dict!{P}(result, opc::DualOpSum{P}, kt::KetSum{P}, i)
    for (o,c) in data(opc), (k,v) in data(kt)
        add_to_dict!(result,
                     setindex(k, blabel(o)[1], i),
                     c'*v*P(klabel(o)[1], k[i]))
    end
    return result
end

##########
# tensor #
##########
tensor{P}(a::OpSum{P}, b::OpSum{P}) = OpSum(P, tensor(data(a), data(b)))
tensor(a::DualOpSum, b::DualOpSum) = tensor(a', b')'
tensor(a::DiracOp, b::DiracOp) = tensor(promote(a,b)...)

###########
# Scaling #
###########
Base.scale!(op::OpSum, c::Number) = (scale!(data(op), c); return op)
Base.scale!(c::Number, op::OpSum) = scale!(op, c)
Base.scale!(opc::DualOpSum, c::Number) = (scale!(opc.op, c'); return opc)
Base.scale!(c::Number, opc::DualOpSum) = scale!(opc, c)

Base.scale{P}(op::OpSum{P}, c::Number) = OpSum(P, scale(data(op), c))
Base.scale(c::Number, op::OpSum) = scale(op, c)
Base.scale(opc::DualOpSum, c::Number) = scale(opc', c')'
Base.scale(c::Number, opc::DualOpSum) = scale(opc, c)

Base.(:*)(c::Number, op::DiracOp) = scale(c, op)
Base.(:*)(op::DiracOp, c::Number) = scale(op, c)
Base.(:/)(op::DiracOp, c::Number) = scale(op, 1/c)

###########
# + and - #
###########
Base.(:+)(a::DiracOp, b::DiracOp) = +(promote(a,b)...)
Base.(:-)(a::DiracOp, b::DiracOp) = a + (-b)

Base.(:-)(op::OpSum) = scale(-1, op)
Base.(:-)(opc::DualOpSum) = ctranspose(-(opc'))

Base.(:+){P,N}(a::OpSum{P,N}, b::OpSum{P,N}) = OpSum(P, data(a) + data(b))
Base.(:-){P,N}(a::OpSum{P,N}, b::OpSum{P,N}) = OpSum(P, data(a) - data(b))

Base.(:+){P,N}(a::DualOpSum{P,N}, b::DualOpSum{P,N}) = ctranspose(a' + b')
Base.(:-){P,N}(a::DualOpSum{P,N}, b::DualOpSum{P,N}) = ctranspose(a' - b')

#################
# Normalization #
#################
Base.norm(op::OpSum) = sqrt(sum(abs2, values(data(op))))
Base.norm(opc::DualOpSum) = norm(opc')

normalize(op::DiracOp) = scale(1/norm(op), op)
normalize!(op::DiracOp) = scale!(1/norm(op), op)

#########
# Trace #
#########
function Base.trace(op::OpSum)
    result = predict_zero(eltype(op))
    for (o,v) in data(op)
        if klabel(o)==blabel(o)
            result += v
        end
    end
    return result
end

Base.trace(opc::DualOpSum) = trace(opc')'

#################
# Partial trace #
#################
ptrace{P}(op::DiracOp{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())

function ptrace_result{P,N}(op::DiracOp{P,N})
    K = eltype(ketlabeltype(op))
    B = eltype(bralabeltype(op))
    return SumDict{OpLabel{N-1,K,B}, eltype(op)}()
end

ptrace{P,N}(op::DiracOp{P,N}, over) = OpSum(P, ptrace_dict!(ptrace_result(op), op, over))

function ptrace_dict!(result, op::OpSum, over)
    for (o,v) in data(op)
        if klabel(o)[over] == blabel(o)[over]
            add_to_dict!(result, except(o, over), v)
        end
    end
    return result
end

function ptrace_dict!(result, opc::DualOpSum, over)
    for (o,v) in data(opc)
        if blabel(o)[over] == klabel(o)[over]
            add_to_dict!(result, except(o', over), v')
        end
    end
    return result
end

#####################
# Partial Transpose #
#####################
function ptrans_result{P,N}(op::DiracOp{P,N})
    T = eltype(promote_type(ketlabeltype(op), bralabeltype(op)))
    return SumDict{OpLabel{N,T,T}, eltype(op)}()
end

function ptrans_result{P}(op::DiracOp{P,1})
    K = eltype(bralabeltype(op))
    B = eltype(ketlabeltype(op))
    return SumDict{OpLabel{N,K,B}, eltype(op)}()
end

ptranspose{P,N}(op::DiracOp{P,N}, over) = OpSum(P, ptrans_dict!(ptrans_result(op), op, over))

function ptrans_dict!(result, op::OpSum, over)
    for (o,v) in data(op)
        add_to_dict!(result, ptranspose(o, over), v)
    end
    return result
end

function ptrans_dict!(result, opc::DualOpSum, over)
    for (o,v) in data(opc)
        add_to_dict!(result, ptranspose(o', over), v')
    end
    return result
end

########################
# Misc. Math Functions #
########################
nfactors{P,N}(::AbsOpSum{P,N}) = N

xsubspace(op::AbsOpSum, x) = filter((k,v)->is_sum_x(k,x), s)
switch(op::AbsOpSum, i, j) = maplabels(l->switch(l,i,j), op)
permute(op::AbsOpSum, perm::Vector) = maplabels(l->permute(l,perm), op)

filternz!(op::AbsOpSum) = (filter!(nzcoeff, data(op)); return op)
filternz{P}(op::OpSum{P}) = OpSum(P, filternz(data(op)))
filternz(opc::DualOpSum) = filternz(opc')'

purity(op::DiracOp) = trace(op^2)
purity(op::DiracOp, i) = purity(ptrace(op,i))

commute(a::DiracOp, b::DiracOp) = (a*b) - (b*a)
anticommute(a::DiracOp, b::DiracOp) = (a*b) + (b*a)

function represent{P}(op::DiracOp{P}, basis)
    T = inner_rettype(op)
    return T[bra(P, i) * op * ket(P, j) for i in basis, j in basis]
end

function represent{P}(op::DiracOp{P}, basis...)
    prodbasis = product(basis...)
    T = inner_rettype(op)
    return T[bra(P, i...) * op * ket(P, j...) for i in prodbasis, j in prodbasis]
end

######################
# Printing Functions #
######################
labelrepr(op::OpSum, o::OpLabel, pad) = "$pad$(op[o]) $(ktstr(klabel(o)))$(brstr(blabel(o)))"
labelrepr(opc::DualOpSum, o::OpLabel, pad) = "$pad$(opc[o']) $(ktstr(blabel(o)))$(brstr(klabel(o)))"

Base.summary(op::DiracOp) = "$(typeof(op)) with $(length(op)) operator(s)"
Base.show(io::IO, op::AbsOpSum) = dirac_show(io, op)
Base.showcompact(io::IO, op::AbsOpSum) = dirac_showcompact(io, op)
Base.repr(op::AbsOpSum) = dirac_repr(op)

export OpSum,
    ptrace,
    ptranspose,
    represent,
    xsubspace,
    nfactors,
    maplabels!,
    mapcoeffs!,
    mapcoeffs,
    maplabels,
    filternz,
    filternz!,
    purity,
    act_on,
    inner_eval
