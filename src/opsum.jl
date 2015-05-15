###################
# OpSum/DualOpSum #
###################
abstract AbsOpSum{P,N,T} <: DiracOp{P,N}

typealias OpDict{N,T} Dict{OpLabel{N},T}

type OpSum{P,N,T} <: AbsOpSum{P,N,T}
    dict::OpDict{N,T}
    OpSum(dict) = new(dict)
    OpSum(dict::OpDict{0}) = error("Cannot construct a 0-factor operator; did you mean to construct a scalar?")
end

OpSum{P,N,T}(::Type{P}, dict::OpDict{N,T}) = OpSum{P,N,T}(dict)
OpSum{P,N,A,B}(kt::Ket{P,N,A}, br::Bra{P,N,B}) = OpSum(P, cons_outer!(OpDict{N,promote_type(A,B)}(), kt, br))

function cons_outer!(result, kt, br)
    for (k,kc) in iter(kt)
        for (b,bc) in iter(br)
            newc = kc * bc'
            if newc != 0
                result[OpLabel(k, b)] = newc
            end
        end
    end
    return result
end

type DualOpSum{P,N,T} <: AbsOpSum{P,N,T}
    op::OpSum{P,N,T}
end

DualOpSum{P,N,T}(op::OpSum{P,N,T}) = DualOpSum{P,N,T}(op)
DualOpSum(items...) = DualOpSum(OpSum(items...))

Base.convert(::Type{OpSum}, opc::DualOpSum) = eager_ctran(opc.op)
Base.convert{P,N,T}(::Type{OpSum{P,N,T}}, opc::DualOpSum{P,N,T}) = convert(OpSum, opc)
Base.promote_rule{O<:OpSum, D<:DualOpSum}(::Type{O}, ::Type{D}) = OpSum

######################
# Accessor functions #
######################
dict(op::OpSum) = op.dict
dict(opc::DualOpSum) = dict(opc.op)

#######################
# Dict-Like Functions #
#######################
Base.eltype{P,N,T}(::AbsOpSum{P,N,T}) = T

Base.copy{P,N,T}(op::OpSum{P,N,T}) = OpSum{P,N,T}(copy(dict(op)))
Base.copy(opc::DualOpSum) = DualOpSum(copy(opc.op))

Base.similar{P}(op::OpSum{P}, d=similar(dict(op))) = OpSum(P, d)
Base.similar{P}(opc::DualOpSum{P}, d=similar(dict(opc))) = DualOpSum(P, d)

Base.(:(==)){P,N}(a::OpSum{P,N}, b::OpSum{P,N}) = dict(filternz(a)) == dict(filternz(b))
Base.(:(==)){P,N}(a::DualOpSum{P,N}, b::DualOpSum{P,N}) = a.op == b.op
Base.(:(==))(a::DiracOp, b::DiracOp) = ==(promote(a,b)...)

Base.hash{P}(op::AbsOpSum{P}) = hash(dict(filternz(op)), hash(P))
Base.hash(op::AbsOpSum, h::Uint64) = hash(hash(op), h)

Base.length(op::AbsOpSum) = length(dict(op))

Base.getindex(op::OpSum, ol::OpLabel) = op.dict[ol]
Base.getindex(op::OpSum, k::StateLabel, b::StateLabel) = op.dict[OpLabel(k,b)]
Base.getindex(opc::DualOpSum, ol::OpLabel) = opc.op[label']'
Base.getindex(opc::DualOpSum, k::StateLabel, b::StateLabel) = opc.op[OpLabel(b,k)]'
Base.getindex(op::AbsOpSum, k, b) = op[StateLabel(k), StateLabel(b)]

Base.setindex!(op::OpSum, c, ol::OpLabel) = (op.dict[ol] = c)
Base.setindex!(op::OpSum, c, k::StateLabel, b::StateLabel) = (op.dict[OpLabel(k,b)] = c)
Base.setindex!(opc::DualOpSum, c, ol::OpLabel) = (opc.op[ol'] = c')
Base.setindex!(opc::DualOpSum, c, k::StateLabel, b::StateLabel) = (opc.op[OpLabel(b,k)] = c')
Base.setindex!(op::AbsOpSum, c, k, b) = setindex!(op, c, StateLabel(k), StateLabel(b))

Base.haskey(op::OpSum, ol::OpLabel) = haskey(dict(op), ol)
Base.haskey(opc::DualOpSum, ol::OpLabel) = haskey(opc.op, ol')
Base.haskey(op::AbsOpSum, k, b) = haskey(op, OpLabel(k, b))

Base.get(op::OpSum, ol::OpLabel, default=predict_zero(eltype(op))) = get(dict(op), ol, default)
Base.get(opc::DualOpSum, ol::OpLabel, default=predict_zero(eltype(opc))) = get(dict(opc), ol, default')'
Base.get(op::AbsOpSum, k, b, default=predict_zero(eltype(op))) = get(op, OpLabel(k, b), default)

Base.delete!(op::OpSum, ol::OpLabel) = (delete!(dict(op), ol); return op)
Base.delete!(opc::DualOpSum, ol::OpLabel) = delete!(opc.op, ol')
Base.delete!(op::AbsOpSum, k, b) = delete!(op, OpLabel(k, b))

########################
# Iteration/Collection #
########################
iter(op::AbsOpSum) = dict(op)

coeffs(op::AbsOpSum) = values(dict(op))
labels(op::AbsOpSum) = keys(dict(op))

Base.collect(op::OpSum) = collect(dict(op))
Base.collect{P,N,T}(opc::DualOpSum{P,N,T}) = collect_pairs!(Array(@compat(Tuple{OpLabel{N}, T}), length(opc)), opc)

function collect_pairs!(result, opc::DualOpSum)
    i = 1
    for (k,v) in iter(opc)
        result[i] = ctpair(k,v)
        i += 1
    end
    return result
end

Base.start(op::AbsOpSum) = start(dict(op))
Base.next(op::OpSum, i) = next(dict(op), i)

function Base.next(opc::DualOpSum, i)
    (k,v), n = next(dict(opc), i)
    return ((k',v'), n)
end

Base.done(op::AbsOpSum, i) = done(dict(op), i)
Base.first(op::AbsOpSum) = next(op, start(op))

#############
# ctranpose #
#############
eager_ctran(op::OpSum) = similar(op, Dict(collect(op')))

Base.ctranspose(op::OpSum) = DualOpSum(op)
Base.ctranspose(opc::DualOpSum) = opc.op

#########
# inner #
#########
function inner{P,N,A,B}(br::Bra{P,N,A}, op::OpSum{P,N,B})
    result = StateDict{N, inner_coefftype(br, op)}()
    return KetSum(P, inner_load!(result, br, op, P))'
end

function inner{P,N,A,B}(op::OpSum{P,N,A}, kt::Ket{P,N,B})
    result = StateDict{N, inner_coefftype(op, kt)}()
    return KetSum(P, inner_load!(result, op, kt, P))
end

function inner{P,N,A,B}(a::OpSum{P,N,A}, b::OpSum{P,N,B})
    result = OpDict{N, inner_coefftype(a, b)}()
    return OpSum(P, inner_load!(result, a, b, P))
end

function inner{P,N,A,B}(a::OpSum{P,N,A}, b::DualOpSum{P,N,B})
    result = OpDict{N, inner_coefftype(a, b)}()
    return OpSum(P, inner_load!(result, a, b, P))
end

function inner{P,N,A,B}(a::DualOpSum{P,N,A}, b::OpSum{P,N,B})
    result = OpDict{N, inner_coefftype(a, b)}()
    return OpSum(P, inner_load!(result, a, b, P))
end

inner(br::Bra, opc::DualOpSum) = inner(opc.op, br')'
inner(opc::DualOpSum, kt::Ket) = inner(kt', opc.op)'
inner(a::DualOpSum, b::DualOpSum) = inner(b.op, a.op)'

function inner_load!{P}(result, a::OpSum, b::OpSum, ::Type{P})
    for (o1,v) in iter(a), (o2,c) in iter(b)
        add_to_dict!(result, 
                     OpLabel(klabel(o1), blabel(o2)),
                     v*c*P(blabel(o1), klabel(o2)))
    end
    return result
end

function inner_load!{P}(result, a::OpSum, b::DualOpSum, ::Type{P})
    for (o1,v) in iter(a), (o2,c) in iter(b)
        add_to_dict!(result, 
                     OpLabel(klabel(o1), klabel(o2)),
                     v*c'*P(blabel(o1), blabel(o2)))
    end
    return result
end

function inner_load!{P}(result, a::DualOpSum, b::OpSum, ::Type{P})
    for (o1,v) in iter(a), (o2,c) in iter(b)
        add_to_dict!(result,
                     OpLabel(blabel(o1), blabel(o2)),
                     v'*c*P(klabel(o1), klabel(o2)))
    end
    return result
end

function inner_load!{K,T,P}(result::Dict{K,T}, br::Bra, op::OpSum, ::Type{P})
    for (o,v) in iter(op)
        add_to_dict!(result, blabel(o), brcoeff(dict(br), P, klabel(o), v, T))
    end
    return result
end

function inner_load!{K,T,P}(result::Dict{K,T}, op::OpSum, kt::Ket, ::Type{P})
    for (o,v) in iter(op)
        add_to_dict!(result, klabel(o), ktcoeff(dict(kt), P, blabel(o), v, T))
    end
    return result
end

function brcoeff{T,P}(brdict, ::Type{P}, klabel, v, ::Type{T})
    coeff = predict_zero(T)
    for (blabel,c) in brdict
        coeff += c' * v * P(klabel, blabel) 
    end
    return coeff'
end

function ktcoeff{T,P}(ktdict, ::Type{P}, blabel, v, ::Type{T})
    coeff = predict_zero(T)
    for (klabel,c) in ktdict
        coeff += c * v * P(klabel, blabel) 
    end
    return coeff
end

Base.(:*)(br::Bra, op::DiracOp) = inner(br,op)
Base.(:*)(op::DiracOp, kt::Ket) = inner(op,kt)
Base.(:*)(a::DiracOp, b::DiracOp) = inner(a,b)

##############
# act/act_on #
##############
act_on(op::AbsOpSum, br::Bra, i) = act_on(op', br', i)'

# clear up ambiguity warnings
act_on{P}(op::OpSum{P,1}, kt::Ket{P,1}, i) = i==1 ? inner(op, kt) : throw(BoundsError())
act_on{P}(opc::DualOpSum{P,1}, kt::Ket{P,1}, i) = i==1 ? inner(opc, kt) : throw(BoundsError())

function act_on{P,N,A,B}(op::OpSum{P,1,A}, kt::Ket{P,N,B}, i)
    result = StateDict{N, inner_coefftype(op, kt)}()
    return KetSum(P, act_on_dict!(result, op, kt, i, P))
end

function act_on{P,N,A,B}(op::DualOpSum{P,1,A}, kt::Ket{P,N,B}, i)
    result = StateDict{N, inner_coefftype(op, kt)}()
    return KetSum(P, act_on_dict!(result, op, kt, i, P))
end

function act_on_dict!{P}(result, op::OpSum, kt::Ket, i, ::Type{P})
    for (o,c) in iter(op), (k,v) in iter(kt)
        add_to_dict!(result, 
                     setindex(k, klabel(o)[1], i),
                     c*v*P(blabel(o)[1], k[i]))
    end
    return result
end

function act_on_dict!{P}(result, op::DualOpSum, kt::Ket, i, ::Type{P})
    for (o,c) in iter(op), (k,v) in iter(kt)
        add_to_dict!(result,
                     setindex(k, blabel(o)[1], i),
                     c'*v*P(klabel(o)[1], k[i]))
    end
    return result
end

##########
# tensor #
##########
tensor{P}(a::OpSum{P}, b::OpSum{P}) = OpSum(P, tensor_merge(dict(a), dict(b)))
tensor(a::DualOpSum, b::DualOpSum) = tensor(a.opc, b.opc)'
tensor(a::DiracOp, b::DiracOp) = tensor(promote(a,b)...)

Base.(:*)(kt::Ket, br::Bra) = tensor(kt,br)

###########
# Scaling #
###########
Base.scale!(op::OpSum, c::Number) = (dscale!(dict(op), c); return op)
Base.scale!(c::Number, op::OpSum) = scale!(op, c)
Base.scale!(opc::DualOpSum, c::Number) = DualOpSum(scale!(opc.op, c'))
Base.scale!(c::Number, opc::DualOpSum) = scale!(opc, c)

Base.scale(op::OpSum, c::Number) = similar(op, dscale(dict(op), c))
Base.scale(c::Number, op::OpSum) = scale(op, c)
Base.scale(opc::DualOpSum, c::Number) = DualOpSum(scale(opc.op, c'))
Base.scale(c::Number, opc::DualOpSum) = scale(opc, c)

Base.(:*)(c::Number, op::DiracOp) = scale(c, op)
Base.(:*)(op::DiracOp, c::Number) = scale(op, c)
Base.(:/)(op::DiracOp, c::Number) = scale(op, 1/c)

###########
# + and - #
###########
Base.(:-)(op::OpSum) = scale(-1, op)
Base.(:-)(opc::DualOpSum) = DualOpSum(-opc.op)

Base.(:+){P,N}(a::OpSum{P,N}, b::OpSum{P,N}) = similar(b, add_merge(dict(a), dict(b)))
Base.(:-){P,N}(a::OpSum{P,N}, b::OpSum{P,N}) = similar(b, sub_merge(dict(a), dict(b)))

Base.(:+){P,N}(a::DualOpSum{P,N}, b::DualOpSum{P,N}) = DualOpSum(a.op + b.op)
Base.(:-){P,N}(a::DualOpSum{P,N}, b::DualOpSum{P,N}) = DualOpSum(a.op - b.op)

Base.(:+)(a::DiracOp, b::DiracOp) = +(promote(a,b)...)
Base.(:-)(a::DiracOp, b::DiracOp) = a + (-b)

#################
# Normalization #
#################
Base.norm(op::OpSum) = sqrt(sum(abs2, coeffs(op)))
Base.norm(opc::DualOpSum) = norm(opc.op)

normalize(op::DiracOp) = scale(1/norm(op), op)
normalize!(op::DiracOp) = scale!(1/norm(op), op)

#########
# Trace #
#########
function Base.trace(op::OpSum{KronDelta})
    result = predict_zero(eltype(op))
    for (o,v) in iter(op)
        if klabel(o)==blabel(o)
            result += v
        end
    end
    return result
end

Base.trace(opc::DualOpSum{KronDelta}) = trace(opc.op)'

#################
# Partial trace #
#################
ptrace{P}(op::DiracOp{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())
ptrace{P,N}(op::DiracOp{P,N}, over) = OpSum(P, ptrace_dict!(OpDict{N-1,eltype(op)}(), op, over))

function ptrace_dict!(result, op::OpSum, over)
    for (o,v) in iter(op)
        if klabel(o)[over] == blabel(o)[over]
            add_to_dict!(result, traceout(o, over), v)
        end
    end
    return result
end

function ptrace_dict!(result, opc::DualOpSum, over)
    for (o,v) in iter(opc)
        if blabel(o)[over] == klabel(o)[over]
            add_to_dict!(result, traceout_dual(o, over), v')
        end
    end
    return result
end

#####################
# Partial Transpose #
#####################
ptranspose{P,N}(op::DiracOp{P,N}, over) = OpSum(P, ptrans_dict!(OpDict{N,eltype(op)}(), op, over))

function ptrans_dict!(result, op::OpSum, over)
    for (o,v) in iter(op)
        add_to_dict!(result, ptranspose(o, over), v)
    end
    return result
end

function ptrans_dict!(result, opc::DualOpSum, over)
    for (o,v) in iter(opc)
        add_to_dict!(result, ptranspose_dual(o, over), v')
    end
    return result
end

########################
# Misc. Math Functions #
########################
nfactors{P,N}(::AbsOpSum{P,N}) = N
xsubspace(op::AbsOpSum, x) = similar(op, filter((k,v)->is_sum_x(k,x), dict(op)))
filternz!(op::AbsOpSum) = (filter!(nzcoeff, dict(op)); return op)
filternz(op::AbsOpSum) = similar(op, filter(nzcoeff, dict(op)))
switch(op::AbsOpSum, i, j) = maplabels(ol->switch(ol,i,j), op)
permute(op::AbsOpSum, perm::Vector) = maplabels(ol->permute(ol,perm), op)

purity(op::DiracOp) = trace(op^2)
purity(op::DiracOp, i) = purity(ptrace(op,i))

commute(a::DiracOp, b::DiracOp) = (a*b) - (b*a)
anticommute(a::DiracOp, b::DiracOp) = (a*b) + (b*a)

inner_eval(f, op::DiracOp) = mapcoeffs(x->inner_eval(f,x),op)

function represent{P}(op::DiracOp{P}, basis)
    T = promote_type(return_type(P), eltype(op))
    return T[bra(P, i) * op * ket(P, j) for i in basis, j in basis]
end

function represent{P}(op::DiracOp{P}, basis...)
    prodbasis = product(basis...)
    T = promote_type(return_type(P), eltype(op))
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
