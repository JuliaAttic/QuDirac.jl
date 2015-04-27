###################
# OpSum/DualOpSum #
###################
abstract AbsOpSum{P,N,T} <: DiracOp{P,N}

typealias OpDict{N,T} Dict{OpLabel{N},T}

type OpSum{P,N,T} <: AbsOpSum{P,N,T}
    ptype::P
    dict::OpDict{N,T}
    OpSum(ptype, dict) = new(ptype, dict)
    OpSum(ptype, dict::OpDict{0}) = error("Cannot construct a 0-factor operator; did you mean to construct a scalar?")
end

OpSum{P,N,T}(ptype::P, dict::OpDict{N,T}) = OpSum{P,N,T}(ptype, dict)
OpSum{P,N,A,B}(kt::Ket{P,N,A}, br::Bra{P,N,B}) = OpSum(ptype(kt), cons_outer!(OpDict{N,promote_type(A,B)}(), kt, br))

function cons_outer!(result, kt, br)
    for (k,kc) in dict(kt)
        for (b,bc) in dict(br)
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
Base.promote_rule{O<:OpSum, D<:DualOpSum}(::Type{O}, ::Type{D}) = OpSum

######################
# Accessor functions #
######################
dict(op::OpSum) = op.dict
dict(opc::DualOpSum) = dict(opc.op)

ptype(op::OpSum) = op.ptype
ptype(opc::DualOpSum) = ptype(opc.op)

#######################
# Dict-Like Functions #
#######################
Base.eltype{P,N,T}(::AbsOpSum{P,N,T}) = T

Base.copy(op::OpSum) = OpSum(ptype(op), copy(dict(op)))
Base.copy(opc::DualOpSum) = DualOpSum(copy(opc.op))

Base.similar(op::OpSum, d=similar(dict(op)); P=ptype(op)) = OpSum(P, d)
Base.similar(opc::DualOpSum, d=similar(dict(opc)); P=ptype(opc)) = DualOpSum(P, d)

Base.(:(==)){P,N}(a::OpSum{P,N}, b::OpSum{P,N}) = ptype(a) == ptype(b) && dict(filternz(a)) == dict(filternz(b))
Base.(:(==)){P,N}(a::DualOpSum{P,N}, b::DualOpSum{P,N}) = a.op == b.op
Base.(:(==))(a::DiracOp, b::DiracOp) = ==(promote(a,b)...)

Base.hash(op::AbsOpSum) = hash(dict(filternz(op)), hash(ptype(op)))
Base.hash(op::AbsOpSum, h::Uint64) = hash(hash(op), h)

Base.length(op::AbsOpSum) = length(dict(op))

Base.getindex(op::OpSum, label::OpLabel) = op.dict[label]
Base.getindex(op::OpSum, k::StateLabel, b::StateLabel) = op.dict[OpLabel(k,b)]
Base.getindex(opc::DualOpSum, label::OpLabel) = opc.op[label']'
Base.getindex(opc::DualOpSum, k::StateLabel, b::StateLabel) = opc.op[OpLabel(b,k)]'
Base.getindex(op::AbsOpSum, k, b) = op[StateLabel(k), StateLabel(b)]

Base.setindex!(op::OpSum, c, label::OpLabel) = (op.dict[label] = c)
Base.setindex!(op::OpSum, c, k::StateLabel, b::StateLabel) = (op.dict[OpLabel(k,b)] = c)
Base.setindex!(opc::DualOpSum, c, label::OpLabel) = (opc.op[label'] = c')
Base.setindex!(opc::DualOpSum, c, k::StateLabel, b::StateLabel) = (opc.op[OpLabel(b,k)] = c')
Base.setindex!(op::AbsOpSum, c, k, b) = setindex!(op, c, StateLabel(k), StateLabel(b))

Base.haskey(op::OpSum, label::OpLabel) = haskey(dict(op), label)
Base.haskey(opc::DualOpSum, label::OpLabel) = haskey(opc.op, label')
Base.haskey(op::AbsOpSum, k, b) = haskey(op, OpLabel(k, b))

Base.get(op::OpSum, label::OpLabel, default=0) = get(dict(op), label, default)
Base.get(opc::DualOpSum, label::OpLabel, default=0) = get(dict(opc), label, default')'
Base.get(op::AbsOpSum, k, b, default=0) = get(op, OpLabel(k, b), default)

Base.delete!(op::OpSum, label::OpLabel) = (delete!(dict(op), label); return op)
Base.delete!(opc::DualOpSum, label::OpLabel) = delete!(opc.op, label')
Base.delete!(op::AbsOpSum, k, b) = delete!(op, OpLabel(k, b))

Base.collect(op::OpSum) = collect(dict(op))
Base.collect{P,N,T}(opc::DualOpSum{P,N,T}) = collect_pairs!(Array(@compat(Tuple{OpLabel{N}, T}), length(opc)), opc)

function collect_pairs!(result, opc::DualOpSum)
    i = 1
    for (k,v) in dict(opc)
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
    prodtype = ptype(op)
    result = StateDict{N, promote_type(A,B,inner_rettype(prodtype))}()
    return Bra(prodtype, inner_load!(result, br, op, prodtype))
end

function inner{P,N,A,B}(op::OpSum{P,N,A}, kt::Ket{P,N,B})
    prodtype = ptype(op)
    result = StateDict{N, promote_type(A,B,inner_rettype(prodtype))}()
    return Ket(prodtype, inner_load!(result, op, kt, prodtype))
end

function inner{P,N,A,B}(a::OpSum{P,N,A}, b::OpSum{P,N,B})
    prodtype = ptype(a)
    result = OpDict{N, promote_type(A,B,inner_rettype(prodtype))}()
    return OpSum(prodtype, inner_load!(result, a, b, prodtype))
end

function inner{P,N,A,B}(a::OpSum{P,N,A}, b::DualOpSum{P,N,B})
    prodtype = ptype(a)
    result = OpDict{N, promote_type(A,B,inner_rettype(prodtype))}()
    return OpSum(prodtype, inner_load!(result, a, b, prodtype))
end

function inner{P,N,A,B}(a::DualOpSum{P,N,A}, b::OpSum{P,N,B})
    prodtype = ptype(a)
    result = OpDict{N, promote_type(A,B,inner_rettype(prodtype))}()
    return OpSum(prodtype, inner_load!(result, a, b, prodtype))
end

inner(br::Bra, opc::DualOpSum) = inner(opc.op, br')'
inner(opc::DualOpSum, kt::Ket) = inner(kt', opc.op)'
inner(a::DualOpSum, b::DualOpSum) = inner(b.op, a.op)'

function inner_load!(result, a::OpSum, b::OpSum, prodtype)
    for (o1,v) in dict(a), (o2,c) in dict(b)
        add_to_dict!(result, 
                     OpLabel(klabel(o1), blabel(o2)),
                     inner_mul(v, c, prodtype, blabel(o1), klabel(o2)))
    end
    return result
end

function inner_load!(result, a::OpSum, b::DualOpSum, prodtype)
    for (o1,v) in dict(a), (o2,c) in dict(b)
        add_to_dict!(result, 
                     OpLabel(klabel(o1), klabel(o2)),
                     inner_mul(v, c', prodtype, blabel(o1), blabel(o2)))
    end
    return result
end

function inner_load!(result, a::DualOpSum, b::OpSum, prodtype)
    for (o1,v) in dict(a), (o2,c) in dict(b)
        add_to_dict!(result,
                     OpLabel(blabel(o1), blabel(o2)),
                     inner_mul(v', c, prodtype, klabel(o1), klabel(o2)))
    end
    return result
end

function inner_load!(result, br::Bra, op::OpSum, prodtype)
    for (o,v) in dict(op)
        add_to_dict!(result, blabel(o), brcoeff(dict(br), prodtype, klabel(o), v))
    end
    return result
end

function inner_load!(result, op::OpSum, kt::Ket, prodtype)
    for (o,v) in dict(op)
        add_to_dict!(result, klabel(o), ktcoeff(dict(kt), prodtype, blabel(o), v))
    end
    return result
end

function brcoeff{K,V}(brdict::Dict{K,V}, prodtype, klabel, v)
    coeff = predict_zero(promote_type(V, typeof(v), inner_rettype(prodtype)))
    for (blabel,c) in brdict
        coeff += inner_mul(c', v, prodtype, klabel, blabel) 
    end
    return coeff'
end

function ktcoeff{K,V}(ktdict::Dict{K,V}, prodtype, blabel, v)
    coeff = predict_zero(promote_type(V, typeof(v), inner_rettype(prodtype)))
    for (klabel,c) in ktdict
        coeff += inner_mul(c, v, prodtype, klabel, blabel)
    end
    return coeff
end

Base.(:*)(br::Bra, op::DiracOp) = inner(br,op)
Base.(:*)(op::DiracOp, kt::Ket) = inner(op,kt)
Base.(:*)(a::DiracOp, b::DiracOp) = inner(a,b)

###################################
# Functional Operator Application #
###################################
immutable DualFunc
    f::Function
end

Base.(:*)(op::Function, kt::Ket) = op(kt)
Base.(:*)(br::Bra, op::Function) = op(br)
Base.(:*)(op::DualFunc, kt::Ket) = (kt' * op.f)'
Base.(:*)(br::Bra, op::DualFunc) = (op.f * br')'

Base.ctranspose(f::Function) = DualFunc(f)
Base.ctranspose(fc::DualFunc) = fc.f

##############
# act/act_on #
##############
act_on(op::AbsOpSum, br::Bra, i) = act_on(op', br', i)'

# clear up ambiguity warnings
act_on{P}(op::OpSum{P,1}, kt::Ket{P,1}, i) = i==1 ? inner(op, kt) : throw(BoundsError())
act_on{P}(opc::DualOpSum{P,1}, kt::Ket{P,1}, i) = i==1 ? inner(opc, kt) : throw(BoundsError())

function act_on{P,N,A,B}(op::OpSum{P,1,A}, kt::Ket{P,N,B}, i)
    prodtype = ptype(op)
    result = StateDict{N, promote_type(A,B,inner_rettype(prodtype))}()
    return Ket(prodtype, act_on_dict!(result, op, kt, i, prodtype))
end

function act_on{P,N,A,B}(op::DualOpSum{P,1,A}, kt::Ket{P,N,B}, i)
    prodtype = ptype(op)
    result = StateDict{N, promote_type(A,B,inner_rettype(prodtype))}()
    return Ket(prodtype, act_on_dict!(result, op, kt, i, prodtype))
end

function act_on_dict!(result, op::OpSum, kt::Ket, i, prodtype)
    for (o,c) in dict(op), (k,v) in dict(kt)
        add_to_dict!(result, 
                     setindex(k, klabel(o)[1], i),
                     inner_mul(c, v, prodtype, blabel(o)[1], k[i]))
    end
    return result
end

function act_on_dict!(result, op::DualOpSum, kt::Ket, i, prodtype)
    for (o,c) in dict(op), (k,v) in dict(kt)
        add_to_dict!(result,
                     setindex(k, blabel(o)[1], i),
                     inner_mul(c', v, prodtype, klabel(o)[1], k[i]))
    end
    return result
end

##########
# tensor #
##########
tensor{P}(a::OpSum{P}, b::OpSum{P}) = OpSum(ptype(a), tensordict(dict(a), dict(b)))
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
Base.norm(op::OpSum) = sqrt(sum(abs2, values(dict(op))))
Base.norm(opc::DualOpSum) = norm(opc.op)

normalize(op::DiracOp) = scale(1/norm(op), op)
normalize!(op::DiracOp) = scale!(1/norm(op), op)

#########
# Trace #
#########
function Base.trace(op::OpSum)
    result = predict_zero(eltype(op))
    for (o,v) in dict(op)
        if klabel(o)==blabel(o)
            result += v
        end
    end
    return result
end

Base.trace(opc::DualOpSum) = trace(opc.op)'

#################
# Partial trace #
#################
ptrace{P}(op::DiracOp{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())
ptrace{P,N}(op::DiracOp{P,N}, over) = OpSum(ptype(op), ptrace_dict!(OpDict{N-1,eltype(op)}(), op, over))

function ptrace_dict!(result, op::OpSum, over)
    for (o,v) in dict(op)
        if klabel(o)[over] == blabel(o)[over]
            add_to_dict!(result, traceout(o, over), v)
        end
    end
    return result
end

function ptrace_dict!(result, opc::DualOpSum, over)
    for (o,v) in dict(opc)
        if blabel(o)[over] == klabel(o)[over]
            add_to_dict!(result, traceout_dual(o, over), v')
        end
    end
    return result
end

#####################
# Partial Transpose #
#####################
ptranspose{P,N}(op::DiracOp{P,N}, over) = OpSum(ptype(op), ptrans_dict!(OpDict{N,eltype(op)}(), op, over))

function ptrans_dict!(result, op::OpSum, over)
    for (o,v) in dict(op)
        add_to_dict!(result, ptranspose(o, over), v)
    end
    return result
end

function ptrans_dict!(result, opc::DualOpSum, over)
    for (o,v) in dict(opc)
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
switch(op::AbsOpSum, i, j) = maplabels(label->switch(label,i,j), op)
permute(op::AbsOpSum, perm::Vector) = maplabels(label->permute(label,perm), op)

purity(op::DiracOp) = trace(op^2)
commute(a::DiracOp, b::DiracOp) = (a*b) - (b*a)
anticommute(a::DiracOp, b::DiracOp) = (a*b) + (b*a)

inner_eval(f, op::DiracOp) = mapcoeffs(x->inner_eval(f,x),op)

function matrep(op::DiracOp, labels)
    T = promote_type(inner_rettype(ptype(op)), eltype(op))
    return T[bra(i) * op * ket(j) for i in labels, j in labels]
end

function matrep(op::DiracOp, labels...)
    iter = product(labels...)
    T = promote_type(inner_rettype(ptype(op)), eltype(op))
    return T[bra(i...) * op * ket(j...) for i in iter, j in iter]
end

function matrep(op::Union(DualFunc, Function), labels)
    return [bra(i) * op * ket(j) for i in labels, j in labels]
end

function matrep(op::Union(DualFunc, Function), labels...)
    iter = Iterators.product(labels...)
    return [bra(i...) * op * ket(j...) for i in iter, j in iter]
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
    matrep,
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
