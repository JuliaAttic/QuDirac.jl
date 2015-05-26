abstract OuterOp{P,N} <: DiracOp{P,N}
abstract AbsOuterSum{P,N,K,B,T} <: OuterOp{P,N}

#####################
# OuterProduct Type #
#####################
immutable OuterProduct{P,N,S,K<:Ket,B<:Bra} <: OuterOp{P,N}
    scalar::S
    kt::K
    br::B
    OuterProduct(scalar::S, kt::Ket{P,N}, br::Bra{P,N}) = new(scalar, kt, br)
end

OuterProduct{P,N,S}(scalar::S, kt::Ket{P,N}, br::Bra{P,N}) = OuterProduct{P,N,S,typeof(kt),typeof(br)}(scalar, kt, br)

#################
# OuterSum Type #
#################
type OuterSum{P,N,K,B,T} <: AbsOuterSum{P,N,K,B,T}
    data::SumDict{OpLabel{N,K,B},T}
end

OuterSum(op::OuterSum) = op
OuterSum{P,N,K,B,T}(::Type{P}, data::SumDict{OpLabel{N,K,B},T}) = OuterSum{P,N,K,B,T}(data)

function OuterSum{P,N}(kt::Ket{P,N}, br::Bra{P,N}, scalar = 1)
    T = promote_type(eltype(kt), eltype(br), typeof(scalar))
    K = labeltype(kt)
    B = labeltype(br)
    result = @compat sizehint!(SumDict{OpLabel{N,K,B},T}(), length(kt) * length(br))
    
    cons_outer!(result, kt, br)
    
    if scalar != 1
        scale!(result, scalar)
    end

    return OuterSum(P, result)
end

function cons_outer!(result::SumDict, kt::KetSum, br::BraSum)
    for (k,kc) in data(kt)
        for (b,bc) in data(br)
            add_to_sum!(result, OpLabel(k, b), kc * bc')
        end
    end
    return result
end

function cons_outer!(result::SumDict, kt::SingleKet, br::BraSum)
    k = label(kt)
    kc = coeff(kt)
    for (b,bc) in data(br)
        add_to_sum!(result, OpLabel(k, b), kc * bc')
    end
    return result
end

function cons_outer!(result::SumDict, kt::KetSum, br::SingleBra)
    b = label(br)
    bc = coeff(br)
    for (k,kc) in data(kt)
        add_to_sum!(result, OpLabel(k, b), kc * bc)
    end
    return result
end

function cons_outer!(result::SumDict, kt::SingleKet, br::SingleBra)
    add_to_sum!(result, OpLabel(label(kt), label(br)), coeff(kt) * coeff(br))
    return result
end

#####################
# DualOuterSum Type #
#####################
type DualOuterSum{P,N,K,B,T} <: AbsOuterSum{P,N,K,B,T}
    op::OuterSum{P,N,B,K,T}
end

##############
# ctranspose #
##############
eager_ctranspose{P}(op::OuterSum{P}) = OuterSum(P, SumDict(collect(op')))

Base.ctranspose(op::OuterSum) = DualOuterSum(op)
Base.ctranspose(opc::DualOuterSum) = opc.op
Base.ctranspose(op::OuterProduct) = OuterProduct(op.scalar', op.br', op.kt')

######################
# Accessor Functions #
######################
data(op::OuterSum) = op.data
data(opc::DualOuterSum) = data(opc')

Base.eltype(op::OuterOp) = eltype(typeof(op))
Base.eltype{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}) = promote_type(S, eltype(K), eltype(B))
Base.eltype{P,N,K,B,T}(::Type{OuterSum{P,N,K,B,T}}) = T
Base.eltype{P,N,K,B,T}(::Type{DualOuterSum{P,N,K,B,T}}) = T

ketlabeltype(op::OuterOp) = ketlabeltype(typeof(op))
ketlabeltype{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}) = labeltype(K)
ketlabeltype{P,N,K,B,T}(::Type{OuterSum{P,N,K,B,T}}) = K
ketlabeltype{P,N,K,B,T}(::Type{DualOuterSum{P,N,K,B,T}}) = K

bralabeltype(op::OuterOp) = bralabeltype(typeof(op))
bralabeltype{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}) = labeltype(B)
bralabeltype{P,N,K,B,T}(::Type{OuterSum{P,N,K,B,T}}) = B
bralabeltype{P,N,K,B,T}(::Type{DualOuterSum{P,N,K,B,T}}) = B

nfactors(op::OuterOp) = nfactors(typeof(op))
nfactors{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}) = N
nfactors{P,N,K,B,T}(::Type{OuterSum{P,N,K,B,T}}) = N
nfactors{P,N,K,B,T}(::Type{DualOuterSum{P,N,K,B,T}}) = N

########################
# Conversion/Promotion #
########################
OuterSum(opc::DualOuterSum) = eager_ctranspose(opc')
OuterSum(op::OuterProduct) = OuterSum(op.kt, op.br, op.scalar)

Base.convert{P,N,K,B,T}(::Type{OuterSum{P,N,K,B,T}}, op::OuterProduct) = convert(OuterSum{P,N,K,B,T}, OuterSum(op))
Base.convert(::Type{OuterSum}, op::OuterProduct) = OuterSum(op)
Base.convert{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}, op::OuterProduct{P,N,S,K,B}) = op

Base.convert{P,N,K,B,T}(::Type{OuterSum{P,N,K,B,T}}, op::OuterSum{P}) = OuterSum(P, convert(SumDict{OpLabel{N,K,B},T}, data(op)))
Base.convert{P,N,K,B,T}(::Type{OuterSum{P,N,K,B,T}}, op::OuterSum{P,N,K,B,T}) = op

Base.convert{P,N,K,B,T}(::Type{OuterSum{P,N,B,K,T}}, opc::DualOuterSum{P,N,K,B,T}) = OuterSum(opc)
Base.convert{P,N,K,B,T}(::Type{OuterSum}, opc::DualOuterSum{P,N,K,B,T}) = OuterSum(opc)
Base.convert{P,N,K,B,T}(::Type{DualOuterSum{P,N,K,B,T}}, opc::DualOuterSum) = DualOuterSum(convert(OuterSum{P,N,B,K,T}, opc'))
Base.convert{P,N,K,B,T}(::Type{DualOuterSum{P,N,K,B,T}}, opc::DualOuterSum{P,N,K,B,T}) = opc

Base.promote_rule{OS<:OuterSum, OP<:OuterProduct}(::Type{OS}, ::Type{OP}) = OS

function Base.promote_rule{P,N,K1,B1,T1,K2,B2,T2}(::Type{OuterSum{P,N,K1,B1,T1}}, ::Type{OuterSum{P,N,K2,B2,T2}})
    return OuterSum{P,N,promote_type(K1,K2),promote_type(B1,B2),promote_type(T1,T2)}
end

Base.promote_rule{P,N,K,B,T}(::Type{OuterSum{P,N,B,K,T}}, ::Type{DualOuterSum{P,N,K,B,T}}) = OuterSum{P,N,B,K,T}

function Base.promote_rule{P,N,K1,B1,T1,K2,B2,T2}(::Type{DualOuterSum{P,N,K1,B1,T1}}, ::Type{DualOuterSum{P,N,K2,B2,T2}})
    return DualOuterSum{P,N,promote_type(K1,K2),promote_type(B1,B2),promote_type(T1,T2)}
end

############################
# Hashing/Equality/Copying #
############################
const outer_hash = hash(OuterSum)

Base.hash{P}(op::OuterSum{P}) = hash(P, hash(data(op), outer_hash))
Base.hash(op::OuterOp) = hash(OuterSum(op))
Base.hash(op::OuterOp, h::Uint64) = hash(hash(op), h)

Base.copy{P}(op::OuterSum{P}) = OuterSum(P, copy(data(op)))
Base.copy(opc::DualOuterSum) = ctranspose(copy(opc'))
Base.copy(op::OuterProduct) = OuterProduct(copy(op.scalar), copy(op.kt), copy(op.br))

Base.(:(==)){P,N}(a::OuterSum{P,N}, b::OuterSum{P,N}) = data(a) == data(b)
Base.(:(==)){P,N}(a::DualOuterSum{P,N}, b::DualOuterSum{P,N}) = a' == b'
Base.(:(==)){P,N}(a::DiracOp{P,N}, b::DiracOp{P,N}) = OuterSum(a) == OuterSum(b)

#######################
# Dict-like Functions #
#######################
Base.length(op::AbsOuterSum) = length(data(op))
Base.length(op::OuterProduct) = length(op.kt)*length(op.br)

Base.getindex(op::OuterProduct, k, b) = op.scalar * op.kt[k] * op.br[b]
Base.getindex(op::OuterProduct, o::OpLabel) = op[klabel(o), blabel(o)]

Base.getindex(op::OuterSum, x::OpLabel) = getindex(data(op), x)
Base.getindex(opc::DualOuterSum, x::OpLabel) = getindex(opc', x')'
Base.getindex(op::AbsOuterSum, k::StateLabel, b::StateLabel) = op[OpLabel(k,b)]
Base.getindex(op::AbsOuterSum, k, b) = op[StateLabel(k),StateLabel(b)]

Base.setindex!(op::OuterSum, x, y::OpLabel) = setindex!(data(op), x, y)
Base.setindex!(opc::DualOuterSum, x, y::OpLabel) = setindex!(opc', x', y')
Base.setindex!(op::AbsOuterSum, x, k::StateLabel, b::StateLabel) = setindex!(op, x, OpLabel(k,b))
Base.setindex!(op::AbsOuterSum, x, k, b) = setindex!(op, x, StateLabel(k), StateLabel(b))

Base.haskey(op::OuterProduct, k, b) = haskey(op.kt, k) && haskey(op.bra, b)
Base.haskey(op::OuterProduct, o::OpLabel) = haskey(op, klabel(o), blabel(o))

Base.haskey(op::OuterSum, x::OpLabel) = haskey(data(op), x)
Base.haskey(opc::DualOuterSum, x::OpLabel) = haskey(opc', x')
Base.haskey(op::AbsOuterSum, k::StateLabel, b::StateLabel) = haskey(op, OpLabel(k,b))
Base.haskey(op::AbsOuterSum, k, b) = haskey(op, StateLabel(k), StateLabel(b))

Base.get(op::OuterProduct, k, b, default=predict_zero(eltype(op))) = haskey(op, k, b) ? op[k,b] : default
Base.get(op::OuterProduct, o::OpLabel, default=predict_zero(eltype(op))) = get(op, klabel(o), blabel(o), default)

Base.get(op::OuterSum, x::OpLabel, default=predict_zero(eltype(op))) = get(data(op), x, default)
Base.get(opc::DualOuterSum, x::OpLabel, default=predict_zero(eltype(opc))) = get(opc', x', default)
Base.get(op::AbsOuterSum, k::StateLabel, b::StateLabel, default=predict_zero(eltype(op))) = get(op, OpLabel(k,b), default)
Base.get(op::AbsOuterSum, k, b, default=predict_zero(eltype(op))) = get(op, StateLabel(k), StateLabel(b), default)

#############
# Iteration #
#############
Base.start(op::AbsOuterSum) = start(data(op))

Base.next(op::OuterSum, i) = next(data(op), i)

function Base.next(opc::DualOuterSum, i)
    (k,v), n = next(data(opc), i)
    return ((k',v'), n)
end

Base.done(op::AbsOuterSum, i) = done(data(op), i)

###########
# collect #
###########
Base.collect(op::OuterSum) = collect(data(op))

function Base.collect(op::Union(OuterProduct, DualOuterSum))
    T = @compat Tuple{OpLabel{nfactors(op),ketlabeltype(op),bralabeltype(op)}, eltype(op)}
    return collect_pairs!(Array(T, length(op)), op)
end

function collect_pairs!(result, op::OuterProduct)
    i = 1
    for (k,kc) in data(op.kt), (b,bc) in data(op.br)
        result[i] = tuple(OpLabel(k, b), op.scalar * kc * bc')
        i += 1
    end
    return result
end

function collect_pairs!(result, opc::DualOuterSum)
    i = 1
    for (k,v) in data(opc)
        result[i] = tuple(k', v')
        i += 1
    end
    return result
end

#########
# inner #
#########
function inner{P,N}(br::Bra{P,N}, op::OuterSum{P,N})
    result = SumDict{StateLabel{N, bralabeltype(op)}, inner_rettype(br, op)}()
    return ctranspose(KetSum(P, inner_load!(result, br, op)))
end

function inner_load!{P}(result::SumDict, br::SingleBra{P}, op::OuterSum{P})
    c = coeff(br)
    b = label(br)
    for (o,v) in data(op)
        add_to_sum!(result, blabel(o), ctranspose(v*c*inner(P, b, klabel(o))))
    end
    return result
end

function inner_load!{P}(result::SumDict, br::BraSum{P}, op::OuterSum{P})
    for (o,v) in data(op)
        k = klabel(o)
        coeff = predict_zero(eltype(result))
        for (b,c) in data(br)
            coeff += c' * v * inner(P, b, k) 
        end
        add_to_sum!(result, blabel(o), ctranspose(coeff))
    end
    return result
end

function inner{P,N,A,B}(op::OuterSum{P,N,A}, kt::Ket{P,N,B})
    result = SumDict{StateLabel{N, ketlabeltype(op)}, inner_rettype(op, kt)}()
    return KetSum(P, inner_load!(result, op, kt))
end

function inner_load!{P}(result::SumDict, op::OuterSum{P}, kt::SingleKet{P})
    c = coeff(kt)
    k = label(kt)
    for (o,v) in data(op)
        add_to_sum!(result, klabel(o), v*c*inner(P, blabel(o), k))
    end
    return result
end

function inner_load!{P}(result::SumDict, op::OuterSum{P}, kt::KetSum{P})
    for (o,v) in data(op)
        b = blabel(o)
        coeff = predict_zero(eltype(result))
        for (k,c) in data(kt)
            coeff += c * v * inner(P, b, k) 
        end
        add_to_sum!(result, klabel(o), coeff)
    end
    return result
end

function inner{P,N,K1,B1,K2,B2}(a::OuterSum{P,N,K1,B1}, b::OuterSum{P,N,K2,B2})
    result = SumDict{OpLabel{N,K1,B2}, inner_rettype(a, b)}()
    return OuterSum(P, inner_load!(result, a, b))
end

function inner_load!{P}(result::SumDict, a::OuterSum{P}, b::OuterSum{P})
    for (o1,v) in data(a), (o2,c) in data(b)
        add_to_sum!(result, 
                    OpLabel(klabel(o1), blabel(o2)),
                    v*c*inner(P, blabel(o1), klabel(o2)))
    end
    return result
end

function inner{P,N,K1,B1,K2,B2}(op::OuterSum{P,N,K1,B1}, opc::DualOuterSum{P,N,K2,B2})
    result = SumDict{OpLabel{N,K1,B2}, inner_rettype(op, opc)}()
    return OuterSum(P, inner_load!(result, op, opc))
end

function inner_load!{P}(result::SumDict, op::OuterSum{P}, opc::DualOuterSum{P})
    for (o,v) in data(op), (oc,c) in data(opc)
        add_to_sum!(result, 
                    OpLabel(klabel(o), klabel(oc)),
                    v*c'*inner(P, blabel(o), blabel(oc)))
    end
    return result
end

function inner{P,N,K1,B1,K2,B2}(opc::DualOuterSum{P,N,K1,B1}, op::OuterSum{P,N,K2,B2})
    result = SumDict{OpLabel{N,K1,B2}, inner_rettype(opc, op)}()
    return OuterSum(P, inner_load!(result, opc, op))
end

function inner_load!{P}(result::SumDict, opc::DualOuterSum{P}, op::OuterSum{P})
    for (oc,v) in data(opc), (o,c) in data(op)
        add_to_sum!(result,
                    OpLabel(blabel(oc), blabel(o)),
                    v'*c*inner(P, klabel(oc), klabel(o)))
    end
    return result
end

inner(br::Bra, opc::DualOuterSum) = inner(opc', br')'
inner(opc::DualOuterSum, kt::Ket) = inner(kt', opc')'
inner(a::DualOuterSum, b::DualOuterSum) = inner(b', a')'

inner(br::Bra, op::OuterProduct) = op.scalar * inner(br, op.kt) * op.br
inner(op::OuterProduct, kt::Ket) = op.scalar * op.kt * inner(op.br, kt)
inner(a::OuterProduct, b::OuterProduct) = OuterProduct(a.scalar * b.scalar * inner(a.br,b.kt), a.kt, b.br)
inner(a::OuterProduct, b::AbsOuterSum) = a.scalar * a.kt * inner(a.br, b)
inner(a::AbsOuterSum, b::OuterProduct) = inner(a, b.kt) * b.br * b.scalar

Base.(:*)(br::Bra, op::DiracOp) = inner(br,op)
Base.(:*)(op::DiracOp, kt::Ket) = inner(op,kt)
Base.(:*)(a::DiracOp, b::DiracOp) = inner(a,b)

inner_eval(f, op::DiracOp) = mapcoeffs(x->inner_eval(f,x),op)

##########
# act_on #
##########
act_on(op::OuterProduct, state::DiracState, i) = act_on(OuterSum(op), state, i)
act_on(op::AbsOuterSum, br::Bra, i) = act_on(op', br', i)'
act_on{P}(op::AbsOuterSum{P,1}, kt::Ket{P,1}, i) = i==1 ? inner(op, kt) : throw(BoundsError())

function act_result{P,N}(op::AbsOuterSum{P,1}, kt::Ket{P,N})
    T = promote_type(ketlabeltype(op), labeltype(kt))
    return SumDict{StateLabel{N,T},inner_rettype(op, kt)}()
end

act_on{P,N}(op::AbsOuterSum{P,1}, kt::Ket{P,N}, i) = KetSum(P, act_on_dict!(act_result(op,kt), op, kt, i))

function act_on_dict!{P}(result, op::OuterSum{P}, kt::SingleKet{P}, i)
    k = label(kt)
    v = coeff(kt)
    k_i = k[i]
    for (o,c) in data(op)
        add_to_sum!(result, 
                    setindex(k, klabel(o)[1], i),
                    c*v*P(blabel(o)[1], k_i))
    end
    return result
end

function act_on_dict!{P}(result, op::OuterSum{P}, kt::KetSum{P}, i)
    for (o,c) in data(op), (k,v) in data(kt)
        add_to_sum!(result, 
                    setindex(k, klabel(o)[1], i),
                    c*v*P(blabel(o)[1], k[i]))
    end
    return result
end

function act_on_dict!{P}(result, opc::DualOuterSum{P}, kt::SingleKet{P}, i)
    k = label(kt)
    v = coeff(kt)
    k_i = k[i]
    for (o,c) in data(opc)
        add_to_sum!(result,
                    setindex(k, blabel(o)[1], i),
                    c'*v*P(klabel(o)[1], k_i))
    end
    return result
end

function act_on_dict!{P}(result, opc::DualOuterSum{P}, kt::KetSum{P}, i)
    for (o,c) in data(opc), (k,v) in data(kt)
        add_to_sum!(result,
                    setindex(k, blabel(o)[1], i),
                    c'*v*P(klabel(o)[1], k[i]))
    end
    return result
end

##########
# tensor #
##########
tensor{P}(a::OuterSum{P}, b::OuterSum{P}) = OuterSum(P, tensor(data(a), data(b)))
tensor(a::DualOuterSum, b::DualOuterSum) = tensor(a', b')'
tensor(a::OuterProduct, b::OuterProduct) = OuterProduct(a.scalar * b.scalar, tensor(a.kt,b.kt), tensor(a.br, b.br))
tensor(a::DiracOp, b::DiracOp) = tensor(OuterSum(a), OuterSum(b))

tensor(kt::Ket, br::Bra) = OuterProduct(1, kt, br)
tensor(br::Bra, kt::Ket) = tensor(kt, br)

Base.(:*)(kt::Ket, br::Bra) = tensor(kt,br)

###########
# Scaling #
###########
Base.scale!(op::OuterSum, c::Number) = (scale!(data(op), c); return op)
Base.scale!(opc::DualOuterSum, c::Number) = (scale!(opc.op, c'); return opc)
Base.scale!(c::Number, op::AbsOuterSum) = scale!(op, c)

Base.scale{P}(op::OuterSum{P}, c::Number) = OuterSum(P, scale(data(op), c))
Base.scale(opc::DualOuterSum, c::Number) = scale(opc', c')'
Base.scale(op::OuterProduct, c::Number) = OuterProduct(op.scalar * c, copy(op.kt), copy(op.br))
Base.scale(c::Number, op::OuterOp) = scale(op, c)

Base.(:*)(c::Number, op::DiracOp) = scale(c, op)
Base.(:*)(op::DiracOp, c::Number) = scale(op, c)
Base.(:/)(op::DiracOp, c::Number) = scale(op, 1/c)

###########
# + and - #
###########
add!{P,N}(a::OuterSum{P,N}, b::OuterSum{P,N}) = (add!(data(a), data(b)); return a)
add!{P,N}(a::DualOuterSum{P,N}, b::DualOuterSum{P,N}) = (add!(data(a), data(b)); return a)
add!{P,N}(a::OuterSum{P,N}, b::OuterOp{P,N}) = add!(a, OuterSum(b))
add!{P,N}(a::DualOuterSum{P,N}, b::OuterOp{P,N}) = add!(a, OuterSum(b)')

sub!{P,N}(a::OuterSum{P,N}, b::OuterSum{P,N}) = (sub!(data(a), data(b)); return a)
sub!{P,N}(a::DualOuterSum{P,N}, b::DualOuterSum{P,N}) = (sub!(data(a), data(b)); return a)
sub!{P,N}(a::OuterSum{P,N}, b::OuterOp{P,N}) = sub!(a, OuterSum(b))
sub!{P,N}(a::DualOuterSum{P,N}, b::OuterOp{P,N}) = sub!(a, OuterSum(b)')

Base.(:+){P,N}(a::OuterSum{P,N}, b::OuterSum{P,N}) = OuterSum(P, data(a) + data(b))
Base.(:+){P,N}(a::DualOuterSum{P,N}, b::DualOuterSum{P,N}) = ctranspose(a' + b')
Base.(:+){P,N}(a::DiracOp{P,N}, b::DiracOp{P,N}) = OuterSum(a) + OuterSum(b)

Base.(:-){P,N}(a::OuterSum{P,N}, b::OuterSum{P,N}) = OuterSum(P, data(a) - data(b))
Base.(:-){P,N}(a::DualOuterSum{P,N}, b::DualOuterSum{P,N}) = ctranspose(a' - b')
Base.(:-){P,N}(a::DiracOp{P,N}, b::DiracOp{P,N}) = OuterSum(a) - OuterSum(b)

Base.(:-)(op::OuterOp) = scale(-1, op)

#################
# Normalization #
#################
Base.norm(op::OuterSum) = sqrt(sum(abs2, values(data(op))))
Base.norm(opc::DualOuterSum) = norm(opc')
function Base.norm(op::OuterProduct)
    result = predict_zero(eltype(op))
    for v in values(data(op.kt)), c in values(data(op.br))
        result += abs2(op.scalar * v * c')
    end
    return sqrt(result)
end

normalize(op::DiracOp) = scale(1/norm(op), op)
normalize!(op::DiracOp) = scale!(1/norm(op), op)

#########
# Trace #
#########
function Base.trace(op::OuterSum)
    result = predict_zero(eltype(op))
    for (o,v) in data(op)
        if klabel(o)==blabel(o)
            result += v
        end
    end
    return result
end

Base.trace(opc::DualOuterSum) = trace(opc')'

function Base.trace(op::OuterProduct)
    result = predict_zero(eltype(op))
    for (k,v) in data(op.kt), (b,c) in data(op.br)
        if b == k
            result += v * c'
        end
    end
    return op.scalar * result
end

#################
# Partial trace #
#################
ptrace{P}(op::DiracOp{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())

function ptrace_result{P,N}(op::DiracOp{P,N})
    K = ketlabeltype(op)
    B = bralabeltype(op)
    return SumDict{OpLabel{N-1,K,B}, eltype(op)}()
end

ptrace{P,N}(op::DiracOp{P,N}, over) = OuterSum(P, ptrace_dict!(ptrace_result(op), op, over))

function ptrace_dict!(result, op::OuterSum, over)
    for (o,v) in data(op)
        if klabel(o)[over] == blabel(o)[over]
            add_to_sum!(result, except(o, over), v)
        end
    end
    return result
end

function ptrace_dict!(result, opc::DualOuterSum, over)
    for (o,v) in data(opc)
        if blabel(o)[over] == klabel(o)[over]
            add_to_sum!(result, except(o', over), v')
        end
    end
    return result
end

function ptrace_dict!(result, op::OuterProduct, over)
    for k in keys(data(op.kt)), b in keys(data(op.br))
        if k[over] == b[over]
            add_to_sum!(result, except(OpLabel(k, b), over), op[k,b])
        end
    end
    return result
end

#####################
# Partial Transpose #
#####################
function ptrans_result{P,N}(op::DiracOp{P,N})
    T = label_promote(ketlabeltype(op), bralabeltype(op))
    return SumDict{OpLabel{N,T,T}, eltype(op)}()
end

function ptrans_result{P}(op::DiracOp{P,1})
    K = bralabeltype(op)
    B = ketlabeltype(op)
    return SumDict{OpLabel{N,K,B}, eltype(op)}()
end

ptranspose{P,N}(op::DiracOp{P,N}, over) = OuterSum(P, ptrans_dict!(ptrans_result(op), op, over))

function ptrans_dict!(result, op::OuterSum, over)
    for (o,v) in data(op)
        add_to_sum!(result, ptranspose(o, over), v)
    end
    return result
end

function ptrans_dict!(result, opc::DualOuterSum, over)
    for (o,v) in data(opc)
        add_to_sum!(result, ptranspose(o', over), v')
    end
    return result
end

function ptrans_dict!(result, op::OuterProduct, over)
    for k in keys(data(op.kt)), b in keys(data(op.br))
        add_to_sum!(result, ptranspose(OpLabel(k, b), over), op[k,b])
    end
    return result
end

########################
# Misc. Math Functions #
########################
xsubspace(op::AbsOuterSum, x) = filter((k,v)->is_sum_x(k,x), op)
xsubspace(op::OuterProduct, x) = xsubspace(OuterSum(op), x)

switch(op::AbsOuterSum, i, j) = maplabels(l->switch(l,i,j), op)
switch(op::OuterProduct, i, j) = switch(OuterSum(op), i, j)

permute(op::AbsOuterSum, perm::Vector) = maplabels(l->permute(l,perm), op)
permute(op::OuterProduct, perm::Vector) = permute(OuterSum(op), perm)

purity(op::DiracOp) = trace(op^2)
purity(op::DiracOp, i) = purity(ptrace(op,i))

commute(a::DiracOp, b::DiracOp) = (a*b) - (b*a)
anticommute(a::DiracOp, b::DiracOp) = (a*b) + (b*a)

represent{P}(op::DiracOp{P}, basis) = [bra(P, i) * op * ket(P, j) for i in basis, j in basis]

function represent{P}(op::DiracOp{P}, basis...)
    prodbasis = product(basis...)
    return [bra(P, i...) * op * ket(P, j...) for i in prodbasis, j in prodbasis]
end

export OuterSum,
    ptrace,
    ptranspose,
    represent,
    xsubspace,
    switch,
    permute,
    nfactors,
    purity,
    act_on,
    inner_eval,
    add!,
    sub!
