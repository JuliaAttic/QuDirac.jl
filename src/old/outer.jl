abstract OuterOp{P} <: DiracOp{P}
abstract AbsOuterSum{P,T,K,B} <: OuterOp{P}

#####################
# OuterProduct Type #
#####################
immutable OuterProduct{P,S,K<:Ket,B<:Bra} <: OuterOp{P}
    coeff::S
    kt::K
    br::B
    function OuterProduct(kt::Ket{P}, br::Bra{P}, coeff::S)
        @assert matching_nfactors(kt, br)
        return new(coeff, kt, br)
    end
end

OuterProduct{P,S}(kt::Ket{P}, br::Bra{P}, coeff::S=true) = OuterProduct{P,S,typeof(kt),typeof(br)}(kt, br, coeff)
OuterProduct{P}(::Type{P}, ol::OuterLabel, coeff=true) = OuterProduct(ket(P, klabel(ol)), bra(P, blabel(ol)), coeff)

#################
# OuterSum Type #
#################
type OuterSum{P,T,K,B} <: AbsOuterSum{P,T,K,B}
    data::SumDict{OuterLabel{K,B},T}
    nfactors::Int
    function OuterSum(data::SumDict{OuterLabel{K,B},T})
        n = isempty(data) ? 0 : nfactors(first(keys(data)))
        return new(data, n)
    end
end

OuterSum(op::OuterSum) = op
OuterSum{P,T,K,B}(::Type{P}, data::SumDict{OuterLabel{K,B},T}) = OuterSum{P,T,K,B}(data)

function OuterSum{P}(kt::Ket{P}, br::Bra{P}, scalar=true)
    @assert matching_nfactors(kt, br)

    T = promote_type(eltype(kt), eltype(br), typeof(scalar))
    K = labeltype(kt)
    B = labeltype(br)
    result = @compat sizehint!(SumDict{OuterLabel{K,B},T}(), length(kt) * length(br))
    
    cons_outer!(result, kt, br)
    
    if scalar != 1
        scale!(result, scalar)
    end

    return OuterSum(P, result)
end

function cons_outer!(result::SumDict, kt::KetSum, br::BraSum)
    for (k,kc) in data(kt)
        for (b,bc) in data(br)
            add_to_sum!(result, OuterLabel(k, b), kc * bc')
        end
    end
    return result
end

function cons_outer!(result::SumDict, kt::SingleKet, br::BraSum)
    k = label(kt)
    kc = coeff(kt)
    for (b,bc) in data(br)
        add_to_sum!(result, OuterLabel(k, b), kc * bc')
    end
    return result
end

function cons_outer!(result::SumDict, kt::KetSum, br::SingleBra)
    b = label(br)
    bc = coeff(br)
    for (k,kc) in data(kt)
        add_to_sum!(result, OuterLabel(k, b), kc * bc)
    end
    return result
end

function cons_outer!(result::SumDict, kt::SingleKet, br::SingleBra)
    add_to_sum!(result, OuterLabel(label(kt), label(br)), coeff(kt) * coeff(br))
    return result
end

#####################
# DualOuterSum Type #
#####################
type DualOuterSum{P,T,K,B} <: AbsOuterSum{P,T,K,B}
    op::OuterSum{P,T,B,K}
end

##############
# ctranspose #
##############
eager_ctranspose{P,T,K,B}(op::OuterSum{P,T,K,B}) = OuterSum(P, load_ct!(SumDict{OuterLabel{B,K},T}(), op))

function load_ct!(result, op)
    for (o,v) in data(op)
        result[o'] = v'
    end
    return result
end

Base.ctranspose(op::OuterSum) = DualOuterSum(op)
Base.ctranspose(opc::DualOuterSum) = opc.op
Base.ctranspose(op::OuterProduct) = OuterProduct(op.br', op.kt', coeff(op)')

######################
# Accessor Functions #
######################
data(op::OuterSum) = op.data
data(opc::DualOuterSum) = data(opc.op)

coeff(op::OuterProduct) = op.coeff

Base.eltype(op::OuterOp) = eltype(typeof(op))
Base.eltype{P,S,K,B}(::Type{OuterProduct{P,S,K,B}}) = promote_type(S, eltype(K), eltype(B))
Base.eltype{P,T,K,B}(::Type{OuterSum{P,T,K,B}}) = T
Base.eltype{P,T,K,B}(::Type{DualOuterSum{P,T,K,B}}) = T

ketlabeltype(op::OuterOp) = ketlabeltype(typeof(op))
ketlabeltype{P,S,K,B}(::Type{OuterProduct{P,S,K,B}}) = labeltype(K)
ketlabeltype{P,T,K,B}(::Type{OuterSum{P,T,K,B}}) = K
ketlabeltype{P,T,K,B}(::Type{DualOuterSum{P,T,K,B}}) = K

bralabeltype(op::OuterOp) = bralabeltype(typeof(op))
bralabeltype{P,S,K,B}(::Type{OuterProduct{P,S,K,B}}) = labeltype(B)
bralabeltype{P,T,K,B}(::Type{OuterSum{P,T,K,B}}) = B
bralabeltype{P,T,K,B}(::Type{DualOuterSum{P,T,K,B}}) = B

nfactors(op::OuterSum) = op.nfactors
nfactors(op::DualOuterSum) = nfactors(op')
nfactors(op::OuterProduct) = min(nfactors(op.kt), nfactors(op.br))

########################
# Conversion/Promotion #
########################
OuterSum(opc::DualOuterSum) = eager_ctranspose(opc')
OuterSum(op::OuterProduct) = OuterSum(op.kt, op.br, coeff(op))

Base.convert(::Type{OuterSum}, op::OuterProduct) = OuterSum(op)
Base.convert(::Type{OuterSum}, opc::DualOuterSum) = OuterSum(opc)
Base.convert{P,T,K,B}(::Type{OuterSum{P,T,B,K}}, opc::DualOuterSum{P,T,K,B}) = OuterSum(opc)
Base.convert{P,T,K,B}(::Type{OuterSum{P,T,K,B}}, op::OuterProduct{P}) = convert(OuterSum{P,T,K,B}, OuterSum(op))
Base.convert{P,T,K,B}(::Type{OuterSum{P,T,K,B}}, op::OuterSum{P}) = OuterSum(P, convert(SumDict{OuterLabel{K,B},T}, data(op)))
Base.convert{P,T,K,B}(::Type{DualOuterSum{P,T,K,B}}, opc::DualOuterSum) = DualOuterSum(convert(OuterSum{P,T,B,K}, opc'))
Base.convert{O<:OuterOp}(::Type{O}, op::O) = op

Base.promote_rule{OS<:OuterSum, OP<:OuterProduct}(::Type{OS}, ::Type{OP}) = OS

function Base.promote_rule{P,T1,K1,B1,T2,K2,B2}(::Type{OuterSum{P,T1,K1,B1}}, ::Type{OuterSum{P,T2,K2,B2}})
    return OuterSum{P,promote_type(K1,K2),promote_type(B1,B2),promote_type(T1,T2)}
end

Base.promote_rule{P,T,K,B}(::Type{OuterSum{P,T,B,K}}, ::Type{DualOuterSum{P,T,K,B}}) = OuterSum{P,T,B,K}

function Base.promote_rule{P,T1,K1,B1,T2,K2,B2}(::Type{DualOuterSum{P,T1,K1,B1}}, ::Type{DualOuterSum{P,T2,K2,B2}})
    return DualOuterSum{P,promote_type(K1,K2),promote_type(B1,B2),promote_type(T1,T2)}
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
Base.copy(op::OuterProduct) = OuterProduct(copy(op.kt), copy(op.br), copy(coeff(op)))

Base.(:(==)){P}(a::OuterSum{P}, b::OuterSum{P}) = data(a) == data(b)
Base.(:(==)){P}(a::DualOuterSum{P}, b::DualOuterSum{P}) = a' == b'
Base.(:(==)){P}(a::DiracOp{P}, b::DiracOp{P}) = OuterSum(a) == OuterSum(b)

#######################
# Dict-like Functions #
#######################
Base.length(op::AbsOuterSum) = length(data(op))
Base.length(op::OuterProduct) = length(op.kt)*length(op.br)

Base.getindex(op::OuterProduct, k, b) = coeff(op) * op.kt[k] * op.br[b]
Base.getindex(op::OuterProduct, o::OuterLabel) = op[klabel(o), blabel(o)]

Base.getindex(op::OuterSum, x::OuterLabel) = getindex(data(op), x)
Base.getindex(opc::DualOuterSum, x::OuterLabel) = getindex(opc', x')'
Base.getindex(op::AbsOuterSum, k::StateLabel, b::StateLabel) = op[OuterLabel(k,b)]
Base.getindex(op::AbsOuterSum, k, b) = op[StateLabel(k),StateLabel(b)]

Base.setindex!(op::OuterSum, x, y::OuterLabel) = setindex!(data(op), x, y)
Base.setindex!(opc::DualOuterSum, x, y::OuterLabel) = setindex!(opc', x', y')
Base.setindex!(op::AbsOuterSum, x, k::StateLabel, b::StateLabel) = setindex!(op, x, OuterLabel(k,b))
Base.setindex!(op::AbsOuterSum, x, k, b) = setindex!(op, x, StateLabel(k), StateLabel(b))

Base.haskey(op::OuterProduct, k, b) = haskey(op.kt, k) && haskey(op.bra, b)
Base.haskey(op::OuterProduct, o::OuterLabel) = haskey(op, klabel(o), blabel(o))

Base.haskey(op::OuterSum, x::OuterLabel) = haskey(data(op), x)
Base.haskey(opc::DualOuterSum, x::OuterLabel) = haskey(opc', x')
Base.haskey(op::AbsOuterSum, k::StateLabel, b::StateLabel) = haskey(op, OuterLabel(k,b))
Base.haskey(op::AbsOuterSum, k, b) = haskey(op, StateLabel(k), StateLabel(b))

Base.get(op::OuterProduct, k, b, default=any_zero(eltype(op))) = haskey(op, k, b) ? op[k,b] : default
Base.get(op::OuterProduct, o::OuterLabel, default=any_zero(eltype(op))) = get(op, klabel(o), blabel(o), default)

Base.get(op::OuterSum, x::OuterLabel, default=any_zero(eltype(op))) = get(data(op), x, default)
Base.get(opc::DualOuterSum, x::OuterLabel, default=any_zero(eltype(opc))) = get(opc', x', default)
Base.get(op::AbsOuterSum, k::StateLabel, b::StateLabel, default=any_zero(eltype(op))) = get(op, OuterLabel(k,b), default)
Base.get(op::AbsOuterSum, k, b, default=any_zero(eltype(op))) = get(op, StateLabel(k), StateLabel(b), default)

#############
# Iteration #
#############
Base.start(op::AbsOuterSum) = start(data(op))

function Base.next{P}(op::OuterSum{P}, i)
    (lab,c), n = next(data(op), i)
    return tuple(OuterProduct(P, lab, c), n)
end

function Base.next(opc::DualOuterSum, i)
    op, n = next(opc', i)
    return tuple(op', n)
end

Base.done(op::AbsOuterSum, i) = done(data(op), i)

Base.collect(op::AbsOuterSum) = [i for i in op]
Base.collect(op::OuterProduct) = [OuterProduct(k,b,coeff(op)) for k in op.kt, b in op.br]

#################
# execute_inner #
#################
inner_ol_kt{P}(::Type{P}, ol::OuterLabel, v, k::StateLabel, c) = SingleKet(P, klabel(ol), v*c*inner(P, blabel(ol), k))
inner_br_ol{P}(::Type{P}, b::StateLabel, c, ol::OuterLabel, v) = ctranspose(SingleKet(P, blabel(ol), v*c*inner(P, b, klabel(ol))))

function inner_ol_ol{P}(::Type{P}, oa::OuterLabel, va, ob::OuterLabel, vb)
    return va * vb * ket(P, klabel(oa)) * bra(P, blabel(oa)) * ket(P, klabel(ob)) * bra(P, blabel(ob))
end

function execute_inner{P<:AbstractInner}(op::OuterSum{P}, kt::SingleKet{P})
    k,c = label(kt), coeff(kt)
    o0,v0 = first(data(op))
    result = inner_ol_kt(P, o0, v0, k, c)
    result -= result

    for (o,v) in data(op)
        result = unsafe_add!(result, inner_ol_kt(P, o, v, k, c))
    end

    return result
end

function execute_inner{P<:AbstractInner}(op::OuterSum{P}, kt::KetSum{P})
    k0,c0 = first(data(kt))
    o0,v0 = first(data(op))
    result = inner_ol_kt(P, o0, v0, k0, c0)
    result -= result

    for (o,v) in data(op)
        b = blabel(o)
        coeff = any_zero(eltype(result))
        for (k,c) in data(kt)
            coeff += c * v * inner(P, b, k) 
        end
        result = unsafe_add!(result, SingleKet(P, klabel(o), coeff))
    end

    return result
end

function execute_inner{P<:AbstractInner}(br::Bra{P}, op::OuterSum{P})
    b,c = label(br), coeff(br)
    o0,v0 = first(data(op))
    result = inner_br_ol(P, b, c, o0, v0)
    result -= result

    for (o,v) in data(op)
        result = unsafe_add!(result, inner_br_ol(P, b, c, o, v))
    end

    return result
end

function execute_inner{P<:AbstractInner}(br::BraSum{P}, op::OuterSum{P})
    b0,c0 = first(data(br))
    o0,v0 = first(data(op))
    result = inner_br_ol(P, b0, c0', o0, v0)
    result -= result

    for (o,v) in data(op)
        k = klabel(o)
        coeff = any_zero(eltype(result))
        for (b,c) in data(br)
            coeff += c' * v * inner(P, b, k) 
        end
        result = unsafe_add!(result, ctranspose(SingleKet(P, blabel(o), coeff')))
    end

    return result
end

function execute_inner{P<:AbstractInner}(a::OuterSum{P}, b::OuterSum{P})
    oa0,va0 = first(data(a))
    ob0,vb0 = first(data(b))
    result = inner_ol_ol(P, oa0, va0, ob0, vb0)
    result -= result

    for (oa,va) in data(a), (ob,vb) in data(b)
        result = unsafe_add!(result, 
                             OuterLabel(klabel(oa), blabel(ob)),
                             va*vb*inner(P, blabel(oa), klabel(ob)))
    end

    return result
end

function execute_inner{P<:AbstractInner}(op::OuterSum{P}, opc::DualOuterSum{P})
    o0,v0 = first(data(op))
    oc0,c0 = first(data(opc))
    result = inner_ol_ol(P, o0, v0, oc0', c0')
    result -= result

    for (o,v) in data(op), (oc,c) in data(opc)
        result = unsafe_add!(result, 
                             OuterLabel(klabel(o), klabel(oc)),
                             v*c'*inner(P, blabel(o), blabel(oc)))
    end
    
    return result
end

function execute_inner{P<:AbstractInner}(opc::DualOuterSum{P}, op::OuterSum{P})
    o0,v0 = first(data(op))
    oc0,c0 = first(data(opc))
    result = inner_ol_ol(P, oc0', c0', o0, v0)
    result -= result

    for (oc,c) in data(opc), (o,v) in data(op)
        result = unsafe_add!(result, 
                             OuterLabel(blabel(oc), blabel(o)),
                             v*c'*inner(P, klabel(oc), klabel(o)))
    end

    return result
end

#####################################################
# execute_inner - optimized for ProvidedInner types #
#####################################################

# we generate these redundant definitions to avoid 
# ambiguity with previous inner() definitions

for T in (:BraSum, :SingleBra)
    @eval begin
        function execute_inner{P<:ProvidedInner}(br::($T){P}, op::OuterSum{P})
            result = SumDict{StateLabel{bralabeltype(op)}, inner_rettype(br, op)}()
            return ctranspose(KetSum(P, inner_load!(result, br, op)))
        end
    end
end

for T in (:KetSum, :SingleKet)
    @eval begin
        function execute_inner{P<:ProvidedInner}(op::OuterSum{P}, kt::($T){P})
            result = SumDict{StateLabel{ketlabeltype(op)}, inner_rettype(op, kt)}()
            return KetSum(P, inner_load!(result, op, kt))
        end
    end
end

function execute_inner{P<:ProvidedInner,T1,K1,B1,T2,K2,B2}(a::OuterSum{P,T1,K1,B1}, b::OuterSum{P,T2,K2,B2})
    result = SumDict{OuterLabel{K1,B2}, inner_rettype(a, b)}()
    return OuterSum(P, inner_load!(result, a, b))
end

function execute_inner{P<:ProvidedInner,T1,K1,B1,T2,K2,B2}(op::OuterSum{P,T1,K1,B1}, opc::DualOuterSum{P,T2,K2,B2})
    result = SumDict{OuterLabel{K1,B2}, inner_rettype(op, opc)}()
    return OuterSum(P, inner_load!(result, op, opc))
end

function execute_inner{P<:ProvidedInner,T1,K1,B1,T2,K2,B2}(opc::DualOuterSum{P,T1,K1,B1}, op::OuterSum{P,T2,K2,B2})
    result = SumDict{OuterLabel{K1,B2}, inner_rettype(opc, op)}()
    return OuterSum(P, inner_load!(result, opc, op))
end

function inner_load!{P<:ProvidedInner}(result::SumDict, br::SingleBra{P}, op::OuterSum{P})
    c = coeff(br)
    b = label(br)
    for (o,v) in data(op)
        add_to_sum!(result, blabel(o), ctranspose(v*c*inner(P, b, klabel(o))))
    end
    return result
end

function inner_load!{P<:ProvidedInner}(result::SumDict, br::BraSum{P}, op::OuterSum{P})
    for (o,v) in data(op)
        k = klabel(o)
        coeff = any_zero(eltype(result))
        for (b,c) in data(br)
            coeff += c' * v * inner(P, b, k) 
        end
        add_to_sum!(result, blabel(o), ctranspose(coeff))
    end
    return result
end

function inner_load!{P<:ProvidedInner}(result::SumDict, op::OuterSum{P}, kt::SingleKet{P})
    c = coeff(kt)
    k = label(kt)
    for (o,v) in data(op)
        add_to_sum!(result, klabel(o), v*c*inner(P, blabel(o), k))
    end
    return result
end

function inner_load!{P<:ProvidedInner}(result::SumDict, op::OuterSum{P}, kt::KetSum{P})
    for (o,v) in data(op)
        b = blabel(o)
        coeff = any_zero(eltype(result))
        for (k,c) in data(kt)
            coeff += c * v * inner(P, b, k) 
        end
        add_to_sum!(result, klabel(o), coeff)
    end
    return result
end

function inner_load!{P<:ProvidedInner}(result::SumDict, a::OuterSum{P}, b::OuterSum{P})
    for (o1,v) in data(a), (o2,c) in data(b)
        add_to_sum!(result, 
                    OuterLabel(klabel(o1), blabel(o2)),
                    v*c*inner(P, blabel(o1), klabel(o2)))
    end
    return result
end

function inner_load!{P<:ProvidedInner}(result::SumDict, op::OuterSum{P}, opc::DualOuterSum{P})
    for (o,v) in data(op), (oc,c) in data(opc)
        add_to_sum!(result, 
                    OuterLabel(klabel(o), klabel(oc)),
                    v*c'*inner(P, blabel(o), blabel(oc)))
    end
    return result
end

function inner_load!{P<:ProvidedInner}(result::SumDict, opc::DualOuterSum{P}, op::OuterSum{P})
    for (oc,v) in data(opc), (o,c) in data(op)
        add_to_sum!(result,
                    OuterLabel(blabel(oc), blabel(o)),
                    v'*c*inner(P, klabel(oc), klabel(o)))
    end
    return result
end

##########################################
# execute_inner - math-based definitions #
##########################################
execute_inner(br::Bra, opc::DualOuterSum) = inner(opc', br')'
execute_inner(opc::DualOuterSum, kt::Ket) = inner(kt', opc')'
execute_inner(a::DualOuterSum, b::DualOuterSum) = inner(b', a')'

execute_inner(br::Bra, op::OuterProduct) = coeff(op) * inner(br, op.kt) * op.br
execute_inner(op::OuterProduct, kt::Ket) = coeff(op) * op.kt * inner(op.br, kt)
execute_inner(a::OuterProduct, b::OuterProduct) = OuterProduct(a.kt, b.br, coeff(a)*coeff(b)*inner(a.br,b.kt))
execute_inner(a::OuterProduct, b::AbsOuterSum) = coeff(a) * a.kt * inner(a.br, b)
execute_inner(a::AbsOuterSum, b::OuterProduct) = inner(a, b.kt) * b.br * coeff(b)

Base.(:*)(br::Bra, op::DiracOp) = inner(br, op)
Base.(:*)(op::DiracOp, kt::Ket) = inner(op, kt)
Base.(:*)(a::DiracOp, b::DiracOp) = inner(a, b)

inner_eval(f, op::DiracOp) = mapcoeffs(x->inner_eval(f,x), op)

##################
# execute_act_on #
##################
execute_act_on(op::DiracOp, br::Bra, i) = act_on(op', br', i)'

function execute_act_on{P}(op::DiracOp, kt::SingleKet{P}, i)
    k = label(kt)
    v = coeff(kt)
    return (v * ket(P, k[1:i-1]...)) * (op * ket(P, k[i])) * ket(P, k[i+1:end]...)
end

execute_act_on(op::DiracOp, kt::KetSum, i) = sum(map(k->act_on(op, k, i), collect(kt)))

execute_act_on{P}(op::OuterProduct{P}, kt::KetSum{P}, i) = act_on(OuterSum(op), kt, i)

function execute_act_on{P}(op::OuterSum{P}, kt::KetSum{P}, i)
    k0,v0 = first(data(kt))
    o0,c0 = first(data(op))
    result = act_on(OuterProduct(P,o0,c0), SingleKet(P,k0,v0), i)
    result -= result

    for (o,c) in data(op), (k,v) in data(kt)
        result = unsafe_add!(result, act_on(OuterProduct(P,o,c), SingleKet(P,k,v), i)) 
    end

    return result
end

function execute_act_on{P}(opc::DualOuterSum{P}, kt::KetSum{P}, i)
    k0,v0 = first(data(kt))
    o0,c0 = first(data(opc))
    result = act_on(OuterProduct(P,o0',c0'), SingleKet(P,k0,v0), i)
    result -= result

    for (o,c) in data(opc), (k,v) in data(kt)
        result = unsafe_add!(result, act_on(OuterProduct(P,o',c'), SingleKet(P,k,v), i)) 
    end

    return result
end

######################################################
# execute_act_on - optimized for ProvidedInner types #
######################################################

# Generate explicit (but redudant) code to resolve method ambiguity
for T in (:KetSum, :SingleKet), OT in (:OuterSum, :DualOuterSum)
    @eval begin
        execute_act_on{P<:ProvidedInner}(op::($OT){P}, kt::($T){P}, i) = KetSum(P, act_on_dict!(act_result(op, kt), op, kt, i))
    end
end

function act_result{P<:ProvidedInner}(op::AbsOuterSum{P}, kt::Ket{P})
    T = promote_type(ketlabeltype(op), labeltype(kt))
    return SumDict{StateLabel{T},inner_rettype(op, kt)}()
end

function act_on_dict!{P<:ProvidedInner}(result, op::OuterSum{P}, kt::SingleKet{P}, i)
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

function act_on_dict!{P<:ProvidedInner}(result, op::OuterSum{P}, kt::KetSum{P}, i)
    for (o,c) in data(op), (k,v) in data(kt)
        add_to_sum!(result, 
                    setindex(k, klabel(o)[1], i),
                    c*v*P(blabel(o)[1], k[i]))
    end
    return result
end

function act_on_dict!{P<:ProvidedInner}(result, opc::DualOuterSum{P}, kt::SingleKet{P}, i)
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

function act_on_dict!{P<:ProvidedInner}(result, opc::DualOuterSum{P}, kt::KetSum{P}, i)
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
tensor(a::OuterProduct, b::OuterProduct) = OuterProduct(tensor(a.kt,b.kt), tensor(a.br, b.br), coeff(a) * coeff(b))
tensor(a::DiracOp, b::DiracOp) = tensor(OuterSum(a), OuterSum(b))

tensor(kt::Ket, br::Bra) = OuterProduct(kt, br)
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
Base.scale(op::OuterProduct, c::Number) = OuterProduct(copy(op.kt), copy(op.br), coeff(op) * c)
Base.scale(c::Number, op::OuterOp) = scale(op, c)

Base.(:*)(c::Number, op::DiracOp) = scale(c, op)
Base.(:*)(op::DiracOp, c::Number) = scale(op, c)
Base.(:/)(op::DiracOp, c::Number) = scale(op, 1/c)

###########
# + and - #
###########
Base.(:+){P}(a::OuterSum{P}, b::OuterSum{P}) = (@assert matching_nfactors(a, b); OuterSum(P, data(a) + data(b)))
Base.(:+){P}(a::DualOuterSum{P}, b::DualOuterSum{P}) = ctranspose(a' + b')
Base.(:+){P}(a::DiracOp{P}, b::DiracOp{P}) = OuterSum(a) + OuterSum(b)

Base.(:-){P}(a::OuterSum{P}, b::OuterSum{P}) = (@assert matching_nfactors(a, b); OuterSum(P, data(a) - data(b)))
Base.(:-){P}(a::DualOuterSum{P}, b::DualOuterSum{P}) = ctranspose(a' - b')
Base.(:-){P}(a::DiracOp{P}, b::DiracOp{P}) = OuterSum(a) - OuterSum(b)

Base.(:-)(op::OuterOp) = scale(-1, op)

# The below methods are unsafe because:
# 1. The first argument is potentially not mutated
# 2. No nfactors check is performed
unsafe_add!(op::OuterSum, ol::OuterLabel, v) = (add!(data(op), SumTerm(ol,v)); op)
unsafe_add!(opc::DualOuterSum, ol::OuterLabel, v) = (unsafe_add!(opc', ol', v'); opc)
unsafe_add!{P}(a::OuterSum{P}, b::OuterSum{P}) = (add!(data(a), data(b)); a)
unsafe_add!{P}(a::DualOuterSum{P}, b::DualOuterSum{P}) = (add!(data(a), data(b)); a)
unsafe_add!{P}(a::OuterSum{P}, b::OuterOp{P}) = unsafe_add!(a, OuterSum(b))
unsafe_add!{P}(a::DualOuterSum{P}, b::OuterOp{P}) = unsafe_add!(a, OuterSum(b)')

unsafe_sub!(op::OuterSum, ol::OuterLabel, v) = (sub!(data(op), SumTerm(ol,v)); op)
unsafe_sub!(opc::DualOuterSum, ol::OuterLabel, v) = (unsafe_sub!(opc', ol', v'); opc)
unsafe_sub!{P}(a::OuterSum{P}, b::OuterSum{P}) = (sub!(data(a), data(b)); a)
unsafe_sub!{P}(a::DualOuterSum{P}, b::DualOuterSum{P}) = (sub!(data(a), data(b)); a)
unsafe_sub!{P}(a::OuterSum{P}, b::OuterOp{P}) = unsafe_sub!(a, OuterSum(b))
unsafe_sub!{P}(a::DualOuterSum{P}, b::OuterOp{P}) = unsafe_sub!(a, OuterSum(b)')

add!{P}(a::AbsOuterSum{P}, b::OuterOp{P}) = (@assert matching_nfactors(a, b); unsafe_add!(a, b))
sub!{P}(a::AbsOuterSum{P}, b::OuterOp{P}) = (@assert matching_nfactors(a, b); unsafe_sub!(a, b))

#################
# Normalization #
#################
Base.norm(op::OuterSum) = sqrt(sum(abs2, values(data(op))))
Base.norm(opc::DualOuterSum) = norm(opc')
function Base.norm(op::OuterProduct)
    result = any_zero(eltype(op))
    for v in values(data(op.kt)), c in values(data(op.br))
        result += abs2(coeff(op) * v * c')
    end
    return sqrt(result)
end

normalize(op::DiracOp) = scale(1/norm(op), op)
normalize!(op::DiracOp) = scale!(1/norm(op), op)

#########
# Trace #
#########
function Base.trace(op::OuterSum)
    result = any_zero(eltype(op))
    for (o,v) in data(op)
        if klabel(o)==blabel(o)
            result += v
        end
    end
    return result
end

Base.trace(opc::DualOuterSum) = trace(opc')'

function Base.trace(op::OuterProduct)
    result = any_zero(eltype(op))
    for (k,v) in data(op.kt), (b,c) in data(op.br)
        if b == k
            result += v * c'
        end
    end
    return coeff(op) * result
end

#################
# Partial trace #
#################
function ptrace_result{P}(op::DiracOp{P})
    K = ketlabeltype(op)
    B = bralabeltype(op)
    return SumDict{OuterLabel{K,B}, eltype(op)}()
end

function ptrace{P}(op::DiracOp{P}, over)
    @assert nfactors(op) > 1 "nfactors(op) must be greater than 1 when using ptrace(op,i); perhaps you meant to use trace(op)?"
    @assert 1 <= over <= nfactors(op) "i is constrained to 1 <= i <= nfactors(op) when using ptrace(op,i)"
    return OuterSum(P, ptrace_dict!(ptrace_result(op), op, over))
end

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
            add_to_sum!(result, except(OuterLabel(k, b), over), op[k,b])
        end
    end
    return result
end

#####################
# Partial Transpose #
#####################
ptranspose{P}(op::OuterOp{P}, over) = OuterSum(P, ptrans_dict(op, over))

ptrans_dict(op::OuterSum, over) = SumDict([ptranspose(o, over) => v for (o,v) in data(op)])
ptrans_dict(op::DualOuterSum, over) = SumDict([ptranspose(o, over) => v for (o,v) in data(op)])
ptrans_dict(op::OuterProduct, over) = SumDict([ptranspose(OuterLabel(k, b), over) => op[k,b] for k in keys(data(op.kt)), b in keys(data(op.br))])

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

represent(op::DiracOp, kets) = [kj' * (op * ki) for ki in kets, kj in kets]

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
