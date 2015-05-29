abstract OuterOp{P,N} <: DiracOp{P,N}
abstract AbsOuterSum{P,N,T,K,B} <: OuterOp{P,N}

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
OuterProduct{P}(::Type{P}, ol::OuterLabel, scalar) = OuterProduct(scalar, ket(P, klabel(ol)), bra(P, blabel(ol)))

#################
# OuterSum Type #
#################
type OuterSum{P,N,T,K,B} <: AbsOuterSum{P,N,T,K,B}
    data::SumDict{OuterLabel{N,K,B},T}
end

OuterSum(op::OuterSum) = op
OuterSum{P,N,T,K,B}(::Type{P}, data::SumDict{OuterLabel{N,K,B},T}) = OuterSum{P,N,T,K,B}(data)

function OuterSum{P,N}(kt::Ket{P,N}, br::Bra{P,N}, scalar = 1)
    T = promote_type(eltype(kt), eltype(br), typeof(scalar))
    K = labeltype(kt)
    B = labeltype(br)
    result = @compat sizehint!(SumDict{OuterLabel{N,K,B},T}(), length(kt) * length(br))
    
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
type DualOuterSum{P,N,T,K,B} <: AbsOuterSum{P,N,T,K,B}
    op::OuterSum{P,N,T,B,K}
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
data(opc::DualOuterSum) = data(opc.op)

Base.eltype(op::OuterOp) = eltype(typeof(op))
Base.eltype{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}) = promote_type(S, eltype(K), eltype(B))
Base.eltype{P,N,T,K,B}(::Type{OuterSum{P,N,T,K,B}}) = T
Base.eltype{P,N,T,K,B}(::Type{DualOuterSum{P,N,T,K,B}}) = T

ketlabeltype(op::OuterOp) = ketlabeltype(typeof(op))
ketlabeltype{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}) = labeltype(K)
ketlabeltype{P,N,T,K,B}(::Type{OuterSum{P,N,T,K,B}}) = K
ketlabeltype{P,N,T,K,B}(::Type{DualOuterSum{P,N,T,K,B}}) = K

bralabeltype(op::OuterOp) = bralabeltype(typeof(op))
bralabeltype{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}) = labeltype(B)
bralabeltype{P,N,T,K,B}(::Type{OuterSum{P,N,T,K,B}}) = B
bralabeltype{P,N,T,K,B}(::Type{DualOuterSum{P,N,T,K,B}}) = B

nfactors(op::OuterOp) = nfactors(typeof(op))
nfactors{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}) = N
nfactors{P,N,T,K,B}(::Type{OuterSum{P,N,T,K,B}}) = N
nfactors{P,N,T,K,B}(::Type{DualOuterSum{P,N,T,K,B}}) = N

########################
# Conversion/Promotion #
########################
OuterSum(opc::DualOuterSum) = eager_ctranspose(opc')
OuterSum(op::OuterProduct) = OuterSum(op.kt, op.br, op.scalar)

Base.convert{P,N,T,K,B}(::Type{OuterSum{P,N,T,K,B}}, op::OuterProduct) = convert(OuterSum{P,N,T,K,B}, OuterSum(op))
Base.convert(::Type{OuterSum}, op::OuterProduct) = OuterSum(op)
Base.convert{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}, op::OuterProduct{P,N,S,K,B}) = op

Base.convert{P,N,T,K,B}(::Type{OuterSum{P,N,T,K,B}}, op::OuterSum{P}) = OuterSum(P, convert(SumDict{OuterLabel{N,K,B},T}, data(op)))
Base.convert{P,N,T,K,B}(::Type{OuterSum{P,N,T,K,B}}, op::OuterSum{P,N,T,K,B}) = op

Base.convert{P,N,T,K,B}(::Type{OuterSum{P,N,T,B,K}}, opc::DualOuterSum{P,N,T,K,B}) = OuterSum(opc)
Base.convert{P,N,T,K,B}(::Type{OuterSum}, opc::DualOuterSum{P,N,T,K,B}) = OuterSum(opc)
Base.convert{P,N,T,K,B}(::Type{DualOuterSum{P,N,T,K,B}}, opc::DualOuterSum) = DualOuterSum(convert(OuterSum{P,N,T,B,K}, opc'))
Base.convert{P,N,T,K,B}(::Type{DualOuterSum{P,N,T,K,B}}, opc::DualOuterSum{P,N,T,K,B}) = opc

Base.promote_rule{OS<:OuterSum, OP<:OuterProduct}(::Type{OS}, ::Type{OP}) = OS

function Base.promote_rule{P,N,T1,K1,B1,T2,K2,B2}(::Type{OuterSum{P,N,T1,K1,B1}}, ::Type{OuterSum{P,N,T2,K2,B2}})
    return OuterSum{P,N,promote_type(K1,K2),promote_type(B1,B2),promote_type(T1,T2)}
end

Base.promote_rule{P,N,T,K,B}(::Type{OuterSum{P,N,T,B,K}}, ::Type{DualOuterSum{P,N,T,K,B}}) = OuterSum{P,N,T,B,K}

function Base.promote_rule{P,N,T1,K1,B1,T2,K2,B2}(::Type{DualOuterSum{P,N,T1,K1,B1}}, ::Type{DualOuterSum{P,N,T2,K2,B2}})
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

Base.get(op::OuterProduct, k, b, default=predict_zero(eltype(op))) = haskey(op, k, b) ? op[k,b] : default
Base.get(op::OuterProduct, o::OuterLabel, default=predict_zero(eltype(op))) = get(op, klabel(o), blabel(o), default)

Base.get(op::OuterSum, x::OuterLabel, default=predict_zero(eltype(op))) = get(data(op), x, default)
Base.get(opc::DualOuterSum, x::OuterLabel, default=predict_zero(eltype(opc))) = get(opc', x', default)
Base.get(op::AbsOuterSum, k::StateLabel, b::StateLabel, default=predict_zero(eltype(op))) = get(op, OuterLabel(k,b), default)
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
    T = @compat Tuple{OuterLabel{nfactors(op),ketlabeltype(op),bralabeltype(op)}, eltype(op)}
    return collect_pairs!(Array(T, length(op)), op)
end

function collect_pairs!(result, op::OuterProduct)
    i = 1
    for (k,kc) in data(op.kt), (b,bc) in data(op.br)
        result[i] = tuple(OuterLabel(k, b), op.scalar * kc * bc')
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
inner_ol_kt{P}(::Type{P}, ol::OuterLabel, v, k::StateLabel, c) = SingleKet(P, klabel(ol), v*c*inner(P, blabel(ol), k))

inner_br_ol{P}(::Type{P}, b::StateLabel, c, ol::OuterLabel, v) = ctranspose(SingleKet(P, blabel(ol), v*c*inner(P, b, klabel(ol))))

function inner_ol_ol{P}(::Type{P}, oa::OuterLabel, va, ob::OuterLabel, vb)
    return va * vb * ket(P, klabel(oa)) * bra(P, blabel(oa)) * ket(P, klabel(ob)) * bra(P, blabel(ob))
end

function inner{P<:AbstractInner,N}(op::OuterSum{P,N}, kt::SingleKet{P,N})
    k,c = label(kt), coeff(kt)
    o0,v0 = first(data(op))
    result = inner_ol_kt(P, o0, v0, k, c)
    result -= result

    for (o,v) in data(op)
        result = unsafe_add!(result, inner_ol_kt(P, o, v, k, c))
    end

    return result
end

function inner{P<:AbstractInner,N}(op::OuterSum{P,N}, kt::KetSum{P,N})
    k0,c0 = first(data(kt))
    o0,v0 = first(data(op))
    result = inner_ol_kt(P, o0, v0, k0, c0)
    result -= result

    for (o,v) in data(op)
        b = blabel(o)
        coeff = predict_zero(eltype(result))
        for (k,c) in data(kt)
            coeff += c * v * inner(P, b, k) 
        end
        result = unsafe_add!(result, SingleKet(P, klabel(o), coeff))
    end

    return result
end

function inner{P<:AbstractInner,N}(br::Bra{P,N}, op::OuterSum{P,N})
    b,c = label(br), coeff(br)
    o0,v0 = first(data(op))
    result = inner_br_ol(P, b, c, o0, v0)
    result -= result

    for (o,v) in data(op)
        result = unsafe_add!(result, inner_br_ol(P, b, c, o, v))
    end

    return result
end

function inner{P<:AbstractInner,N}(br::BraSum{P,N}, op::OuterSum{P,N})
    b0,c0 = first(data(br))
    o0,v0 = first(data(op))
    result = inner_br_ol(P, b0, c0', o0, v0)
    result -= result

    for (o,v) in data(op)
        k = klabel(o)
        coeff = predict_zero(eltype(result))
        for (b,c) in data(br)
            coeff += c' * v * inner(P, b, k) 
        end
        result = unsafe_add!(result, ctranspose(SingleKet(P, blabel(o), coeff')))
    end

    return result
end

function inner{P<:AbstractInner,N}(a::OuterSum{P,N}, b::OuterSum{P,N})
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

function inner{P<:AbstractInner,N}(op::OuterSum{P,N}, opc::DualOuterSum{P,N})
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

function inner{P<:AbstractInner,N}(opc::DualOuterSum{P,N}, op::OuterSum{P,N})
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

#############################################
# inner - optimized for ProvidedInner types #
#############################################

# we generate these redundant definitions to avoid 
# ambiguity with previous inner() definitions

for T in (:BraSum, :SingleBra)
    @eval begin
        function inner{P<:ProvidedInner,N}(br::($T){P,N}, op::OuterSum{P,N})
            result = SumDict{StateLabel{N, bralabeltype(op)}, inner_rettype(br, op)}()
            return ctranspose(KetSum(P, inner_load!(result, br, op)))
        end
    end
end

for T in (:KetSum, :SingleKet)
    @eval begin
        function inner{P<:ProvidedInner,N}(op::OuterSum{P,N}, kt::($T){P,N})
            result = SumDict{StateLabel{N, ketlabeltype(op)}, inner_rettype(op, kt)}()
            return KetSum(P, inner_load!(result, op, kt))
        end
    end
end

function inner{P<:ProvidedInner,N,T1,K1,B1,T2,K2,B2}(a::OuterSum{P,N,T1,K1,B1}, b::OuterSum{P,N,T2,K2,B2})
    result = SumDict{OuterLabel{N,K1,B2}, inner_rettype(a, b)}()
    return OuterSum(P, inner_load!(result, a, b))
end

function inner{P<:ProvidedInner,N,T1,K1,B1,T2,K2,B2}(op::OuterSum{P,N,T1,K1,B1}, opc::DualOuterSum{P,N,T2,K2,B2})
    result = SumDict{OuterLabel{N,K1,B2}, inner_rettype(op, opc)}()
    return OuterSum(P, inner_load!(result, op, opc))
end

function inner{P<:ProvidedInner,N,T1,K1,B1,T2,K2,B2}(opc::DualOuterSum{P,N,T1,K1,B1}, op::OuterSum{P,N,T2,K2,B2})
    result = SumDict{OuterLabel{N,K1,B2}, inner_rettype(opc, op)}()
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
        coeff = predict_zero(eltype(result))
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
        coeff = predict_zero(eltype(result))
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

##############################
# inner - simple definitions #
##############################
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
act_on(op::DiracOp, br::Bra, i) = act_on(op', br', i)'
act_on{P}(op::DiracOp{P,1}, kt::Ket{P,1}, i) = i==1 ? inner(op, kt) : throw(BoundsError())

function act_on{P,N}(op::DiracOp{P,1}, kt::SingleKet{P,N}, i)
    k = label(kt)
    v = coeff(kt)
    return (v * ket(P, k[1:i-1]...)) * (op * ket(P, k[i])) * ket(P, k[i+1:end]...)
end

act_on{P,N}(op::OuterProduct{P,1}, kt::KetSum{P,N}, i) = act_on(OuterSum(op), kt, i)

function act_on{P,N}(op::OuterSum{P,1}, kt::KetSum{P,N}, i)
    k0,v0 = first(data(kt))
    o0,c0 = first(data(op))
    result = act_on(OuterProduct(P,o0,c0), SingleKet(P,k0,v0), i)
    result -= result

    for (o,c) in data(op), (k,v) in data(kt)
        result = unsafe_add!(result, act_on(OuterProduct(P,o,c), SingleKet(P,k,v), i)) 
    end

    return result
end

function act_on{P,N}(opc::DualOuterSum{P,1}, kt::KetSum{P,N}, i)
    k0,v0 = first(data(kt))
    o0,c0 = first(data(opc))
    result = act_on(OuterProduct(P,o0',c0'), SingleKet(P,k0,v0), i)
    result -= result

    for (o,c) in data(opc), (k,v) in data(kt)
        result = unsafe_add!(result, act_on(OuterProduct(P,o',c'), SingleKet(P,k,v), i)) 
    end

    return result
end

##############################################
# act_on - optimized for ProvidedInner types #
##############################################

# Generate explicit (but redudant) code to resolve method ambiguity
for T in (:KetSum, :SingleKet), OT in (:OuterSum, :DualOuterSum)
    @eval begin
        act_on{P<:ProvidedInner,N}(op::($OT){P,1}, kt::($T){P,N}, i) = KetSum(P, act_on_dict!(act_result(op, kt), op, kt, i))
    end
end

function act_result{P<:ProvidedInner,N}(op::AbsOuterSum{P,1}, kt::Ket{P,N})
    T = promote_type(ketlabeltype(op), labeltype(kt))
    return SumDict{StateLabel{N,T},inner_rettype(op, kt)}()
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

unsafe_add!{P}(op::OuterSum{P}, ol::OuterLabel, v) = OuterSum(P, data(op) + SumTerm(ol,v))
unsafe_add!(opc::DualOuterSum, ol::OuterLabel, v) = unsafe_add!(opc', ol', v')'

unsafe_add!{P,N,T,K,B}(op::OuterSum{P,N,T,K,B}, ol::OuterLabel{N,K,B}, v::T) = (add!(data(op), SumTerm(ol,v)); return op)
unsafe_add!{P,N,T,K,B}(opc::DualOuterSum{P,N,T,K,B}, ol::OuterLabel{N,K,B}, v::T) = (add!(data(opc), SumTerm(ol',v')); return opc)

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
    return SumDict{OuterLabel{N-1,K,B}, eltype(op)}()
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
            add_to_sum!(result, except(OuterLabel(k, b), over), op[k,b])
        end
    end
    return result
end

#####################
# Partial Transpose #
#####################
function ptrans_result{P,N}(op::DiracOp{P,N})
    T = label_promote(ketlabeltype(op), bralabeltype(op))
    return SumDict{OuterLabel{N,T,T}, eltype(op)}()
end

function ptrans_result{P}(op::DiracOp{P,1})
    K = bralabeltype(op)
    B = ketlabeltype(op)
    return SumDict{OuterLabel{N,K,B}, eltype(op)}()
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
        add_to_sum!(result, ptranspose(OuterLabel(k, b), over), op[k,b])
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
