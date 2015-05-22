################
# OuterProduct #
################
immutable OuterProduct{P,N,S,K<:Ket,B<:Bra} <: DiracOp{P,N}
    scalar::S
    kt::K
    br::B
    OuterProduct(scalar::S, kt::Ket{P,N}, br::Bra{P,N}) = new(scalar, kt, br)
end

OuterProduct{P,N,S}(scalar::S, kt::Ket{P,N}, br::Bra{P,N}) = OuterProduct{P,N,S,typeof(kt),typeof(br)}(scalar, kt, br)

Base.copy(op::OuterProduct) = OuterProduct(copy(op.scalar), copy(op.kt), copy(op.br))
Base.eltype(op::OuterProduct) = promote_type(typeof(op.scalar), eltype(op.kt), eltype(op.br))
ketlabeltype(op::OuterProduct) = labeltype(op.kt)
bralabeltype(op::OuterProduct) = labeltype(op.br)

OpSum(op::OuterProduct) = OpSum(op.kt, op.br, op.scalar)

Base.convert{P,N,K,B,T}(::Type{OpSum{P,N,K,B,T}}, op::OuterProduct) = convert(OpSum{P,N,K,B,T}, OpSum(op))
Base.convert(::Type{OpSum}, op::OuterProduct) = OpSum(op)
Base.convert{P,N,S,K,B}(::Type{OuterProduct{P,N,S,K,B}}, op::OuterProduct{P,N,S,K,B}) = op

Base.promote_rule{OS<:OpSum, OP<:OuterProduct}(::Type{OS}, ::Type{OP}) = OS

Base.(:(==)){P,N}(a::OuterProduct{P,N}, b::OuterProduct{P,N}) = OpSum(a) == OpSum(b) 
Base.hash(op::OuterProduct) = hash(OpSum(op))
Base.hash(op::OuterProduct, h::Uint64) = hash(hash(op), h)

Base.length(op::OuterProduct) = length(op.kt)*length(op.br)

Base.getindex(op::OuterProduct, k, b) = op.scalar * op.kt[k] * op.br[b]
Base.getindex(op::OuterProduct, o::OpLabel) = op[klabel(o), blabel(o)]

Base.haskey(op::OuterProduct, k, b) = hasket(op, k) && hasbra(op, b)
Base.haskey(op::OuterProduct, o::OpLabel) = haskey(op, klabel(o), blabel(o))

Base.get(op::OuterProduct, k, b, default=predict_zero(eltype(op))) = haskey(op, k, b) ? op[k,b] : default
Base.get(op::OuterProduct, o::OpLabel, default=predict_zero(eltype(op))) = get(op, klabel(o), blabel(o), default)

function Base.collect{P,N}(op::OuterProduct{P,N})
    T = @compat Tuple{OpLabel{N,ketlabeltype(op),bralabeltype(op)}, eltype(op)}
    return collect_pairs!(Array(T, length(op)), op)
end

function collect_pairs!(result, op::OuterProduct)
    i = 1
    for (k,kc) in data(op.kt)
        for (b,bc) in data(op.br)
            result[i] = (OpLabel(k, b), op.scalar * kc * bc')
            i += 1
        end
    end
    return result
end

##############
# ctranspose #
##############
Base.ctranspose(op::OuterProduct) = OuterProduct(op.scalar', op.br', op.kt')

#########
# inner #
#########
inner(br::Bra, op::OuterProduct) = op.scalar * inner(br, op.kt) * op.br
inner(op::OuterProduct, kt::Ket) = op.scalar * op.kt * inner(op.br, kt)
inner(a::OuterProduct, b::OuterProduct) = OuterProduct(a.scalar * b.scalar * inner(a.br,b.kt), a.kt, b.br)
inner(a::OuterProduct, b::AbsOpSum) = a.scalar * a.kt * inner(a.br, b)
inner(a::AbsOpSum, b::OuterProduct) = inner(a, b.kt) * b.br * b.scalar

##########
# act_on #
##########
act_on(op::OuterProduct, state::DiracState, i) = act_on(OpSum(op), state, i)

##########
# tensor #
##########
tensor(kt::Ket, br::Bra) = OuterProduct(1, kt, br)
tensor(br::Bra, kt::Ket) = tensor(kt, br)
tensor(a::OuterProduct, b::OuterProduct) = OuterProduct(a.scalar * b.scalar, tensor(a.kt,b.kt), tensor(a.br, b.br))

Base.(:*)(kt::Ket, br::Bra) = tensor(kt,br)

###########
# Scaling #
###########
Base.scale(c::Number, op::OuterProduct) = OuterProduct(c * op.scalar, copy(op.kt), copy(op.br))
Base.scale(op::OuterProduct, c::Number) = OuterProduct(op.scalar * c, copy(op.kt), copy(op.br))

###########
# + and - #
###########
Base.(:-)(op::OuterProduct) = scale(-1, op)
Base.(:+)(a::OuterProduct, b::OuterProduct) = OpSum(a) + OpSum(b)

#################
# Normalization #
#################
function Base.norm(op::OuterProduct)
    result = predict_zero(eltype(op))
    for v in values(data(op.kt)), c in values(data(op.br))
        result += abs2(op.scalar * v * c')
    end
    return sqrt(result)
end

#########
# Trace #
#########
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
# Partial Trace #
#################
function ptrace_dict!(result, op::OuterProduct, over)
    for k in keys(data(op.kt)), b in keys(data(op.br))
        if k[over] == b[over]
            add_to_sum!(result, except(OpLabel(k, b), over), op[k,b])
        end
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
nfactors{P,N}(::OuterProduct{P,N}) = N
xsubspace(op::OuterProduct,x) = xsubspace(OpSum(op), x)
filternz(op::OuterProduct) = filternz(OpSum(op))
switch(op::OuterProduct, i, j) = switch(OpSum(op), i, j)
permute(op::OuterProduct, perm::Vector) = permute(OpSum(op), perm)

######################
# Printing Functions #
######################
labelrepr(op::OuterProduct, k, b, pad) = "$pad$(op[k,b]) $(ktstr(k))$(brstr(b))"

Base.summary{P,N}(op::OuterProduct{P,N}) = "OuterProduct{$P,$N,$(eltype(op))} with $(length(op)) operator(s)"

function Base.show(io::IO, op::OuterProduct)
    print(io, summary(op)*":")
    pad = "  "
    maxlen = 10
    i = 1
    for k in keys(data(op.kt)), b in keys(data(op.br))
        if i > maxlen
            break
        else
            i += 1
            println(io)
            print(io, labelrepr(op, k, b, pad))
        end
    end
    if length(op) > maxlen
        println(io)
        print(io, "$pad$vdots")
    end
end

export ptrace,
    purity,
    filternz,
    xsubspace,
    nfactors
