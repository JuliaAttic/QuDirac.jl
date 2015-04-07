################
# OuterProduct #
################
type OuterProduct{P,N} <: DiracOp{P,N}
    scalar::Number
    kt::Ket{P,N}
    br::Bra{P,N}
end


fact(op::OuterProduct) = fact(op.kt)
ptype(op::OuterProduct) = ptype(op.kt)

Base.copy(op::OuterProduct) = OuterProduct(copy(op.scalar), copy(op.kt), copy(op.br))

Base.convert(::Type{GenericOp}, op::OuterProduct) = scale!(op.scalar, GenericOp(op.kt, op.br))
Base.convert{P}(::Type{GenericOp{P}}, op::OuterProduct{P}) = convert(GenericOp, op)
Base.convert{P,N}(::Type{GenericOp{P,N}}, op::OuterProduct{P,N}) = convert(GenericOp, op)

Base.promote_rule(::Type{GenericOp}, ::Type{OuterProduct}) = GenericOp
Base.promote_rule{P}(::Type{GenericOp{P}}, ::Type{OuterProduct{P}}) = GenericOp{P}
Base.promote_rule{P,N}(::Type{GenericOp{P,N}}, ::Type{OuterProduct{P,N}}) = GenericOp{P,N}

#######################
# Dict-Like Functions #
#######################
Base.(:(==)){P}(a::OuterProduct{P}, b::OuterProduct{P}) = a.scalar == b.scalar && a.kt == b.kt && a.br == b.br

Base.hash(op::OuterProduct) = hash(op.scalar, hash(op.kt, hash(op.br)))
Base.length(op::OuterProduct) = length(op.kt)*length(op.br)

Base.getindex(op::OuterProduct, k::Array, b::Array) = op.scalar * op.kt[k] * op.br[b]
Base.getindex(op::OuterProduct, label::OpLabel) = op[ktlabel(label), brlabel(label)]
Base.getindex(op::OuterProduct, k, b) = op[[k],[b]]

# would be great if the below worked with normal indexing
# notation (e.g. op[k,:]) but the Colon is apparently
# special and doesn't dispatch directly to getindex...
# Base.getindex(op::OuterProduct, k, ::Colon) = (op.scalar * op.kt[k]) * op.br
# Base.getindex(op::OuterProduct, ::Colon, b) = (op.scalar * op.br[b]) * op.kt
# Base.getindex(op::OuterProduct, ::Colon, ::Colon) = convert(GenericOp, op)

getbra(op::OuterProduct, k::Array) = (op.scalar * get(op.kt,k)) * op.br
getket(op::OuterProduct, b::Array) = (op.scalar * get(op.br,b)) * op.kt

Base.haskey(op::OuterProduct, k::Array, b::Array) = hasket(op,k) && hasbra(op, b)
Base.haskey(op::OuterProduct, label::OpLabel) = haskey(op, ktlabel(label), brlabel(label))
hasket(op::OuterProduct, label::Array) = haskey(op.kt, label)
hasbra(op::OuterProduct, label::Array) = haskey(op.br, label)

Base.get(op::OuterProduct, k::Array, b::Array, default=0) = haskey(op, k, b) ? op[k,b] : default
Base.get(op::OuterProduct, k, b, default=0) = get(op, collect(k), collect(b), default)
Base.get(op::OuterProduct, label::OpLabel, default=0) = get(op, ktlabel(label), brlabel(label), default)

label_from_pair(pair) = OpLabel(pair[1],pair[2])
labels(op::OuterProduct) = imap(label_from_pair, product(labels(op.kt), labels(op.br)))
QuBase.coeffs(op::OuterProduct) = imap(v->op.scalar*prod(v), product(coeffs(op.kt), coeffs(op.br)))

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
inner(a::OuterProduct, b::GeneralOp) = a.scalar * a.kt * inner(a.br, b)
inner(a::GeneralOp, b::OuterProduct) = inner(a, b.kt) * b.br * b.scalar

##########
# act_on #
##########
act_on(op::OuterProduct, kt::Ket, i) = act_on(convert(GenericOp, op), kt, i)

##########
# tensor #
##########
QuBase.tensor(kt::Ket, br::Bra) = OuterProduct(1, kt, br)
QuBase.tensor(br::Bra, kt::Ket) = tensor(kt, br)
QuBase.tensor(a::OuterProduct, b::OuterProduct) = OuterProduct(a.scalar * b.scalar, tensor(a.kt,b.kt), tensor(a.br, b.br))

###########
# Scaling #
###########
Base.scale!(c::Number, op::OuterProduct) = (op.scalar = c*op.scalar; return op)
Base.scale!(op::OuterProduct, c::Number) = (op.scalar = op.scalar*c; return op)

Base.scale(c::Number, op::OuterProduct) = scale!(c,copy(op))
Base.scale(op::OuterProduct, c::Number) = scale!(copy(op),c)

###########
# + and - #
###########
Base.(:-)(op::OuterProduct) = scale(-1, op)
Base.(:+)(a::OuterProduct, b::OuterProduct) = convert(GenericOp, a) + convert(GenericOp, b)

#################
# Normalization #
#################
function Base.norm(op::OuterProduct)
    result = 0
    for v in values(dict(op.kt)), c in values(dict(op.br))
        result += abs2(op.scalar * v * c')
    end
    return sqrt(result)
end

#########
# Trace #
#########
function Base.trace(op::OuterProduct{Orthonormal})
    result = 0
    for k in labels(op.kt), b in labels(op.br)
        if b==k
            result += op[k,b]
        end
    end
    return result
end

function Base.trace(op::OuterProduct)
    result = 0
    prodtype = ptype(op)
    for i in labels(op.kt), (k,v) in dict(op.kt), (b,c) in dict(op.br)
        result += v*c'*inner_rule(prodtype, i, k) * inner_rule(prodtype, b, i)
    end
    return op.scalar * result
end

#################
# Partial Trace #
#################
ptrace{P}(op::OuterProduct{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())
ptrace(op::OuterProduct, over) = ptrace_proj(op, over)

function ptrace_proj(op::OuterProduct{Orthonormal}, over)
    result = OpDict()
    for k in labels(op.kt), b in labels(op.br)
        if k[over] == b[over]
            add_to_dict!(result,
                         OpLabel(except(k, over), except(b, over)),
                         op[k,b])
        end
    end
    return GenericOp(result, ptype(op), decr(fact(op)))
end

function ptrace_proj{P,N}(op::OuterProduct{P,N}, over)
    result = OpDict()
    prodtype = ptype(op)
    for i in labels(op.kt), (k,v) in dict(op.kt), (b,c) in dict(op.br)
        add_to_dict!(result,
                     OpLabel(except(k, over), except(b, over)),
                     op.scalar*v*c'*inner_rule(prodtype, i[over], k[over])
                     *inner_rule(prodtype, b[over], i[over]))
    end
    return GenericOp(result, ptype(op), decr(fact(op)))
end

########################
# Misc. Math Functions #
########################
nfactors{P,N}(::OuterProduct{P,N}) = N

xsubspace(op::OuterProduct,x) = xsubspace(convert(GenericOp, op), x)
filternz(op::OuterProduct) = filternz(convert(GenericOp, op))
purity(op::OuterProduct) = trace(op^2)

######################
# Printing Functions #
######################
labelrepr(op::OuterProduct, k, b, pad) = "$pad$(op[k,b]) $(ktstr(k))$(brstr(b))"

function Base.show(io::IO, op::OuterProduct)
    print(io, summary(op)*":")
    pad = "  "
    maxlen = 4
    for k in take(keys(dict(op.kt)), maxlen),
        b in take(keys(dict(op.br)), maxlen)
        println(io)
        print(io, labelrepr(op, k, b, pad))
    end
    if length(op) > maxlen^2
        println(io)
        print(io, "$pad$vdots")
    end
end

export getbra,
    getket,
    hasket,
    hasbra,
    ptrace,
    purity,
    labels,
    filternz,
    xsubspace,
    nfactors
