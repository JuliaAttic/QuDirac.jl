################
# OuterProduct #
################
type OuterProduct{P,N,S,K<:Ket,B<:Bra} <: DiracOp{P,N}
    scalar::S
    kt::K
    br::B
    OuterProduct(scalar::S, kt::Ket{P,N}, br::Bra{P,N}) = new(scalar, kt, br)
end

OuterProduct{P,N,S}(scalar::S, kt::Ket{P,N}, br::Bra{P,N}) = OuterProduct{P,N,S,typeof(kt),typeof(br)}(scalar, kt, br)

ptype(op::OuterProduct) = ptype(op.kt)

Base.copy(op::OuterProduct) = OuterProduct(copy(op.scalar), copy(op.kt), copy(op.br))

Base.convert(::Type{GenericOp}, op::OuterProduct) = scale!(op.scalar, GenericOp(op.kt, op.br))
Base.promote_rule{G<:GenericOp, O<:OuterProduct}(::Type{G}, ::Type{O}) = GenericOp

#######################
# Dict-Like Functions #
#######################
Base.eltype(op::OuterProduct) = promote_type(typeof(op.scalar), eltype(op.kt), eltype(op.br))

# these equality/hash functions are pretty inefficient...
Base.(:(==)){P,N}(a::OuterProduct{P,N}, b::OuterProduct{P,N}) = convert(GenericOp, a) == convert(GenericOp, b) 
Base.hash(op::OuterProduct) = hash(convert(GenericOp, op))
Base.hash(op::OuterProduct, h::Uint64) = hash(hash(op), h)

Base.length(op::OuterProduct) = length(op.kt)*length(op.br)

Base.getindex(op::OuterProduct, k, b) = op.scalar * op.kt[k] * op.br[b]
Base.getindex(op::OuterProduct, o::OuterLabel) = op[klabel(o), blabel(o)]

# would be great if the below worked with normal indexing
# notation (e.g. op[k,:]) but the Colon is apparently
# special and doesn't dispatch directly to getindex...
# Base.getindex(op::OuterProduct, k, ::Colon) = (op.scalar * op.kt[k]) * op.br
# Base.getindex(op::OuterProduct, ::Colon, b) = (op.scalar * op.br[b]) * op.kt
# Base.getindex(op::OuterProduct, ::Colon, ::Colon) = convert(GenericOp, op)

getbra(op::OuterProduct, k) = (op.scalar * get(op.kt,k)) * op.br
getket(op::OuterProduct, b) = (op.scalar * get(op.br,b)) * op.kt

Base.haskey(op::OuterProduct, k, b) = hasket(op, k) && hasbra(op, b)
Base.haskey(op::OuterProduct, o::OuterLabel) = haskey(op, klabel(o), blabel(o))

Base.get(op::OuterProduct, k, b, default=0) = haskey(op, k, b) ? op[k,b] : default
Base.get(op::OuterProduct, o::OuterLabel, default=0) = get(op, klabel(o), blabel(o), default)

Base.collect{P,N}(op::OuterProduct{P,N}) = collect_pairs!(Array((OuterLabel{N}, eltype(op)), length(op)), op)

function collect_pairs!(result, op::OuterProduct)
    i = 1
    for (k,kc) in dict(op.kt)
        for (b,bc) in dict(op.br)
            result[i] = (OuterLabel(k, b), op.scalar * kc * bc')
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
inner(a::OuterProduct, b::GeneralOp) = a.scalar * a.kt * inner(a.br, b)
inner(a::GeneralOp, b::OuterProduct) = inner(a, b.kt) * b.br * b.scalar

##########
# act_on #
##########
act_on(op::OuterProduct, kt::Ket, i) = act_on(convert(GenericOp, op), kt, i)

##########
# tensor #
##########
tensor(kt::Ket, br::Bra) = OuterProduct(1, kt, br)
tensor(br::Bra, kt::Ket) = tensor(kt, br)
tensor(a::OuterProduct, b::OuterProduct) = OuterProduct(a.scalar * b.scalar, tensor(a.kt,b.kt), tensor(a.br, b.br))

###########
# Scaling #
###########
Base.scale!(c::Number, op::OuterProduct) = (op.scalar = c*op.scalar; return op)
Base.scale!(op::OuterProduct, c::Number) = (op.scalar = op.scalar*c; return op)

Base.scale(c::Number, op::OuterProduct) = OuterProduct(c * op.scalar, copy(op.kt), copy(op.br))
Base.scale(op::OuterProduct, c::Number) = OuterProduct(op.scalar * c, copy(op.kt), copy(op.br))

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
function Base.trace(op::OuterProduct)
    result = 0
    for (k,v) in dict(op.kt), (b,c) in dict(op.br)
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
    for k in keys(dict(op.kt)), b in keys(dict(op.br))
        if k[over] == b[over]
            add_to_dict!(result, OuterLabel(except(k, over), except(b, over)), op[k,b])
        end
    end
    return result
end

########################
# Misc. Math Functions #
########################
nfactors{P,N}(::OuterProduct{P,N}) = N
xsubspace(op::OuterProduct,x) = xsubspace(convert(GenericOp, op), x)
filternz(op::OuterProduct) = filternz(convert(GenericOp, op))
switch(op::OuterProduct, i, j) = switch(convert(GenericOp, op), i, j)
permute(op::OuterProduct, perm::Vector) = permute(convert(GenericOp, op), perm)

######################
# Printing Functions #
######################
labelrepr(op::OuterProduct, k, b, pad) = "$pad$(op[k,b]) $(ktstr(k))$(brstr(b))"

Base.summary(op::OuterProduct) = "OuterProduct with $(length(op)) operator(s); $(typeof(op.kt)) * $(typeof(op.br))"

function Base.show(io::IO, op::OuterProduct)
    print(io, summary(op)*":")
    pad = "  "
    maxlen = 10
    i = 1
    for k in keys(dict(op.kt)), b in keys(dict(op.br))
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

export getbra,
    getket,
    hasket,
    hasbra,
    ptrace,
    purity,
    filternz,
    xsubspace,
    nfactors
