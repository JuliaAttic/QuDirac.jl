####################
# GenericOp/DualOp #
####################
abstract GeneralOp{P,N,T} <: DiracOp{P,N,T}

typealias OpDict{N,T} Dict{OuterLabel{N},T}

type GenericOp{P,N,T} <: GeneralOp{P,N,T}
    ptype::P
    dict::OpDict{N,T}
    GenericOp(ptype, dict) = new(ptype, dict)
    GenericOp(ptype, dict::OpDict{0}) = error("Cannot construct a 0-factor operator; did you mean to construct a scalar?")
end

GenericOp{P,N,T}(ptype::P, dict::OpDict{N,T}) = GenericOp{P,N,T}(ptype, dict)
GenericOp{P,N,A,B}(kt::Ket{P,N,A}, br::Bra{P,N,B}) = cons_outer!(OpDict{N,promote_type(A,B)}(), kt, br)

type DualOp{P,N,T} <: GeneralOp{P,N,T}
    op::GenericOp{P,N,T}
end

DualOp{P,N,T}(op::GenericOp{P,N,T}) = DualOp{P,N,T}(op)
DualOp(items...) = DualOp(GenericOp(items...))

Base.convert(::Type{GenericOp}, opc::DualOp) = eager_ctran(opc.op)
Base.promote_rule{G<:GenericOp, D<:DualOp}(::Type{G}, ::Type{D}) = GenericOp

################
# Constructors #
################
function cons_outer!(result, kt, br)
    for (k,kc) in dict(kt)
        for (b,bc) in dict(br)
            newc = kc * bc'
            if newc != 0
                result[OuterLabel(k, b)] = newc
            end
        end
    end
    return GenericOp(ptype(kt), result)
end

func_op(f, kt::Ket) = cons_func_op(f, collect(keys(dict(kt))), kt)

function cons_func_op(f, keys, kt)
    pairs = map(f, keys)
    result = OpDict{eltype(pairs)...}()
    return GenericOp(ptype(kt), load_func_dict!(result, keys, pairs))
end

function load_func_dict!(result, keys, pairs)
    for j=1:length(pairs)
        setfunckv!(result, keys[j], pairs[j][1], pairs[j][2])
    end
    return result
end

function setfunckv!(result, label, c, new_label)
    if c != 0
        result[OuterLabel(new_label, label)] = c
    end
    return c
end

func_op(f, kt::Ket, i) = cons_func_op(f, collect(keys(dict(kt))), kt, i)

function cons_func_op{P,N}(f, keys, kt::Ket{P,N}, i)
    pairs = map(label->f(label[i]), keys)
    result = OpDict{N, eltype(pairs)[1]}()
    return GenericOp(ptype(kt), load_func_dict!(result, keys, pairs, i))
end

function load_func_dict!(result, keys, pairs, i)
    for j=1:length(pairs)
        setfunckv!(result, keys[j], pairs[j][1], pairs[j][2], i)
    end
    return result
end

function setfunckv!(result, label, c, new_i, i)
    if c != 0
        result[OuterLabel(setindex(label, new_i, i), label)] = c
    end
    return c
end

######################
# Accessor functions #
######################
dict(op::GenericOp) = op.dict
dict(opc::DualOp) = dict(opc.op)

ptype(op::GenericOp) = op.ptype
ptype(opc::DualOp) = ptype(opc.op)

#######################
# Dict-Like Functions #
#######################
Base.eltype{P,N,T}(::GeneralOp{P,N,T}) = T

Base.copy(op::GenericOp) = GenericOp(ptype(op), copy(dict(op)))
Base.copy(opc::DualOp) = DualOp(copy(opc.op))

Base.similar(op::GenericOp, d=similar(dict(kt)); P=ptype(op)) = GenericOp(P, d)
Base.similar(opc::DualOp, d=similar(dict(br)); P=ptype(opc)) = DualOp(P, d)

Base.(:(==)){P,N}(a::GenericOp{P,N}, b::GenericOp{P,N}) = dict(a) == dict(b)
Base.(:(==)){P,N}(a::DualOp{P,N}, b::DualOp{P,N}) = a.op == b.op
Base.(:(==))(a::DiracOp, b::DiracOp) = ==(promote(a,b)...)

Base.hash(op::GeneralOp) = hash(dict(op), hash(typeof(op)))

Base.length(op::GeneralOp) = length(dict(op))

Base.getindex(op::GenericOp, label::OuterLabel) = op.dict[label]
Base.getindex(op::GenericOp, k::StateLabel, b::StateLabel) = op.dict[OuterLabel(k,b)]
Base.getindex(opc::DualOp, label::OuterLabel) = opc.op[reverse(label)]'
Base.getindex(opc::DualOp, k::StateLabel, b::StateLabel) = opc.op[OuterLabel(b,k)]'
Base.getindex(op::GeneralOp, k, b) = op[StateLabel(k), StateLabel(b)]

Base.setindex!(op::GenericOp, c, label::OuterLabel) = (op.dict[label] = c)
Base.setindex!(op::GenericOp, c, k::StateLabel, b::StateLabel) = (op.dict[OuterLabel(k,b)] = c)
Base.setindex!(opc::DualOp, c, label::OuterLabel) = (opc.op[reverse(label)] = c')
Base.setindex!(opc::DualOp, c, k::StateLabel, b::StateLabel) = (opc.op[OuterLabel(b,k)] = c')
Base.setindex!(op::GeneralOp, c, k, b) = setindex!(op, c, StateLabel(k), StateLabel(b))

Base.haskey(op::GenericOp, label::OuterLabel) = haskey(dict(op), label)
Base.haskey(opc::DualOp, label::OuterLabel) = haskey(opc.op, reverse(label))
Base.haskey(op::GeneralOp, k, b) = haskey(op, OuterLabel(k, b))

Base.get(op::GenericOp, label::OuterLabel, default=0) = get(dict(op), label, default)
Base.get(opc::DualOp, label::OuterLabel, default=0) = get(dict(opc), label, default')'
Base.get(op::GeneralOp, k, b, default=0) = get(op, OuterLabel(k, b), default)

Base.delete!(op::GenericOp, label::OuterLabel) = (delete!(dict(op), label); return op)
Base.delete!(opc::DualOp, label::OuterLabel) = delete!(opc.op, reverse(label))
Base.delete!(op::GeneralOp, k, b) = delete!(op, OuterLabel(k, b))

labels(op::GenericOp) = keys(dict(op))
labels(opc::DualOp) = imap(reverse, labels(opc.op))
QuBase.coeffs(op::GenericOp) = values(dict(op))
QuBase.coeffs(opc::DualOp) = imap(conj, coeffs(opc.op))

#############
# ctranpose #
#############
eager_ctran(op::GenericOp) = map(ctpair, op)

Base.ctranspose(op::GenericOp) = DualOp(op)
Base.ctranspose(opc::DualOp) = opc.op

#########
# inner #
#########
function inner{P,N,A,B}(br::Bra{P,N,A}, op::GenericOp{P,N,B})
    prodtype = ptype(op)
    result = StateDict{N, promote_type(A,B,inner_type(prodtype,N))}()
    return Bra(prodtype, inner_load!(result, br, op, prodtype))
end

function inner{P,N,A,B}(op::GenericOp{P,N,A}, kt::Ket{P,N,B})
    prodtype = ptype(op)
    result = StateDict{N, promote_type(A,B,inner_type(prodtype,N))}()
    return Ket(prodtype, inner_load!(result, op, kt, prodtype))
end

function inner{P,N,A,B}(a::GenericOp{P,N,A}, b::GenericOp{P,N,B})
    prodtype = ptype(a)
    result = OpDict{N, promote_type(A,B,inner_type(prodtype,N))}()
    return GenericOp(prodtype, inner_load!(result, a, b, prodtype))
end

function inner{P,N,A,B}(a::GenericOp{P,N,A}, b::DualOp{P,N,B})
    prodtype = ptype(a)
    result = OpDict{N, promote_type(A,B,inner_type(prodtype,N))}()
    return GenericOp(prodtype, inner_load!(result, a, b, prodtype))
end

function inner{P,N,A,B}(a::DualOp{P,N,A}, b::GenericOp{P,N,B})
    prodtype = ptype(a)
    result = OpDict{N, promote_type(A,B,inner_type(prodtype,N))}()
    return GenericOp(prodtype, inner_load!(result, a, b, prodtype))
end

inner(br::Bra, opc::DualOp) = inner(opc.op, br')'
inner(opc::DualOp, kt::Ket) = inner(kt', opc.op)'
inner(a::DualOp, b::DualOp) = inner(a.op, b.op)'

function inner_load!(result, a::GenericOp, b::GenericOp, prodtype)
    for (o1,v) in dict(a), (o2,c) in dict(b)
        add_to_dict!(result, 
                     OuterLabel(klabel(o1), blabel(o2)),
                     v*c*inner_rule(prodtype, blabel(o1), klabel(o2)))
    end
    return result
end

function inner_load!(result, a::GenericOp, b::DualOp, prodtype)
    for (o1,v) in dict(a), (o2,c) in dict(b)
        add_to_dict!(result, 
                     OuterLabel(klabel(o1), klabel(o2)),
                     v*c'*inner_rule(prodtype, blabel(o1), blabel(o2)))
    end
    return result
end

function inner_load!(result, a::DualOp, b::GenericOp, prodtype)
    for (o1,v) in dict(a), (o2,c) in dict(b)
        add_to_dict!(result,
                     OuterLabel(blabel(o1), blabel(o2)),
                     v'*c*inner_rule(prodtype, klabel(o1), klabel(o2)))
    end
    return result
end

function inner_load!(result, br::Bra, op::GenericOp, prodtype)
    for (o,v) in dict(op)
        add_to_dict!(result, blabel(o), brcoeff(dict(br), prodtype, klabel(o), v))
    end
    return result
end

function inner_load!(result, op::GenericOp, kt::Ket, prodtype)
    for (o,v) in dict(op)
        add_to_dict!(result, klabel(o), ktcoeff(dict(kt), prodtype, blabel(o), v))
    end
    return result
end

function brcoeff(brdict, prodtype, klabel, v)
    coeff = 0
    for (blabel,c) in brdict
        coeff += c'*v*inner_rule(prodtype, klabel, blabel) 
    end
    return coeff'
end

function ktcoeff(ktdict, prodtype, blabel, v)
    coeff = 0
    for (klabel,c) in ktdict
        coeff += c*v*inner_rule(prodtype, klabel, blabel) 
    end
    return coeff
end

Base.(:*)(br::Bra, op::DiracOp) = inner(br,op)
Base.(:*)(op::DiracOp, kt::Ket) = inner(op,kt)
Base.(:*)(a::DiracOp, b::DiracOp) = inner(a,b)

##########
# act_on #
##########
act_on(op::GeneralOp, kt::Ket, i) = error("inner(op::DiracOp,k::Ket,i) is only defined when nfactors(op) == 1")

# clear up ambiguity warnings
act_on{P}(op::GenericOp{P,1}, kt::Ket{P,1}, i) = invoke(act_on, (GeneralOp{P,1}, Ket{P,1}, Any), op, kt, i)
act_on{P}(op::DualOp{P,1}, kt::Ket{P,1}, i) = invoke(act_on, (GeneralOp{P,1}, Ket{P,1}, Any), op, kt, i)

function act_on{P}(op::GeneralOp{P,1}, kt::Ket{P,1}, i)
    if i==1
        return inner(op, kt)
    else
        throw(BoundsError())
    end
end

function act_on{P,N,A,B}(op::GenericOp{P,1,A}, kt::Ket{P,N,B}, i)
    prodtype = ptype(op)
    result = StateDict{N, promote_type(A,B,inner_type(prodtype,N))}()
    return Ket(prodtype, act_on_dict!(result, op, kt, i, prodtype))
end

function act_on{P,N,A,B}(op::DualOp{P,1,A}, kt::Ket{P,N,B}, i)
    prodtype = ptype(op)
    result = StateDict{N, promote_type(A,B,inner_type(prodtype,N))}()
    return Ket(prodtype, act_on_dict!(result, op, kt, i, prodtype))
end

function act_on_dict!(result, op::GenericOp, kt::Ket, i, prodtype)
    for (o,c) in dict(op), (k,v) in dict(kt)
        add_to_dict!(result, 
                     setindex(k, klabel(o)[1], i),
                     c * v * inner_rule(prodtype, blabel(o)[1], k[i]))
    end
    return result
end

function act_on_dict!(result, op::DualOp, kt::Ket, i, prodtype)
    for (o,c) in dict(op), (k,v) in dict(kt)
        add_to_dict!(result,
                     setindex(k, blabel(o)[1], i),
                     c' * v * inner_rule(prodtype, klabel(o)[1], k[i]))
    end
    return result
end

##########
# tensor #
##########
QuBase.tensor{P,N,M,A,B}(a::GenericOp{P,N,A}, b::GenericOp{P,M,B}) = GenericOp(ptype(a), tensordict!(OpDict{N+M,promote_type(A,B)}(), dict(a), dict(b)))
QuBase.tensor(a::DualOp, b::DualOp) = tensor(a.opc, b.opc)'
QuBase.tensor(a::DiracOp, b::DiracOp) = tensor(promote(a,b)...)

Base.(:*)(kt::Ket, br::Bra) = tensor(kt,br)

###########
# Scaling #
###########
Base.scale!(op::GenericOp, c::Number) = (dscale!(dict(op), c); return op)
Base.scale!(c::Number, op::GenericOp) = scale!(op, c)
Base.scale!(opc::DualOp, c::Number) = DualOp(scale!(opc.op, c'))
Base.scale!(c::Number, opc::DualOp) = scale!(opc, c)

Base.scale(op::GenericOp, c::Number) = similar(op, dscale(dict(op), c))
Base.scale(c::Number, op::GenericOp) = scale(op, c)
Base.scale(opc::DualOp, c::Number) = DualOp(scale(opc.op, c'))
Base.scale(c::Number, opc::DualOp) = scale(opc, c)

Base.(:*)(c::Number, op::DiracOp) = scale(c, op)
Base.(:*)(op::DiracOp, c::Number) = scale(op, c)
Base.(:/)(op::DiracOp, c::Number) = scale(op, 1/c)

###########
# + and - #
###########
Base.(:-)(op::GenericOp) = mapcoeffs(-, op)
Base.(:-)(opc::DualOp) = DualOp(-opc.op)

Base.(:+){P,N}(a::GenericOp{P,N}, b::GenericOp{P,N}) = similar(b, add_merge(dict(a), dict(b)))
Base.(:-){P,N}(a::GenericOp{P,N}, b::GenericOp{P,N}) = similar(b, sub_merge(dict(a), dict(b)))

Base.(:+){P,N}(a::DualOp{P,N}, b::DualOp{P,N}) = DualOp(a.op + b.op)
Base.(:-){P,N}(a::DualOp{P,N}, b::DualOp{P,N}) = DualOp(a.op - b.op)

Base.(:+)(a::DiracOp, b::DiracOp) = +(promote(a,b)...)
Base.(:-)(a::DiracOp, b::DiracOp) = a + (-b)

#################
# Normalization #
#################
Base.norm(op::GenericOp) = sqrt(sum(abs2, values(dict(op))))
Base.norm(opc::DualOp) = norm(opc.op)

QuBase.normalize(op::DiracOp) = scale(1/norm(op), op)
QuBase.normalize!(op::DiracOp) = scale!(1/norm(op), op)

#########
# Trace #
#########
function Base.trace(op::GenericOp{KroneckerDelta})
    result = 0
    for o in labels(op)
        if klabel(o)==blabel(o)
            result += op[o]
        end
    end
    return result
end

function Base.trace(op::GenericOp)
    result = 0
    prodtype = ptype(op)
    for i in distinct(imap(klabel, labels(op))), (o,v) in dict(op)
        result += v * inner_rule(prodtype, i, klabel(o)) * inner_rule(prodtype, blabel(o), i)
    end
    return result
end

Base.trace(opc::DualOp) = trace(opc.op)'

#################
# Partial trace #
#################
ptrace{P}(op::DiracOp{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())
ptrace{N}(op::DiracOp{KroneckerDelta,N}, over) = GenericOp(ptype(op), ortho_ptrace!(OpDict{N-1,eltype(op)}(), op, over))
ptrace{P,N}(op::DiracOp{P,N}, over) = GenericOp(ptype(op), general_ptrace!(OpDict{N-1,promote_type(eltype(op),inner_type(ptype(op),N))}(), op, over))

function ortho_ptrace!(result, op::GeneralOp, over)
    for o in labels(op)
        if klabel(o)[over] == blabel(o)[over]
            add_to_dict!(result, OuterLabel(except(klabel(o), over), except(blabel(o), over)), op[o])
        end
    end
    return result
end

function general_ptrace!(result, op::GeneralOp, over)
    prodtype = ptype(op)
    for i in distinct(imap(klabel, labels(op))), (o,v) in dict(op)
        add_to_dict!(result, 
                     OuterLabel(except(klabel(o), over), except(blabel(o), over)),
                     v * inner_rule(prodtype, i[over], klabel(o)[over]) * inner_rule(prodtype, blabel(o)[over], i[over]))
    end
    return result
end

function general_ptrace!(result, opc::DualOp, over)
    prodtype = ptype(opc)
    for i in distinct(imap(klabel, labels(opc))), (o,v) in dict(opc)
        add_to_dict!(result, 
                     OuterLabel(except(blabel(o), over), except(klabel(o), over)),
                     v' * inner_rule(prodtype, i[over], blabel(o)[over]) * inner_rule(prodtype, klabel(o)[over], i[over]))
    end
    return result
end

########################
# Misc. Math Functions #
########################
nfactors{P,N}(::GeneralOp{P,N}) = N
xsubspace(op::GeneralOp, x) = similar(op, filter((k,v)->is_sum_x(k,x), dict(op)))
filternz!(op::GeneralOp) = (filter!(nzcoeff, dict(op)); return op)
filternz(op::GeneralOp) = similar(op, filter(nzcoeff, dict(op)))

purity(op::GeneralOp) = trace(op^2)

inner_eval(f::Function, op::DiracOp) = mapcoeffs(x->inner_eval(f,x),op)
inner_eval(op::DiracOp) = mapcoeffs(inner_eval,op)

######################
# Printing Functions #
######################
labelrepr(op::GenericOp, o::OuterLabel, pad) = "$pad$(op[o]) $(ktstr(klabel(o)))$(brstr(blabel(o)))"
labelrepr(opc::DualOp, o::OuterLabel, pad) = "$pad$(opc[reverse(o)]) $(ktstr(blabel(o)))$(brstr(klabel(o)))"

Base.summary(op::DiracOp) = "$(typeof(op)) with $(length(op)) operator(s)"
Base.show(io::IO, op::GeneralOp) = dirac_show(io, op)
Base.showcompact(io::IO, op::GeneralOp) = dirac_showcompact(io, op)
Base.repr(op::GeneralOp) = dirac_repr(op)

export GenericOp,
    ptrace,
    xsubspace,
    nfactors,
    maplabels!,
    mapcoeffs!,
    mapcoeffs,
    maplabels,
    filternz,
    filternz!,
    purity,
    labels,
    act_on,
    inner_eval,
    func_op
