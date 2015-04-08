####################
# GenericOp/DualOp #
####################
abstract GeneralOp{P,N} <: DiracOp{P,N}

typealias OpDict{N} Dict{OuterLabel{N},Number}

type GenericOp{P,N} <: GeneralOp{P,N}
    ptype::P
    dict::OpDict{N}
    GenericOp(ptype, dict) = new(ptype, dict)
    GenericOp(ptype, dict::OpDict{0}) = error("Cannot construct a 0-factor operator; did you mean to construct a scalar?")
end

GenericOp{P,N}(ptype::P, dict::OpDict{N}) = GenericOp{P,N}(ptype, dict)

type DualOp{P,N} <: GeneralOp{P,N}
    op::GenericOp{P,N}
end

DualOp{P,N}(op::GenericOp{P,N}) = DualOp{P,N}(op)
DualOp(items...) = DualOp(GenericOp(items...))

Base.convert(::Type{GenericOp}, opc::DualOp) = eager_ctran(opc.op)
Base.convert{P}(::Type{GenericOp{P}}, opc::DualOp{P}) = convert(GenericOp, opc)
Base.convert{P,N}(::Type{GenericOp{P,N}}, opc::DualOp{P,N}) = convert(GenericOp, opc)

Base.promote_rule(::Type{GenericOp}, ::Type{DualOp}) = GenericOp
Base.promote_rule{P}(::Type{GenericOp{P}}, ::Type{DualOp{P}}) = GenericOp{P}
Base.promote_rule{P,N}(::Type{GenericOp{P,N}}, ::Type{DualOp{P,N}}) = GenericOp{P,N}

################
# Constructors #
################
function GenericOp{P,N}(kt::Ket{P,N}, br::Bra{P,N})
    result = OpDict{N}()
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

function func_op{P,N}(f::Function, kt::Ket{P,N})
    result = OpDict{N}()
    for label in labels(kt)
        (c, new_label) = f(label)
        if c != 0
            result[OuterLabel(new_label, label)] = c
        end
    end
    return GenericOp(ptype(kt), result)
end

function func_op{P,N}(f::Function, kt::Ket{P,N}, i)
    result = OpDict{N}()
    for label in labels(kt)
        (c, new_i) = f(label[i])
        if c != 0
            result[OuterLabel(setindex(label, new_i, i), label)] = c
        end
    end
    return GenericOp(ptype(kt), result)
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

Base.get(op::GenericOp, label::OuterLabel, default=0) = get(dict(op), label, default)
Base.get(opc::DualOp, label::OuterLabel, default=0) = haskey(opc, label) ? opc[label] : default

Base.delete!(op::GenericOp, label::OuterLabel) = (delete!(dict(op), label); return op)
Base.delete!(opc::DualOp, label::OuterLabel) = delete!(opc.op, reverse(label))

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
function inner{P,N}(br::Bra{P,N}, op::GenericOp{P,N})
    result = StateDict{N}()
    prodtype = ptype(op)
    for (o,v) in dict(op)
        coeff = 0
        for (b,c) in dict(br)
            coeff += c'*v*inner_rule(prodtype, klabel(o), b) 
        end
        add_to_dict!(result, blabel(o), coeff')
    end
    return Bra(prodtype, result)
end

function inner{P,N}(op::GenericOp{P,N}, kt::Ket{P,N})
    result = StateDict{N}()
    prodtype = ptype(op)
    for (o,c) in dict(op)
        coeff = 0
        for (k,v) in dict(kt)
            coeff += c*v*inner_rule(prodtype, blabel(o), k) 
        end
        add_to_dict!(result, klabel(o), coeff)
    end
    return Ket(prodtype, result)
end

function inner{P,N}(a::GenericOp{P,N}, b::GenericOp{P,N})
    result = OpDict{N}()
    prodtype = ptype(a)
    for (o1,v) in dict(a), (o2,c) in dict(b)
        add_to_dict!(result,
                     OuterLabel(klabel(o1), blabel(o2)),
                     v*c*inner_rule(prodtype, blabel(o1), klabel(o2)))
    end
    return GenericOp(prodtype, result)
end

function inner{P,N}(a::GenericOp{P,N}, b::DualOp{P,N})
    result = OpDict{N}()
    prodtype = ptype(a)
    for (o1,v) in dict(a), (o2,c) in dict(b)
        add_to_dict!(result,
                     OuterLabel(klabel(o1), klabel(o2)),
                     v*c'*inner_rule(prodtype, blabel(o1), blabel(o2)))
    end
    return GenericOp(prodtype, result)
end

function inner{P,N}(a::DualOp{P,N}, b::GenericOp{P,N})
    result = OpDict{N}()
    prodtype = ptype(a)
    for (o1,v) in dict(a), (o2,c) in dict(b)
        add_to_dict!(result,
                     OuterLabel(blabel(o1), blabel(o2)),
                     v'*c*inner_rule(prodtype, klabel(o1), klabel(o2)))
    end
    return GenericOp(prodtype, result)
end

inner(br::Bra, opc::DualOp) = inner(opc.op, br')'
inner(opc::DualOp, kt::Ket) = inner(kt', opc.op)'
inner(a::DualOp, b::DualOp) = inner(a.op, b.op)'

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

function act_on{P,N}(op::GenericOp{P,1}, kt::Ket{P,N}, i)
    result = StateDict{N}()
    prodtype = ptype(op)
    for (o,c) in dict(op), (k,v) in dict(kt)
        add_to_dict!(result,
                     setindex(k, klabel(o)[1], i),
                     c * v * inner_rule(prodtype, blabel(o)[1], k[i]))
    end
    return Ket(prodtype, result)
end

function act_on{P,N}(op::DualOp{P,1}, kt::Ket{P,N}, i)
    result = StateDict{N}()
    prodtype = ptype(op)
    for (o,c) in dict(op), (k,v) in dict(kt)
        add_to_dict!(result,
                     setindex(k, blabel(o)[1], i),
                     c' * v * inner_rule(prodtype, klabel(o)[1], k[i]))
    end
    return Ket(prodtype, result)
end

##########
# tensor #
##########
QuBase.tensor{P,A,B}(a::GenericOp{P,A}, b::GenericOp{P,B}) = GenericOp(ptype(a), tensordict!(OpDict{A+B}(), dict(a), dict(b)))
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
function Base.trace(op::GenericOp{Orthonormal})
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
ptrace{P}(op::GeneralOp{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())
ptrace(op::GeneralOp, over) = ptrace_op(op, over)

function ptrace_op{N}(op::GeneralOp{Orthonormal,N}, over)
    result = OpDict{N-1}()
    for o in labels(op)
        if klabel(o)[over] == blabel(o)[over]
            add_to_dict!(result, traced_label(op,o,over), op[o])
        end
    end
    return GenericOp(ptype(op), result)
end

function ptrace_op{P,N}(op::GenericOp{P,N}, over)
    result = OpDict{N-1}()
    prodtype = ptype(op)
    for i in distinct(imap(klabel, labels(op))), (o,v) in dict(op)
        add_to_dict!(result, 
                     OuterLabel(except(klabel(o), over), except(blabel(o), over)), 
                     v * inner_rule(prodtype, i[over], klabel(o)[over]) 
                     * inner_rule(prodtype, blabel(o)[over], i[over]))
    end
    return GenericOp(prodtype, result)
end

function ptrace_op{P,N}(opc::DualOp{P,N}, over)
    result = OpDict{N-1}()
    prodtype = ptype(opc)
    for i in distinct(imap(klabel, labels(opc))), (o,v) in dict(opc)
        add_to_dict!(result, 
                     OuterLabel(except(blabel(o), over), except(klabel(o), over)), 
                     v' * inner_rule(prodtype, i[over], blabel(o)[over]) 
                     * inner_rule(prodtype, klabel(o)[over], i[over]))
    end
    return GenericOp(prodtype, result)
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
labelrepr(op::GenericOp, label, pad) = "$pad$(op[label]) $(ktstr(first(label)))$(brstr(second(label)))"
labelrepr(opc::DualOp, label, pad) = "$pad$(opc[reverse(label)]) $(ktstr(second(label)))$(brstr(first(label)))"

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
