####################
# GenericOp/DualOp #
####################
abstract GeneralOp{P,N} <: DiracOp{P,N}

typealias OpDict{N} Dict{(NTuple{N},NTuple{N}),Number}

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
                result[k,b] = newc
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
            result[new_label, label] = c
        end
    end
    return GenericOp(ptype(kt), result)
end

function func_op{P,N}(f::Function, kt::Ket{P,N}, i)
    result = OpDict{N}()
    for label in labels(kt)
        (c, new_i) = f(label[i])
        if c != 0
            result[placeat(label, new_i, i), label] = c
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

Base.similar(op::GenericOp, d::OpDict=similar(dict(kt)); P=ptype(op)) = GenericOp(P, d)
Base.similar(opc::DualOp, d::OpDict=similar(dict(br)); P=ptype(opc)) = DualOp(P, d)

Base.(:(==)){P,N}(a::GenericOp{P,N}, b::GenericOp{P,N}) = dict(a) == dict(b)
Base.(:(==)){P,N}(a::DualOp{P,N}, b::DualOp{P,N}) = a.op == b.op
Base.(:(==))(a::DiracOp, b::DiracOp) = ==(promote(a,b)...)

Base.hash(op::GeneralOp) = hash(dict(op), hash(typeof(op)))

Base.length(op::GeneralOp) = length(dict(op))

Base.getindex(op::GenericOp, label::(Tuple,Tuple)) = op.dict[label]
Base.getindex(op::GenericOp, k::Tuple, b::Tuple) = op.dict[k,b]
Base.getindex(opc::DualOp, label::(Tuple,Tuple)) = opc.op[reverse(label)]'
Base.getindex(opc::DualOp, k::Tuple, b::Tuple) = opc.op[b,k]'
Base.getindex(op::GeneralOp, k, b) = op[tuple(k), tuple(b)]

Base.setindex!(op::GenericOp, c, label::(Tuple,Tuple)) = (op.dict[label] = c)
Base.setindex!(op::GenericOp, c, k::Tuple, b::Tuple) = (op.dict[k,b] = c)
Base.setindex!(opc::DualOp, c, label::(Tuple,Tuple)) = (opc.op[reverse(label)] = c')
Base.setindex!(opc::DualOp, c, k::Tuple, b::Tuple) = (opc.op[b,k] = c')
Base.setindex!(op::GeneralOp, c, k, b) = setindex!(op, c, tuple(k), tuple(b))

Base.haskey(op::GenericOp, label::(Tuple,Tuple)) = haskey(dict(op), label)
Base.haskey(opc::DualOp, label::(Tuple,Tuple)) = haskey(opc.op, reverse(label))
Base.haskey(op::GeneralOp, k::Tuple, b::Tuple) = haskey(op, (k,b))

Base.get(op::GenericOp, label::(Tuple,Tuple), default=0) = get(dict(op), label, default)
Base.get(opc::DualOp, label::(Tuple,Tuple), default=0) = haskey(opc, label) ? opc[label] : default

Base.delete!(op::GenericOp, label::(Tuple,Tuple)) = (delete!(dict(op), label); return op)
Base.delete!(opc::DualOp, label::(Tuple,Tuple)) = delete!(opc.op, reverse(label))
Base.delete!(op::GeneralOp, k::Tuple, b::Tuple) = delete!(op, (k,b))
Base.delete!(op::GeneralOp, k, b) = delete!(op, tuple(k), tuple(b))

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
    for (label,v) in dict(op)
        coeff = 0
        for (b,c) in dict(br)
            coeff += c'*v*inner_rule(prodtype, first(label), b) 
        end
        add_to_dict!(result, second(label), coeff')
    end
    return Bra(prodtype, result)
end

function inner{P,N}(op::GenericOp{P,N}, kt::Ket{P,N})
    result = StateDict{N}()
    prodtype = ptype(op)
    for (label,c) in dict(op)
        coeff = 0
        for (k,v) in dict(kt)
            coeff += c*v*inner_rule(prodtype, second(label), k) 
        end
        add_to_dict!(result, first(label), coeff)
    end
    return Ket(prodtype, result)
end

function inner{P,N}(a::GenericOp{P,N}, b::GenericOp{P,N})
    result = OpDict{N}()
    prodtype = ptype(a)
    for (label1,v) in dict(a), (label2,c) in dict(b)
        add_to_dict!(result,
                     (first(label1), second(label2)),
                     v*c*inner_rule(prodtype, second(label1), first(label2)))
    end
    return GenericOp(prodtype, result)
end

function inner{P,N}(a::GenericOp{P,N}, b::DualOp{P,N})
    result = OpDict{N}()
    prodtype = ptype(a)
    for (label1,v) in dict(a), (label2,c) in dict(b)
        add_to_dict!(result,
                     (first(label1), first(label2)),
                     v*c'*inner_rule(prodtype, second(label1), second(label2)))
    end
    return GenericOp(prodtype, result)
end

function inner{P,N}(a::DualOp{P,N}, b::GenericOp{P,N})
    result = OpDict{N}()
    prodtype = ptype(a)
    for (label1,v) in dict(a), (label2,c) in dict(b)
        add_to_dict!(result,
                     (second(label1), second(label2)),
                     v'*c*inner_rule(prodtype, first(label1), first(label2)))
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
    for (label, c) in dict(op), (k,v) in dict(kt)
        add_to_dict!(result,
                     placeat(k, first(label)[1], i),
                     c * v * inner_rule(prodtype, first(second(label)), k[i]))
    end
    return Ket(prodtype, result)
end

function act_on{P,N}(op::DualOp{P,1}, kt::Ket{P,N}, i)
    result = StateDict{N}()
    prodtype = ptype(op)
    for (label, c) in dict(op), (k,v) in dict(kt)
        add_to_dict!(result,
                     placeat(k, second(label)[1], i),
                     c' * v * inner_rule(prodtype, first(first(label)), k[i]))
    end
    return Ket(prodtype, result)
end

##########
# tensor #
##########
QuBase.tensor{P,A,B}(a::GenericOp{P,A}, b::GenericOp{P,B}) = GenericOp(ptype(a), tensorop!(OpDict{A+B}(), dict(a), dict(b)))
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
    for label in labels(op)
        if first(label)==second(label)
            result += op[label]
        end
    end
    return result
end

function Base.trace(op::GenericOp)
    result = 0
    prodtype = ptype(op)
    for i in distinct(imap(first, labels(op))), (label,v) in dict(op)
        result += v * inner_rule(prodtype, i, first(label)) * inner_rule(prodtype, second(label), i)
    end
    return result
end

Base.trace(opc::DualOp) = trace(opc.op)'

#################
# Partial trace #
#################
ptrace{P}(op::GeneralOp{P,1}, over) = over == 1 ? trace(op) : throw(BoundsError())
ptrace(op::GeneralOp, over) = ptrace_op(op, over)

traced_label(::GenericOp, label, over) = (except(first(label), over), except(second(label), over))
traced_label(::DualOp, label, over) = (except(second(label), over), except(first(label), over))

traced_coeff(op::GenericOp, i, label, v, over) = v * inner_rule(ptype(op), i[over], first(label)[over]) * inner_rule(ptype(op), second(label)[over], i[over])
traced_coeff(op::DualOp, i, label, v, over) = v' * inner_rule(ptype(op), i[over], second(label)[over]) * inner_rule(ptype(op), first(label)[over], i[over])

function ptrace_op{N}(op::GeneralOp{Orthonormal,N}, over)
    result = OpDict{N-1}()
    for label in labels(op)
        if first(label)[over] == second(label)[over]
            add_to_dict!(result, traced_label(op,label,over), op[label])
        end
    end
    return GenericOp(ptype(op), result)
end

function ptrace_op{P,N}(op::GeneralOp{P,N}, over)
    result = OpDict{N-1}()
    for i in distinct(imap(first, labels(op))), (label,v) in dict(op)
        add_to_dict!(result, traced_label(op,label,over), traced_coeff(op,i,label,v,over))
    end
    return GenericOp(ptype(op), result)
end

########################
# Misc. Math Functions #
########################
nfactors{P,N}(::GeneralOp{P,N}) = N
xsubspace(op::GeneralOp, x) = similar(op, filter((k,v)->isx(k,x), dict(op)))
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
