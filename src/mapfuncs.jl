#############
# Filtering #
#############
Base.filter!(f::Function, kt::Ket) = (filter!(f, dict(kt)); return kt)
Base.filter!(f::Function, br::Bra) = (filter!((k,v)->f(k,v'), dict(br)); return br)
Base.filter!(f::Function, op::GenericOp) = (filter!(f, dict(op)); return op)
Base.filter!(f::Function, opc::DualOp) = (filter!((k,v)->f(reverse(k),v'), dict(opc)); return opc)

Base.filter(f::Function, kt::Ket) = similar(kt, filter(f, dict(kt)))
Base.filter(f::Function, br::Bra) = similar(br, filter((k,v)->f(k,v'), dict(br)))
Base.filter(f::Function, op::GenericOp) = similar(op, filter(f, dict(op)))
Base.filter(f::Function, opc::DualOp) = similar(opc, filter((k,v)->f(reverse(k),v'), dict(opc)))
Base.filter(f::Function, op::OuterProduct) = filter(f, convert(GenericOp, op))

########################
# mapcoeffs!/mapcoeffs #
########################
mapcoeffs!(f::Function, k::Ket) = (mapvals!(f, dict(k)); return k)
mapcoeffs!(f::Function, b::Bra) = (mapvals!(v->f(v')', dict(b)); return b)
mapcoeffs!(f::Function, op::GenericOp) = (mapvals!(f, dict(op)); return op)
mapcoeffs!(f::Function, opc::DualOp) = (mapvals!(v->f(v')', dict(opc)); return opc)

mapcoeffs(f::Function, kt::Ket) = similar(kt, mapvals(f, dict(kt)))
mapcoeffs(f::Function, br::Bra) = similar(br, mapvals(v->f(v')', dict(br)))
mapcoeffs(f::Function, op::GenericOp) = similar(op, mapvals(f, dict(op)))
mapcoeffs(f::Function, opc::DualOp) = similar(opc, mapvals(v->f(v')', dict(opc)))
mapcoeffs(f::Function, op::OuterProduct) = mapcoeffs(f, convert(GenericOp, op))

########################
# maplabels!/maplabels #
########################
function maplabels!(f::Function, obj::AbstractDirac)
    for (label,v) in dict(obj)
        replace_label!(dict(obj), label, f(label), v)
    end
    return obj
end

maplabels!(f::Function, ::OuterProduct) = error("cannot mutate OuterProduct; perhaps you meant to convert it to GenericOp before calling maplabels!")

function maplabels!(f::Function, opc::DualOp)
    for (label,v) in dict(opc)
        replace_label!(dict(opc), label, reverse(f(reverse(label))), v)
    end
    return opc
end

function maplabels(f::Function, s::DiracState)
    nguess = nguess_label(f, dict(s))
    result = StateDict{nguess}()
    load_labels!(f, result, s)
    return similar(s, result)
end

function maplabels(f::Function, op::GeneralOp)
    nguess = nguess_label(f, dict(op))
    result = OpDict{nguess}()
    load_labels!(f, result, op)
    return GenericOp(ptype(op), result)
end

maplabels(f::Function, op::OuterProduct) = maplabels(f, convert(GenericOp, op))

#######
# map #
#######
function Base.map(f::Function, s::DiracState)
    nguess = nguess_pair(f, dict(s))
    result = mapload!(f, StateDict{nguess}(), s) 
    return similar(s, result)
end

function Base.map(f::Function, op::GeneralOp)
    nguess = nguess_pair(f, dict(op))
    result = mapload!(f, OpDict{nguess}(), op) 
    return GenericOp(ptype(op), result)
end

Base.map(f::Function, op::OuterProduct) = map(f, convert(GenericOp, op))

################
# Helper funcs #
################
function load_labels!(f, result, obj)
    for (label, v) in dict(obj)
        result[f(label)] = v
    end
    return result
end

function load_labels!(f, result, opc::DualOp)
    for (label,v) in dict(obj)
        result[f(reverse(label))] = v'
    end
    return result
end

function replace_label!(d, old_label, new_label, v)
    delete!(d, old_label)  
    d[new_label] = v
    return d
end

nguess_pair(f, d::StateDict) = length(f(first(d)...)[1])
nguess_pair(f, d::OpDict) = length(f(first(d)...)[1][1])
nguess_label(f, d::StateDict) = length(f(first(keys(d))))
nguess_label(f, d::OpDict) = length(f(first(keys(d))[1]))

function mapload!(f, result, obj::AbstractDirac)
    for (label,v) in dict(obj)
        (new_label, new_v) = f(label, v)
        result[new_label] = new_v
    end
    return result
end

function mapload!(f, result, br::Bra)
    for (label,v) in dict(br)
        (new_label, new_v) = f(label, v')
        result[new_label] = new_v'
    end
    return result
end

function mapload!(f, result, opc::DualOp)
    for (label,v) in dict(op)
        (new_label, new_v) = f(reverse(label), v')
        result[new_label] = new_v
    end
    return result
end

export maplabels!,
    mapcoeffs!,
    maplabels,
    mapcoeffs