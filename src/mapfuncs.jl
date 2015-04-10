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
maplabels!(f::Function, obj::AbstractDirac) = mapkeys!(f, dict(obj))
maplabels!(f::Function, opc::DualOp) = mapkeys!(label->reverse(f(reverse(label))), dict(opc))
maplabels!(f::Function, ::OuterProduct) = error("cannot mutate OuterProduct; perhaps you meant to convert it to GenericOp before calling maplabels!")

maplabels(f::Function, s::DiracState) = similar(s, mapkeys(f, dict(s)))
maplabels(f::Function, op::GenericOp) = similar(op, mapkeys(f, dict(op)))
maplabels(f::Function, opc::DualOp) = GenericOp(ptype(opc), mapkeys(label->f(reverse(label)), dict(op)))
maplabels(f::Function, op::OuterProduct) = maplabels(f, convert(GenericOp, op))

#######
# map #
#######
Base.map(f::Function, obj::AbstractDirac) = similar(obj, mapkv(kv->f(kv[1], kv[2]), dict(obj)))
Base.map(f::Function, br::Bra) = similar(br, mapkv(kv->br_tup(f(kv[1], kv[2]')), dict(d)))
Base.map(f::Function, opc::DualOp) = GenericOp(ptype(opc), mapkv(kv->f(reverse(kv[1]), kv[2]'), dict(d)))
Base.map(f::Function, op::OuterProduct) = map(f, convert(GenericOp, op))

br_tup(tup) = (tup[1], tup[2]')

export maplabels!,
    mapcoeffs!,
    maplabels,
    mapcoeffs