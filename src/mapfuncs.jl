##############################
# Mapping functions on Dicts #
##############################
function mapvals!(f, d)
    for (k,v) in d
        d[k] = f(v)
    end
    return d
end
mapvals(f, d) = Dict(zip(collect(keys(d)), map(f, collect(values(d)))))

mapkeys(f, d) = Dict(zip(map(f, collect(keys(d))), collect(values(d))))
mapkv(f, d) = Dict(map(f, collect(d)))

#############
# Filtering #
#############
Base.filter!(f::Function, kt::Ket) = (filter!(f, dict(kt)); return kt)
Base.filter!(f::Function, br::Bra) = (filter!((k,v)->f(k,v'), dict(br)); return br)
Base.filter!(f::Function, op::OpSum) = (filter!(f, dict(op)); return op)
Base.filter!(f::Function, opc::DualOpSum) = (filter!((k,v)->f(k',v'), dict(opc)); return opc)

Base.filter(f::Function, kt::Ket) = similar(kt, filter(f, dict(kt)))
Base.filter(f::Function, br::Bra) = similar(br, filter((k,v)->f(k,v'), dict(br)))
Base.filter(f::Function, op::OpSum) = similar(op, filter(f, dict(op)))
Base.filter(f::Function, opc::DualOpSum) = similar(opc, filter((k,v)->f(k',v'), dict(opc)))
Base.filter(f::Function, op::OuterProduct) = filter(f, convert(OpSum, op))

########################
# mapcoeffs!/mapcoeffs #
########################
mapcoeffs!(f::Function, k::Ket) = (mapvals!(f, dict(k)); return k)
mapcoeffs!(f::Function, b::Bra) = (mapvals!(v->f(v')', dict(b)); return b)
mapcoeffs!(f::Function, op::OpSum) = (mapvals!(f, dict(op)); return op)
mapcoeffs!(f::Function, opc::DualOpSum) = (mapvals!(v->f(v')', dict(opc)); return opc)

mapcoeffs(f::Function, kt::Ket) = similar(kt, mapvals(f, dict(kt)))
mapcoeffs(f::Function, br::Bra) = similar(br, mapvals(v->f(v')', dict(br)))
mapcoeffs(f::Function, op::OpSum) = similar(op, mapvals(f, dict(op)))
mapcoeffs(f::Function, opc::DualOpSum) = similar(opc, mapvals(v->f(v')', dict(opc)))
mapcoeffs(f::Function, op::OuterProduct) = mapcoeffs(f, convert(OpSum, op))

########################
# maplabels!/maplabels #
########################
maplabels(f::Function, s::DiracState) = similar(s, mapkeys(f, dict(s)))
maplabels(f::Function, op::OpSum) = similar(op, mapkeys(f, dict(op)))
maplabels(f::Function, opc::DualOpSum) = OpSum(ptype(opc), mapkeys(label->f(label'), dict(op)))
maplabels(f::Function, op::OuterProduct) = maplabels(f, convert(OpSum, op))

#######
# map #
#######
Base.map(f::Function, obj::AbstractDirac) = similar(obj, mapkv(kv->f(kv[1], kv[2]), dict(obj)))
Base.map(f::Function, br::Bra) = similar(br, mapkv(kv->br_tup(f(kv[1], kv[2]')), dict(d)))
Base.map(f::Function, opc::DualOpSum) = OpSum(ptype(opc), mapkv(kv->f(kv[1]', kv[2]'), dict(d)))
Base.map(f::Function, op::OuterProduct) = map(f, convert(OpSum, op))

br_tup(tup) = (tup[1], tup[2]')

export maplabels!,
    mapcoeffs!,
    maplabels,
    mapcoeffs
