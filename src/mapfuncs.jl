##############################
# Mapping functions on Dicts #
##############################
function mapvals!(f, d)
    for (k,v) in d
        d[k] = f(v)
    end
    return d
end
mapvals(f, d) = Dict(collect(keys(d)), map(f, collect(values(d))))

mapkeys(f, d) = Dict(map(f, collect(keys(d))), collect(values(d)))
mapkv(f, d) = Dict(map(f, collect(d)))

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
mapcoeffs!(f::Union(Function,DataType), k::Ket) = (mapvals!(f, dict(k)); return k)
mapcoeffs!(f::Union(Function,DataType), b::Bra) = (mapvals!(v->f(v')', dict(b)); return b)
mapcoeffs!(f::Union(Function,DataType), op::GenericOp) = (mapvals!(f, dict(op)); return op)
mapcoeffs!(f::Union(Function,DataType), opc::DualOp) = (mapvals!(v->f(v')', dict(opc)); return opc)

mapcoeffs(f::Union(Function,DataType), kt::Ket) = similar(kt, mapvals(f, dict(kt)))
mapcoeffs(f::Union(Function,DataType), br::Bra) = similar(br, mapvals(v->f(v')', dict(br)))
mapcoeffs(f::Union(Function,DataType), op::GenericOp) = similar(op, mapvals(f, dict(op)))
mapcoeffs(f::Union(Function,DataType), opc::DualOp) = similar(opc, mapvals(v->f(v')', dict(opc)))
mapcoeffs(f::Union(Function,DataType), op::OuterProduct) = mapcoeffs(f, convert(GenericOp, op))

########################
# maplabels!/maplabels #
########################
maplabels(f::Union(Function,DataType), s::DiracState) = similar(s, mapkeys(f, dict(s)))
maplabels(f::Union(Function,DataType), op::GenericOp) = similar(op, mapkeys(f, dict(op)))
maplabels(f::Union(Function,DataType), opc::DualOp) = GenericOp(ptype(opc), mapkeys(label->f(reverse(label)), dict(op)))
maplabels(f::Union(Function,DataType), op::OuterProduct) = maplabels(f, convert(GenericOp, op))

#######
# map #
#######
Base.map(f::Union(Function,DataType), obj::AbstractDirac) = similar(obj, mapkv(kv->f(kv[1], kv[2]), dict(obj)))
Base.map(f::Union(Function,DataType), br::Bra) = similar(br, mapkv(kv->br_tup(f(kv[1], kv[2]')), dict(d)))
Base.map(f::Union(Function,DataType), opc::DualOp) = GenericOp(ptype(opc), mapkv(kv->f(reverse(kv[1]), kv[2]'), dict(d)))
Base.map(f::Union(Function,DataType), op::OuterProduct) = map(f, convert(GenericOp, op))

br_tup(tup) = (tup[1], tup[2]')

export maplabels!,
    mapcoeffs!,
    maplabels,
    mapcoeffs