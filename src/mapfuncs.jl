#############
# Filtering #
#############
Base.filter!(f::Function, kt::Ket) = (filter!(f, data(kt)); return kt)
Base.filter!(f::Function, br::Bra) = (filter!((k,v)->f(k,v'), data(br)); return br)
Base.filter!(f::Function, op::OpSum) = (filter!(f, data(op)); return op)
Base.filter!(f::Function, opc::DualOpSum) = (filter!((k,v)->f(k',v'), data(opc)); return opc)

Base.filter{P}(f::Function, kt::Ket{P}) = Ket(P, filter(f, data(kt)))
Base.filter{P}(f::Function, op::OpSum{P}) = OpSum(P, filter(f, data(op)))

Base.filter(f::Function, br::Bra) = Bra(filter((k,v)->f(k,v'), br.kt))
Base.filter(f::Function, opc::DualOpSum) = DualOpSum(filter((k,v)->f(k',v'), opc.op))

Base.filter(f::Function, op::OuterProduct) = filter(f, OpSum(op))

########################
# mapcoeffs!/mapcoeffs #
########################
mapcoeffs!(f::Union(Function,DataType), kt::Ket) = (mapvals!(f, data(kt)); return kt)
mapcoeffs!(f::Union(Function,DataType), op::OpSum) = (mapvals!(f, data(op)); return op)

mapcoeffs!(f::Union(Function,DataType), br::Bra) = (mapvals!(v->f(v')', data(br)); return br)
mapcoeffs!(f::Union(Function,DataType), opc::DualOpSum) = (mapvals!(v->f(v')', data(opc)); return opc)

mapcoeffs{P}(f::Union(Function,DataType), kt::Ket{P}) = Ket(P, mapvals(f, data(kt)))
mapcoeffs{P}(f::Union(Function,DataType), op::OpSum{P}) = OpSum(P, mapvals(f, data(op)))

mapcoeffs(f::Union(Function,DataType), br::Bra) = Bra(mapvals(v->f(v')', br.kt))
mapcoeffs(f::Union(Function,DataType), opc::DualOpSum) = DualOpSum(mapvals(v->f(v')', opc.op))

mapcoeffs(f::Union(Function,DataType), op::OuterProduct) = mapcoeffs(f, OpSum(op))

########################
# maplabels!/maplabels #
########################
maplabels{P}(f::Union(Function,DataType), kt::Ket{P}) = Ket(P, mapkeys(f, data(kt)))
maplabels{P}(f::Union(Function,DataType), op::OpSum{P}) = OpSum(P, mapkeys(f, data(op)))

maplabels(f::Union(Function,DataType), br::Bra) = Bra(maplabels(f, br.kt))
maplabels(f::Union(Function,DataType), opc::DualOpSum) = DualOpSum(maplabels(l->f(l')', opc.op))

maplabels(f::Union(Function,DataType), op::OuterProduct) = maplabels(f, OpSum(op))

#######
# map #
#######
Base.map{P}(f::Union(Function,DataType), kt::Ket{P}) = Ket(P, map(f, data(kt)))
Base.map{P}(f::Union(Function,DataType), op::OpSum{P}) = OpSum(P, map(f, data(op)))

Base.map(f::Union(Function,DataType), br::Bra) = Bra(map((k,v)->br_tup(f(k,v')), br.kt))
Base.map{P}(f::Union(Function,DataType), opc::DualOpSum{P}) = OpSum(P, map((k,v)->f(k', v'), data(opc)))

Base.map(f::Union(Function,DataType), op::OuterProduct) = map(f, OpSum(op))

br_tup(tup) = (tup[1], tup[2]')

export maplabels!,
    mapcoeffs!,
    maplabels,
    mapcoeffs