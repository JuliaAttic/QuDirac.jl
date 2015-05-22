#############
# Filtering #
#############
Base.filter!(f::Function, kt::KetSum) = (filter!(f, data(kt)); return kt)
Base.filter!(f::Function, br::BraSum) = (filter!((k,v)->f(k,v'), data(br)); return br)
Base.filter!(f::Function, op::OpSum) = (filter!(f, data(op)); return op)
Base.filter!(f::Function, opc::DualOpSum) = (filter!((k,v)->f(k',v'), data(opc)); return opc)

Base.filter{P}(f::Function, kt::Ket{P}) = make_kt(P, filter(f, data(kt)))
Base.filter{P}(f::Function, op::OpSum{P}) = OpSum(P, filter(f, data(op)))

Base.filter(f::Function, br::Bra) = ctranspose(filter((k,v)->f(k,v'), br'))
Base.filter(f::Function, opc::DualOpSum) = ctranspose(filter((k,v)->f(k',v'), opc'))

Base.filter(f::Function, op::OuterProduct) = filter(f, OpSum(op))

########################
# mapcoeffs!/mapcoeffs #
########################
mapcoeffs!(f::Union(Function,DataType), kt::KetSum) = (mapvals!(f, data(kt)); return kt)
mapcoeffs!(f::Union(Function,DataType), op::OpSum) = (mapvals!(f, data(op)); return op)

mapcoeffs!(f::Union(Function,DataType), br::BraSum) = (mapvals!(v->f(v')', data(br)); return br)
mapcoeffs!(f::Union(Function,DataType), opc::DualOpSum) = (mapvals!(v->f(v')', data(opc)); return opc)

mapcoeffs{P}(f::Union(Function,DataType), kt::Ket{P}) = make_kt(P, mapvals(f, data(kt)))
mapcoeffs{P}(f::Union(Function,DataType), op::OpSum{P}) = OpSum(P, mapvals(f, data(op)))

mapcoeffs(f::Union(Function,DataType), br::Bra) = ctranspose(mapvals(v->f(v')', br'))
mapcoeffs(f::Union(Function,DataType), opc::DualOpSum) = ctranspose(mapvals(v->f(v')', opc'))

mapcoeffs(f::Union(Function,DataType), op::OuterProduct) = mapcoeffs(f, OpSum(op))

########################
# maplabels!/maplabels #
########################
maplabels{P}(f::Union(Function,DataType), kt::Ket{P}) = make_kt(P, mapkeys(f, data(kt)))
maplabels{P}(f::Union(Function,DataType), op::OpSum{P}) = OpSum(P, mapkeys(f, data(op)))

maplabels(f::Union(Function,DataType), br::Bra) = ctranspose(maplabels(f, br'))
maplabels(f::Union(Function,DataType), opc::DualOpSum) = ctranspose(maplabels(l->f(l')', opc'))

maplabels(f::Union(Function,DataType), op::OuterProduct) = maplabels(f, OpSum(op))

#######
# map #
#######
Base.map{P}(f::Union(Function,DataType), kt::Ket{P}) = make_kt(P, map(f, data(kt)))
Base.map{P}(f::Union(Function,DataType), op::OpSum{P}) = OpSum(P, map(f, data(op)))

Base.map(f::Union(Function,DataType), br::Bra) = ctranspose(map((k,v)->br_tup(f(k,v')), br'))
Base.map{P}(f::Union(Function,DataType), opc::DualOpSum{P}) = OpSum(P, map((k,v)->f(k', v'), data(opc)))

Base.map(f::Union(Function,DataType), op::OuterProduct) = map(f, OpSum(op))

br_tup(tup) = (tup[1], tup[2]')

export maplabels!,
    mapcoeffs!,
    maplabels,
    mapcoeffs