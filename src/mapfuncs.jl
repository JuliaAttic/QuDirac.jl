#############
# Filtering #
#############
Base.filter!(f::Function, kt::KetSum) = (filter!(f, data(kt)); return kt)
Base.filter!(f::Function, br::BraSum) = (filter!((k,v)->f(k,v'), data(br)); return br)
Base.filter!(f::Function, op::OuterSum) = (filter!(f, data(op)); return op)
Base.filter!(f::Function, opc::DualOuterSum) = (filter!((k,v)->f(k',v'), data(opc)); return opc)

Base.filter{P}(f::Function, kt::Ket{P}) = make_kt(P, filter(f, data(kt)))
Base.filter{P}(f::Function, op::OuterSum{P}) = OuterSum(P, filter(f, data(op)))

Base.filter(f::Function, br::Bra) = ctranspose(filter((k,v)->f(k,v'), br'))
Base.filter(f::Function, opc::DualOuterSum) = ctranspose(filter((k,v)->f(k',v'), opc'))

Base.filter(f::Function, op::OuterProduct) = filter(f, OuterSum(op))

########################
# mapcoeffs!/mapcoeffs #
########################
mapcoeffs!(f::Union(Function,DataType), kt::KetSum) = (mapvals!(f, data(kt)); return kt)
mapcoeffs!(f::Union(Function,DataType), op::OuterSum) = (mapvals!(f, data(op)); return op)

mapcoeffs!(f::Union(Function,DataType), br::BraSum) = (mapvals!(v->f(v')', data(br)); return br)
mapcoeffs!(f::Union(Function,DataType), opc::DualOuterSum) = (mapvals!(v->f(v')', data(opc)); return opc)

mapcoeffs{P}(f::Union(Function,DataType), kt::Ket{P}) = make_kt(P, mapvals(f, data(kt)))
mapcoeffs{P}(f::Union(Function,DataType), op::OuterSum{P}) = OuterSum(P, mapvals(f, data(op)))

mapcoeffs(f::Union(Function,DataType), br::Bra) = ctranspose(mapvals(v->f(v')', br'))
mapcoeffs(f::Union(Function,DataType), opc::DualOuterSum) = ctranspose(mapvals(v->f(v')', opc'))

mapcoeffs(f::Union(Function,DataType), op::OuterProduct) = mapcoeffs(f, OuterSum(op))

########################
# maplabels!/maplabels #
########################
maplabels{P}(f::Union(Function,DataType), kt::Ket{P}) = make_kt(P, mapkeys(f, data(kt)))
maplabels{P}(f::Union(Function,DataType), op::OuterSum{P}) = OuterSum(P, mapkeys(f, data(op)))

maplabels(f::Union(Function,DataType), br::Bra) = ctranspose(maplabels(f, br'))
maplabels(f::Union(Function,DataType), opc::DualOuterSum) = ctranspose(maplabels(l->f(l')', opc'))

maplabels(f::Union(Function,DataType), op::OuterProduct) = maplabels(f, OuterSum(op))

#######
# map #
#######
Base.map{P}(f::Union(Function,DataType), kt::Ket{P}) = make_kt(P, map(f, data(kt)))
Base.map{P}(f::Union(Function,DataType), op::OuterSum{P}) = OuterSum(P, map(f, data(op)))

Base.map(f::Union(Function,DataType), br::Bra) = ctranspose(map((k,v)->br_tup(f(k,v')), br'))
Base.map{P}(f::Union(Function,DataType), opc::DualOuterSum{P}) = OuterSum(P, map((k,v)->f(k', v'), data(opc)))

Base.map(f::Union(Function,DataType), op::OuterProduct) = map(f, OuterSum(op))

br_tup(tup) = (tup[1], tup[2]')

export maplabels!,
    mapcoeffs!,
    maplabels,
    mapcoeffs