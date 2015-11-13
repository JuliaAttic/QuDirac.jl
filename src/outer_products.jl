#######################
# Outer Product Types #
#######################

# Abstract Types #
#----------------#
abstract AbstractOuterProduct{P,K,B,T} <: AbstractDirac{P}

# Outer Product Types #
#---------------------#
immutable BasisOuterProduct{P,K,B,T} <: AbstractOuterProduct{P,K,B,T}
    data::LabelTerm{OuterLabel{K,B},T}
end

BasisOuterProduct{P,K,B,T}(::Type{P}, data::LabelTerm{OuterLabel{K,B},T}) = BasisOuterProduct{P,K,B,T}(data)
BasisOuterProduct{P}(::Type{P}, label::StateLabel, coeff) = BasisOuterProduct(P, LabelTerm(label, coeff))

type OuterProductSum{P,K,B,T} <: AbstractOuterProduct{P,K,B,T}
    data::LabelSum{OuterLabel{K,B},T}
end

######################
# Property Functions #
######################
data(bop::BasisOuterProduct) = bop.data
data(ops::OuterProductSum) = ops.data

coeff(bop::BasisOuterProduct) = coeff(data(bop))
label(bop::BasisOuterProduct) = label(data(bop))

Base.eltype(ops::OuterProductSum) = eltype(typeof(ops))
Base.eltype{P,K,B,T}(::Type{OuterProductSum{P,K,B,T}}) = BasisOuterProduct{P,K,B,T}

coefftype(op::AbstractOuterProduct) = coefftype(typeof(op))
labeltype(op::AbstractOuterProduct) = labeltype(typeof(op))
nfactors(op::AbstractOuterProduct) = nfactors(typeof(op))

for S in (:AbstractOuterProduct, :BasisOuterProduct, :OuterProductSum)
    @eval begin
        innertype{P,K,B,T}(::Type{($S){P,K,B,T}}) = P
        labeltype{P,K,B,T}(::Type{($S){P,K,B,T}}) = (K, B)
        coefftype{P,K,B,T}(::Type{($S){P,K,B,T}}) = T
        nfactors{P,K,B,T}(::Type{($S){P,K,B,T}}) = nfactors(K)
    end
end
