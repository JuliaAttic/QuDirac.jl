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

immutable OuterProductSum{P,K,B,T} <: AbstractOuterProduct{P,K,B,T}
    data::LabelSum{OuterLabel{K,B},T}
end
