#######################
# Outer Product Types #
#######################

# Abstract Types #
#----------------#
abstract OuterProduct{P,K,B,T} <: DiracOperator{P}

# Outer Product Types #
#---------------------#
immutable BasisOuterProduct{P,K,B,T} <: OuterProduct{P,K,B,T}
    data::LabelTerm{OuterLabel{K,B},T}
end

immutable OuterProductSum{P,K,B,T} <: OuterProduct{P,K,B,T}
    data::LabelSum{OuterLabel{K,B},T}
end