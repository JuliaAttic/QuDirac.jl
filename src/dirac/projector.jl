#############
# Projector #
#############
    type Projector{S} <: AbstractOperator{S}
        ket::Ket{S}
        bra::Bra{S}
    end
