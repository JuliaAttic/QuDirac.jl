# This file presents an idea for implementing
# Operators, which are essentially functions on
# States that can undergo operator arithmetic:
#
#   @defop " A | x,y,z > = √y * | x,y+1,z > "
#   @defop " < x,y,z | A = √(y-1) * < x,y-1,z | "
#
#   B = A'
#   C = 1/√2 * (A + B)
#   D = tensor(C, C)
#   result = D * | 1, 2, 3 >

###########
# OpLabel #
###########
immutable OpLabel{P<:InnerProduct,K,B} <: AbstractLabel
    onket # onket(::StateLabel) -> AbstractKet{P,K}
    onbra # onbra(::StateLabel) -> AbstractBra{P,B}
    function OpLabel(onket, onbra)
        @assert length(K.parameters) == length(B.parameters)
        return new(onket, onbra)
    end
end

@generated function tensor{P,K1,B1,K2,B2}(a::OpLabel{P,K1,B1}, b::OpLabel{P,K2,B2})
    K = Tuple{K1.parameters..., K2.parameters...}
    B = Tuple{B1.parameters..., B2.parameters...}
    i = length(K1.parameters)
    quote
        function onket{$K}(lbl::StateLabel{$K})
            albl, blbl = split(lbl, Val{$i})
            akt = a.onket(albl)
            bkt = b.onket(blbl)
            return tensor(akt, bkt)
        end

        function onbra{$B}(lbl::StateLabel{$B})
            albl, blbl = split(lbl, Val{$i})
            abr = a.onbra(albl)
            bbr = b.onbra(blbl)
            return tensor(abr, bbr)
        end

        return OpLabel{P,$K,$B}(onket, onbra)
    end
end

##################
# Operator Types #
##################

# Abstract Types #
#----------------#
abstract AbstractOperator{P,L,T} <: AbstractDirac{P}

# Operator Types #
#----------------#
immutable BasisOperator{P,L,T} <: AbstractOperator{P,L,T}
    data::LabelTerm{OpLabel{L},T}
end

type OperatorSum{P,L,T} <: AbstractOperator{P,L,T}
    data::LabelSum{OpLabel{L},T}
end
