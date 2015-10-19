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

#################
# FunctionLabel #
#################
immutable FunctionLabel{P<:InnerProduct,L} <: AbstractLabel
    onket::Function # onket(::StateLabel) -> AbstractKet{P,L}
    onbra::Function # onbra(::StateLabel) -> AbstractKet{P,L}
end

@generated function tensor{P,A,B}(a::FunctionLabel{P,A}, b::FunctionLabel{P,B})
    aparams, bparams = A.parameters, B.parameters
    L = Tuple{aparams..., bparams...}
    i = length(aparams)
    quote
        function onket{$L}(lbl::StateLabel{$L})
            albl, blbl = split(lbl, Val{$i})
            akt = a.onket(albl)
            bkt = b.onket(blbl)
            return tensor(akt, bkt)
        end

        function onbra{$L}(lbl::StateLabel{$L})
            albl, blbl = split(lbl, Val{$i})
            abr = a.onbra(albl)
            bbr = b.onbra(blbl)
            return tensor(abr, bbr)
        end

        return FunctionLabel{P,$L}(onket, onbra)
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
    data::LabelTerm{FunctionLabel{L},T}
end

immutable OperatorSum{P,L,T} <: AbstractOperator{P,L,T}
    data::LabelSum{FunctionLabel{L},T}
end
