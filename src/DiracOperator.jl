# This file presents an idea for implementing
# DiracOperators, which are essentially functions on
# DiracStates that can undergo operator arithmetic:
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
immutable FunctionLabel{P,L} <: DiracLabel
    # func has signatures:
    #   func(::StateLabel{L}, ::Type{BasisKet{P,L}}) -> AbstractKet{P,L}
    #   func(::StateLabel{L}, ::Type{BasisBra{P,L}}) -> AbstractBra{P,L}
    func::Function
end

@generated function tensor{P,A,B}(a::FunctionLabel{P,A}, b::FunctionLabel{P,B})
    aparams, bparams = A.parameters, B.parameters
    L = Tuple{aparams..., bparams...}
    i = length(aparams)
    quote
        @inline function f{P,$L}(lbl::StateLabel{$L}, ::Type{BasisKet{P,$L}})
            albl, blbl = split(lbl, Val{$i})
            akt = a.func(albl, BasisKet{P,A})
            bkt = b.func(blbl, BasisKet{P,B})
            return tensor(akt, bkt)
        end

        @inline function f{P,$L}(lbl::StateLabel{$L}, ::Type{BasisBra{P,$L}})
            albl, blbl = split(lbl, Val{$i})
            abr = a.func(albl, BasisBra{P,A})
            bbr = b.func(blbl, BasisBra{P,B})
            return tensor(abr, bbr)
        end

        return FunctionLabel{P,$L}(f)
    end
end

#######################
# DiracOperator Types #
#######################

# Abstract Types #
#----------------#
abstract DiracOperator{P,L,T} <: AbstractDirac{P}

# Operator Types #
#----------------#
immutable BasisDiracOperator{P,L,T} <: DiracOperator{P,L,T}
    data::LabelTerm{FunctionLabel{L},T}
end

immutable DiracOperatorSum{P,L,T} <: DiracOperator{P,L,T}
    data::LabelSum{FunctionLabel{L},T}
end
