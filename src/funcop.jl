abstract FuncOp{P} <: DiracOp{P}
abstract FuncOpTerm{P} <: FuncOp{P}

Base.one(::FuncOp) = 1
Base.zero(::FuncOp) = 0

apply_op_def() = error("`apply_op_def` takes arguments")

#############
# FuncOpDef #
#############
abstract FuncOpDef{P} <: FuncOpTerm{P}

execute_inner(op::FuncOpDef, kt::Ket) = apply_op_def(op, kt)
execute_inner(br::Bra, op::FuncOpDef) = apply_op_def(op, br)

Base.show(io::IO, op::FuncOpDef) = print(io, split(string(typeof(op)), '_')[1])

#############
# FuncOpExp #
#############
immutable FuncOpExp{P,F<:FuncOpDef} <: FuncOpTerm{P}
    op::F
    n::Int
    FuncOpExp(op::FuncOpDef{P}, n::Int) = new(op, n)
end

FuncOpExp{P}(op::FuncOpDef{P}, n) = FuncOpExp{P,typeof(op)}(op, n)

nfactors(ope::FuncOpExp) = nfactors(ope.op)

function apply_op_exp(ope::FuncOpExp, state)
    result = state
    for i=1:ope.n
        result = ope.op * result
    end
    return result
end

execute_inner(ope::FuncOpExp, kt::Ket) = apply_op_exp(ope, kt)
execute_inner(br::Bra, ope::FuncOpExp) = apply_op_exp(ope, br)

execute_inner{F<:FuncOpDef}(a::FuncOpExp{F}, b::FuncOpExp{F}) = FuncOpExp(a.op, a.n + b.n)
execute_inner{F<:FuncOpDef}(op::F, ope::FuncOpExp{F}) = FuncOpExp(ope.op, ope.n + 1)
execute_inner{F<:FuncOpDef}(ope::FuncOpExp{F}, op::F) = inner(op,ope)
execute_inner{F<:FuncOpDef}(a::F, b::F) = FuncOpExp(a,2)

Base.(:^)(op::FuncOpDef, n::Int) = FuncOpExp(op,n)
Base.(:^)(ope::FuncOpExp, n::Int) = FuncOpExp(ope.op, ope.n * n)

Base.show(io::IO, ope::FuncOpExp) = print(io, "$(ope.op)^$(ope.n)")

# ###############
# # FuncOpChain #
# ###############
# immutable FuncOpChain{P} <: FuncOpTerm{P}
#     ops::Vector{FuncOp{P}} # This is wrong
# end

# function apply_op_chain(opch::FuncOpChain, kt::Ket)
#     result = kt
#     ops = opch.ops
#     for i=1:length(ops)
#         result = ops[i] * result
#     end
#     return result
# end

# function apply_op_chain(opch::FuncOpChain, br::Bra)
#     result = br
#     ops = opch.ops
#     for i=reverse(1:length(ops))
#         result = result * ops[i]
#     end
#     return result
# end

# execute_inner(a::FuncOpChain, b::FuncOpChain) = FuncOpChain(vcat(a.ops, b.ops))
# execute_inner(opch::FuncOpChain, op::FuncOp) = FuncOpChain(vcat(opch.ops, op))
# execute_inner(op::FuncOp, opch::FuncOpChain) = FuncOpChain(vcat(op, opch.ops))
# execute_inner(a::FuncOp, b::FuncOp) = FuncOpChain(vcat(a, b))

# execute_inner(opch::FuncOpChain, kt::Ket) = apply_op_chain(opch, kt)
# execute_inner(br::Bra, opch::FuncOpChain) = apply_op_chain(opch, br)

# Base.show(io::IO, op::FuncOpChain) = print(io, "("*join(op.ops, "*")*")")

################
# ScaledFuncOp #
################
immutable ScaledFuncOp{P,F<:FuncOp, N<:Number} <: FuncOp{P}
    term::SumTerm{F,N}
    ScaledFuncOp(op::FuncOp{P}, c::N) = new(SumTerm(op, c))
end

ScaledFuncOp{P,N}(op::FuncOp{P}, c::N) = ScaledFuncOp{P,typeof(op), N}(SumTerm(op, c))

coeff(sop::ScaledFuncOp) = val(sop.term)
oper(sop::ScaledFuncOp) = key(sop.term)

nfactors(sop::ScaledFuncOp) = nfactors(oper(sop))

execute_inner(sop::ScaledFuncOp, kt::Ket) = coeff(sop) * inner(oper(sop), kt)
execute_inner(br::Bra, sop::ScaledFuncOp) = coeff(sop) * inner(br, oper(sop))

execute_inner(a::ScaledFuncOp, b::ScaledFuncOp) = ScaledFuncOp(inner(oper(a), oper(b)), coeff(a) * coeff(b))

# generate redundant code to resolve ambiguity warnings
# for T in (:FuncOpChain, :FuncOp)
for T in (:FuncOp,)
    @eval begin
        execute_inner(op::($T), sop::ScaledFuncOp) = ScaledFuncOp(inner(op, oper(sop)), coeff(sop))
        execute_inner(sop::ScaledFuncOp, op::($T)) = ScaledFuncOp(inner(oper(sop), op), coeff(sop))
    end
end

Base.scale(op::FuncOp, c::Number) = ScaledFuncOp(op, c)
Base.scale(c::Number, op::FuncOp) = scale(op, c) 
Base.scale(sop::ScaledFuncOp, c::Number) = ScaledFuncOp(oper(sop), c * coeff(sop))
Base.scale(c::Number, sop::ScaledFuncOp) = scale(sop, c) 

Base.ctranspose(sop::ScaledFuncOp) = ScaledFuncOp(oper(sop)', coeff(sop)')

Base.show(io::IO, op::ScaledFuncOp) = print(io, "($(coeff(op)) * $(oper(op)))")

# #############
# # FuncOpSum #
# #############
# type FuncOpSum{F<:FuncOpTerm, N<:Number} <: FuncOp
#     data::SumDict{F,N}    
# end

##############
# DualFuncOp #
##############
immutable DualFuncOp{P,F<:FuncOp} <: FuncOp{P}
    op::F
    DualFuncOp(op::FuncOp{P}) = new(op)
end

DualFuncOp{P}(op::FuncOp{P}) = DualFuncOp{P,typeof(op)}(op)

nfactors(opc::DualFuncOp) = nfactors(opc.op)

Base.ctranspose(op::FuncOp) = DualFuncOp(op)
Base.ctranspose(opc::DualFuncOp) = opc.op

execute_inner(a::DualFuncOp, b::DualFuncOp) = inner(b',a')'
execute_inner(opc::DualFuncOp, kt::Ket) = inner(kt',opc')'
execute_inner(br::Bra, opc::DualFuncOp) = inner(opc',br')'

###########
# @def_op #
###########
macro def_op(str)
    result_expr = def_op_expr(OpDefStr(str))
    return esc(result_expr)
end 

function def_op_expr(ods::OpDefStr)
    odex = OpDefExpr(ods)

    lhs_type = odex.lhs_type
    single_lhs_type = symbol("Single"*string(lhs_type))

    name = odex.op_sym
    T = symbol(string(name) * "Def")
    T_sym = Expr(:quote, T)

    on_label = symbol("_" * string(T) * "_on_" * string(lhs_type) * "_label")
    on_pair = symbol("_" * string(T) * "_on_" * string(lhs_type) * "_pair")

    coeff_sym = lhs_type == :Ket ? :c : :(c')
    
    if isa(odex.label_args, Expr)
        len_args = length(odex.label_args.args)
        on_label_args_def = quote 
            $on_label($(odex.label_args.args...)) = $(odex.rhs)
        end
    else # isa(odex.label_args, Symbol)
        len_args = 1
        on_label_args_def = quote 
            $on_label($(odex.label_args)) = $(odex.rhs)
        end
    end

    result = quote
        
        local ProdType = innertype(ket())

        if isdefined($T_sym)
            @assert QuDirac.nfactors($name) == $len_args "Cannot redefine "*string($T)*" on states with different number of factors"
            @assert $T <: QuDirac.FuncOpDef{ProdType} "Cannot redefine functional operator type "*string($T)*" with new inner product type $(ProdType)"
        else
            immutable $T <: QuDirac.FuncOpDef{ProdType} end
            QuDirac.nfactors(::$T) = $len_args
            const $name = $T()
        end

        $on_label_args_def

        $on_label(label::StateLabel) = $on_label(label...)

        function $on_pair(pair::Tuple)
          label, c = pair # c and coeff_sym are RELATED; see coeff_sym def above
          return $(coeff_sym) * $on_label(label) 
        end 

        QuDirac.apply_op_def(::$T, state::$(lhs_type)) = isempty(state) ? state : sum($on_pair, QuDirac.data(state))
        QuDirac.apply_op_def(::$T, state::QuDirac.$(single_lhs_type)) = QuDirac.coeff(state) * $on_label(QuDirac.label(state))
    
    end

    return result
end

export represent, @def_op