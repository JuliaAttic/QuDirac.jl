abstract FuncOp <: DiracOp
abstract FuncOpTerm <: FuncOp

Base.one(::FuncOp) = 1
Base.zero(::FuncOp) = 0

#############
# FuncOpDef #
#############
abstract FuncOpDef <: FuncOpTerm

inner{F<:FuncOpDef}(::F, kt::Ket) = F(kt)
inner{F<:FuncOpDef}(br::Bra, ::F) = F(br)

Base.show(io::IO, op::FuncOpDef) = print(io, split(string(typeof(op)), '_')[1])

#############
# FuncOpExp #
#############
immutable FuncOpExp{F<:FuncOpDef} <: FuncOpTerm
    op::F
    n::Int
end

function apply_op_exp{F<:FuncOpDef}(ope::FuncOpExp{F}, state)
    result = state
    for i=1:ope.n
        result = F(result)
    end
    return result
end

inner(ope::FuncOpExp, kt::Ket) = apply_op_exp(ope, kt)
inner(br::Bra, ope::FuncOpExp) = apply_op_exp(ope, br)

inner{F<:FuncOpDef}(a::FuncOpExp{F}, b::FuncOpExp{F}) = FuncOpExp(a.op, a.n + b.n)
inner{F<:FuncOpDef}(op::F, ope::FuncOpExp{F}) = FuncOpExp(ope.op, ope.n + 1)
inner{F<:FuncOpDef}(ope::FuncOpExp{F}, op::F) = inner(op,ope)
inner{F<:FuncOpDef}(a::F, b::F) = FuncOpExp(a,2)

Base.(:^)(op::FuncOpDef, n::Int) = FuncOpExp(op,n)
Base.(:^)(ope::FuncOpExp, n::Int) = FuncOpExp(ope.op, ope.n * n)

Base.show(io::IO, ope::FuncOpExp) = print(io, "$(ope.op)^$(ope.n)")

###############
# FuncOpChain #
###############
immutable FuncOpChain{F<:FuncOpTerm} <: FuncOpTerm
    ops::Vector{F}
end

function apply_op_chain(opch::FuncOpChain, kt::Ket)
    result = kt
    ops = opch.ops
    for i=1:length(ops)
        result = ops[i] * result
    end
    return result
end

function apply_op_chain(opch::FuncOpChain, br::Bra)
    result = br
    ops = opch.ops
    for i=reverse(1:length(ops))
        result = result * ops[i]
    end
    return result
end

inner(a::FuncOpChain, b::FuncOpChain) = FuncOpChain(vcat(a.ops, b.ops))
inner(opch::FuncOpChain, op::FuncOp) = FuncOpChain(vcat(opch.ops, op))
inner(op::FuncOp, opch::FuncOpChain) = FuncOpChain(vcat(op, opch.ops))
inner(a::FuncOp, b::FuncOp) = FuncOpChain(vcat(a, b))

inner(opch::FuncOpChain, kt::Ket) = apply_op_chain(opch, kt)
inner(br::Bra, opch::FuncOpChain) = apply_op_chain(opch, br)

Base.show(io::IO, op::FuncOpChain) = print(io, "("*join(op.ops, "*")*")")

################
# ScaledFuncOp #
################
immutable ScaledFuncOp{F<:FuncOp, N<:Number} <: FuncOp
    term::SumTerm{F,N}
end

ScaledFuncOp(op, c) = ScaledFuncOp(SumTerm(op, c))

coeff(sop::ScaledFuncOp) = val(sop.term)
oper(sop::ScaledFuncOp) = key(sop.term)

inner(sop::ScaledFuncOp, kt::Ket) = coeff(sop) * inner(oper(sop), kt)
inner(br::Bra, sop::ScaledFuncOp) = coeff(sop) * inner(br, oper(sop))

inner(a::ScaledFuncOp, b::ScaledFuncOp) = ScaledFuncOp(inner(oper(a), oper(b)), coeff(a) * coeff(b))

# generate redundant code to resolve ambiguity warnings
for T in (:FuncOpChain, :FuncOp)
    @eval begin
        inner(op::($T), sop::ScaledFuncOp) = ScaledFuncOp(inner(op, oper(sop)), coeff(sop))
        inner(sop::ScaledFuncOp, op::($T)) = ScaledFuncOp(inner(oper(sop), op), coeff(sop))
    end
end

Base.scale(op::FuncOp, c::Number) = ScaledFuncOp(op, c)
Base.scale(c::Number, op::FuncOp) = scale(op, c) 
Base.scale(sop::ScaledFuncOp, c::Number) = ScaledFuncOp(oper(sop), c * coeff(sop))
Base.scale(c::Number, sop::ScaledFuncOp) = scale(sop, c) 

Base.ctranspose(sop::ScaledFuncOp) = ScaledFuncOp(oper(sop)', coeff(sop)')

Base.show(io::IO, op::ScaledFuncOp) = print(io, "($(coeff(op)) * $(oper(op)))")

#############
# FuncOpSum #
#############
type FuncOpSum{F<:FuncOpTerm, N<:Number} <: FuncOp
    data::SumDict{F,N}    
end

##############
# DualFuncOp #
##############
immutable DualFuncOp{F<:FuncOp} <: FuncOp
    op::F
end

Base.ctranspose(op::FuncOp) = DualFuncOp(op)
Base.ctranspose(opc::DualFuncOp) = opc.op

inner(a::DualFuncOp, b::DualFuncOp) = inner(b',a')'
inner(opc::DualFuncOp, kt::Ket) = inner(kt',opc')'
inner(br::Bra, opc::DualFuncOp) = inner(opc',br')'

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

    if isa(odex.label_args, Expr)
        label_args = Any[odex.label_args.args...]
    else # isa(odex.label_args, Symbol)
        label_args = [odex.label_args]
    end

    args_ex = parse(join(label_args,","))
    len_args = length(label_args)

    on_label = symbol("_" * string(T) * "_on_" * string(lhs_type) * "_label")
    on_pair = symbol("_" * string(T) * "_on_" * string(lhs_type) * "_pair")

    coeff_sym = lhs_type == :Ket ? :c : :(c')

    result = quote
        if ! isdefined($T_sym)
            immutable $T <: QuDirac.FuncOpDef end
        end

        function $on_label(label::StateLabel{$len_args})
            $args_ex = label
            $(odex.rhs)
        end

        function $on_pair(pair::Tuple)
          label, c = pair # c and coeff_sym are RELATED; see coeff_sym def above
          return $(coeff_sym) * $on_label(label) 
        end 

        $T{P}(state::$(lhs_type){P,$len_args}) = isempty(state) ? state : sum($on_pair, QuDirac.data(state))
        $T{P}(state::QuDirac.$(single_lhs_type){P,$len_args}) = QuDirac.coeff(state) * $on_label(QuDirac.label(state))

        const $name = $T()
    end

    return result
end

export represent, @def_op