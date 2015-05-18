abstract FuncOp <: DiracOp
abstract FuncOpDef <: FuncOp

Base.one(::FuncOp) = 1
Base.zero(::FuncOp) = 0

apply_op_def() = error("apply_op_def requires arguments")

inner(op::FuncOpDef, kt::Ket) = apply_op_def(op, kt)
inner(br::Bra, op::FuncOpDef) = apply_op_def(op, br)

Base.copy(op::FuncOpDef) = op

##############
# DualFuncOp #
##############
immutable DualFuncOp{T<:FuncOp} <: FuncOp
    op::T
end

Base.ctranspose(op::FuncOp) = DualFuncOp(op)
Base.ctranspose(f::DualFuncOp) = f.op

inner(a::DualFuncOp, b::DualFuncOp) = inner(b',a')'
inner(opc::DualFuncOp, kt::Ket) = inner(kt',opc')'
inner(br::Bra, opc::DualFuncOp) = inner(opc',br')'

#############
# FuncOpExp #
#############
immutable FuncOpExp{T<:FuncOpDef} <: FuncOp
    op::T
    n::Int
end

function apply_op_exp{T<:FuncOpDef}(ope::FuncOpExp{T}, state)
    result = state
    for i=1:f.n
        result = apply_op_def(ope.op, result)
    end
    return result
end

inner{T<:FuncOpDef}(ope::FuncOpExp{T}, kt::Ket) = apply_op_exp(ope, kt)
inner{T<:FuncOpDef}(br::Bra, ope::FuncOpExp{T}) = apply_op_exp(ope, br)

inner{T<:FuncOpDef}(a::FuncOpExp{T}, b::FuncOpExp{T}) = FuncOpExp(a.op, a.n + b.n)
inner{T<:FuncOpDef}(op::T, ope::FuncOpExp{T}) = FuncOpExp(ope.op, ope.n + 1)
inner{T<:FuncOpDef}(ope::FuncOpExp{T}, op::T) = inner(op,ope)
inner{T<:FuncOpDef}(a::T, b::T) = FuncOpExp(a,2)

Base.(:^)(op::FuncOpDef, n::Int) = FuncOpExp(op,n)
Base.(:^)(ope::FuncOpExp, n::Int) = FuncOpExp(ope.op, ope.n * n)

###############
# FuncOpChain #
###############
immutable FuncOpChain{T<:FuncOp} <: FuncOp
    ops::Vector{T}
end

function apply_op_chain(ops::FuncOpChain, kt::Ket)
    result = kt
    for i=1:length(f.ops)
        result = f.ops[i] * result
    end
    return result
end

function apply_op_chain(ops::FuncOpChain, br::Bra)
    result = br
    for i=reverse(1:length(f.ops))
        result = result * f.ops[i]
    end
    return result
end

inner(a::FuncOpChain, b::FuncOpChain) = FuncOpChain(vcat(a.ops, b.ops))
inner(opch::FuncOpChain, op::FuncOp) = FuncOpChain(vcat(opch.ops, op))
inner(op::FuncOp, opch::FuncOpChain) = FuncOpChain(vcat(op, opch.ops))
inner(a::FuncOp, b::FuncOp) = FuncOpChain(vcat(a, b))

inner(opch::FuncOpChain, kt::Ket) = apply_op_chain(opch, kt)
inner(br::Bra, opch::FuncOpChain) = apply_op_chain(opch, br)

function represent(op::FuncOp, basis)
    return [bra(i) * op * ket(j) for i in basis, j in basis]
end

function represent(op::FuncOp, basis...)
    prodbasis = product(basis...)
    return [bra(i...) * op * ket(j...) for i in prodbasis, j in prodbasis]
end

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
    T = symbol(string(name) * "FuncOpDef")
    T_sym = Expr(:quote, T)

    if isa(odex.label_args, Expr)
        label_args = {odex.label_args.args...}
    else # isa(odex.label_args, Symbol)
        label_args = {odex.label_args}
    end

    args_ex = parse(join(label_args,","))
    len_args = length(label_args)

    on_label = symbol("_" * string(T) * "_on_" * string(lhs_type) * "_label")
    on_pair = symbol("_" * string(T) * "_on_" * string(lhs_type) * "_pair")

    on_label_def = quote 
        function $on_label(label::StateLabel{$len_args})
            $args_ex = label
            $(odex.rhs)
        end
    end

    coeff_sym = lhs_type == :Ket ? :c : :(c')

    result = quote
        if ! isdefined($T_sym)
            immutable $T <: QuDirac.FuncOpDef end
        end

        $on_label_def

        function $on_pair(pair::Tuple)
          label, c = pair # c and coeff_sym are RELATED; see coeff_sym def above
          return $(coeff_sym) * $on_label(label) 
        end  

        QuDirac.apply_op_def(::$T, state::QuDirac.$(single_lhs_type)) = QuDirac.coeff(state) * $on_label(QuDirac.label(state))
        QuDirac.apply_op_def(::$T, state::$(lhs_type)) = sum($on_pair, QuDirac.data(state))

        const $name = $T()
    end

    return result
end

export represent, @def_op