abstract FuncOp <: DiracOp

############################
# Functional Operator Math #
############################
abstract DualFunc{T<:FuncOp}

Base.(:*){T<:FuncOp}(::Type{T}, kt::Ket) = T(kt)
Base.(:*){T<:FuncOp}(br::Bra, ::Type{T}) = T(br)
Base.(:*){T<:FuncOp}(::Type{DualFunc{T}}, kt::Ket) = T(kt')'
Base.(:*){T<:FuncOp}(br::Bra, ::Type{DualFunc{T}}) = T(br')'

Base.ctranspose{T<:FuncOp}(::Type{T}) = DualFunc{T}
Base.ctranspose{T<:FuncOp}(::Type{DualFunc{T}}) = T

function represent{T<:FuncOp}(op::Union(Type{DualFunc{T}}, Type{T}), basis)
    return [bra(i) * op * ket(j) for i in basis, j in basis]
end

function represent{T<:FuncOp}(op::Union(Type{DualFunc{T}}, Type{T}), basis...)
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

    T = odex.op_sym
    T_sym = Expr(:quote, T)

    on_args = symbol("_" * string(T) * "_on_" * string(lhs_type) * "_args")
    on_label = symbol("_" * string(T) * "_on_" * string(lhs_type) * "_label")
    on_pair = symbol("_" * string(T) * "_on_" * string(lhs_type) * "_pair")

    if isa(odex.label_args, Expr)
        on_label_def = quote 
            $on_args($(odex.label_args.args...)) = $(odex.rhs)
            $on_label(label::StateLabel) = $on_args(label...)
        end
    else # isa(label_expr, Symbol)
        on_label_def = quote 
            $on_args($(odex.label_args)) = $(odex.rhs)
            $on_label(label::StateLabel) = $on_args(first(label))
        end
    end

    coeff_sym = lhs_type == :Ket ? :c : :(c')

    result = quote
        if ! isdefined($T_sym)
            immutable $T <: QuDirac.FuncOp end
            
        end

        $on_label_def

        function $on_pair(pair::Tuple)
          label, c = pair # c and coeff_sym are RELATED; see coeff_sym def above
          return $(coeff_sym) * $on_label(label) 
        end  

        $T(state::QuDirac.$(single_lhs_type)) = QuDirac.coeff(state) * $on_label(QuDirac.label(state))
        $T(state::$(lhs_type)) = sum($on_pair, QuDirac.dict(state))
    end

    return result
end

export represent, @def_op