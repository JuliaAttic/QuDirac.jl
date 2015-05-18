abstract FuncOp <: DiracOp
abstract FuncOpDef <: FuncOp

Base.one{T<:FuncOp}(::Type{T}) = 1
Base.one(::FuncOp) = 1
Base.zero{T<:FuncOp}(::Type{T}) = 0
Base.zero(::FuncOp) = 0

Base.(:*){T<:FuncOpDef}(::Union(Type{T},T), kt::Ket) = T(kt)
Base.(:*){T<:FuncOpDef}(br::Bra, ::Union(Type{T},T)) = T(br)

############################
# Functional Operator Math #
############################
immutable DualFunc{T<:FuncOp} <: FuncOp
    op::T
end

Base.ctranspose{T<:FuncOpDef}(::Type{T}) = DualFunc(T())
Base.ctranspose(op::FuncOp) = DualFunc(op)
Base.ctranspose(f::DualFunc) = f.op

Base.(:*)(a::DualFunc, b::DualFunc) = (b' * a')'
Base.(:*)(f::DualFunc, kt::Ket) = (kt' * f')'
Base.(:*)(br::Bra, f::DualFunc) = (f' * br')'

immutable FuncOpExp{T<:FuncOpDef} <: FuncOp
    op::Type{T}
    n::Integer
end

apply_op{T<:FuncOpDef}(::Type{T}, state, n) = n == 0 ? state : apply_op(T, T(state), n-1)
Base.(:*){T<:FuncOpDef}(f::FuncOpExp{T}, kt::Ket) = apply_op(T, kt, f.n)
Base.(:*){T<:FuncOpDef}(br::Bra, f::FuncOpExp{T}) = apply_op(T, br, f.n)

Base.(:*){T<:FuncOpDef}(::Type{T}, ::Type{T}) = FuncOpExp(T, 2)
Base.(:*){T<:FuncOpDef}(::Type{T}, ::T) = T * T
Base.(:*){T<:FuncOpDef}(::T, ::Type{T}) = T * T
Base.(:*){T<:FuncOpDef}(::T, ::T) = T * T

Base.(:*){T<:FuncOpDef}(a::FuncOpExp{T}, b::FuncOpExp{T}) = FuncOpExp(T, a.n + b.n)
Base.(:*){T<:FuncOpDef}(::Type{T}, f::FuncOpExp{T}) = FuncOpExp(T, f.n+1)
Base.(:*){T<:FuncOpDef}(::T, f::FuncOpExp{T}) = T * f

Base.(:*){T<:FuncOpDef}(f::FuncOpExp{T}, ::Type{T}) = FuncOpExp(T, f.n+1)
Base.(:*){T<:FuncOpDef}(f::FuncOpExp{T}, ::T) = f * T

immutable FuncOpChain{T<:FuncOp} <: FuncOp
    ops::Vector{T}
end

apply_op(ops, kt::Ket) = length(ops) == 1 ? first(ops) * kt : apply_op(ops[1:end-1], last(ops) * kt)
apply_op(ops, br::Bra) = length(ops) == 1 ? br * last(ops) : apply_op(ops[2:end], br * first(ops))

Base.(:*)(a::FuncOpChain, b::FuncOpChain) = FuncOpChain(vcat(a.ops, b.ops))
Base.(:*)(a::FuncOpChain, b::FuncOp) = FuncOpChain(vcat(a.ops, b))
Base.(:*)(a::FuncOp, b::FuncOpChain) = FuncOpChain(vcat(a, b.ops))
Base.(:*)(a::FuncOp, b::FuncOp) = FuncOpChain(vcat(a, b))

Base.(:*){A<:FuncOpDef, B<:FuncOpDef}(::Type{A}, ::Type{B}) = A() * B()
Base.(:*){A<:FuncOpDef}(::Type{A}, f::FuncOp) = A() * f
Base.(:*){B<:FuncOpDef}(f::FuncOp, ::Type{B}) = f * B()

Base.(:*)(f::FuncOpChain, kt::Ket) = apply_op(f.ops, kt)
Base.(:*)(br::Bra, f::FuncOpChain) = apply_op(f.ops, br)

function represent{T<:FuncOpDef}(::Union(FuncOp, Type{T}), basis)
    return [bra(i) * T * ket(j) for i in basis, j in basis]
end

function represent{T<:FuncOpDef}(::Union(FuncOp, Type{T}), basis...)
    prodbasis = product(basis...)
    return [bra(i...) * T * ket(j...) for i in prodbasis, j in prodbasis]
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

        $T(state::QuDirac.$(single_lhs_type)) = QuDirac.coeff(state) * $on_label(QuDirac.label(state))
        $T(state::$(lhs_type)) = sum($on_pair, QuDirac.data(state))
    end

    return result
end

export represent, @def_op