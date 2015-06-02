##########
# @d_str #
##########
const ktpat = r"\|.*?\>"
const brpat = r"\<.*?\|"
const inpat = r"\<[^\|]*\|[^\|]*\>"

ktrep(str) = "ket("*str[2:end-1]*")"
brrep(str) = "bra("*str[2:end-1]*")"

function inrep(str)
    i = split(str, '|')
    return "(bra("i[1][2:end]")*ket("*i[2][1:end-1]*"))"
end

function prune_dirac(str)
    return replace(replace(replace(str, inpat, inrep), brpat, brrep), ktpat, ktrep)
end

macro d_str(str)
    return esc(parse(prune_dirac(replace(strip(str), '\n', ';'))))
end

###################################
# Definition String Parsing Utils #
###################################
immutable OpDefStr
    op_name::AbstractString
    label_args::AbstractString
    lhs_type::AbstractString
    rhs::AbstractString
end

function OpDefStr(str::AbstractString)
    str = rm_whspace(str)
    left, right = @compat(split(str, '=', limit=2))

    if ismatch(ktpat, left)
        lhs_type = "Ket"
    elseif ismatch(brpat, left)
        lhs_type = "Bra"
    else
        error("Couldn't detect whether operator is acting on Bra or Ket...")
    end

    if lhs_type == "Ket"
        op_name, label_args = leftside_ket(left)
    else # type_str == "Bra"
        op_name, label_args = leftside_bra(left)
    end

    return OpDefStr(op_name, label_args, lhs_type, prune_dirac(right))
end

immutable OpDefExpr
    op_sym::Symbol
    label_args::Union(Symbol, Expr)
    lhs_type::Symbol
    rhs::Union(Symbol, Expr)
end

OpDefExpr(ods::OpDefStr) = OpDefExpr(symbol(ods.op_name), parse(ods.label_args), symbol(ods.lhs_type), parse(ods.rhs))

rm_whspace(str) = join(split(str, r"\s"))

function leftside_ket(str)
    op_name, label_args = @compat(split(str, '|', limit=2))
    if last(label_args) == '>'
        return op_name, label_args[1:end-1]
    else
        error("Left side of definition string is malformed; couldn't find terminating '>' for the Ket.")
    end
end

function leftside_bra(str)
    label_args, op_name = @compat(split(str, '|', limit=2))
    if first(label_args) == '<'
        return op_name, label_args[2:end]
    else
        error("Left side of definition string is malformed; couldn't find initial '<' for the Bra.")
    end
end

###########
# @rep_op #
###########
macro rep_op(str, bases...)
    result_expr = rep_op_expr(OpDefStr(str), build_prod_basis(bases...))
    return esc(result_expr)
end 

build_prod_basis(basis) = basis
build_prod_basis(first, bases...) = :(Iterators.product($first, $(bases...)))

function rep_op_expr(ods::OpDefStr, basis)
    odex = OpDefExpr(ods)
    if odex.lhs_type == :Ket
        return gen_ket_repr_expr(odex, basis)
    else
        return gen_bra_repr_expr(odex, basis)
    end
    return ex
end

function gen_ket_repr_expr(odex::OpDefExpr, basis)
    if isa(odex.label_args, Expr)
        ex = quote 
            local func_on_args($(odex.label_args.args...)) = $(odex.rhs) * bra($(odex.label_args.args...))
            $(odex.op_sym) = sum(args->func_on_args(args...), $basis)
        end
    else
        ex = quote 
            $(odex.op_sym) = sum($(odex.label_args) -> $(odex.rhs) * bra($(odex.label_args)), $basis)
        end
    end
    return ex
end

function gen_bra_repr_expr(odex::OpDefExpr, basis)
    if isa(odex.label_args, Expr)
        ex = quote 
            local func_on_args($(odex.label_args.args...)) = ket($(odex.label_args.args...)) * $(odex.rhs)
            $(odex.op_sym) = sum(args->func_on_args(args...), $basis)
        end
    else
        ex = quote 
            $(odex.op_sym) = sum($(odex.label_args) -> ket($(odex.label_args)) * $(odex.rhs), $basis)
        end
    end
    return ex
end

export @d_str,
       @rep_op