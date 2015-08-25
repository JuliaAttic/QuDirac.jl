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

function replace_dirac(str)
    return replace(replace(replace(str, inpat, inrep), brpat, brrep), ktpat, ktrep)
end

function prune_dstr(str)
    return string("(", replace_dirac(replace(strip(str), r"(\n|\r)", "; ")), ")")
end

macro d_str(str)
    result = prune_dstr(str)
    return esc(parse(result))
end

# if v"0.3-" <= VERSION < v"0.4-":
macro d_mstr(str)
    result = prune_dstr(str)
    return esc(parse(result))
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

    return OpDefStr(op_name, label_args, lhs_type, replace_dirac(right))
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

export @d_str, @d_mstr