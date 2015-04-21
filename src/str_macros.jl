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
    return esc(parse(prune_dirac(str)))
end

macro d_mstr(str)
    return esc(quote
        local s;
        for s in filter(x -> ! isempty(x), split(strip($str), '\n'))
            eval(parse(QuDirac.prune_dirac(s)))
        end
    end)
end

############
# @repr_op #
############

# Accepted grammar is:
# @repr_op " OP_NAME | LABEL_ARG > = KET " basis_iterable

macro repr_op(str, basis)
    op_name, label_arg, ket_str = parse_defstr(str)
    result_expr = op_expr(op_name, label_arg, ket_str, basis)
    return esc(result_expr)
end 

rm_whspace(str) = join(split(str, r"\s"))

function sep_leftside(str)
    op_name, label_arg = @compat split(str, '|'; limit=2)
    if last(label_arg) == '>'
        return op_name, label_arg[1:end-1]
    else
        error("Left side of definition string is malformed; couldn't find terminating '>' for the Ket.")
    end
end

function parse_defstr(str)
    str = rm_whspace(str)
    left, right = @compat split(str, '='; limit=2)
    op_name, label_arg = sep_leftside(left)
    return (op_name, label_arg, prune_dirac(right))
end

function op_expr(op_name, label_arg, ket_str, basis)
    op_sym = symbol(op_name)
    label_expr = parse(label_arg)
    ket_expr = parse(ket_str)
    if isa(label_expr, Expr)
        ex = quote 
            local f($(label_expr.args...)) = $ket_expr * bra($(label_expr.args...))
            $op_sym = sum(args->f(args...), $basis)
        end
    else
        ex = quote 
            $op_sym = sum($label_expr -> $ket_expr * bra($label_expr), $basis)
        end
    end
    return ex
end

############
# @def_op #
############

# Accepted grammar is:
# @def_op " OP_NAME | LABEL_ARG > = KET "

macro def_op(str)
    op_name, label_arg, ket_str = parse_defstr(str)
    result_expr = func_op_expr(op_name, label_arg, ket_str)
    return esc(result_expr)
end 

function func_op_expr(op_name, label_arg, ket_str)
    op_sym = symbol(op_name)
    label_expr = parse(label_arg)
    ket_expr = parse(ket_str)
    if isa(label_expr, Expr)
        ex = quote 
            $(op_sym)($(label_expr.args...)) = $ket_expr
            $(op_sym)(s::StateLabel) = $(op_sym)(s...)
        end
    else
        ex = quote 
            $(op_sym)($label_expr) = $ket_expr
            $(op_sym)(s::StateLabel) = $(op_sym)(s...)
        end
    end
    return ex
end

export @d_str,
       @d_mstr,
       @repr_op,
       @def_op