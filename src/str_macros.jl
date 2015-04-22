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

###########
# @rep_op #
###########

macro rep_op(str, basis)
    op_name, label_args, def_str = parse_defstr(str)
    result_expr = rep_op_expr(op_name, label_args, def_str, basis)
    return esc(result_expr)
end 

rm_whspace(str) = join(split(str, r"\s"))

function leftside_ket(str)
    op_name, label_args = @compat(split(str, '|'; limit=2))
    if last(label_args) == '>'
        return op_name, label_args[1:end-1]
    else
        error("Left side of definition string is malformed; couldn't find terminating '>' for the Ket.")
    end
end

function leftside_bra(str)
    label_args, op_name = @compat(split(str, '|'; limit=2))
    if first(label_args) == '<'
        return op_name, label_args[2:end]
    else
        error("Left side of definition string is malformed; couldn't find initial '<' for the Bra.")
    end
end

function parse_defstr(str)
    str = rm_whspace(str)
    left, right = @compat(split(str, '='; limit=2))

    if ismatch(ktpat, left)
        type_str = "Ket"
    elseif ismatch(brpat, left)
        type_str = "Bra"
    else
        error("Couldn't detect whether operator is acting on Bra or Ket...")
    end

    if type_str == "Ket"
        op_name, label_arg = leftside_ket(left)
    else # type_str == "Bra"
        op_name, label_arg = leftside_bra(left)
    end

    return (op_name, label_arg, prune_dirac(right), type_str)
end

function rep_op_expr(op_name, label_arg, ket_str, basis)
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

###########
# @def_op #
###########

macro def_op(str)
    op_name, label_args, def_str, type_str = parse_defstr(str)
    result_expr = def_op_expr(op_name, label_args, def_str, type_str)
    return esc(result_expr)
end 

function def_op_expr(op_name, label_args, def_str, type_str)
    op_sym = symbol(op_name)
    type_sym = symbol(type_str)
    label_expr = parse(label_args)
    def_expr = parse(def_str)

    func_label = symbol("_" * op_name * "_on_" * type_str * "_label")
    func_pair = symbol("_" * op_name * "_on_" * type_str * "_pair")

    if isa(label_expr, Expr)
        func_label_def = quote 
            $(func_label)($(label_expr.args...)) = $def_expr 
            $(func_label)(label::StateLabel) = $(func_label)(label...)
        end
    else # isa(label_expr, Symbol)
        func_label_def = quote 
            $(func_label)($label_expr) = $def_expr 
            $(func_label)(label::StateLabel) = $(func_label)(first(label))
        end
    end

    result = quote
        $func_label_def

        function $(func_pair)(pair::Tuple)
          label, c = pair
          return c * $(func_label)(label)
        end  

        $(op_sym)(state::$(type_sym)) = sum($(func_pair), QuDirac.dict(state))
    end

    return result
end

export @d_str,
       @d_mstr,
       @rep_op,
       @def_op