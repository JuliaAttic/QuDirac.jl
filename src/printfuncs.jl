######################
# Printing Functions #
######################
function dirac_show(io::IO, dirac::AbstractDirac)
    print(io, summary(dirac)*":")
    pad = "  "
    maxlen = 10
    i = 1
    for j in keys(data(dirac))
        if i > maxlen
            break
        else
            i += 1
            println(io)
            print(io, labelrepr(dirac, j, pad))
        end
    end
    if length(dirac) > maxlen
        println(io)
        print(io, "$pad$vdots")
    end
end

function dirac_showcompact(io::IO, dirac::AbstractDirac)
    pad = " "
    maxlen = 4
    i = 1
    for j in keys(data(dirac))
        if i > maxlen
            break
        else
            print(io, labelrepr(dirac, j, pad))
            i += 1
            if i <= length(dirac)
                print(io, "$pad+")
            end
        end
    end
    if length(dirac) > maxlen
        print(io, "$pad...")
    end
end

function dirac_repr(dirac::AbstractDirac)
    tempio = IOBuffer()
    showcompact(tempio, dirac)
    return takebuf_string(tempio)
end

labelstr(s::StateLabel) = join(map(repr, s.label), ',')
ktstr(sl) = "| $(labelstr(sl)) $rang"
brstr(sl) = "$lang $(labelstr(sl)) |"

labelrepr(kt::Ket, sl, pad) = "$pad$(kt[sl]) $(ktstr(sl))"
labelrepr(br::Bra, sl, pad) = "$pad$(br[sl]) $(brstr(sl))"
labelrepr(op::OuterProduct, k, b, pad) = "$pad$(op[k,b]) $(ktstr(k))$(brstr(b))"
labelrepr(op::OuterSum, o::OpLabel, pad) = "$pad$(op[o]) $(ktstr(klabel(o)))$(brstr(blabel(o)))"
labelrepr(opc::DualOuterSum, o::OpLabel, pad) = "$pad$(opc[o']) $(ktstr(blabel(o)))$(brstr(klabel(o)))"

Base.repr(s::StateLabel) = repr(typeof(s)) * "(" * labelstr(s) * ")"
Base.repr(o::OpLabel) = repr(typeof(o)) * "(" * ktstr(o.k) * "," * brstr(o.b) * ")"
Base.repr(s::DiracState) = dirac_repr(s)
Base.repr(op::AbsOuterSum) = dirac_repr(op)

Base.summary(s::DiracState) = "$(typeof(s)) with $(length(s)) state(s)"
Base.summary{P,N}(op::OuterProduct{P,N}) = "OuterProduct{$P,$N,$(eltype(op))} with $(length(op)) operator(s)"
Base.summary(op::DiracOp) = "$(typeof(op)) with $(length(op)) operator(s)"

Base.show(io::IO, s::StateLabel) = print(io, repr(s))
Base.show(io::IO, o::OpLabel) = print(io, repr(o))
Base.show(io::IO, s::DiracState) = dirac_show(io, s)
Base.show(io::IO, op::AbsOuterSum) = dirac_show(io, op)
function Base.show(io::IO, op::OuterProduct)
    print(io, summary(op)*":")
    pad = "  "
    maxlen = 10
    i = 1
    for k in keys(data(op.kt)), b in keys(data(op.br))
        if i > maxlen
            break
        else
            i += 1
            println(io)
            print(io, labelrepr(op, k, b, pad))
        end
    end
    if length(op) > maxlen
        println(io)
        print(io, "$pad$vdots")
    end
end

Base.showcompact(io::IO, s::DiracState) = dirac_showcompact(io, s)
Base.showcompact(io::IO, op::AbsOuterSum) = dirac_showcompact(io, op)

# Placing here for now - may be useful for
# subscripting label factors in the future
digit_subscript(i::Int) = Base.REPLCompletions.latex_completions("\\_$i",3)[2][1][1]
subscript(i::Int) = reduce(*, reverse(map(digit_subscript, digits(i))))
