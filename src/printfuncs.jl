######################
# Printing Functions #
######################
function dirac_show(io::IO, dirac)
    print(io, summary(dirac)*":")
    pad = "  "
    maxlen = 10
    i = 1
    for label in keys(dict(dirac))
        if i > maxlen
            break
        else
            i += 1
            println(io)
            print(io, labelrepr(dirac, label, pad))
        end
    end
    if length(dirac) > maxlen
        println(io)
        print(io, "$pad$vdots")
    end
end

function dirac_showcompact(io::IO, dirac)
    pad = " "
    maxlen = 4
    i = 1
    for label in keys(dict(dirac))
        if i > maxlen
            break
        else
            print(io, labelrepr(dirac, label, pad))
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

function dirac_repr(dirac)
    tempio = IOBuffer()
    showcompact(tempio, dirac)
    return takebuf_string(tempio)
end

# Placing here for now - may be useful for
# subscripting label factors in the future
digit_subscript(i::Int) = Base.REPLCompletions.latex_completions("\\_$i",3)[2][1][1]
subscript(i::Int) = reduce(*, reverse(map(digit_subscript, digits(i))))
