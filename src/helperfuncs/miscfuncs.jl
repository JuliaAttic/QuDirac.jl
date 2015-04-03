####################
# Helper Functions #
####################
isx(label::OpLabel, x) = sum(ktlabel(label))==sum(brlabel(label))==x
isx(k::Array, x) = sum(k) == x
ctpair(k,v) = (reverse(k), v')
nzcoeff(k,v) = v!=0
second(t) = t[2]
except(arr, i) = deleteat!(copy(arr), i)
function switch!(arr, i, j)
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp
    return arr
end
switch(arr, i, j) = switch!(copy(arr), i, j)
permute(arr, perm) = permute!(copy(arr), perm)
placeat(arr, x, y) = setindex!(copy(arr), x, y)

######################
# Printing Functions #
######################
function dirac_show(io::IO, dirac)
    print(io, summary(dirac)*":")
    pad = "  "
    maxlen = 16
    for label in take(keys(dict(dirac)), maxlen)
        println(io)
        print(io, labelrepr(dirac, label, pad))
    end
    if length(dirac) > maxlen
        println(io)
        print(io, "$pad$vdots")
    end
end

function dirac_showcompact(io::IO, dirac)
    pad = " "
    maxlen = 4
    for label in take(keys(dict(dirac)), maxlen)
        print(io, labelrepr(dirac, label, pad))
        i += 1
        if i <= length(dirac)
            print(io, "$pad+")
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
