#############
# Filtering #
#############
function Base.filter(f, s::LabelSum)
    result = similar(s)
    for t in s
        if f(t)
            add!(result, t)
        end
    end
    return result
end

function Base.filter!(f, s::LabelSum)
    for t in s
        if !(f(t))
            delete!(s, label(t))
        end
    end
    return s
end

hasnzcoeff{L,T}(t::LabelTerm{L,T}) = coeff(t) != zero(T)

filternz(s::LabelSum) = filter(hasnzcoeff, s)
filternz!(s::LabelSum) = filter!(hasnzcoeff, s)

###########
# Mapping #
###########
function Base.map(f, s::LabelSum)
    result = f(first(s))
    result -= result
    for t in s
        add!(result, f(t))
    end
    return result
end

###########
# replace #
###########
replace_label(f, label) = f(label)
replace_label(f, t::LabelTerm) = coeff(t) * replace_label(f, label(t))
replace_label(f, s::LabelSum) = replace(f, s)

function Base.replace(f, s::LabelSum)
    result = zero(coefftype(s))
    for t in s
        result += replace_label(f, t)
    end
    return result
end
