##################
# Multiplication #
##################

# AbstractLabelSum * Number #
#--------------------------#
function Base.scale!(s::LabelSum, x::Number)
    for t in s
        s[label(t)] = coeff(t)*x
    end
    return s
end

Base.scale!(x::Number, s::LabelSum) = scale!(s, x)

Base.scale(s::LabelSum, x::Number) = LabelSum([label(t) => coeff(t)*x for t in s])
Base.scale(t::LabelTerm, x::Number) = LabelTerm(label(t), coeff(t)*x)
Base.scale(x::Number, s::AbstractLabelSum) = scale(s, x)

Base.(:*)(s::AbstractLabelSum, x::Number) = scale(s, x)
Base.(:*)(x::Number, s::AbstractLabelSum) = scale(x, s)

Base.(:-)(s::LabelSum) = scale(s, -one(coefftype(s)))
Base.(:-)(t::LabelTerm) = LabelTerm(label(t), -coeff(t))

Base.(:+)(s::AbstractLabelSum) = s

# AbstractLabelSum * AbstractLabelSum #
#-----------------------------------#
Base.(:*)(a::LabelTerm, b::LabelTerm) = LabelTerm(label(a)*label(b), coeff(a)*coeff(b))

function Base.(:*)(t::LabelTerm, s::LabelSum)
    lbl, c = label(t), coeff(t)
    return LabelSum([lbl*label(x) => c*coeff(x) for x in s])
end

function Base.(:*)(s::LabelSum, t::LabelTerm)
    lbl, c = label(t), coeff(t)
    return LabelSum([label(x)*lbl => coeff(x)*c for x in s])
end

Base.(:*)(a::LabelSum, b::LabelSum) = LabelSum([label(x)*label(y) => coeff(x)*coeff(y) for x in a, y in b])

############
# Division #
############
Base.(:/)(s::LabelSum, x::Number) = LabelSum([label(t) => coeff(t)/x for t in s])
Base.(:/)(t::LabelTerm, x::Number) = LabelTerm(label(t), coeff(t)/x)

############
# Addition #
############
add_result{K1,V1,K2,V2}(::AbstractLabelSum{K1,V1}, ::AbstractLabelSum{K2,V2}) = LabelSum{promote_type(K1,K2), promote_type(V1,V2)}()

function add!(f, result::LabelSum, itr)
    for i in itr
        add!(result, f(i))
    end
    return result
end

function add!(result::LabelSum, other::LabelTerm)
    lbl, c = label(other), coeff(other)
    result[lbl] = get(result, lbl, zero(c)) + c
    return result
end

function add!(result::LabelSum, other::LabelSum)
    for t in other
        add!(result, t)
    end
    return result
end

Base.(:+)(a::LabelSum, b::AbstractLabelSum) = add!(merge!(add_result(a, b), a), b)
Base.(:+)(a::LabelTerm, b::LabelSum) = b + a

function Base.(:+)(a::LabelTerm, b::LabelTerm)
    result = add_result(a, b)
    add!(result, a)
    add!(result, b)
    return result
end

###############
# Subtraction #
###############
sub!(result::LabelSum, other::LabelSum) = add!(-, result, other)
sub!(result::LabelSum, other::LabelTerm) = add!(result, -other)

Base.(:-)(a::LabelSum, b::AbstractLabelSum) = sub!(merge!(add_result(a, b), a), b)

function Base.(:-)(a::LabelTerm, b::LabelSum)
    result = scale!(merge!(add_result(a, b), b), -one(coefftype(b)))
    return add!(result, term)
end

function Base.(:-)(a::LabelTerm, b::LabelTerm)
    result = add_result(a, b)
    add!(result, a)
    sub!(result, b)
    return result
end
