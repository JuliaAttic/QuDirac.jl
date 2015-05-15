##########
# KVPair #
##########
# TODO: copy, hash, convert, promotion, 
#       ==, keys, values, iteration
type KVPair{K,V}
    key::K
    val::V
end

key(pair::KVPair) = pair.key
val(pair::KVPair) = pair.val

###########
# Scaling #
###########
function scale_vals!(d::Dict, c)
    for k in keys(d)
        d[k] *= c
    end
    return result
end

scale_vals!(pair::KVPair, c) = (pair.val *= c; return pair)

scale_result{K,V,T}(d::Dict{K,V}, ::T) = @compat(sizehint!(Dict{K, promote_type(T,V)}(), length(d)))

scale_vals(dict::Dict, c) = dscale!(merge!(scale_result(dict,c), dict), c)
scale_vals(pair::KVPair, c) = KVPair(key(pair), val(pair) * c)

##########
# Tensor #
##########
function tensor_merge!(result::Dict, a::Dict, b::Dict)
    for (k,v) in a
        for (l,c) in b
            result[tensor(k,l)] = v*c
        end
    end
    return result
end

function tensor_merge!(result::Dict, dict::Dict, pair::KVPair)
    k0,v0 = key(pair), val(pair)
    for (k,v) in dict
        result[tensor(k,k0)] = v*v0
    end
    return result
end

function tensor_merge!(result::Dict, pair::KVPair, dict::Dict)
    k0,v0 = key(pair), val(pair)
    for (k,v) in dict
        result[tensor(k0,k)] = v*v0
    end
    return result
end

tensor_result{A,B,T,V}(a::Dict{A,T}, b::Dict{B,V}) = @compat(sizehint!(Dict{tensor_type(A,B), promote_type(T,V)}(), length(a) * length(b)))
tensor_result{K,V,L,C}(d::Dict{K,V}, ::KVPair{L,C}) = @compat(sizehint!(Dict{tensor_type(K,L), promote_type(V,C)}(), length(d)))

tensor_merge(a::Dict, b::Dict) = tensor_merge!(tensor_result(a, b), a, b)
tensor_merge(dict::Dict, pair::KVPair) = tensor_merge!(tensor_result(dict, pair), dict, pair)
tensor_merge(pair::KVPair, dict::Dict) = tensor_merge!(tensor_result(dict, pair), pair, dict)
tensor_merge(a::KVPair, b::KVPair) = DiracTerm(tensor(key(a), key(b)), val(a) * val(b))

############
# Addition #
############
function add_to_dict!(dict, k, v)
    if v != 0
        dict[k] = get(dict, k, 0) + v
    end
    return dict
end

function add_merge!(result::Dict, d::Dict)
    for (k,v) in d
        add_to_dict!(result, k, v)
    end
    return result
end

add_result{A,B,T,V}(a::Dict{A,T}, b::Dict{B,V}) = @compat(sizehint!(Dict{promote_type(A,B), promote_type(T,V)}(), max(length(a), length(b))))
add_result{K,V,L,C}(d::Dict{K,V}, ::KVPair{L,C}) = @compat(sizehint!(Dict{promote_type(K,L), promote_type(V,C)}(), length(d)))
add_result{K,V,L,C}(::KVPair{K,V}, ::KVPair{L,C}) = Dict{promote_type(K,L), promote_type(V,C)}()

add_merge(a::Dict, b::Dict) = add_merge!(merge!(add_result(a,b), a), b)
add_merge(dict::Dict, pair::KVPair) = add_to_dict!(merge!(add_result(dict,pair), dict), key(pair), val(pair))
add_merge(pair::KVPair, dict::Dict) = add_merge(dict, pair)
function add_merge(a::KVPair, b::KVPair)
    result = add_result(a,b)
    add_to_dict!(result, a)
    add_to_dict!(result, b)
    return result
end

###############
# Subtraction #
###############
function sub_merge!(result::Dict, d::Dict)
    for (k,v) in d
        add_to_dict!(result, k, -v)
    end
    return result
end

sub_merge(a::Dict, b::Dict) = sub_merge!(merge!(add_result(a,b), a), b)
sub_merge(pair::KVPair, dict::Dict) = add_to_dict!(scale_vals!(merge!(add_result(dict,pair), dict),-1), key(pair), val(pair))
sub_merge(dict::Dict, pair::KVPair) = add_to_dict!(merge!(add_result(dict,pair), dict), key(pair), -val(pair))
function sub_merge(a::KVPair, b::KVPair)
    result = add_result(a,b)
    add_to_dict!(result, a)
    add_to_dict!(result, -b)
    return result
end

