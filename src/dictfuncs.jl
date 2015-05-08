function tensor_dicts!(result, a, b)
    for (k,v) in a
        for (l,c) in b
            result[tensor(k,l)] = v*c
        end
    end
    return result
end

function tensor_single!(result, dict::Dict, state::SingleKet)
    c = state.coeff
    l = state.label
    for (k,v) in dict
        result[tensor(k,l)] = v*c
    end
    return result
end

function tensor_single!(result, state::SingleKet, dict::Dict)
    c = state.coeff
    l = state.label
    for (k,v) in dict
        result[tensor(l,k)] = v*c
    end
    return result
end

tensor_result{N,M,T,V}(::StateDict{N,T}, ::StateDict{M,V}) = StateDict{N+M, promote_type(T,V)}()
tensor_result{N,M,T,V}(::OpDict{N,T}, ::OpDict{M,V}) = OpDict{N+M, promote_type(T,V)}()
tensor_result{N,M,T,V}(::OpDict{N,T}, ::OpDict{M,V}) = OpDict{N+M, promote_type(T,V)}()
tensor_result{P,N,M,T,V}(::StateDict{N,T}, ::SingleKet{P,M,V}) = StateDict{N+M, promote_type(T,V)}()

tensor_merge(a::Dict, b::Dict) = tensor_dicts!(tensor_result(a, b), a, b)
tensor_merge(a::Dict, b::Dict) = tensor_dicts!(tensor_result(a, b), a, b)
tensor_merge(dict::Dict, state::SingleKet) = tensor_single!(tensor_result(dict, state), dict, state)
tensor_merge(state::SingleKet, dict::Dict) = tensor_single!(tensor_result(dict, state), state, dict)

function add_merge!(result, other)
    for (k,v) in other
        add_to_dict!(result, k, v)
    end
    return result
end

add_merge!(result, label, coeff) = add_to_dict!(result, label, coeff)

function sub_merge!(result, other)
    for (k,v) in other
        add_to_dict!(result, k, -v)
    end
    return result
end

sub_merge!(result, label, coeff) = add_to_dict!(result, label, -coeff)

merge_result{K,T,V}(::Dict{K,T}, ::Dict{K,V}) = Dict{K, promote_type(T,V)}()
merge_result{K,T,V}(::Dict{K,T}, ::V) = Dict{K, promote_type(T,V)}()

add_merge(a::Dict, b::Dict) = add_merge!(merge(merge_result(a,b), a), b)
function add_merge(dict::Dict, state::SingleKet)
    result = merge(merge_result(dict, state.coeff), dict)
    return add_merge!(result, state.label, state.coeff)
end

sub_merge(a::Dict, b::Dict) = sub_merge!(merge(merge_result(a,b), a), b)
function sub_merge(dict::Dict, state::SingleKet)
    result = merge(merge_result(dict, state.coeff), dict)
    sub_merge!(result, state.label, state.coeff)
end

function dscale!(result, d, c)
    for (k,v) in d
        result[k] = c*v
    end
    return result
end
dscale!(d, c) = dscale!(d, d, c)

scale_result{K,V,T}(::Dict{K,V}, ::T) = Dict{K, promote_type(T,V)}()
dscale(d, c) = dscale!(scale_result(d,c), d, c)

function add_to_dict!(dict, label, c)
    if c != 0
        dict[label] = get(dict, label, 0) + c
    end
    return dict
end
