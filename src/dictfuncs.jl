function tensordict!(result, a, b)
    for (k,v) in a
        for (l,c) in b
            result[tensor(k,l)] = v*c
        end
    end
    return result
end

tensor_result{N,M,T,V}(::StateDict{N,T}, ::StateDict{M,V}) = StateDict{N+M, promote_type(T,V)}()
tensor_result{N,M,T,V}(::OpDict{N,T}, ::OpDict{M,V}) = OpDict{N+M, promote_type(T,V)}()
tensordict(a, b) = tensordict!(tensor_result(a, b), a, b)

function add_merge!(result, other)
    for (k,v) in other
        add_to_dict!(result, k, v)
    end
    return result
end

function sub_merge!(result, other)
    for (k,v) in other
        add_to_dict!(result, k, -v)
    end
    return result
end

merge_result{K,T,V}(::Dict{K,T}, ::Dict{K,V}) = Dict{K, promote_type(T,V)}()
add_merge(a, b) = add_merge!(merge(merge_result(a,b), a), b)
sub_merge(a, b) = sub_merge!(merge(merge_result(a,b), a), b)

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
