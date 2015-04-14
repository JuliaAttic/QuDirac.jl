function tensordict!(result, a, b)
    for (k,v) in a
        for (l,c) in b
            result[tensor(k,l)] = v*c
        end
    end
    return result
end

function add_merge!(result, other)
    for (k,v) in other
        add_to_dict!(result, k, v)
    end
    return result
end

function add_merge{K,A,B}(a::Dict{K,A}, b::Dict{K,B})
    return add_merge!(merge(Dict{K,promote_type(A,B)}(), a), b)
end

function sub_merge!(result, other)
    for (k,v) in other
        add_to_dict!(result, k, -v)
    end
    return result
end

function sub_merge{K,A,B}(a::Dict{K,A}, b::Dict{K,B})
    return sub_merge!(merge(Dict{K,promote_type(A,B)}(), a), b)
end

function dscale!(result, d, c)
    for (k,v) in d
        result[k] = c*v
    end
    return result
end

dscale!(d, c) = dscale!(d, d, c)
dscale{K,V,T}(d::Dict{K,V}, c::T) = dscale!(Dict{K,promote_type(V,T)}(), d, c)

function add_to_dict!(dict, label, c)
    if c != 0
        dict[label] = get(dict, label, 0) + c
    end
    return dict
end
