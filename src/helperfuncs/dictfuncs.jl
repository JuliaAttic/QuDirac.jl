function singlet_dict(label, c)
    o = ObjectIdDict()
    o[label] = c
    return o
end

function mergecart!(f::Function, result, d)
    for pairs in product(d...)
        (k,v) = f(pairs)
        result[k] = v
    end
    return result
end

function mergef!(f::Function, d, others...)
    for other in others
        for (k,v) in other
            if haskey(d, k)
                d[k] = f(d[k], v)
            else   
                d[k] = v
            end
        end
    end
    return d
end

function mergef(f::Function, d, others...)
    return mergef!(f, similar(d), d, others...)
end

castvals(f::Function, a::ObjectIdDict, b::ObjectIdDict) = mergef(f, a, b)
castvals(f::Function, d::ObjectIdDict, c) = mapvals(v->f(c,v), d)
castvals(f::Function, c, d::ObjectIdDict) = mapvals(v->f(c,v), d)

function mapkv!(f::Function, result, d)
    for (k,v) in d
        (k0,v0) = f(k,v)
        result[k0] = v0
    end
    return result
end

mapkv(f::Function, d) = mapkv!(f, similar(d), d)

function mapvals!(f::Function, result, d)
    for (k,v) in d
        result[k] = f(v)
    end
    return result
end

mapvals(f::Function, d) = mapvals!(f, similar(d), d)

function mapkeys!(f::Function, result, d)
    for (k,v) in d
        result[f(k)] = v
    end
    return result
end

mapkeys(f::Function, d) = mapkeys(f, similar(d), d)