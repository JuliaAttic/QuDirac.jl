single_dict(dict, label, c) = (dict[label] = c; return dict)

function mergecart!(f::Function, result, d::Tuple)
    for pairs in product(d...)
        (k,v) = f(pairs)
        result[k] = v
    end
    return result
end

mergecart!(f::Function, result, d...) = mergecart!(f, result, d)

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

castvals(f::Function, d::Associative, c::Number) = mapvals(v->f(c,v), d)
castvals(f::Function, c::Number, d::Associative) = mapvals(v->f(c,v), d)

castvals!(f::Function, d::Associative, c::Number) = mapvals!(v->f(c,v), d)
castvals!(f::Function, c::Number, d::Associative) = mapvals!(v->f(c,v), d)

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

mapvals!(f::Function, d) = mapvals!(f,d,d)
mapvals(f::Function, d) = mapvals!(f, similar(d), d)

function mapkeys!(f::Function, result, d)
    for (k,v) in d
        result[f(k)] = v
    end
    return result
end

mapkeys(f::Function, d) = mapkeys!(f, similar(d), d)