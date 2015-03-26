#########
# mapkv #
#########
function mapkv!(f::Function, result, d)
    for (k,v) in d
        (k0,v0) = f(k,v)
        result[k0] = v0
    end
    return result
end

function mapkv!(f::Function, d)
    for (k,v) in d
        delete!(d,k)
        (k0,v0) = f(k,v)
        d[k0] = v0
    end
    return d
end

mapkv(f::Function, d) = mapkv!(f, similar(d), d)

###########
# mapvals #
###########
function mapvals!(f::Function, result, d)
    for (k,v) in d
        result[k] = f(v)
    end
    return result
end

mapvals!(f::Function, d) = mapvals!(f,d,d)
mapvals(f::Function, d) = mapvals!(f, similar(d), d)

###########
# mapkeys #
###########
function mapkeys!(f::Function, result, d)
    for (k,v) in d
        result[f(k)] = v
    end
    return result
end

function mapkeys!(f::Function, d)
    for (k,v) in d
        delete!(d,k)
        d[f(k)] = v
    end
    return d
end

mapkeys(f::Function, d) = mapkeys!(f, similar(d), d)

##########
# tensor #
##########
function tensorstate!(result, d::Tuple)
    for pairs in product(d...)
        result[vcat(map(first, pairs)...)] = prod(second, pairs)
    end
    return result
end

function tensorop!(result, d::Tuple)
    for pairs in product(d...)
        labels = map(first, pairs)
        result[OpLabel(vcat(map(ktlabel, labels)...), vcat(map(brlabel, labels)...))] = prod(second, pairs)
    end
    return result
end

#########
# merge #
#########
function addmerge!(d, others...)
    for other in others
        for (k,v) in other
            d[k] = v + get(d,k,0)
        end
    end
    return d
end

function addmerge(d, others...)
    return addmerge!(similar(d), d, others...)
end

##########
# dscale #
##########
function dscale!(result, d, c)
    for (k,v) in d
        result[k] = c*v
    end
    return result
end

dscale!(d, c) = dscale!(d, d, c)
dscale(d, c) = dscale!(similar(d), d, c)

########
# misc #
########
function single_dict(dict, label, c)
    dict[label] = c
    return dict
end

function add_to_dict!(dict, label, c)
    dict[label] = get(dict, label, 0)+c
    return dict
end
