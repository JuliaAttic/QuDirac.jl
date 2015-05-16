abstract SumAssoc{K,V} <: Associative{K,V}

Base.eltype{K,V}(::SumAssoc{K,V}) = V
labeltype{K,V}(::SumAssoc{K,V}) = K

###########
# SumTerm #
###########
type SumTerm{K,V} <: SumAssoc{K,V}
    key::K
    val::V
end

@compat SumTerm{K,V}(tup::Tuple{K,V}) = SumTerm{K,V}(tup[1], tup[2])

key(term::SumTerm) = term.key
val(term::SumTerm) = term.val

Base.eltype{K,V}(::Type{SumTerm{K,V}}) = V
labeltype{K,V}(::Type{SumTerm{K,V}}) = K

Base.copy(term::SumTerm) = SumTerm(copy(key(term)), copy(val(term)))

Base.promote_rule{K1,V1,K2,V2}(::Type{SumTerm{K1,V1}}, ::Type{SumTerm{K2,V2}}) = SumTerm{promote_type(K1,K2), promote_type(V1,V2)}
Base.convert{K,V}(::Type{SumTerm{K,V}}, term::SumTerm) = SumTerm(convert(K,key(term)), convert(V,val(term)))

Base.keys(term::SumTerm) = tuple(key(term))
Base.values(term::SumTerm) = tuple(values(term))

Base.length(term::SumTerm) = 1

Base.getindex(term::SumTerm, x) = haskey(term, x) ? val(term) : throw(KeyError(x))
Base.setindex!(term::SumTerm, x, y) = haskey(term, y) ? term.val = x : throw(KeyError(y))

Base.start(term::SumTerm) = false
Base.next(term::SumTerm, i) = (tuple(key(term), val(term)), true)
Base.done(term::SumTerm, i) = i
Base.collect(term::SumTerm) = [first(term)]

Base.haskey(term::SumTerm, k) = k == key(term)
Base.get(term::SumTerm, k, default=predict_zero(eltype(term))) = haskey(term, k) ? val(term) : default

###########
# SumTerm #
###########
type SumDict{K,V} <: SumAssoc{K,V}
    data::Dict{K,V}
    SumDict(data::Dict{K,V}) =new(data)
    SumDict(args...) = SumDict{K,V}(Dict{K,V}(args...))
end

SumDict{K,V}(data::Dict{K,V}) = SumDict{K,V}(data)
SumDict(term::SumTerm) = @compat SumDict(Dict(key(term) => val(term)))
SumDict(args...) = SumDict(Dict(args...))

Base.eltype{K,V}(::Type{SumDict{K,V}}) = V
labeltype{K,V}(::Type{SumDict{K,V}}) = K

Base.sizehint(dict::SumDict, len) = @compat (sizehint!(dict.data, len); return dict)
Base.copy(dict::SumDict) = SumDict(copy(dict.data))

Base.merge!(result::SumDict, dict::SumDict) = (merge!(result.data, dict.data); return result)

Base.hash(s::SumAssoc) = hash(filternz(s).data)
Base.hash(s::SumAssoc, h::Uint64) = hash(hash(s), h)
Base.(:(==))(a::SumAssoc, b::SumAssoc) = filternz(a).data == filternz(b).data

Base.promote_rule{K1,V1,K2,V2}(::Type{SumDict{K1,V1}}, ::Type{SumDict{K2,V2}}) = SumDict{promote_type(K1,K2), promote_type(V1,V2)}
Base.promote_rule{K1,V1,K2,V2}(::Type{SumTerm{K1,V1}}, ::Type{SumDict{K2,V2}}) = SumDict{promote_type(K1,K2), promote_type(V1,V2)}
Base.convert{K,V}(::Type{SumDict{K,V}}, term::SumTerm) = SumDict(convert(SumTerm{K,V}, term))
Base.convert{K,V}(::Type{SumDict{K,V}}, dict::SumDict) = SumDict(convert(Dict{K,V}, dict.data))

Base.keys(dict::SumDict) = keys(dict.data)
Base.values(dict::SumDict) = values(dict.data)

Base.length(dict::SumDict) = length(dict.data)

Base.getindex(dict::SumDict, x) = dict.data[x]
Base.setindex!(dict::SumDict, x, y) = haskey(dict, y) ? setindex!(dict.data, x, y) : throw(KeyError(y))

Base.start(dict::SumDict) = start(dict.data)
Base.next(dict::SumDict, i) = next(dict.data, i)
Base.done(dict::SumDict, i) = done(dict.data, i)
Base.collect(dict::SumDict) = collect(dict.data)

Base.haskey(dict::SumDict, k) = haskey(dict.data, k)
Base.get(dict::SumDict, k, default=predict_zero(eltype(dict))) = get(dict.data, k, default)

#############
# Filtering #
#############
filter(f::Union(Function,DataType), term::SumTerm) = filter(f, SumDict(term))
filter(f::Union(Function,DataType), dict::SumDict) = SumDict(filter(f, dict.data))
filter!(f::Union(Function,DataType), dict::SumDict) = (filter!(f, dict.data); return dict)

nzcoeff(k,v) = v != 0

filternz(s::SumAssoc) = filter(nzcoeff, d)
filternz!(dict::SumDict) = filter!(nzcoeff, dict)

###########
# Mapping #
###########
Base.map(f::Union(Function,DataType), term::SumTerm) = SumTerm(f(key(term), val(term)))
Base.map(f::Union(Function,DataType), dict::SumDict) = SumDict(map(f, collect(dict)))

function mapvals!(f, dict::SumDict)
    for (k,v) in dict
        dict.data[k] = f(v)
    end
    return dict
end
mapvals!(f, term::SumTerm) = (term.val = f(val(term)); return term)

mapvals(f, term::SumTerm) = SumTerm(key(term), f(val(term)))
mapvals(f, d::SumDict) = SumDict(zip(collect(keys(d)), map(f, collect(values(d)))))
mapkeys(f, term::SumTerm) = SumTerm(f(key(term)), val(term))
mapkeys(f, d::SumDict) = SumDict(zip(map(f, collect(keys(d))), collect(values(d))))

###########
# Scaling #
###########
scale_result{K,V,T}(d::SumDict{K,V}, ::T) = @compat sizehint!(SumDict{K, promote_type(T,V)}(), length(d))

function Base.scale!(dict::SumDict, c::Number)
    for k in keys(dict)
        dict.data[k] *= c
    end
    return dict
end

Base.scale!(term::SumTerm, c::Number) = (term.val *= c; return term)
Base.scale!(c::Number, sa::SumAssoc) = scale!(sa, c)

Base.scale(dict::SumDict, c::Number) = Base.scale!(merge!(scale_result(dict,c), dict), c)
Base.scale(term::SumTerm, c::Number) = SumTerm(key(term), val(term) * c)
Base.scale(c::Number, sa::SumAssoc) = scale(sa, c)

Base.(:*)(a::SumAssoc, b::SumAssoc) = mul_sums(a, b)
Base.(:*)(c::Number, s::SumAssoc) = scale(c, s)
Base.(:*)(s::SumAssoc, c::Number) = scale(s, c)
Base.(:/)(s::SumAssoc, c::Number) = scale(s, 1/c)
Base.(:-)(s::SumAssoc) = scale(s, -1)

#######
# Mul #
#######
function mul_merge!(result::SumDict, a::SumDict, b::SumDict)
    for (k,v) in a
        for (l,c) in b
            result.data[k*l] = v*c
        end
    end
    return result
end

function mul_merge!(result::SumDict, dict::SumDict, term::SumTerm)
    k0,v0 = key(term), val(term)
    for (k,v) in dict
        result.data[k*k0] = v*v0
    end
    return result
end

function mul_merge!(result::SumDict, term::SumTerm, dict::SumDict)
    k0,v0 = key(term), val(term)
    for (k,v) in dict
        result.data[k0*k] = v*v0
    end
    return result
end

mul_result{A,B,T,V}(a::SumDict{A,T}, b::SumDict{B,V}) = @compat sizehint!(SumDict{tensor_type(A,B), promote_type(T,V)}(), length(a) * length(b))
mul_result{K,V,L,C}(d::SumDict{K,V}, ::SumTerm{L,C}) = @compat sizehint!(SumDict{tensor_type(K,L), promote_type(V,C)}(), length(d))

mul_sums(a::SumDict, b::SumDict) = mul_merge!(mul_result(a, b), a, b)
mul_sums(dict::SumDict, term::SumTerm) = mul_merge!(mul_result(dict, term), dict, term)
mul_sums(term::SumTerm, dict::SumDict) = mul_merge!(mul_result(dict, term), term, dict)
mul_sums(a::SumTerm, b::SumTerm) = SumTerm(key(a) * key(b), val(a) * val(b))

##########
# Tensor #
##########
function tensor_merge!(result::SumDict, a::SumDict, b::SumDict)
    for (k,v) in a
        for (l,c) in b
            result.data[tensor(k,l)] = v*c
        end
    end
    return result
end

function tensor_merge!(result::SumDict, dict::SumDict, term::SumTerm)
    k0,v0 = key(term), val(term)
    for (k,v) in dict
        result.data[tensor(k,k0)] = v*v0
    end
    return result
end

function tensor_merge!(result::SumDict, term::SumTerm, dict::SumDict)
    k0,v0 = key(term), val(term)
    for (k,v) in dict
        result.data[tensor(k0,k)] = v*v0
    end
    return result
end

tensor_result{A,B,T,V}(a::SumDict{A,T}, b::SumDict{B,V}) = @compat sizehint!(SumDict{tensor_type(A,B), promote_type(T,V)}(), length(a) * length(b))
tensor_result{K,V,L,C}(d::SumDict{K,V}, ::SumTerm{L,C}) = @compat sizehint!(SumDict{tensor_type(K,L), promote_type(V,C)}(), length(d))

tensor(a::SumDict, b::SumDict) = tensor_merge!(tensor_result(a, b), a, b)
tensor(dict::SumDict, term::SumTerm) = tensor_merge!(tensor_result(dict, term), dict, term)
tensor(term::SumTerm, dict::SumDict) = tensor_merge!(tensor_result(dict, term), term, dict)
tensor(a::SumTerm, b::SumTerm) = SumTerm(tensor(key(a), key(b)), val(a) * val(b))

############
# Addition #
############
function add_to_dict!(dict::SumDict, k, v)
    if v != 0
        dict.data[k] = get(dict, k) + v
    end
    return dict
end

function add_merge!(result::SumDict, dict::SumDict)
    for (k,v) in dict
        add_to_dict!(result, k, v)
    end
    return result
end

add_result{A,B,T,V}(a::SumDict{A,T}, b::SumDict{B,V}) = @compat sizehint!(SumDict{promote_type(A,B), promote_type(T,V)}(), max(length(a), length(b)))
add_result{K,V,L,C}(d::SumDict{K,V}, ::SumTerm{L,C}) = @compat sizehint!(SumDict{promote_type(K,L), promote_type(V,C)}(), length(d))
add_result{K,V,L,C}(::SumTerm{K,V}, ::SumTerm{L,C}) = SumDict{promote_type(K,L), promote_type(V,C)}()

Base.(:+)(a::SumDict, b::SumDict) = add_merge!(merge!(add_result(a,b), a), b)
Base.(:+)(dict::SumDict, term::SumTerm) = add_to_dict!(merge!(add_result(dict,term), dict), key(term), val(term))
Base.(:+)(term::SumTerm, dict::SumDict) = dict + term
function Base.(:+)(a::SumTerm, b::SumTerm)
    result = add_result(a,b)
    add_to_dict!(result, key(a), val(a))
    add_to_dict!(result, key(b), val(b))
    return result
end

###############
# Subtraction #
###############
function sub_merge!(result::SumDict, d::SumDict)
    for (k,v) in d
        add_to_dict!(result, k, -v)
    end
    return result
end

Base.(:-)(a::SumDict, b::SumDict) = sub_merge!(merge!(add_result(a,b), a), b)
Base.(:-)(term::SumTerm, dict::SumDict) = add_to_dict!(scale!(merge!(add_result(dict,term), dict), -1), key(term), val(term))
Base.(:-)(dict::SumDict, term::SumTerm) = add_to_dict!(merge!(add_result(dict,term), dict), key(term), -val(term))
function Base.(:-)(a::SumTerm, b::SumTerm)
    result = add_result(a,b)
    add_to_dict!(result, key(a), val(a))
    add_to_dict!(result, key(b), -val(b))
    return result
end

export SumTerm, SumDict, filternz, filternz!, mapvals, mapkeys
