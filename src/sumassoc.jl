abstract SumAssoc{K,V} <: Associative{K,V}

###################
# SumTerm/SumDict #
###################
immutable SumTerm{K,V} <: SumAssoc{K,V}
    key::K
    val::V
end

@compat SumTerm{K,V}(tup::Tuple{K,V}) = SumTerm{K,V}(tup[1], tup[2])

type SumDict{K,V} <: SumAssoc{K,V}
    data::Dict{K,V}
    SumDict(data::Dict{K,V}) = new(data)
    SumDict(args...) = SumDict{K,V}(Dict{K,V}(args...))
end

SumDict{K,V}(data::Dict{K,V}) = SumDict{K,V}(data)
SumDict(term::SumTerm) = @compat SumDict(Dict(key(term) => val(term)))
SumDict(args...) = SumDict(Dict(args...))

####################
# Access Functions #
####################
key(term::SumTerm) = term.key
val(term::SumTerm) = term.val

Base.eltype{K,V}(::SumAssoc{K,V}) = V
Base.eltype{K,V}(::Type{SumTerm{K,V}}) = V
Base.eltype{K,V}(::Type{SumDict{K,V}}) = V

labeltype{K,V}(::SumAssoc{K,V}) = K
labeltype{K,V}(::Type{SumTerm{K,V}}) = K
labeltype{K,V}(::Type{SumDict{K,V}}) = K

########################
# Conversion/Promotion #
########################
Base.promote_rule{K1,V1,K2,V2}(::Type{SumTerm{K1,V1}}, ::Type{SumTerm{K2,V2}}) = SumTerm{promote_type(K1,K2), promote_type(V1,V2)}
Base.promote_rule{K1,V1,K2,V2}(::Type{SumDict{K1,V1}}, ::Type{SumDict{K2,V2}}) = SumDict{promote_type(K1,K2), promote_type(V1,V2)}
Base.promote_rule{K1,V1,K2,V2}(::Type{SumTerm{K1,V1}}, ::Type{SumDict{K2,V2}}) = SumDict{promote_type(K1,K2), promote_type(V1,V2)}

Base.convert{K,V}(::Type{SumTerm{K,V}}, term::SumTerm) = SumTerm(convert(K,key(term)), convert(V,val(term)))
Base.convert{K,V}(::Type{SumTerm{K,V}}, term::SumTerm{K,V}) = term
Base.convert{K,V}(::Type{SumDict{K,V}}, term::SumTerm) = SumDict(convert(SumTerm{K,V}, term))
Base.convert{K,V}(::Type{SumDict{K,V}}, dict::SumDict) = SumDict(convert(Dict{K,V}, dict.data))
Base.convert{K,V}(::Type{SumDict{K,V}}, dict::SumDict{K,V}) = dict

############################
# Hashing/Equality/Copying #
############################
Base.copy(term::SumTerm) = SumTerm(key(term), val(term))
Base.copy(dict::SumDict) = SumDict(copy(dict.data))

Base.similar(term::SumTerm) = similar(SumDict(term))
Base.similar(dict::SumDict) = SumDict(similar(dict.data))

Base.hash(s::SumAssoc) = hash(filternz(s).data)
Base.hash(s::SumAssoc, h::Uint64) = hash(hash(s), h)

Base.(:(==))(a::SumAssoc, b::SumAssoc) = filternz(a).data == filternz(b).data

#########################
# Associative Functions #
#########################
Base.merge!(result::SumDict, dict::SumDict) = (merge!(result.data, dict.data); return result)

Base.length(term::SumTerm) = 1
Base.length(dict::SumDict) = length(dict.data)

Base.getindex(term::SumTerm, x) = haskey(term, x) ? val(term) : throw(KeyError(x))
Base.getindex(dict::SumDict, x) = dict.data[x]

Base.setindex!(dict::SumDict, x, y) = setindex!(dict.data, x, y)

Base.haskey(term::SumTerm, k) = k == key(term)
Base.haskey(dict::SumDict, k) = haskey(dict.data, k)

Base.get(term::SumTerm, k, default=predict_zero(eltype(term))) = haskey(term, k) ? val(term) : default
Base.get(dict::SumDict, k, default=predict_zero(eltype(dict))) = get(dict.data, k, default)

Base.isempty(::SumTerm) = false
Base.isempty(dict::SumDict) = isempty(dict.data)
Base.empty!(dict::SumDict) = (empty!(dict.data); return dict)

if VERSION >= v"0.4-"
    Base.sizehint!(dict::SumDict, len) = (sizehint!(dict.data, len); return dict)
else
    Base.sizehint(dict::SumDict, len) = @compat (sizehint!(dict.data, len); return dict)
end

########################
# Iteration/Collection #
########################
Base.keys(term::SumTerm) = tuple(key(term))
Base.keys(dict::SumDict) = keys(dict.data)

Base.values(term::SumTerm) = tuple(values(term))
Base.values(dict::SumDict) = values(dict.data)

Base.start(term::SumTerm) = false
Base.next(term::SumTerm, i) = tuple(tuple(key(term), val(term)), true)
Base.done(term::SumTerm, i) = i

Base.start(dict::SumDict) = start(dict.data)
Base.next(dict::SumDict, i) = next(dict.data, i)
Base.done(dict::SumDict, i) = done(dict.data, i)

Base.collect(term::SumTerm) = [first(term)]
Base.collect(dict::SumDict) = collect(dict.data)

#############
# Filtering #
#############
Base.filter(f::Function, term::SumTerm) = filter(f, SumDict(term))
Base.filter(f::Function, dict::SumDict) = SumDict(filter(f, dict.data))
Base.filter!(f::Function, dict::SumDict) = (filter!(f, dict.data); return dict)

nzcoeff(k,v) = v != 0

filternz(s::SumAssoc) = filter(nzcoeff, s)
filternz!(dict::SumDict) = filter!(nzcoeff, dict)

###########
# Mapping #
###########
Base.map(f::Union(Function,DataType), term::SumTerm) = SumTerm(f(key(term), val(term)))
Base.map(f::Union(Function,DataType), dict::SumDict) = SumDict(map(kv -> f(kv[1], kv[2]), collect(dict)))

function mapvals!(f, dict::SumDict)
    for (k,v) in dict
        dict.data[k] = f(v)
    end
    return dict
end

mapvals(f, term::SumTerm) = SumTerm(key(term), f(val(term)))
mapvals(f, d::SumDict) = SumDict(zip(collect(keys(d)), map(f, collect(values(d)))))
mapkeys(f, term::SumTerm) = SumTerm(f(key(term)), val(term))
mapkeys(f, d::SumDict) = SumDict(zip(map(f, collect(keys(d))), collect(values(d))))

##########
# Tensor #
##########
tensor(a::SumDict, b::SumDict) = SumDict([tensor(ka,kb) => va*vb for (ka,va) in a, (kb,vb) in b])
tensor(dict::SumDict, term::SumTerm) = SumDict([tensor(k,key(term)) => v*val(term) for (k,v) in dict])
tensor(term::SumTerm, dict::SumDict) = SumDict([tensor(key(term),k) => val(term)*v for (k,v) in dict])
tensor(a::SumTerm, b::SumTerm) = SumTerm(tensor(key(a), key(b)), val(a) * val(b))

###########
# Scaling #
###########
function Base.scale!(dict::SumDict, c::Number)
    for k in keys(dict)
        dict.data[k] *= c
    end
    return dict
end

Base.scale!(c::Number, dict::SumDict) = scale!(dict, c)

Base.scale(dict::SumDict, c::Number) = SumDict([k => c * v for (k,v) in dict])
Base.scale(term::SumTerm, c::Number) = SumTerm(key(term), val(term) * c)
Base.scale(c::Number, s::SumAssoc) = scale(s, c)

Base.(:*)(a::SumAssoc, b::SumAssoc) = tensor(a, b)
Base.(:*)(c::Number, s::SumAssoc) = scale(c, s)
Base.(:*)(s::SumAssoc, c::Number) = scale(s, c)
Base.(:/)(s::SumAssoc, c::Number) = scale(s, 1/c)
Base.(:-)(s::SumAssoc) = scale(s, -1)

############
# Addition #
############
function add_to_sum!(dict::SumDict, k, v)
    if v != 0
        dict.data[k] = get(dict.data, k, 0) + v
    end
    return dict
end

function add!(result::SumDict, dict::SumDict)
    for (k,v) in dict
        add_to_sum!(result, k, v)
    end
    return result
end

add!(result::SumDict, term::SumTerm) = add_to_sum!(result, key(term), val(term))

add_result{A,B,T,V}(a::SumDict{A,T}, b::SumDict{B,V}) = SumDict{promote_type(A,B), promote_type(T,V)}()
add_result{K,V,L,C}(d::SumDict{K,V}, ::SumTerm{L,C}) = SumDict{promote_type(K,L), promote_type(V,C)}()
add_result{K,V,L,C}(::SumTerm{K,V}, ::SumTerm{L,C}) = SumDict{promote_type(K,L), promote_type(V,C)}()

Base.(:+)(a::SumDict, b::SumAssoc) = add!(merge!(add_result(a,b), a), b)
Base.(:+)(term::SumTerm, dict::SumDict) = dict + term
function Base.(:+)(a::SumTerm, b::SumTerm)
    result = add_result(a,b)
    add!(result, a)
    add!(result, b)
    return result
end

###############
# Subtraction #
###############
function sub!(result::SumDict, d::SumDict)
    for (k,v) in d
        add_to_sum!(result, k, -v)
    end
    return result
end

sub!(result::SumDict, term::SumTerm) = add_to_sum!(result, key(term), -val(term))

Base.(:-)(a::SumDict, b::SumAssoc) = sub!(merge!(add_result(a,b), a), b)

function Base.(:-)(term::SumTerm, dict::SumDict)
    result = scale!(merge!(add_result(dict,term), dict), -1)
    return add!(result, term)
end

function Base.(:-)(a::SumTerm, b::SumTerm)
    result = add_result(a,b)
    add!(result, a)
    sub!(result, b)
    return result
end
