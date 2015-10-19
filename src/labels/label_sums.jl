abstract AbstractLabelSum{L,T}

any_zero{T}(::Type{T}) = zero(T)
any_zero(::Type{Any}) = false

#########
# Types #
#########
immutable LabelTerm{L,T} <: AbstractLabelSum{L,T}
    label::L
    coeff::T
end

LabelTerm(pair::Pair) = LabelTerm(pair.first, pair.second)

term(label) = LabelTerm(label, 1)

type LabelSum{L,T} <: AbstractLabelSum{L,T}
    data::Dict{L,T}
    LabelSum(data::Dict{L,T}) = new(data)
    LabelSum(args...) = LabelSum{L,T}(Dict{L,T}(args...))
end

LabelSum{L,T}(data::Dict{L,T}) = LabelSum{L,T}(data)
LabelSum(args::LabelTerm...) = LabelSum([label(i)=>coeff(i) for i in args])
LabelSum(args...) = LabelSum(Dict(args...))

######################
# Property Functions #
######################
label(t::LabelTerm) = t.label
coeff(t::LabelTerm) = t.coeff

data(s::LabelSum) = s.data

Base.eltype{L,T}(::Type{LabelSum{L,T}}) = LabelTerm{L,T}

labeltype{L,T}(::AbstractLabelSum{L,T}) = L
labeltype{L,T}(::Type{LabelTerm{L,T}}) = L
labeltype{L,T}(::Type{LabelSum{L,T}}) = L

coefftype{L,T}(::AbstractLabelSum{L,T}) = T
coefftype{L,T}(::Type{LabelTerm{L,T}}) = T
coefftype{L,T}(::Type{LabelSum{L,T}}) = T

########################
# Conversion/Promotion #
########################
Base.promote_rule{L1,L2,T1,T2}(::Type{LabelTerm{L1,T1}}, ::Type{LabelTerm{L2,T2}}) = LabelTerm{promote_type(L1,L2),promote_type(T1,T2)}
Base.promote_rule{L1,L2,T1,T2}(::Type{LabelTerm{L1,T1}}, ::Type{LabelSum{L2,T2}}) = LabelSum{promote_type(L1,L2),promote_type(T1,T2)}
Base.promote_rule{L1,L2,T1,T2}(::Type{LabelSum{L1,T1}}, ::Type{LabelSum{L2,T2}}) = LabelSum{promote_type(L1,L2),promote_type(T1,T2)}

Base.convert{L,T}(::Type{LabelTerm{L,T}}, t::LabelTerm) = LabelTerm(convert(L, label(t)), convert(T, coeff(t)))
Base.convert{L,T}(::Type{LabelTerm{L,T}}, t::LabelTerm{L,T}) = t
Base.convert{L,T}(::Type{LabelSum{L,T}}, t::LabelTerm) = LabelSum(label(t) => coeff(t))
Base.convert{L,T}(::Type{LabelSum{L,T}}, s::LabelSum) = LabelSum(convert(Dict{L,T}, data(s)))
Base.convert{L,T}(::Type{LabelSum{L,T}}, s::LabelSum{L,T}) = s

############################
# Hashing/Equality/Copying #
############################
Base.copy(t::LabelTerm) = LabelTerm(copy(label(t)), copy(coeff(t)))
Base.copy(s::LabelSum) = LabelSum(copy(data(s)))

const labelsum_hash = hash(AbstractLabelSum)

Base.hash(s::LabelSum) = hash(data(s), labelsum_hash)
Base.hash(t::LabelTerm) = hash(LabelSum(t))
Base.hash(s::AbstractLabelSum, h::UInt64) = hash(hash(s), h)

Base.(:(==))(a::LabelSum, b::LabelSum) = data(a) == data(b)
Base.(:(==))(a::LabelTerm, b::LabelTerm) = coeff(a) == coeff(b) && label(a) == label(b)
Base.(:(==))(t::LabelTerm, s::LabelSum) = s == t
Base.(:(==))(s::LabelSum, t::LabelTerm) = length(s) == 1 && first(s) == t

####################################
# Associative Functions on LabelSum #
####################################
Base.length(s::LabelSum) = length(data(s))

Base.similar(s::LabelSum) = LabelSum(similar(data(s)))

Base.get(s::LabelSum, k, default=any_zero(coefftype(s))) = get(data(s), k, default)

Base.getindex(s::LabelSum, k) = getindex(data(s), k)

Base.setindex!(s::LabelSum, v, k) = setindex!(data(s), v, k)

Base.isempty(s::LabelSum) = isempty(data(s))

Base.empty!(s::LabelSum) = (empty!(data(s)); return s)

Base.sizehint!(s::LabelSum, len) = (sizehint!(data(s), len); return s)

Base.delete!(s::LabelSum, k) = (delete!(data(s), k); return s)

Base.merge!(result::LabelSum, other::LabelSum) = (merge!(data(result), data(other)); return result)
Base.merge!(result::LabelSum, other::LabelTerm) = (setindex!(result, coeff(other), label(other)); return result)

Base.merge(a::LabelSum, b::LabelSum) = LabelSum(merge(data(a), data(b)))
Base.merge(s::LabelSum, t::LabelTerm) = merge(s, LabelSum(t))
Base.merge(t::LabelTerm, s::LabelSum) = merge(s, t)

haslabel(s::LabelSum, k) = haskey(data(s), k)

########################
# Iteration/Collection #
########################
labels(s::LabelSum) = keys(data(s))
coeffs(s::LabelSum) = values(data(s))

Base.start(s::LabelSum) = start(data(s))

function Base.next(s::LabelSum, i)
    item, state = next(data(s), i)
    return LabelTerm(item), state
end

Base.done(s::LabelSum, i) = done(data(s), i)

function Base.collect(s::LabelSum)
    arr = Vector{eltype(s)}(length(s))
    i = 1
    for t in s
        arr[i] = t
        i += 1
    end
    return arr
end

###################
# Pretty Printing #
###################
Base.showcompact(io::IO, t::LabelTerm) = print(io, "$(coeff(t)) * $(repr(label(t)))")

function Base.show(io::IO, t::LabelTerm)
    println(io, typeof(t), ": ")
    showcompact(io, t)
end

function show_loop(io::IO, s::LabelSum, cutoff::Int)
    i = 1
    for t in s
        if i == length(s)
            showcompact(io, t)
        elseif i > cutoff
            print(io, "...")
            break
        else
            showcompact(io, t)
            print(io, " + ")
            i += 1
        end
    end
end

function Base.show(io::IO, s::LabelSum)
    println(io, typeof(s), ": ")
    show_loop(io, s, 15)
end

Base.showcompact(io::IO, s::LabelSum) = show_loop(io, s, 3)
