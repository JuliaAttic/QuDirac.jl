##############
# StateLabel #
##############
immutable StateLabel{N,T}
    label::Vector{T}
    hash::Uint64
    StateLabel(label::Vector{T}, hash::Uint64) = new(label, hash)
    StateLabel(label::Vector{T}) = StateLabel{N,T}(label, hash(label))
end

StateLabel(s::StateLabel) = s
StateLabel(label::NTuple{0}) = error("Cannot have a 0-factor StateLabel")
StateLabel{N,T}(label::NTuple{N,T}) = StateLabel{N,T}(collect(label))
StateLabel{N}(label::NTuple{N}) = StateLabel{N,Any}(collect(label))
StateLabel(i...) = StateLabel(i)

nfactors{N,T}(::Type{StateLabel{N,T}}) = N
Base.eltype{N,T}(::StateLabel{N,T}) = T
Base.eltype{N,T}(::Type{StateLabel{N,T}}) = T
Base.eltype(s::StateLabel) = eltype(typeof(s))
Base.getindex(s::StateLabel, i) = s.label[i]
Base.getindex(s::StateLabel, arr::AbstractArray) = s.label[arr]

except{N,T}(s::StateLabel{N,T}, i) = StateLabel{N-1,T}(deleteat!(copy(s.label), i))
setindex{N,T}(s::StateLabel{N,T}, x, y) = StateLabel(s.label[1:y-1]..., x, s.label[y+1:end]...)
setindex{N,T}(s::StateLabel{N,T}, x::T, y) = StateLabel{N,T}(setindex!(copy(s.label), x, y))
switch{N,T}(s::StateLabel{N,T}, i, j) =  StateLabel{N,T}(switch!(copy(s.label), i, j))
permute{N,T}(label::StateLabel{N,T}, perm::Vector) = StateLabel{N,T}(label[perm])

Base.copy{N,T}(s::StateLabel{N,T}) = StateLabel{N,T}(s.label, s.hash)
Base.hash(s::StateLabel) = s.hash
Base.hash(s::StateLabel, h::Uint64) = hash(hash(s), h)

Base.(:(==)){N}(a::StateLabel{N},b::StateLabel{N}) = hash(a) == hash(b)

Base.start(s::StateLabel) = start(s.label)
Base.next(s::StateLabel, i) = next(s.label, i)
Base.done(s::StateLabel, i) = done(s.label, i)

Base.first(s::StateLabel) = first(s.label)
Base.last(s::StateLabel) = last(s.label)
Base.endof(s::StateLabel) = endof(s.label)

Base.length{N}(::StateLabel{N}) = N

is_sum_x(s::StateLabel, x) = sum(s) == x

tensor(a::StateLabel, b::StateLabel) = StateLabel(a.label..., b.label...)
tensor{N,M,T}(a::StateLabel{N,T}, b::StateLabel{M,T}) = StateLabel{N+M,T}(vcat(a.label, b.label))

Base.(:*)(a::StateLabel, b::StateLabel) = tensor(a,b)

tensor_type{N,M,A,B}(::Type{StateLabel{N,A}}, ::Type{StateLabel{M,B}}) = StateLabel{N+M,Any}
tensor_type{N,M,T}(::Type{StateLabel{N,T}}, ::Type{StateLabel{M,T}}) = StateLabel{N+M,T}

label_promote{L1,L2}(::Type{L1}, ::Type{L2}) = Any
label_promote{L}(::Type{L}, ::Type{L}) = L

Base.promote_type{N,A,B}(::Type{StateLabel{N,A}}, ::Type{StateLabel{N,B}}) = StateLabel{N,label_promote(A,B)}

Base.convert{N,T}(::Type{StateLabel{N,T}}, s::StateLabel{N}) = StateLabel{N,T}(convert(Vector{T}, s.label))
Base.convert{N,T}(::Type{StateLabel{N,T}}, s::StateLabel{N,T}) = s

###########
# OuterLabel #
###########
immutable OuterLabel{N,K,B}
    k::StateLabel{N,K}
    b::StateLabel{N,B}
end

OuterLabel(o::OuterLabel) = o
OuterLabel(::StateLabel, ::StateLabel) = error("OuterLabel can only be constructed if both StateLabels have the same number of factors")
OuterLabel{N,K,B}(k::StateLabel{N,K}, b::StateLabel{N,B}) = OuterLabel{N,K,B}(k, b)
OuterLabel(k, b) = OuterLabel(StateLabel(k), StateLabel(b))

klabel(o::OuterLabel) = o.k
blabel(o::OuterLabel) = o.b

Base.copy(o::OuterLabel) = OuterLabel(o.k, o.b)
Base.hash(o::OuterLabel) = hash(o.k, hash(o.b))
Base.hash(o::OuterLabel, h::Uint64) = hash(hash(o), h)

Base.(:(==)){N}(a::OuterLabel{N}, b::OuterLabel{N}) = a.k == b.k && a.b == b.b

Base.length{N}(::OuterLabel{N}) = N

Base.ctranspose(o::OuterLabel) = OuterLabel(o.b, o.k)
tensor(o1::OuterLabel, o2::OuterLabel) = OuterLabel(tensor(o1.k, o2.k), tensor(o1.b, o2.b))

is_sum_x(o::OuterLabel, x) = sum(klabel(o))==sum(blabel(o))==x

ptranspose(o::OuterLabel, i) = OuterLabel(setindex(o.k, o.b[i], i), setindex(o.b, o.k[i], i))

permute(o::OuterLabel, perm::Vector) = OuterLabel(permute(o.k, perm), permute(o.b, perm))
except(o::OuterLabel, i) = OuterLabel(except(o.k, i), except(o.b, i))
switch(o::OuterLabel, i, j) = OuterLabel(switch(o.k, i, j), switch(o.b, i, j))

tensor_type{N,M,K1,K2,B1,B2}(::Type{OuterLabel{N,K1,B1}}, ::Type{OuterLabel{M,K2,B2}}) = OuterLabel{N+M,Any,Any}
tensor_type{N,M,K,B}(::Type{OuterLabel{N,K,B}}, ::Type{OuterLabel{M,K,B}}) = OuterLabel{N+M,K,B}
tensor_type{N,M,K1,K2,B}(::Type{OuterLabel{N,K1,B}}, ::Type{OuterLabel{M,K2,B}}) = OuterLabel{N+M,Any,B}
tensor_type{N,M,K,B1,B2}(::Type{OuterLabel{N,K,B1}}, ::Type{OuterLabel{M,K,B2}}) = OuterLabel{N+M,K,Any}

Base.promote_type{N,K1,K2,B1,B2}(::Type{OuterLabel{N,K1,B1}}, ::Type{OuterLabel{N,K2,B2}}) = OuterLabel{N,label_promote(K1,K2),label_promote(B1,B2)}

Base.convert{N,K,B}(::Type{OuterLabel{N,K,B}}, o::OuterLabel{N}) = OuterLabel{N,K,B}(convert(StateLabel{N,K}, o.k), convert(StateLabel{N,B}, o.b))
Base.convert{N,K,B}(::Type{OuterLabel{N,K,B}}, o::OuterLabel{N,K,B}) = o

####################
# Helper Functions #
####################
function switch!(arr, i, j)
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp
    return arr
end

export StateLabel,
    OuterLabel,
    klabel,
    blabel
