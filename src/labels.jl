##############
# StateLabel #
##############
immutable StateLabel{T}
    label::Vector{T}
    hash::Uint64
    StateLabel(label::Vector{T}) = new(label, hash(label))
end

StateLabel(s::StateLabel) = s
StateLabel{T}(v::Vector{T}) = StateLabel{T}(v)
StateLabel(t::Tuple) = StateLabel(t...)
StateLabel(i...) = StateLabel(vcat(i...))

Base.length(s::StateLabel) = length(s.label)
Base.eltype{T}(::Type{StateLabel{T}}) = T
Base.eltype(s::StateLabel) = eltype(typeof(s))
Base.getindex(s::StateLabel, i) = s.label[i]
Base.getindex(s::StateLabel, arr::AbstractArray) = s.label[arr]

nfactors(s::StateLabel) = length(s)
except{T}(s::StateLabel{T}, i) = StateLabel(deleteat!(copy(s.label), i))
setindex{T}(s::StateLabel{T}, x, y) = StateLabel(vcat(s.label[1:y-1], x, s.label[y+1:end]))
setindex{T}(s::StateLabel{T}, x::T, y) = StateLabel(setindex!(copy(s.label), x, y))
switch{T}(s::StateLabel{T}, i, j) =  StateLabel(switch!(copy(s.label), i, j))
permute{T}(label::StateLabel{T}, perm::Vector) = StateLabel(label[perm])

Base.copy(s::StateLabel) = s
Base.hash(s::StateLabel) = s.hash
Base.(:(==))(a::StateLabel, b::StateLabel) = hash(a) == hash(b)

Base.start(s::StateLabel) = start(s.label)
Base.next(s::StateLabel, i) = next(s.label, i)
Base.done(s::StateLabel, i) = done(s.label, i)

Base.first(s::StateLabel) = first(s.label)
Base.last(s::StateLabel) = last(s.label)
Base.endof(s::StateLabel) = endof(s.label)

is_sum_x(s::StateLabel, x) = sum(s) == x

tensor(a::StateLabel, b::StateLabel) = StateLabel(vcat(a.label, b.label))

Base.(:*)(a::StateLabel, b::StateLabel) = tensor(a,b)

Base.promote_type{A,B}(::Type{StateLabel{A}}, ::Type{StateLabel{B}}) = StateLabel{promote_type(A,B)}

Base.convert{T}(::Type{StateLabel{T}}, s::StateLabel) = StateLabel(convert(Vector{T}, s.label))
Base.convert{T}(::Type{StateLabel{T}}, s::StateLabel{T}) = s

##############
# OuterLabel #
##############
immutable OuterLabel{K,B}
    k::StateLabel{K}
    b::StateLabel{B}
    function OuterLabel(k::StateLabel{K}, b::StateLabel{B})
        @assert nfactors(k) == nfactors(b)
        return new(k, b)
    end
end

OuterLabel(o::OuterLabel) = o
OuterLabel{K,B}(k::StateLabel{K}, b::StateLabel{B}) = OuterLabel{K,B}(k,b)
OuterLabel(k, b) = OuterLabel(StateLabel(k), StateLabel(b))

klabel(o::OuterLabel) = o.k
blabel(o::OuterLabel) = o.b

Base.copy(o::OuterLabel) = o
Base.hash(o::OuterLabel) = hash(hash(o.k), hash(o.b))
Base.(:(==))(a::OuterLabel, b::OuterLabel) = a.k == b.k && a.b == b.b

nfactors(o::OuterLabel) = nfactors(o.k)

Base.ctranspose(o::OuterLabel) = OuterLabel(o.b, o.k)
tensor(o1::OuterLabel, o2::OuterLabel) = OuterLabel(tensor(o1.k, o2.k), tensor(o1.b, o2.b))

is_sum_x(o::OuterLabel, x) = sum(klabel(o))==sum(blabel(o))==x

ptranspose(o::OuterLabel, i) = OuterLabel(setindex(o.k, o.b[i], i), setindex(o.b, o.k[i], i))

permute(o::OuterLabel, perm::Vector) = OuterLabel(permute(o.k, perm), permute(o.b, perm))
except(o::OuterLabel, i) = OuterLabel(except(o.k, i), except(o.b, i))
switch(o::OuterLabel, i, j) = OuterLabel(switch(o.k, i, j), switch(o.b, i, j))

Base.promote_type{K1,K2,B1,B2}(::Type{OuterLabel{K1,B1}}, ::Type{OuterLabel{K2,B2}}) = OuterLabel{promote_type(K1,K2),promote_type(B1,B2)}

Base.convert{K,B}(::Type{OuterLabel{K,B}}, o::OuterLabel) = OuterLabel{K,B}(convert(StateLabel{K}, o.k), convert(StateLabel{B}, o.b))
Base.convert{K,B}(::Type{OuterLabel{K,B}}, o::OuterLabel{K,B}) = o

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
