##############
# StateLabel #
##############
immutable StateLabel{N,T}
    label::Vector{T}
    StateLabel(label::Vector{T}) = new(label)
end

StateLabel(s::StateLabel) = s
StateLabel(label::NTuple{0}) = error("Cannot have a 0-factor StateLabel")
StateLabel{N,T}(label::NTuple{N,T}) = StateLabel{N,T}(collect(label))
StateLabel{N}(label::NTuple{N}) = StateLabel{N,Any}(collect(label))
StateLabel(i...) = StateLabel(i)

nfactors{N,T}(::Type{StateLabel{N,T}}) = N
Base.getindex(s::StateLabel, i) = s.label[i]
Base.getindex(s::StateLabel, arr::AbstractArray) = s.label[arr]

except{N,T}(s::StateLabel{N,T}, i) = StateLabel{N-1,T}(deleteat!(copy(s.label), i))
setindex{N,T}(s::StateLabel{N,T}, x, y) = StateLabel{N,T}(setindex!(copy(s.label), i))
switch{N,T}(s::StateLabel{N,T}, i, j) =  StateLabel{N,T}(switch!(copy(s.label), i, j))
permute{N,T}(label::StateLabel{N,T}, perm::Vector) = StateLabel{N,T}(label[perm])

Base.copy{N,T}(s::StateLabel{N,T}) = StateLabel{N,T}(s.label)
Base.hash(s::StateLabel) = hash(s.label)
Base.hash(s::StateLabel, h::Uint64) = hash(hash(s), h)

Base.(:(==)){N}(a::StateLabel{N},b::StateLabel{N}) = a.label == b.label

Base.start(s::StateLabel) = start(s.label)
Base.next(s::StateLabel, i) = next(s.label, i)
Base.done(s::StateLabel, i) = done(s.label, i)

Base.first(s::StateLabel) = first(s.label)
Base.last(s::StateLabel) = last(s.label)
Base.endof(s::StateLabel) = endof(s.label)

Base.length{N}(::StateLabel{N}) = N

is_sum_x(s::StateLabel, x) = sum(s) == x

tensor{N,M,A,B}(a::StateLabel{N,A}, b::StateLabel{M,B}) = StateLabel{N+M,promote_type(A,B)}(vcat(a.label, b.label))
Base.(:*)(a::StateLabel, b::StateLabel) = tensor(a,b)

labelstr(s::StateLabel) = join(map(repr, s.label), ',')
Base.repr(s::StateLabel) = repr(typeof(s)) * "(" * labelstr(s) * ")"
Base.show(io::IO, s::StateLabel) = print(io, repr(s))

tensor_type{N,M,A,B}(::Type{StateLabel{N,A}}, ::Type{StateLabel{M,B}}) = StateLabel{N+M,promote_type(A,B)}
Base.promote_type{N,A,B}(::Type{StateLabel{N,A}}, ::Type{StateLabel{N,B}}) = StateLabel{N,promote_type(A,B)}
Base.convert{N,T}(::Type{StateLabel{N,T}}, s::StateLabel{N}) = StateLabel{N,T}(convert(Vector{T}, s.label))

###########
# OpLabel #
###########
immutable OpLabel{N,K,B}
    k::StateLabel{N,K}
    b::StateLabel{N,B}
end

OpLabel(op::OpLabel) = op
OpLabel(::StateLabel, ::StateLabel) = error("OpLabel can only be constructed if both StateLabels have the same number of factors")
OpLabel{N,K,B}(k::StateLabel{N,K}, b::StateLabel{N,B}) = OpLabel{N,K,B}(k, b)
OpLabel(k, b) = OpLabel(StateLabel(k), StateLabel(b))

klabel(o::OpLabel) = o.k
blabel(o::OpLabel) = o.b

Base.copy(o::OpLabel) = OpLabel(o.k, o.b)
Base.hash(o::OpLabel) = hash(o.k, hash(o.b))
Base.hash(o::OpLabel, h::Uint64) = hash(hash(o), h)

Base.(:(==)){N}(o1::OpLabel{N}, o2::OpLabel{N}) = o1.k == o2.k && o1.b == o2.b

Base.length{N}(::OpLabel{N}) = N

Base.ctranspose(o::OpLabel) = OpLabel(o.b, o.k)
tensor(o1::OpLabel, o2::OpLabel) = OpLabel(tensor(o1.k, o2.k), tensor(o1.b, o2.b))

is_sum_x(o::OpLabel, x) = sum(klabel(o))==sum(blabel(o))==x

ptranspose{N}(k::StateLabel{N}, b::StateLabel{N}, i) = OpLabel(setindex(k, b[i], i), setindex(b, k[i], i))
ptranspose(o::OpLabel, i) = ptranspose(o.k, o.b, i)
ptranspose_dual(o::OpLabel, i) = ptranspose(o.b, o.k, i)

permute(o::OpLabel, perm::Vector) = OpLabel(permute(o.k, perm), permute(o.b, perm))
except(o::OpLabel, i) = OpLabel(except(o.k, i), except(o.b, i))
switch(o::OpLabel, i, j) = OpLabel(switch(o.k, i, j), switch(o.b, i, j))

Base.repr(o::OpLabel) = repr(typeof(o)) * "(" * ktstr(o.k) * "," * brstr(o.b) * ")"
Base.show(io::IO, o::OpLabel) = print(io, repr(o))

tensor_type{N,M,K1,K2,B1,B2}(::Type{OpLabel{N,K1,B1}}, ::Type{OpLabel{M,K2,B2}}) = OpLabel{N+M,promote_type(K1,K2), promote_type(B1,B2)}
Base.promote_type{N,K1,K2,B1,B2}(::Type{OpLabel{N,K1,B1}}, ::Type{OpLabel{N,K2,B2}}) = OpLabel{N,promote_type(K1,K2), promote_type(B1,B2)}
Base.convert{N,K,B}(::Type{OpLabel{N,K,B}}, s::OpLabel{N}) = OpLabel{N,K,B}(convert(StateLabel{N,K}, o.k), convert(StateLabel{N,B}, o.b))

####################
# Helper Functions #
####################
ctpair(k,v) = (k', v')

function switch!(arr, i, j)
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp
    return arr
end

export StateLabel,
    OpLabel,
    klabel,
    blabel
