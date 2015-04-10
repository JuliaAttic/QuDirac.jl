##############
# StateLabel #
##############
immutable StateLabel{N}
    label::NTuple{N}
    hash::Uint64
    StateLabel(label::NTuple{N}) = new(label, hash(label))
    StateLabel(label::NTuple{N}, h::Uint64) = new(label, h)
end

StateLabel(label::StateLabel) = label
StateLabel{N}(label::NTuple{N}) = StateLabel{N}(label)
StateLabel(i...) = StateLabel(i)

Base.getindex(s::StateLabel, i) = s.label[i]
Base.getindex(s::StateLabel, arr::AbstractArray) = s.label[arr]

Base.copy{N}(s::StateLabel{N}) = StateLabel{N}(s.label, s.hash)
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

Base.map(f::Union(Function,DataType), s::StateLabel) = StateLabel(map(f, s.label))

QuBase.tensor(s::StateLabel...) = StateLabel(apply(tuple, s...))

labelstr(s::StateLabel) = join(map(repr, s.label), ',')
Base.repr(s::StateLabel) = repr(typeof(s)) * "(" * labelstr(s) * ")"
Base.show(io::IO, s::StateLabel) = print(io, repr(s))

##############
# OuterLabel #
##############
immutable OuterLabel{N}
    k::StateLabel{N}
    b::StateLabel{N}
    hash::Uint64
    OuterLabel(k::StateLabel{N}, b::StateLabel{N}) = new(k, b, hash(k, hash(b)))
    OuterLabel(k::StateLabel{N}, b::StateLabel{N}, h::Uint64) = new(k, b, h)
end

OuterLabel(::StateLabel, ::StateLabel) = error("OuterLabel can only be constructed if both StateLabels have the same number of factors")
OuterLabel{N}(k::StateLabel{N}, b::StateLabel{N}) = OuterLabel{N}(k, b)
OuterLabel(k, b) = OuterLabel(StateLabel(k), StateLabel(b))

klabel(o::OuterLabel) = o.k
blabel(o::OuterLabel) = o.b

Base.copy{N}(o::OuterLabel{N}) = OuterLabel{N}(o.k, o.b, o.hash)
Base.hash(o::OuterLabel) = o.hash
Base.hash(o::OuterLabel, h::Uint64) = hash(hash(o), h)

Base.(:(==)){N}(a::OuterLabel{N}, b::OuterLabel{N}) = hash(a) == hash(b)

Base.length{N}(::OuterLabel{N}) = N

Base.reverse(o::OuterLabel) = OuterLabel(o.b, o.k)
QuBase.tensor(o1::OuterLabel, o2::OuterLabel) = OuterLabel(tensor(o1.k, o2.k), tensor(o1.b, o2.b))

Base.repr(o::OuterLabel) = repr(typeof(o)) * "(" * ktstr(o.k) * "," * brstr(o.b) * ")"
Base.show(io::IO, o::OuterLabel) = print(io, repr(o))

####################
# Helper Functions #
####################
is_sum_x(o::OuterLabel, x) = sum(klabel(o))==sum(blabel(o))==x
is_sum_x(s::StateLabel, x) = sum(s) == x

ctpair(k,v) = (reverse(k), v')
nzcoeff(k,v) = v!=0
second(t) = t[2]
except(label::StateLabel, i) = StateLabel(label[1:i-1]..., label[i+1:end]...)
setindex(label::StateLabel, x, y) = StateLabel(label[1:y-1]..., x, label[y+1:end]...)
permute{N}(label::StateLabel{N}, perm::Vector) = StateLabel{N}(label[perm])

function switch!(arr, i, j)
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp
    return arr
end
switch{N}(label::StateLabel{N}, i, j) =  StateLabel{N}(label[switch!([1:N], i, j)])

export StateLabel,
    OuterLabel,
    klabel,
    blabel

