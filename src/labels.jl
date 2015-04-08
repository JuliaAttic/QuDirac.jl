##############
# StateLabel #
##############
type StateLabel{N}
    label::NTuple{N}
end

StateLabel(i...) = StateLabel(i)

Base.getindex(s::StateLabel, i) = s.label[i]
Base.getindex(s::StateLabel, arr::AbstractArray) = StateLabel(s.label[arr])

Base.copy(s::StateLabel) = StateLabel(s.label)

Base.start(s::StateLabel) = start(s.label)
Base.next(s::StateLabel, i) = next(s.label, i)
Base.done(s::StateLabel, i) = done(s.label, i)

Base.first(s::StateLabel) = first(s.label)
Base.last(s::StateLabel) = last(s.label)
Base.endof(s::StateLabel) = endof(s.label)

Base.length(s::StateLabel) = length(s.label)

QuBase.tensor(s::StateLabel...) = StateLabel(apply(tuple, s...))
bob(a::StateLabel, b::StateLabel) = StateLabel(a..., b...)

##############
# OuterLabel #
##############
type OuterLabel{N}
    ktlabel::StateLabel{N}
    brlabel::StateLabel{N}
end

OuterLabel(k, b) = OuterLabel(StateLabel(k), StateLabel(b))

brlabel(o::OuterLabel) = i.brlabel
ktlabel(o::OuterLabel) = i.ktlabel

##############
# InnerLabel #
##############
type InnerLabel{P<:AbstractInner,N} <: DiracScalar
    ptype::P
    brlabel::StateLabel{N}
    ktlabel::StateLabel{N}
end

InnerLabel(ptype, b, k) = InnerLabel(ptype, StateLabel(k), StateLabel(br))

brlabel(i::InnerLabel) = i.brlabel
ktlabel(i::InnerLabel) = i.ktlabel

####################
# Helper Functions #
####################
isx(label::(Tuple,Tuple), x) = sum(label[1])==sum(label[2])==x
isx(label::Tuple, x) = sum(label) == x
ctpair(k,v) = (reverse(k), v')
nzcoeff(k,v) = v!=0
second(t) = t[2]

except(tup, i) = tuple(tup[1:i-1]..., tup[i+1:end]...)
function switch!(arr, i, j)
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp
    return arr
end

switch(tup, i, j) = tup[switch!([1:length(tup)], i, j)]
permute(tup, perm::Vector) = tup[perm]
placeat(tup, x, y) = tuple(tup[1:y-1]..., x, tup[y+1:end]...)
join_tup(a, b) = tuple(a..., b...)

