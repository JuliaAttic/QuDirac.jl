##############
# DiracLabel #
##############
# Types L <: DiracLabel should define:
#   nfactors(::Type{L})
#   hash(::L)
#   tensor(a::L, b::L)
#   conversion/promotion

abstract DiracLabel

nfactors(label::DiracLabel) = nfactors(typeof(label))

Base.(:*)(a::DiracLabel, b::DiracLabel) = tensor(a, b)

##############
# StateLabel #
##############
immutable StateLabel{L} <: DiracLabel
    items::L
    hsh::Uint64
    StateLabel(items::Tuple) = new(items, hash(items))
end

# Constructors #
#--------------#
StateLabel(items::Tuple) = StateLabel{typeof(items)}(items)
StateLabel(items...) = StateLabel(items)

# Generic functions #
#-------------------#
Base.copy(label::StateLabel) = label
Base.hash(label::StateLabel) = label.hsh
Base.(:(==))(a::StateLabel, b::StateLabel) = hash(a) == hash(b)

# generated so that nfactors is a compile-time constant based on L
@generated nfactors{L}(::Type{StateLabel{L}}) = length(L.parameters)

# Conversion/Promotion #
#----------------------#
Base.convert{L}(::Type{StateLabel{L}}, label::StateLabel) = StateLabel(convert(L, label.items))
Base.convert{L}(::Type{StateLabel{L}}, label::StateLabel{L}) = label

# Container-like functions #
#--------------------------#
Base.getindex{i}(label::StateLabel, ::Type{Val{i}}) = label.items[i]

Base.length(label::StateLabel) = nfactors(label)

Base.start(label::StateLabel) = start(label.items)
Base.next(label::StateLabel, i) = next(label.items, i)
Base.done(label::StateLabel, i) = done(label.items, i)

# Label Transformations #
#-----------------------#
@generated function except{i}(label::StateLabel, ::Type{Val{i}})
    items = Any[:(label.items[$x]) for x in 1:nfactors(label)]
    deleteat!(items, i)
    ex = Expr(:tuple, items...)    
    return :(StateLabel($ex))
end

@generated function setindex{i}(label::StateLabel, y, ::Type{Val{i}})
    items = Any[:(label.items[$x]) for x in 1:nfactors(label)]
    items[i] = :y
    ex = Expr(:tuple, items...)    
    return :(StateLabel($ex))
end

@generated function switch{i,j}(label::StateLabel, ::Type{Val{i}}, ::Type{Val{j}})
    items = Any[:(label.items[$x]) for x in 1:nfactors(label)]
    tmp = items[i]
    items[i] = items[j]
    items[j] = tmp
    ex = Expr(:tuple, items...)    
    return :(StateLabel($ex))
end

@generated function permute{T<:Tuple}(label::StateLabel, ::Type{T})
    items = Any[:(label.items[$x]) for x in T.parameters]
    ex = Expr(:tuple, items...)    
    return :(StateLabel($ex))
end

@generated function tensor(a::StateLabel, b::StateLabel)
    as = Any[:(a.items[$x]) for x in 1:nfactors(a)]
    bs = Any[:(b.items[$x]) for x in 1:nfactors(b)]
    ex = Expr(:tuple, as..., bs...)
    return :(StateLabel($ex))
end

##############
# OuterLabel #
##############
immutable OuterLabel{K,B} <: DiracLabel
    k::StateLabel{K}
    b::StateLabel{B}
    hsh::Uint64
    function OuterLabel(k, b)
        @assert nfactors(K) == nfactors(B)
        return new(k, b, hash(hash(k), hash(b)))
    end
end

# Constructors #
#--------------#
OuterLabel{K,B}(k::StateLabel{K}, b::StateLabel{B}) = OuterLabel{K,B}(k, b)
OuterLabel(k, b) = OuterLabel(StateLabel(k), StateLabel(b))

# Generic functions #
#-------------------#
Base.copy(o::OuterLabel) = o
Base.hash(o::OuterLabel) = o.hsh
Base.(:(==))(a::OuterLabel, b::OuterLabel) = hash(a) == hash(b)

nfactors{K,B}(::Type{OuterLabel{K,B}}) = length(K.parameters)
ketlabel(o::OuterLabel) = o.k
bralabel(o::OuterLabel) = o.b

# Conversion/Promotion #
#----------------------#
Base.convert{K,B}(::Type{OuterLabel{K,B}}, o::OuterLabel) = OuterLabel(convert(StateLabel{K}, ketlabel(o)), convert(StateLabel{B}, bralabel(o)))
Base.convert{K,B}(::Type{OuterLabel{K,B}}, o::OuterLabel{K,B}) = o

# Label Transformations #
#-----------------------#
tensor(a::OuterLabel, b::OuterLabel) = OuterLabel(tensor(ketlabel(a), ketlabel(b)), tensor(bralabel(a), bralabel(b)))

Base.ctranspose(o::OuterLabel) = OuterLabel(bralabel(o), ketlabel(o))

function ptranspose{i}(o::OuterLabel, V::Type{Val{i}})
    return OuterLabel(setindex(ketlabel(o), bralabel(o)[V], V), setindex(bralabel(o), ketlabel(o)[V], V))
end

permute{T<:Tuple}(o::OuterLabel, P::Type{T}) = OuterLabel(permute(ketlabel(o), P), permute(bralabel(o), P))
except{i}(o::OuterLabel, V::Type{Val{i}}) = OuterLabel(except(ketlabel(o), V), except(bralabel(o), V))
switch{i,j}(o::OuterLabel, A::Type{Val{i}}, B::Type{Val{j}}) = OuterLabel(switch(ketlabel(o), A, B), switch(bralabel(o), A, B))

export StateLabel,
       OuterLabel,
       ketlabel,
       bralabel,
       tensor,
       switch,
       except,
       setindex,
       permute
