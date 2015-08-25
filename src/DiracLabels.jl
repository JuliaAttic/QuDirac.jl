##############
# DiracLabel #
##############
# Types L <: DiracLabel should define:
#   nfactors(::Type{L})
#   hash(::L)
#   tensor(a::L, b::L)
#   conversion/promotion

abstract DiracLabel{N}

nfactors{N}(::DiracLabel{N}) = N

Base.(:*)(a::DiracLabel, b::DiracLabel) = tensor(a, b)

##############
# StateLabel #
##############
immutable StateLabel{N,T} <: DiracLabel{N}
    items::T
    hsh::Uint64
    StateLabel(items::NTuple{N}) = new(items, hash(items))
end

# Constructors #
#--------------#
StateLabel{N}(items::NTuple{N}) = StateLabel{N,typeof(items)}(items)
StateLabel(items...) = StateLabel(items)

# Generic functions #
#-------------------#
Base.copy(label::StateLabel) = label
Base.hash(label::StateLabel) = label.hsh
Base.(:(==))(a::StateLabel, b::StateLabel) = hash(a) == hash(b)

nfactors{N,T}(::Type{StateLabel{N,T}}) = N

# Conversion/Promotion #
#----------------------#
Base.promote_type{N,A,B}(::Type{StateLabel{N,A}}, ::Type{StateLabel{N,B}}) = StateLabel{N,promote_type(A,B)}

Base.convert{N,T}(::Type{StateLabel{N,T}}, label::StateLabel{N}) = StateLabel(convert(T, label.items))
Base.convert{N,T}(::Type{StateLabel{N,T}}, label::StateLabel{N,T}) = label

# Container-like functions #
#--------------------------#
Base.getindex{i}(label::StateLabel, ::Type{Val{i}}) = label.items[i]

Base.length(label::StateLabel) = nfactors(label)

Base.start(label::StateLabel) = start(label.items)
Base.next(label::StateLabel, i) = next(label.items, i)
Base.done(label::StateLabel, i) = done(label.items, i)

Base.first(label::StateLabel) = first(label.items)
Base.last(label::StateLabel) = last(label.items)
Base.endof(label::StateLabel) = Val{nfactors(label)}

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
immutable OuterLabel{N,K,B} <: DiracLabel{N}
    k::StateLabel{N,K}
    b::StateLabel{N,B}
    hsh::Uint64
    OuterLabel(k,b) = new(k, b, hash(hash(k), hash(b)))
end

# Constructors #
#--------------#
OuterLabel{N,K,B}(k::StateLabel{N,K}, b::StateLabel{N,B}) = OuterLabel{N,K,B}(k, b)
OuterLabel{N,M}(k::StateLabel{N}, b::StateLabel{M}) = error("OuterLabel must obey nfactors(bralabel(o)) == nfactors(ketlabel(o))")
OuterLabel(k, b) = OuterLabel(StateLabel(k), StateLabel(b))


# Generic functions #
#-------------------#
Base.copy(o::OuterLabel) = o
Base.hash(o::OuterLabel) = o.hsh
Base.(:(==))(a::OuterLabel, b::OuterLabel) = hash(a) == hash(b)

nfactors{N,K,B}(::Type{OuterLabel{N,K,B}}) = N
ketlabel(o::OuterLabel) = o.k
bralabel(o::OuterLabel) = o.b

# Conversion/Promotion #
#----------------------#
Base.promote_type{N,K1,K2,B1,B2}(::Type{OuterLabel{K1,B1}}, ::Type{OuterLabel{K2,B2}}) = OuterLabel{promote_type(K1,K2),promote_type(B1,B2)}

Base.convert{N,K,B}(::Type{OuterLabel{N,K,B}}, o::OuterLabel) = OuterLabel{N,K,B}(convert(StateLabel{K}, ketlabel(o)), convert(StateLabel{B}, bralabel(o)))
Base.convert{N,K,B}(::Type{OuterLabel{N,K,B}}, o::OuterLabel{N,K,B}) = o

# Label Transformations #
#-----------------------#
tensor(a::OuterLabel, b::OuterLabel) = OuterLabel(tensor(ketlabel(a), ketlabel(b)), tensor(bralabel(a), bralabel(b)))

Base.ctranspose(o::OuterLabel) = OuterLabel(bralabel(o), ketlabel(o))

function ptranspose{i}(o::OuterLabel, V::Type{Val{i}})
    return OuterLabel(setindex(ketlabel(o), bralabel(o)[V], V), setindex(bralabel(o), ketlabel(o)[V], V))
end

permute{T<:Tuple}(o::OuterLabel, P::Type{T}) = OuterLabel(permute(ketlabel(o), P), permute(bralabel(o), P))
except{i}(o::OuterLabel, i) = OuterLabel(except(ketlabel(o), i), except(bralabel(o), i))
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
