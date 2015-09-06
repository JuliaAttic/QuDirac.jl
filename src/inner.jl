#########
# inner #
#########
macro newinner(name)
    return quote 
        immutable $name <: AbstractInner end
    end
end

@generated function inner{P<:AbstractInner}(::Type{P}, b::StateLabel, k::StateLabel)
    @assert nfactors(b) == nfactors(k)
    return quote
        @inbounds result = P(b[1], k[1])
        @inbounds for i=2:nfactors(b)
            result *= P(b[i], k[i])
        end
        return result
    end
end

#############
# KronDelta #
#############
@newinner KronDelta
KronDelta(b, k) = b == k

rettype(::Type{KronDelta}) = Bool

# we can hack in this optimization for KronDelta
@generated function inner(::Type{KronDelta}, b::StateLabel, k::StateLabel) 
    @assert nfactors(b) == nfactors(k)
    return :(b == k)
end


export @newinner, KronDelta
