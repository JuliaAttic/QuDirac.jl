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

    ex = quote
        @inbounds result = $P(b[Val{1}], k[Val{1}])
    end

    for i in 2:nfactors(b)
        ex = quote 
            $ex
            @inbounds result *= $P(b[Val{$i}], k[Val{$i}])
        end
    end

    return quote
        $ex
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
