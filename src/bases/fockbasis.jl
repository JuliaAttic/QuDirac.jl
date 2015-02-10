#############
# FockBasis #
#############
    # A FockBasis uses precomputed values to efficiently 
    # generate labels for given indices in the basis, or
    # vice versa (an index for a given state in the basis).
    # For example:
    #
    #   julia> f=FockBasis(2,2,2)
    #    FockBasis{AbstractStructure,3}(0:2,0:2,0:2)    
    #   
    #   julia> labelvec(f)
    #    8-element Array{(Any...,),1}:
    #    (0,0,0)
    #    (1,0,0)
    #    (0,1,0)
    #    (1,1,0)
    #    (0,0,1)
    #    (1,0,1)
    #    (0,1,1)
    #    (1,1,1)
    #
    #   julia> f[6]
    #    (1,0,1)
    #
    #   julia> f[(1,0,1)]
    #    6
    # 
    # Because the labels are generated rather than actually 
    # stored, one can represent very large bases without 
    # any storage overhead:
    #
    #   julia> f = FockBasis(221,135,31,42,321,3)
    #   FockBasis{AbstractStructure,6}(0:221,0:135,0:31,0:42,0:321,0:3)
    #   
    #   julia> length(f)
    #   53508919296
    #   
    #   julia> last(f)
    #   (221,135,31,42,321,3)
    #   
    #   julia> f[34234134]
    #   (179,119,13,35,0,0)
    #   
    #   julia> f[(179,119,13,35,0,0)]
    #   34234134
    # 
    # Arbitrary numeric ranges are supported for labels:
    #   
    #   julia> f=FockBasis(0.0:0.1:0.2, 4:7)
    #   FockBasis{AbstractStructure,2}(0.0:0.1:0.2,4:7)     
    #       
    #   julia> collect(f)
    #   12-element Array{(Any...,),1}:
    #    (0.0,4)
    #    (0.1,4)
    #    (0.2,4)
    #    (0.0,5)
    #    (0.1,5)
    #    (0.2,5)
    #    (0.0,6)
    #    (0.1,6)
    #    (0.2,6)
    #    (0.0,7)
    #    (0.1,7)
    #    (0.2,7)        
    #
    #   julia> f[f[(0.0, 7)]] == (0.0, 7)
    #   true
    
    immutable FockBasis{S<:AbstractStructure,N} <: AbstractFiniteBasis{S}
        ranges::NTuple{N,Range}
        denoms::NTuple{N,Float64}
        FockBasis(ranges, denoms, ::Type{BypassFlag}) = new(ranges, denoms)
        FockBasis(::()) = error("")
        FockBasis() = error("")
        function FockBasis(ranges::NTuple{N,Range})
            # reverse is done to match cartesianmap order
            return FockBasis{S,N}(ranges, precompute_denoms(reverse(map(length,ranges))), BypassFlag) 
        end

    end
    
    FockBasis{N,S<:AbstractStructure}(::Type{S}, lens::NTuple{N,Range}) = FockBasis{S,N}(lens)
    FockBasis{S<:AbstractStructure}(::Type{S}, lens::Tuple) = FockBasis(S, map(torange, lens))
    FockBasis(lens::Tuple) = FockBasis(AbstractStructure, lens)
    FockBasis(lens...) = FockBasis(AbstractStructure, lens)

    Base.convert{A,B,N}(::Type{FockBasis{A,N}}, f::FockBasis{B,N}) = FockBasis{A,N}(f.ranges, f.denoms, BypassFlag)
    Base.convert{A,B,N}(::Type{FiniteBasis{A,N}}, f::FockBasis{B,N}) = FiniteBasis{A,N}(size(f))

    Base.copy{S,N}(f::FockBasis{S,N}) = FockBasis{S,N}(copy(f.ranges), copy(f.denoms), BypassFlag)

    ####################
    # Helper Functions #
    ####################
    # This function precomputes the 
    # denominators for each factor of 
    # the cartesian product.
    #
    # This site offers a thorough 
    # explanation of this method:
    # http://phrogz.net/lazy-cartesian-product
    #
    function precompute_denoms(lens)
        # storing as Floats avoids number precision issues for 
        # outrageously large bases
        total_divisor = prod(float,lens) 
        function get_denom(i)
            total_divisor = div(total_divisor, i)
            return max(1.0, total_divisor)
        end
        # reverse is done to match cartesianmap order
        return reverse(map(get_denom, lens))
    end

    torange(n::Number) = zero(eltype(n)):(n)
    torange(r::Range) = r

    ######################
    # Property Functions #
    ######################
    Base.size(f::FockBasis) = map(length, ranges(f))
    Base.size(f::FockBasis, i) = length(ranges(f, i))
    Base.length(f::FockBasis) = prod(length, ranges(f))
    Base.ndims(f::FockBasis) = nfactors(f)

    QuBase.structure{S}(::Type{FockBasis{S}}) = S
    QuBase.structure{S,N}(::Type{FockBasis{S,N}}) = S
    QuBase.nfactors{S,N}(::FockBasis{S,N}) = N
    QuBase.checkcoeffs(coeffs, dim, f::FockBasis) = size(coeffs, dim) == length(f)

    ranges(f::FockBasis) = f.ranges
    ranges(f::FockBasis, i) = f.ranges[i]

    ######################
    # Accessor Functions #
    ######################
    ind_value(n, range, denom, modulus) = range[(div(n, denom) % modulus)+1]
    tuple_at_ind{S,N}(f::FockBasis{S,N}, i) = ntuple(N, x->ind_value(i-1, ranges(f,x), f.denoms[x], size(f,x)))
    pos_in_range(r::Range, i) = i in r ? (i-first(r))/step(r) : throw(BoundsError())
    
    getpos(f::FockBasis, label) = int(sum(map(*, map(pos_in_range, ranges(f), label), f.denoms))+1)

    Base.in(label, f::FockBasis) = reduce(&, map(in, label, ranges(f)))

    Base.getindex(f::FockBasis, i) = tuple_at_ind(f, i)
    Base.getindex(f::FockBasis, t::Tuple) = getpos(f, t)
    Base.getindex(f::FockBasis, arr::AbstractArray) = [f[i] for i in arr]

    ######################
    # Iterator Functions #
    ######################
    Base.start(::FockBasis) = 1
    Base.done(f::FockBasis, state) = length(f) == state-1
    Base.next(f::FockBasis, state) = f[state], state+1
    Base.endof(f::FockBasis) = length(f)
    Base.last(f::FockBasis) = f[length(f)]
    Base.first(f::FockBasis) = f[1]
    Base.collect(f::FockBasis) = f[1:end]

    ##########################
    # Mathematical Functions #
    ##########################
    QuBase.tensor{S}(a::FockBasis{S}, b::FockBasis{S}) = FockBasis(S, tensor_tup(a.ranges, b.ranges))

    ######################
    # Printing Functions #
    ######################
    Base.repr(f::FockBasis) = "$(typeof(f))$(ranges(f))"
    Base.show(io::IO, f::FockBasis) = print(io, repr(f))

export FockBasis
