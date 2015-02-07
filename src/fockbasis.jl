import Base: getindex,
    length,
    size,
    ndims,
    start,
    done,
    next,
    endof,
    last,
    first,
    collect,
    ctranspose,
    repr,
    show,
    in,
    copy

#############
# FockBasis #
#############
    # A FockBasis uses precomputed values to efficiently 
    # generate StateLabels for given indices in the basis, or
    # vice versa (an index for a given state in the basis).
    # For example:
    #
    #   julia> f=FockBasis(2,2,2)
    #    FockBasis{AbstractStructure,3}(2,2,2)    
    #   
    #   julia> labelvec(f)
    #    8-element Array{StateLabel{3},1}:
    #    StateLabel(0,0,0)
    #    StateLabel(1,0,0)
    #    StateLabel(0,1,0)
    #    StateLabel(1,1,0)
    #    StateLabel(0,0,1)
    #    StateLabel(1,0,1)
    #    StateLabel(0,1,1)
    #    StateLabel(1,1,1)
    #
    #   julia> f[6]
    #    StateLabel(1,0,1)
    #
    #   julia> getpos(f,StateLabel(1,0,1))
    #    6
    # 
    # Because the labels are generated rather than actually 
    # stored, one can represent very large bases without 
    # any storage overhead:
    #
    #   julia> f = FockBasis(221,135,31,42,321,3)
    #    FockBasis{AbstractStructure,6}(221,135,31,42,321,3)
    #
    #   julia> length(f)
    #    37407898710
    #
    #   julia> last(f)
    #    StateLabel(220,134,30,41,320,2)
    #
    #   julia> f[34234134]
    #    StateLabel(128,60,0,37,0,0)
    #
    #   julia> getpos(f, StateLabel(128,60,0,37,0,0))
    #    34234134
    # 
    # Arbitrary numeric ranges are supported for labels:
    #   
    #   julia> f=FockBasis(0.0:0.1:0.2, 4:7)
    #   FockBasis{AbstractStructure,2}(0.0:0.1:0.2,4:7)     
    #       
    #   julia> collect(f)
    #   12-element Array{StateLabel{2},1}:
    #    StateLabel(0.0,4)
    #    StateLabel(0.1,4)
    #    StateLabel(0.2,4)
    #    StateLabel(0.0,5)
    #    StateLabel(0.1,5)
    #    StateLabel(0.2,5)
    #    StateLabel(0.0,6)
    #    StateLabel(0.1,6)
    #    StateLabel(0.2,6)
    #    StateLabel(0.0,7)
    #    StateLabel(0.1,7)
    #    StateLabel(0.2,7)        
    #
    #   julia> f[getpos(f,(0.0, 7))] == StateLabel(0.0, 7)
    #   true
    
    immutable FockBasis{S<:AbstractStructure,N} <: AbstractLabelBasis{S,N}
        ranges::NTuple{N, Range}
        denoms::NTuple{N, Float64}
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

    convert{A,B,N}(::Type{FockBasis{A,N}}, f::FockBasis{B,N}) = FockBasis{A,N}(f.ranges, f.denoms, BypassFlag)
    convert{A,B,N}(::Type{FiniteBasis{A,N}}, f::FockBasis{B,N}) = FiniteBasis{A,N}(size(f))

    copy{S,N}(f::FockBasis{S,N}) = FockBasis{S,N}(copy(f.ranges), copy(f.denoms), BypassFlag)

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
    structure{S}(::Type{FockBasis{S}}) = S
    structure{S,N}(::Type{FockBasis{S,N}}) = S

    labelvec(f::FockBasis) = collect(f)
    ranges(f::FockBasis) = f.ranges
    ranges(f::FockBasis, i) = f.ranges[i]

    size(f::FockBasis) = map(length, ranges(f))
    size(f::FockBasis, i) = length(ranges(f, i))
    length(f::FockBasis) = prod(length, ranges(f))
    nfactors{S,N}(::FockBasis{S,N}) = N
    ndims(f::FockBasis) = nfactors(f)

    samelabels(a::FockBasis, b::FockBasis) = ranges(a) == ranges(b)
    checkcoeffs(coeffs, dim, f::FockBasis) = size(coeffs, dim) == length(f)

    ######################
    # Accessor Functions #
    ######################
    ind_value(n, range, denom, modulus) = range[(div(n, denom) % modulus)+1]
    tuple_at_ind(f::FockBasis, i) = ntuple(nfactors(f), x->ind_value(i-1, ranges(f,x), f.denoms[x], size(f,x)))
    pos_in_range(r::Range, i) = i in r ? (i-first(r))/step(r) : throw(BoundsError())

    in(label, f::FockBasis) = reduce(&, map(in, label, ranges(f)))
    
    getpos(f::FockBasis, s::AbstractState) = getpos(f, label(s))
    getpos(f::FockBasis, label) = int(sum(map(*, map(pos_in_range, ranges(f), label), f.denoms)))+1

    getindex{S,N}(f::FockBasis{S,N}, i) = StateLabel{N}(tuple_at_ind(f, i))
    getindex(f::FockBasis, arr::AbstractArray) = [f[i] for i in arr]
    getindex(f::FockBasis, t::Tuple) = getpos(f, t)

    ######################
    # Iterator Functions #
    ######################
    start(::FockBasis) = 1
    done(f::FockBasis, state) = length(f) == state-1
    next(f::FockBasis, state) = f[state], state+1
    endof(f::FockBasis) = length(f)
    last(f::FockBasis) = f[length(f)]
    first(f::FockBasis) = f[1]
    collect(f::FockBasis) = f[1:end]

    ##########################
    # Mathematical Functions #
    ##########################
    tensor{S,A,B}(a::FockBasis{S,A}, b::FockBasis{S,B}) = FockBasis{S,A+B}(tuple(a.ranges..., b.ranges...))

    ######################
    # Printing Functions #
    ######################
    repr(f::FockBasis) = "$(typeof(f))$(ranges(f))"
    show(io::IO, f::FockBasis) = print(io, repr(f))

export FockBasis,
    structure,
    tensor,
    nfactors,
    samelabels,
    labelvec