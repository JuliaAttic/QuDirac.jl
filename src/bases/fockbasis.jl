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
    in

#############
# FockBasis #
#############
    # A FockBasis uses precomputed values to efficiently 
    # generate StateLabels for given indices in the basis, or
    # vice versa (an index for a given state in the basis).
    # For example:
    #
    #   julia> f=FockBasis(2,2,2)
    #    FockBasis{AbstractStructure}(2,2,2)    
    #   
    #   julia> labelvec(f)
    #    8-element Array{StateLabel,1}:
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
    #   julia> f[StateLabel(1,0,1)]
    #    6
    # 
    # Because the labels are generated rather than actually 
    # stored, one can represent very large bases without 
    # any storage overhead:
    #
    #   julia> f = FockBasis(221,135,31,42,321,3)
    #    FockBasis{AbstractStructure}(221,135,31,42,321,3)
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
    #   julia> f[StateLabel(128,60,0,37,0,0)]
    #    34234134
    # 
    # Arbitrary numeric ranges are supported for labels:
    #   
    #   julia> f=FockBasis(0.0:0.1:0.2, 4:7)
    #   FockBasis{AbstractStructure}(0.0:0.1:0.2,4:7)     
    #       
    #   julia> collect(f)
    #   12-element Array{StateLabel,1}:
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
    #   julia> f[f[(0.0, 7)]] == StateLabel(0.0, 7)
    #   true
    
    immutable FockBasis{S<:AbstractStructure} <: AbstractLabelBasis{S}
        ranges::(Range...)
        denoms::(Float64...)
        FockBasis(ranges, denoms, ::Type{BypassFlag}) = new(ranges, denoms)

        FockBasis(::()) = error("")
        FockBasis() = error("")

        function FockBasis(ranges::(Range...))
            # reverse is done to match cartesianmap order
            return FockBasis{S}(ranges, precompute_denoms(reverse(map(length,ranges))), BypassFlag) 
        end

        FockBasis(lens::(Int...)) = FockBasis{S}(map(x->0:x, lens))
        FockBasis(lens::Int...) = FockBasis{S}(lens)
        FockBasis(lens::Range...) = FockBasis{S}(lens)
    end

    FockBasis(lens...) = FockBasis{AbstractStructure}(lens)

    convert{S}(::Type{FockBasis{S}}, f::FockBasis) = FockBasis{S}(f.ranges, f.denoms, BypassFlag)
    convert{S}(::Type{FiniteBasis{S}}, f::FockBasis) = FiniteBasis{S}(size(f))

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
        total_divisor = prod(map(float,lens)) 
        function get_denom(i)
            total_divisor = div(total_divisor, i)
            return max(1.0, total_divisor)
        end
        # reverse is done to match cartesianmap order
        return reverse(map(get_denom, lens))
    end

    ######################
    # Property Functions #
    ######################
    @defstructure FockBasis

    labelvec(f::FockBasis) = collect(f)
    ranges(f::FockBasis) = f.ranges
    ranges(f::FockBasis, i) = f.ranges[i]

    size(f::FockBasis) = map(length, ranges(f))
    size(f::FockBasis, i) = length(ranges(f, i))
    length(f::FockBasis) = prod(size(f))
    ndims(f::FockBasis) = length(ranges(f))
    nfactors(f::FockBasis) = ndims(f)

    samelabels(a::FockBasis, b::FockBasis) = ranges(a) == ranges(b)

    checkcoeffs(coeffs::AbstractArray, dim::Int, f::FockBasis) = size(coeffs, dim) = length(f)

    ######################
    # Accessor Functions #
    ######################
    ind_value(n, range, denom, modulus) = range[(div(n, denom) % modulus)+1]
    tuple_at_ind(f::FockBasis, i) = ntuple(ndims(f), x->ind_value(i-1, ranges(f,x), f.denoms[x], size(f,x)))
    pos_in_range(r::Range, i) = i in r ? (i-first(r))/step(r) : throw(BoundsError())

    in(label, f::FockBasis) = reduce(&, map(in, label, ranges(f)))
    getpos(f::FockBasis, s::AbstractState) = getpos(f, label(s))
    getpos(f::FockBasis, label) = int(sum(map(*, map(pos_in_range, ranges(f), label), f.denoms)))+1

    getindex(f::FockBasis, i) = StateLabel(tuple_at_ind(f, i))
  
    getindex(f::FockBasis, s::AbstractState) = getpos(f, s) 
    getindex(f::FockBasis, label::StateLabel) = getpos(f, label) 
    getindex(f::FockBasis, label::Tuple) = getpos(f, label)

    getindex(f::FockBasis, arr::AbstractArray) = [f[i] for i in arr]

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
    tensor(a::FockBasis, b::FockBasis) = FockBasis(tuple(a.ranges..., b.ranges...))

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