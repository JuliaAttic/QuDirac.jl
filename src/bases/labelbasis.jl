import Base:
    convert,
    copy,
    size,
    length,
    in,
    getindex,
    filter,
    map,
    ctranspose,
    vcat,
    setdiff,
    summary,
    show,
    start,
    done,
    endof,
    next,
    last,
    first,
    ==,
    hash,
    hcat,
    vcat

####################
# Helper Functions #
####################

    # XORing hashes is commutative; using this instead
    # of the normal hash function to mix hashes would allow 
    # for constant time checking of whether two LabelBasis objects have 
    # the same labels regardless of label order. 
    # This is not currently in use because, in most cases, I
    # think the user will only care if the labels are both
    # the same and in order...
    # 
    # combine_hashes(hash1, hash2) = hash1 $ hash2 

    nfactors(labels) = length(first(labels))
    function samelength(labels)
        n = nfactors(labels)
        return all(map(x->length(x)==n, labels))
    end

##############
# LabelBasis #
##############
    # A `LabelBasis` is a basis type for arbitrary storage of 
    # `StateLabel`s in a manner that supports the behavior 
    # an AbstractLabelBasis (see list of expected methods in
    # bases.jl). It also supports tensor products.
    # 
    # The point of this type is that it allows the representation
    # of bases that do NOT assume a fock-like structure.
    # It also supports arbitrary labeling, which is nice.
    # 
    # This type will be the backbone for implementing
    # bra-ket additive arthimetic involving DiracArrays -
    # it allows additive expressions like this:
    #   
    #   julia> 1/sqrt(2) * (ket('a','b') + ket(1,2))
    #   2-element DiracVector{Ket, AbstractStructure, Float64, LabelBasis{AbstractBasis}}:
    #    0.7071067811865475 | 'a','b' ⟩
    #    0.7071067811865475 | 1,2 ⟩
    # 
    # Obviously, this is a silly state - it's merely a good 
    # demonstration of the generality enabled by the `LabelBasis` 
    # type.
    #
    # Note that DiracArrays are not actually implemented yet.
    # The above example is from QuDirac. This functionality 
    # will be present here soon.

    immutable LabelBasis{S<:AbstractStructure} <: AbstractLabelBasis{S}
        labels::Vector{StateLabel}      # stores StateLabels
        labelmap::Dict{StateLabel, Int} # used for basis[state] -> index functionality 
        labels_hash::Uint64             # used for constant-time samelabels operation
        
        function LabelBasis(labels::Vector{StateLabel}, 
                            labelmap::Dict{StateLabel, Int}, 
                            labels_hash::Uint64)
            return new(labels, labelmap, labels_hash)
        end 

        function LabelBasis(labels::Vector{StateLabel}, ::Type{BypassFlag})
            labelmap = Dict{StateLabel,Int}()
            sizehint(labelmap, length(labels))
            
            labelmap[labels[1]] = 1
            labels_hash = hash(labels[1])
            for i=2:length(labels)
                labels_hash = hash(labels_hash, hash(labels[i]))
                labelmap[labels[i]] = i
            end
            
            if length(labels)==length(labelmap) # ensures uniqueness and correct lengths
                return LabelBasis{S}(labels, labelmap, labels_hash)
            else
                error("Basis states are malformed w.r.t. their index mapping; perhaps the labels are not unique?")
            end
        end

        function LabelBasis(labels::Vector{StateLabel})
            if samelength(labels)
                return LabelBasis{S}(labels, BypassFlag)
            else
                error("input labels must be of uniform length")
            end
        end

        LabelBasis(labels) = LabelBasis{S}(labelvec(labels))
        LabelBasis(arrs::AbstractArray...) = LabelBasis{S}(cart_prod(map(labelvec, arrs)))
        LabelBasis(labels...) = LabelBasis{S}(labelvec(labels...))
    end
    
    LabelBasis(args...) = LabelBasis{AbstractStructure}(args...)

    convert{S}(::Type{LabelBasis{S}}, b::LabelBasis) = LabelBasis{S}(b.labels, b.labelmap, b.labels_hash)
    copy{S}(b::LabelBasis{S}) = LabelBasis{S}(copy(b.labels), copy(b.labelmap), copy(b.labels_hash))

    ######################
    # Property Functions #
    ######################
    @defstructure LabelBasis

    size(basis::LabelBasis) = size(basis.labels)
    length(basis::LabelBasis) = length(basis.labels)

    labelvec(basis::LabelBasis) = basis.labels

    samelabels(a::LabelBasis, b::LabelBasis) = a.labels_hash == b.labels_hash
    =={S}(a::LabelBasis{S}, b::LabelBasis{S}) = samelabels(a,b)
    =={A,B}(a::LabelBasis{A}, b::LabelBasis{B}) = false
    hash{S}(basis::LabelBasis{S}) = hash(basis.labels_hash, hash(S))
    checkcoeffs(coeffs::AbstractArray, dim::Int, basis::LabelBasis) = size(coeffs, dim) == length(basis) 

    ########################
    # Array-like Functions #
    ########################
    getpos(basis::LabelBasis, label::StateLabel) = basis.labelmap[label]
    getpos(basis::LabelBasis, s::AbstractState) = getpos(basis, label(s)) 
    getpos(basis::LabelBasis, label::Tuple) = getpos(basis, StateLabel(label))
    
    in(label::StateLabel, basis::LabelBasis) = haskey(basis.labelmap, label)

    getindex(basis::LabelBasis, s::AbstractState) = getpos(basis, s) 
    getindex(basis::LabelBasis, label::StateLabel) = getpos(basis, label) 
    getindex(basis::LabelBasis, label::Tuple) = getpos(basis, label)
    
    getindex(basis::LabelBasis, i) = basis.labels[i]

    #####################
    # Joining Functions #
    #####################
    function append_label!(map::Dict, label::StateLabel)
        map[label] = length(map)+1
        return map
    end

    function append{S}(b::LabelBasis{S}, label::StateLabel, ::Type{BypassFlag})
        if nfactors(b) == length(label)
            return LabelBasis{S}(vcat(b.labels, label), 
                        append_label!(copy(b.labelmap), label), 
                        hash(b.labels_hash, hash(label)))
        else
            error("input labels not of uniform length")
        end
    end

    append(basis::LabelBasis, label::StateLabel) = label in basis ? basis : append(basis, label, BypassFlag)

    function append{S}(a::LabelBasis{S}, b::LabelBasis{S}) 
        if nfactors(a) == nfactors(b)
            labelmap = copy(a.labelmap)
            labels_hash = a.labels_hash
            for label in b.labels
                if ! in(label,a)
                    append_label!(labelmap, label)
                    labels_hash = hash(labels_hash, hash(label))
                end
            end
            return LabelBasis{S}(unique(vcat(a.labels, b.labels)), labelmap, labels_hash)
        else
            error("input labels not of uniform length")
        end
    end

    function setdiff{S}(a::AbstractLabelBasis{S}, b::AbstractLabelBasis{S})
        return LabelBasis{S}(setdiff(labelvec(a), labelvec(b)))
    end

    function hcat_method(a::Vector{StateLabel}, b::Vector{StateLabel})
        if length(a)==length(b)
            return StateLabel[combine(a[i], b[i]) for i=1:length(a)]
        else
            error("Could not take direct product of bases of differing length")
        end
    end

    hcat_method(labels::(Vector{StateLabel}...,)) = reduce(dir_prod, labels)

    function hcat{S}(bases::AbstractLabelBasis{S}...)
        return LabelBasis{S}(hcat_method(map(labelvec, bases)), BypassFlag)
    end

    function cart_prod(labels::(Vector{StateLabel}...,))
        lens = map(length, labels)
        arr = Array(StateLabel, prod(lens))
        index = 1

        function set_ind!(inds...)
            arr[index] = combine(map(getindex, labels, inds))
            index += 1
        end
        
        cartesianmap(set_ind!, lens)
        return arr
    end


    function tensor(bases::AbstractLabelBasis...)
        return LabelBasis{S}(cart_prod(map(labelvec, bases)), BypassFlag) 
    end

    function factorize{S}(basis::LabelBasis{S})
        n = nfactors(basis)
        sets = ntuple(n, x->OrderedSet{StateLabel}()) # use OrderedSet here
        for i=1:length(basis)
            for j=1:n
                push!(sets[j], StateLabel(basis[i][j]))
            end
        end
        return map(labels->LabelBasis{S}(labelvec(labels), BypassFlag), sets)
    end

    ######################
    # Iterator Functions #
    ######################
    start(::LabelBasis) = 1
    done(b::LabelBasis, state) = length(b) == state-1
    next(b::LabelBasis, state) = b[state], state+1
    endof(b::LabelBasis) = length(b)
    last(b::LabelBasis) = b[length(b)]
    first(b::LabelBasis) = b[1]
    collect(b::LabelBasis) = b[1:end]

    ######################
    # Printing Functions #
    ######################
    summary(basis::LabelBasis) = "$(length(basis))-element $(typeof(basis))"
    function show(io::IO, basis::LabelBasis)
        print(io, "$(summary(basis)):")
        
        if length(basis) > 25
            
            println(io)
            for i=1:10
                println(io, " ($(labelstr(basis[i])))")
            end

            print(io, " $vdots") # vdots is repo-wide const
            
            start_ind = length(basis)-10
        else
            start_ind = 1
        end

        for i=start_ind:length(basis)
            println(io)
            print(io, " ($(labelstr(basis[i])))")
        end
    end

export LabelBasis,
    structure,
    nfactors,
    factorize,
    labelvec,
    getpos,
    samelabels,
    append