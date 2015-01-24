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
    hcat

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
    #   KetVector in LabelBasis{AbstractStructure,2} with 2 Complex{Float64} entries:
    #    0.7071067811865475 + 0.0im | 'a','b' ⟩
    #    0.7071067811865475 + 0.0im | 1,2 ⟩
    # 
    # Obviously, this is a silly state - it's merely a good 
    # demonstration of the generality enabled by the `LabelBasis` 
    # type.

    immutable LabelBasis{S<:AbstractStructure,N} <: AbstractLabelBasis{S,N}
        labels::Vector{StateLabel{N}}      # stores StateLabels
        labelmap::Dict{StateLabel{N}, Int} # used for getpos(basis,label) -> index functionality 
        labels_hash::Uint64                # used for constant-time samelabels operation
        
        function LabelBasis(labels::Vector{StateLabel{N}}, 
                            labelmap::Dict{StateLabel{N}, Int}, 
                            labels_hash::Uint64)
            return new(labels, labelmap, labels_hash)
        end 

        function LabelBasis(labels::Vector{StateLabel{N}})
            labelmap = Dict{StateLabel{N},Int}()
            sizehint(labelmap, length(labels))
            
            labelmap[labels[1]] = 1
            labels_hash = hash(labels[1])
            for i=2:length(labels)
                labels_hash = hash(labels_hash, hash(labels[i]))
                labelmap[labels[i]] = i
            end
            
            if length(labels)==length(labelmap) # ensures uniqueness and correct lengths
                return LabelBasis{S,N}(labels, labelmap, labels_hash)
            else
                error("Basis states are malformed w.r.t. their index mapping; perhaps the labels are not unique?")
            end
        end
    end
    
    LabelBasis{S<:AbstractStructure, N}(::Type{S}, labels::Vector{StateLabel{N}}) = LabelBasis{S,N}(labels)
    LabelBasis{S<:AbstractStructure}(::Type{S}, labels) = LabelBasis(S, labelvec(labels))
    LabelBasis{S<:AbstractStructure}(::Type{S}, arrs::AbstractArray...) = LabelBasis(S, cart_prod(map(labelvec, arrs)))
    LabelBasis{S<:AbstractStructure}(::Type{S}, labels...) = LabelBasis(S, labelvec(labels...))
    LabelBasis(args...) = LabelBasis(AbstractStructure, args...)

    convert{A,B,N}(::Type{LabelBasis{A,N}}, b::LabelBasis{B,N}) = LabelBasis{A,N}(b.labels, b.labelmap, b.labels_hash)
    convert{A,B,N}(::Type{LabelBasis{A,N}}, b::AbstractLabelBasis{B,N}) = LabelBasis(A, labelvec(b))

    copy{S,N}(b::LabelBasis{S,N}) = LabelBasis{S,N}(copy(b.labels), copy(b.labelmap), copy(b.labels_hash))

    ######################
    # Property Functions #
    ######################
    structure{S}(::Type{LabelBasis{S}}) = S
    structure{S,N}(::Type{LabelBasis{S,N}}) = S

    size(basis::LabelBasis) = size(basis.labels)
    length(basis::LabelBasis) = length(basis.labels)
    nfactors{S,N}(::LabelBasis{S,N}) = N

    labelvec(basis::LabelBasis) = basis.labels

    samelabels(a::LabelBasis, b::LabelBasis) = a.labels_hash == b.labels_hash
    =={A,B,C,D}(::LabelBasis{A,B}, ::LabelBasis{C,D}) = false
    =={S,N}(a::LabelBasis{S,N}, b::LabelBasis{S,N}) = samelabels(a,b)

    hash{S}(basis::LabelBasis{S}) = hash(basis.labels_hash, hash(S))
    checkcoeffs(coeffs, dim, basis::LabelBasis) = size(coeffs, dim) == length(basis) 

    ########################
    # Array-like Functions #
    ########################
    getpos(basis::LabelBasis, label::StateLabel) = basis.labelmap[label]
    getpos(basis::LabelBasis, label::Tuple) = getpos(basis, StateLabel(label))
    
    in(label::StateLabel, basis::LabelBasis) = haskey(basis.labelmap, label)

    getindex(basis::LabelBasis, i) = basis.labels[i]

    #####################
    # Joining Functions #
    #####################
    append{S,N}(basis::LabelBasis{S,N}, label::StateLabel{N}) = label in basis ? basis : append_label(basis, label)
    append{S,N}(a::LabelBasis{S,N}, b::AbstractLabelBasis{S,N}) = append_labelvec(a, labelvec(b))
    append{S,N}(a::AbstractLabelBasis{S,N}, b) = append(convert(LabelBasis{S,N}, a), b)

    function setdiff{S,N}(a::LabelBasis{S,N}, b::LabelBasis{S,N})
        return LabelBasis{S}(setdiff(labelvec(a), labelvec(b)))
    end

    function hcat{S}(bases::AbstractLabelBasis{S}...)
        return LabelBasis(S, hcat_labels(map(labelvec, bases)))
    end

    function factorize{S,N}(basis::LabelBasis{S,N})
        n = nfactors(basis)
        sets = ntuple(n, x->OrderedSet{StateLabel,N}())
        for i=1:length(basis)
            for j=1:n
                push!(sets[j], StateLabel{1}(basis[i][j]))
            end
        end
        return map(labels->LabelBasis{S,1}(labelvec(labels), BypassFlag), sets)
    end

    # Tensor product of labels
    function tensor{S}(bases::AbstractLabelBasis{S}...)
        return LabelBasis(S, cart_prod(map(labelvec, bases))) 
    end
    function tensor{S}(basis::AbstractLabelBasis{S}, label::StateLabel)
        return LabelBasis(S, map(s->combine(s, label), labelvec(basis))) 
    end
    function tensor{S}(label::StateLabel, basis::AbstractLabelBasis{S})
        return LabelBasis(S, map(s->combine(label, s), labelvec(basis))) 
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

    ####################
    # Helper Functions #
    ####################
    function append_label!(map::Dict, label::StateLabel)
        map[label] = length(map)+1
        return map
    end
    
    function append_label{S,N}(b::LabelBasis{S,N}, label::StateLabel{N})
        return LabelBasis{S,N}(vcat(b.labels, label), 
                    append_label!(copy(b.labelmap), label), 
                    hash(b.labels_hash, hash(label)))
    end

    function append_labelvec{S,N}(a::LabelBasis{S,N}, b::Vector{StateLabel{N}})
        labelmap = copy(a.labelmap)
        labels_hash = a.labels_hash
        for label in b
            if ! in(label,a)
                append_label!(labelmap, label)
                labels_hash = hash(labels_hash, hash(label))
            end
        end
        return LabelBasis{S}(unique(vcat(a.labels, b)), labelmap, labels_hash)
    end

    function hcat_labels{A,B}(a::Vector{StateLabel{A}}, b::Vector{StateLabel{B}})
        if length(a)==length(b)
            return StateLabel{A+B}[combine(a[i], b[i]) for i=1:length(a)]
        else
            error("Could not take direct product of bases of differing length")
        end
    end

    hcat_labels(labels::(Vector...,)) = reduce(hcat_labels, labels)

    function cart_prod(labels::(Vector...,))
        lens = map(length, labels)
        N = sum(map(nfactors, map(eltype, labels)))
        arr = Array(StateLabel{N}, prod(lens))
        index = 1

        function set_ind!(inds...)
            arr[index] = combine(map(getindex, labels, inds))
            index += 1
        end
        
        cartesianmap(set_ind!, lens)
        return arr
    end

export LabelBasis,
    structure,
    nfactors,
    factorize,
    labelvec,
    getpos,
    samelabels,
    append,
    tensor