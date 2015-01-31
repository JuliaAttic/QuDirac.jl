import Base: getindex,
    length,
    convert,
    repr,
    show,
    start,
    done,
    next,
    endof,
    last,
    first,
    collect,
    map

##############
# StateLabel #
##############
    immutable StateLabel{N}
        label::NTuple{N}
    end

    StateLabel(s::StateLabel) = StateLabel(gettuple(s))
    StateLabel(s::AbstractState) = label(s)
    StateLabel(label...) = StateLabel(label)

    convert(::Type{StateLabel}, t::Tuple) = StateLabel(t)

    ###############################
    # Accessor/Property Functions #
    ###############################
    label(s::StateLabel) = s
    gettuple(s::StateLabel) = s.label
    getindex(s::StateLabel, i) = getindex(s.label, i)
    nfactors{N}(::StateLabel{N}) = N
    nfactors{N}(::Type{StateLabel{N}}) = N
    length{N}(::StateLabel{N}) = N

    #####################
    # Joining Functions #
    #####################
    tupletensor(a, b) = tuple(a..., b...)
    binarytensor(a::StateLabel, b::StateLabel) = StateLabel(tupletensor(gettuple(a), gettuple(b))) 
    tensor(s::(StateLabel...)) = reduce(binarytensor, s)
    tensor(s::StateLabel...) = tensor(s)

    separate(s::StateLabel) = map(StateLabel, gettuple(s))

    ######################
    # Printing Functions #
    ######################
    labelstr(label::Tuple) = strip(repr(label)[2:end-1], ',')
    labelstr(s::StateLabel) = "$(labelstr(s.label))"
    repr(s::StateLabel) = "StateLabel($(labelstr(s)))"
    show(io::IO, s::StateLabel) = print(io, repr(s))

    #################################
    # Iterator/Array-like Functions #
    #################################
    start(::StateLabel) = 1
    done(s::StateLabel, state) = length(s) == state-1
    next(s::StateLabel, state) = s[state], state+1
    endof(s::StateLabel) = length(s)
    last(s::StateLabel) = s[length(s)]
    first(s::StateLabel) = s[1]
    collect(s::StateLabel) = s[1:end]

    map(f::Union(Function,DataType), s::StateLabel) = StateLabel(map(f, gettuple(s)))
    map(f, s::StateLabel) = StateLabel(map(f, gettuple(s)))

    permute(s::StateLabel, p) = StateLabel(permute!(collect(gettuple(s)), p)...)

    function switch(s::StateLabel, i, j)
        v = collect(gettuple(s)) 
        tmp = v[i]
        v[i] = v[j]
        v[j] = tmp
        return StateLabel(v...)
    end

    labelvec(labels::Array) = map(StateLabel, labels)
    labelvec(labels::AbstractArray) = collect(map(StateLabel, labels))
    labelvec(labels...) = labelvec(collect(labels))
    labelvec{S<:StateLabel}(labels::AbstractArray{S}) = collect(labels)
    labelvec{S<:StateLabel}(labels::Set{S}) = collect(labels)
