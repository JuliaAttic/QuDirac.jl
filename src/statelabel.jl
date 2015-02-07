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
    type StateLabel{T<:Tuple}
        label::T
    end

    typealias NLabel{N,T} StateLabel{NTuple{N,T}}

    StateLabel(s::StateLabel) = StateLabel(gettuple(s))
    StateLabel(s::AbstractState) = label(s)
    StateLabel(label...) = StateLabel(label)

    convert(::Type{StateLabel}, t::Tuple) = StateLabel(t)

    ###############################
    # Accessor/Property Functions #
    ###############################
    label(s::StateLabel) = s
    gettuple(s::StateLabel) = s.label
    nfactors{N}(::NLabel{N}) = N

    #####################
    # Joining Functions #
    #####################
    tensor(s::StateLabel...) = apply(StateLabel, s...)
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
    length(s::StateLabel) = length(s.label)
    getindex(s::StateLabel, i) = getindex(s.label, i)

    start(s::StateLabel) = start(s.label)
    done(s::StateLabel, state) = done(s.label, state)
    next(s::StateLabel, state) = next(s.label, state)
    endof(s::StateLabel) = endof(s.label)
    last(s::StateLabel) = last(s.label)
    first(s::StateLabel) = first(s.label)
    collect(s::StateLabel) = collect(s.label)

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
