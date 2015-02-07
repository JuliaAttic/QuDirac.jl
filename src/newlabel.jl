import Base: getindex,
    length,
    start,
    done,
    next,
    endof,
    last,
    first



type StateLabel
    label::Tuple
    StateLabel(t::Tuple) = new(t)
    StateLabel(i...) = new(i)    
end

# type StateLabel
#     label::Tuple
#     StateLabel(v::Vector{T}) = new(v)
#     StateLabel(i::T...) = StateLabel{T}(vcat(i...))    
# end

getindex(s::StateLabel, i) = getindex(s.label, i)
length(s::StateLabel) = length(s.label)

start(s::StateLabel) = start(s.label)
done(s::StateLabel, state) = done(s.label, state)
next(s::StateLabel, state) = next(s.label, state)
endof(s::StateLabel) = endof(s.label)
last(s::StateLabel) = last(s.label)
first(s::StateLabel) = first(s.label)
collect(s::StateLabel) = collect(s.label)

# in order from best to worst
tuptensor(t::Tuple...) = apply(tuple, t...)
tensor(s::StateLabel...) = apply(StateLabel, s...)
stensor(s::StateLabel...) = apply(tuple, s...)
