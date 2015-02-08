nfactors{N}(::NTuple{N}) = N
nfactors{N}(::Type{NTuple{N}}) = N

second(t::Tuple) = t[2]

#####################
# Joining Functions #
#####################
tensor_tup(items) = apply(tuple, items...)
separate(t::Tuple) = map(tuple, t)

######################
# Printing Functions #
######################
labelstr(label::Tuple) = strip(repr(label)[2:end-1], ',')

############################
# Combinatorical Functions #
############################
permute(t::Tuple, p) = tuple(permute!(collect(gettuple(s)), p)...)

function switch(t::Tuple, i, j)
    v = collect(t) 
    tmp = v[i]
    v[i] = v[j]
    v[j] = tmp
    return tuple(v...)
end
