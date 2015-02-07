

nfactors{N}(::NTuple{N}) = N
nfactors{N}(::Type{NTuple{N}}) = N

#####################
# Joining Functions #
#####################
tupletensor(items) = apply(tuple, items...) 
tupletensor(items...) = apply(tuple, items...) 

tensor(t::Tuple...) = tupletensor(t)
separate(t::Tuple) = map(tuple, t)

######################
# Printing Functions #
######################
labelstr(label::Tuple) = strip(repr(label)[2:end-1], ',')

#################################
# Iterator/Array-like Functions #
#################################
permute(t::Tuple, p) = tuple(permute!(collect(gettuple(s)), p)...)

function switch(t::Tuple, i, j)
    v = collect(t) 
    tmp = v[i]
    v[i] = v[j]
    v[j] = tmp
    return tuple(v...)
end
