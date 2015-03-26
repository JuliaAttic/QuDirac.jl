###########
# Ket/Bra #
###########
    typealias StateDict Dict{Vector{Any},Number}

    type Ket{P,N} <: AbstractState{P,N}
        dict::StateDict
        fact::Factors{N}
        Ket(dict,fact) = new(dict,fact)
        Ket(dict,::Factors{0}) = error("Cannot construct a 0-factor state; did you mean to construct a scalar?")
    end

    Ket{P,N}(::Type{P}, dict::StateDict, fact::Factors{N}) = Ket{P,N}(dict,fact)

    ket{P<:AbstractInner,N}(::Type{P}, label::NTuple{N}) = Ket(P,single_dict(StateDict(), collect(label), 1), Factors{N}())
    ket{P<:AbstractInner}(::Type{P}, items...) = ket(P,items)
    ket(items...) = ket(DEFAULT_INNER, items)

    type Bra{P,N} <: AbstractState{P,N}
        kt::Ket{P,N}
        Bra(kt::Ket{P,N}) = new(kt)
        Bra(items...) = new(Ket{P,N}(items...))
    end

    Bra{P,N}(kt::Ket{P,N}) = Bra{P,N}(kt)
    Bra(items...) = Bra(Ket(items...))
    bra(items...) = Bra(ket(items...))

    dict(k::Ket) = k.dict
    dict(b::Bra) = dict(b.kt)

    fact(k::Ket) = k.fact
    fact(b::Bra) = fact(b.kt)

################
# Constructors #
################
    Base.copy(s::AbstractState) = typeof(s)(copy(dict(s)), fact(s))
    Base.similar(s::AbstractState, d::StateDict=similar(dict(s))) = typeof(s)(d, fact(s))

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){P,N}(a::Ket{P,N}, b::Ket{P,N}) = dict(a) == dict(b)
    Base.(:(==)){P,N}(a::Bra{P,N}, b::Bra{P,N}) = a.kt == b.kt
    Base.hash(s::AbstractState) = hash(dict(s), hash(typeof(s)))

    Base.length(s::AbstractState) = length(dict(s))

    Base.getindex(k::Ket, label::Array) = dict(k)[label]
    Base.getindex(b::Bra, label::Array) = b.kt[label]'
    Base.getindex(s::AbstractState, i...) = s[collect(i)]

    _setindex!(k::Ket, c, label::Array) = setindex!(dict(k), c, label)
    _setindex!(b::Bra, c, label::Array) = setindex!(b.kt, c', label)

    function Base.setindex!{P,N}(s::AbstractState{P,N}, c, label::Array)
        if length(label) == N
            return _setindex!(s, c, label)
        else
            throw(BoundsError())
        end
    end

    Base.setindex!{P,N}(::AbstractState{P,N}, c, ::Tuple) =  throw(BoundsError())
    Base.setindex!{P,N}(s::AbstractState{P,N}, label::NTuple{N}) =  _setindex!(s, c, collect(label))
    Base.setindex!{P,N}(s::AbstractState{P,N}, c, i...) = setindex!(s,c,i)

    Base.haskey(s::AbstractState, label::Array) = haskey(dict(s), label)
    Base.get(k::Ket, label::Array, default) = get(dict(k), label, default)
    Base.get(b::Bra, label::Array, default) = haskey(b, label) ? b[label] : default

    Base.delete!(s::AbstractState, label::Array) = (delete!(dict(s), label); return s)

    labels(s::AbstractState) = keys(dict(s))
    coeffs(kt::Ket) = values(dict(kt))
    coeffs(br::Bra) = imap(ctranspose, coeffs(br.kt))

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!(f::Function, kt::Ket) = (filter!(f, dict(kt)); return kt)
    Base.filter!(f::Function, br::Bra) = (filter!((k,v)->f(k,v'), dict(br)); return br)

    Base.filter(f::Function, kt::Ket) = similar(kt, filter(f, dict(kt)))
    Base.filter(f::Function, br::Bra) = similar(br, filter((k,v)->f(k,v'), dict(br)))

    labelcheck(label::Array, N) = length(label) == N ? label : throw(BoundsError())
    ktprune(kv,N) = length(kv[1]) == N ? kv : throw(BoundsError())
    brprune(kv,N) = (labelcheck(kv[1], N), kv[2]')

    Base.map{P,N}(f::Function, kt::Ket{P,N}) = similar(kt, mapkv((k,v)->ktprune(f(k,v),N), dict(kt)))
    Base.map{P,N}(f::Function, br::Bra{P,N}) = similar(br, mapkv((k,v)->brprune(f(k,v'),N), dict(br)))

    mapcoeffs!(f::Function, k::Ket) = (mapvals!(f, dict(k)); return k)
    mapcoeffs!(f::Function, b::Bra) = (mapvals!(v->f(v')', dict(b)); return b)
    mapcoeffs(f::Function, kt::Ket) = similar(kt, mapvals(f, dict(kt)))
    mapcoeffs(f::Function, br::Bra) = similar(br, mapvals(v->f(v')', dict(br)))

    maplabels!{P,N}(f::Function, s::AbstractState{P,N}) = (mapkeys!(label->labelcheck(f(label), N), dict(s)); return s)
    maplabels{P,N}(f::Function, s::AbstractState{P,N}) = similar(s, mapkeys(label->labelcheck(f(label), N), dict(s)))

    function wavefunc(f::Function, kt::Ket)
        return (args...) -> sum(pair->pair[2]*f(pair[1])(args...), dict(kt))
    end

##########################
# Mathematical Functions #
##########################
    nfactors{P,N}(::AbstractState{P,N}) = N

    inner(br::Bra, kt::Ket) = error("inner(b::Bra,k::Ket) is only defined when nfactors(b) == nfactors(k)")
    inner(br::Bra, kt::Ket, i) = error("inner(b::Bra,k::Ket,i) is only defined when nfactors(b) == 1")

    function inner{A,B,N}(br::Bra{A,N}, kt::Ket{B,N})
        result = 0
        for (b,c) in dict(br), (k,v) in dict(kt)
            result += c'*v*inner_eval(A,B,b,k)
        end
        return result  
    end

    function inner{A,B}(br::Bra{A,1}, kt::Ket{B,1}, i)
        if i==1
            return inner(br, kt)
        else
            throw(BoundsError())
        end
    end

    function inner{A,B,N}(br::Bra{A,1}, kt::Ket{B,N}, i)
        result = StateDict()
        P = typejoin(A,B)
        for (b,c) in dict(br), (k,v) in dict(kt)
            add_to_dict!(result, 
                         except(k,i),
                         c'*v*inner_rule(P,b[1],k[i]))
        end
        return Ket(P, result, Factors{N-1}())
    end 

    function ortho_inner{A<:Orthonormal,B<:Orthonormal}(a::AbstractState{A}, b::AbstractState{B})
        result = 0
        for label in keys(dict(b))
            if haskey(a, label)
                result += a[label]*b[label]*inner_eval(A,B,label,label)
            end
        end
        return result
    end

    function inner{A<:Orthonormal,B<:Orthonormal,N}(br::Bra{A,N}, kt::Ket{B,N})
        if length(br) < length(kt)
            return ortho_inner(kt, br)
        else
            return ortho_inner(br, kt)
        end
    end

    Base.scale!(k::Ket, c::Number) = (dscale!(dict(k), c); return k)
    Base.scale!(c::Number, k::Ket) = scale!(k,c)
    Base.scale!(b::Bra, c::Number) = Bra(scale!(b.kt, c'))
    Base.scale!(c::Number, b::Bra) = scale!(b,c)

    Base.scale(k::Ket, c::Number) = similar(k, dscale(dict(k), c))
    Base.scale(c::Number, k::Ket) = scale(k,c)
    Base.scale(b::Bra, c::Number) = Bra(scale(b.kt, c'))
    Base.scale(c::Number, b::Bra) = scale(b,c)
    
    Base.(:+){P,N}(a::Ket{P,N}, b::Ket{P,N}) = similar(b, addmerge(dict(a), dict(b)))
    Base.(:-){P,N}(a::Ket{P,N}, b::Ket{P,N}) = a + (-b)
    Base.(:-){P,N}(kt::Ket{P,N}) = scale(-1,kt)

    Base.(:+)(a::Bra, b::Bra) = Bra(a.kt+b.kt)
    Base.(:-)(a::Bra, b::Bra) = Bra(a.kt-b.kt)
    Base.(:-)(br::Bra) = Bra(-br.kt)

    Base.(:*)(br::Bra, kt::Ket) = inner(br,kt)
    Base.(:*)(a::Ket, b::Ket) = tensor(a,b)
    Base.(:*)(a::Bra, b::Bra) = tensor(a,b)

    Base.(:*)(c::Number, s::AbstractState) = scale(c, s)
    Base.(:*)(s::AbstractState, c::Number) = scale(s, c)
    Base.(:/)(s::AbstractState, c::Number) = scale(s, 1/c)

    Base.ctranspose(k::Ket) = Bra(k)
    Base.ctranspose(b::Bra) = b.kt
    Base.norm(s::AbstractState) = sqrt(sum(sqr, values(dict(s))))

    QuBase.tensor{P}(a::Ket{P}, b::Ket{P}) = Ket(P,tensorstate!(StateDict(), dict(a), dict(b)), fact(a)+fact(b))
    QuBase.tensor(a::Bra, b::Bra) = tensor(a.kt, b.kt)'

    QuBase.normalize(s::AbstractState) = (1/norm(s))*s
    QuBase.normalize!(s::AbstractState) = scale!(1/norm(s), s)

    xsubspace(s::AbstractState, x) = similar(s, filter((k,v)->isx(k,x), dict(s)))
    switch(s::AbstractState, i, j) = similar(s, mapkeys(label->switch(label,i,j), dict(s)))
    permute(s::AbstractState, perm) = similar(s, mapkeys(label->permute(label,perm), dict(s)))
    switch!(s::AbstractState, i, j) = (mapkeys!(label->switch(label,i,j), dict(s)); return s)
    Base.permute!(s::AbstractState, perm::AbstractVector) = (mapkeys!(label->permute(label,perm), dict(s)); return s)

    filternz!(s::AbstractState) = (filter!(nzcoeff, dict(s)); return s)
    filternz(s::AbstractState) = similar(s, filter(nzcoeff, dict(s)))

    # should always be pure, of course,
    # but makes a good sanity check function
    purity(kt::Ket) = purity(kt*kt')
    purity(br::Bra) = purity(br.kt)

    queval(f, s::AbstractState) = mapcoeffs(x->queval(f,x),s)

######################
# Printing Functions #
######################
    labelstr(label) = join(map(repr, label), ',')
    ktstr(label) = "| $(labelstr(label)) $rang"
    brstr(label) = "$lang $(labelstr(label)) |"

    labelrepr(kt::Ket, label, pad) = "$pad$(kt[label]) $(ktstr(label))"
    labelrepr(br::Bra, label, pad) = "$pad$(br[label]) $(brstr(label))"

    Base.summary(s::AbstractState) = "$(typeof(s)) with $(length(s)) state(s)"
    Base.show(io::IO, s::AbstractState) = dirac_show(io, s)
    Base.showcompact(io::IO, s::AbstractState) = dirac_showcompact(io, s)
    Base.repr(s::AbstractState) = dirac_repr(s)

export ket,
    bra,
    nfactors,
    maplabels!,
    mapcoeffs!,
    maplabels,
    mapcoeffs,
    xsubspace,
    switch,
    permute,
    switch!,
    permute!,
    filternz!,
    filternz,
    purity,
    wavefunc,
    labels,
    coeffs,
    inner