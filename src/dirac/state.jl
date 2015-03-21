###########
# Ket/Bra #
###########
    DEFAULT_INNER() = Orthonormal

    typealias StateDict Dict{Vector{Any},Number}

    type Ket{P<:AbstractInner} <: AbstractState{P}
        dict::StateDict
        Ket() = new(StateDict())
        Ket(dict::Dict) = new(dict)
        Ket(label::Vector) = new(single_dict(StateDict(), label, 1))
        Ket(items...) = Ket{P}(collect(items))
    end

    Ket(items...) = Ket{DEFAULT_INNER()}(items...)

    type Bra{P} <: AbstractState{P}
        ket::Ket{P}
        Bra(items...) = new(Ket{P}(items...))
        Bra(ket::Ket{P}) = new(ket)
    end

    Bra{P}(ket::Ket{P}) = Bra{P}(ket)
    Bra(items...) = Bra(Ket(items...))

    dict(k::Ket) = k.dict
    dict(b::Bra) = dict(b.ket)

################
# Constructors #
################
    Base.copy(s::AbstractState) = typeof(s)(copy(dict(s)))
    Base.similar(s::AbstractState) = typeof(s)(similar(dict(s)))

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){P}(a::Ket{P}, b::Ket{P}) = dict(a) == dict(b)
    Base.(:(==)){P}(a::Bra{P}, b::Bra{P}) = a.ket == b.ket
    Base.hash(s::AbstractState) = hash(dict(s), hash(typeof(s)))

    Base.length(s::AbstractState) = length(dict(s))

    Base.getindex(k::Ket, label::Array) = dict(k)[label]
    Base.getindex(b::Bra, label::Array) = b.ket[label]'
    Base.getindex(s::AbstractState, i...) = s[collect(i)]

    Base.setindex!(k::Ket, c, label::Array) = setindex!(dict(k), c, label)
    Base.setindex!(b::Bra, c, label::Array) = setindex!(b.ket, c', label)
    Base.setindex!(s::AbstractState, c, i...) = setindex!(s, c, collect(i))

    Base.haskey(s::AbstractState, label::Array) = haskey(dict(s), label)
    Base.get(k::Ket, label::Array, default) = get(dict(k), label, default)
    Base.get(b::Bra, label::Array, default) = haskey(b, label) ? b[label] : default

    Base.delete!(s::AbstractState, label::Array) = (delete!(dict(s), label); return s)

    labels(s::AbstractState) = keys(dict(s))
    coeffs(ket::Ket) = values(dict(ket))
    coeffs(bra::Bra) = imap(ctranspose, coeffs(bra.ket))

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!(f::Function, ket::Ket) = (filter!(f, dict(ket)); return ket)
    Base.filter!(f::Function, bra::Bra) = (filter!((k,v)->f(k,v'), bra.ket); return bra)

    Base.filter{P}(f::Function, ket::Ket{P}) = Ket{P}(filter(f, dict(ket)))
    Base.filter(f::Function, bra::Bra) = Bra(filter((k,v)->f(k,v'), bra.ket))

    Base.map{P}(f::Function, ket::Ket{P}) = Ket{P}(mapkv(f, dict(ket)))

    # By mutating an existing Bra instance, coefficients are
    # properly conjugated when they're both accessed *and* set
    Base.map(f::Function, bra::Bra) = mapkv!((k,v)->f(k,v'), similar(bra), bra.ket)

    mapcoeffs!(f::Function, k::Ket) = (mapvals!(f, dict(k)); return k)
    mapcoeffs!(f::Function, b::Bra) = (mapvals!(f, b, dict(b)); return b)
    mapcoeffs{P}(f::Function, ket::Ket{P}) = Ket{P}(mapvals(f, dict(ket)))
    mapcoeffs(f::Function, bra::Bra) = mapvals!(v->f(v'), similar(bra), bra.ket)

    maplabels!(f::Function, s::AbstractState) = (mapkeys!(f, dict(s)); return s)
    maplabels(f::Function, s::AbstractState) = typeof(s)(mapkeys(f, dict(s)))

    function wavefunc(f::Function, ket::Ket)
        return (args...) -> sum(pair->pair[2]*f(pair[1])(args...), dict(ket))
    end

##########################
# Mathematical Functions #
##########################
    nfactor_guess(s::AbstractState) = length(first(labels(s)))
    is_single(s::AbstractState) = nfactor_guess(s) == 1

    function inner{A,B}(bra::Bra{A}, ket::Ket{B}, i...)
        result = 0
        for (b,c) in dict(bra), (k,v) in dict(ket)
            result += c'*v*inner_eval(A,B,b,k,i...)
        end
        return result  
    end

    function inner{A,B}(bra::Bra{A}, ket::Ket{B}, i)
        if is_single(ket)
            return inner(bra, ket)
        else
            result = StateCoeffs()
            for (b,c) in dict(bra), (k,v) in dict(ket)
                add_to_dict!(result, 
                             except(k,i), 
                             c'*v*inner_eval(A,B,b,k,i,i))
            end
            return Ket{B}(result)
        end
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

    function inner{A<:Orthonormal,B<:Orthonormal}(bra::Bra{A}, ket::Ket{B})
        if length(bra) < length(ket)
            return ortho_inner(ket, bra)
        else
            return ortho_inner(bra, ket)
        end
    end

    Base.scale!(c::Number, k::Ket) = (castvals!(*, c, dict(k)); return k)
    Base.scale!(k::Ket, c::Number) = (castvals!(*, dict(k), c); return k)
    Base.scale!(c::Number, b::Bra) = Bra(scale!(c', b.ket))
    Base.scale!(b::Bra, c::Number) = Bra(scale!(b.ket, c'))

    Base.scale(c::Number, k::Ket) = typeof(k)(castvals(*, c, dict(k)))
    Base.scale(k::Ket, c::Number) = typeof(k)(castvals(*, dict(k), c))
    Base.scale(c::Number, b::Bra) = Bra(scale(c', b.ket))
    Base.scale(b::Bra, c::Number) = Bra(scale(b.ket, c'))

    Base.(:+){P}(a::Ket{P}, b::Ket{P}) = Ket{P}(mergef(+, dict(a), dict(b)))
    Base.(:-){P}(a::Ket{P}, b::Ket{P}) = a + (-b)
    Base.(:-){P}(ket::Ket{P}) = mapcoeffs(-, ket)

    Base.(:+)(a::Bra, b::Bra) = Bra(a.ket+b.ket)
    Base.(:-)(a::Bra, b::Bra) = Bra(a.ket-b.ket)
    Base.(:-)(bra::Bra) = Bra(-bra.ket)

    Base.(:*)(bra::Bra, ket::Ket) = inner(bra,ket)
    Base.(:*)(a::Ket, b::Ket) = tensor(a,b)
    Base.(:*)(a::Bra, b::Bra) = tensor(a,b)

    Base.(:*)(c::Number, s::AbstractState) = scale(c, s)
    Base.(:*)(s::AbstractState, c::Number) = scale(s, c)
    Base.(:/)(s::AbstractState, c::Number) = scale(s, 1/c)

    Base.ctranspose(k::Ket) = Bra(k)
    Base.ctranspose(b::Bra) = b.ket
    Base.norm(s::AbstractState) = sqrt(sum(v->v^2, values(dict(s))))

    QuBase.tensor{P}(kets::Ket{P}...) = Ket{P}(mergecart!(tensor_state, StateDict(), map(dict, kets)))
    QuBase.tensor{P}(bras::Bra{P}...) = Bra(tensor(map(ctranspose, bras)...))

    QuBase.normalize(s::AbstractState) = (1/norm(s))*s
    QuBase.normalize!(s::AbstractState) = scale!(1/norm(s), s)

    xsubspace(s::AbstractState, x) = filter((k,v)->sum(k)==x, s)
    switch(s::AbstractState, i, j) = maplabels(label->switch(label,i,j), s)
    permute(s::AbstractState, perm) = maplabels(label->permute(label,perm), s)
    switch!(s::AbstractState, i, j) = maplabels!(label->switch!(label,i,j), s)
    Base.permute!(s::AbstractState, perm::AbstractVector) = maplabels!(label->permute!(label,perm), s)

    filternz!(s::AbstractState) = filter!((k, v) -> v != 0, s)
    filternz(s::AbstractState) = filter((k, v) -> v != 0, s)

    # should always be pure, of course,
    # but makes a good sanity check function
    purity(ket::Ket) = purity(ket*ket')
    purity(bra::Bra) = purity(bra.ket)

######################
# Printing Functions #
######################
    labelstr(label) = join(map(repr, label), ',')
    ketstr(label) = "| $(labelstr(label)) $rang"
    brastr(label) = "$lang $(labelstr(label)) |"

    labelrepr(ket::Ket, label, pad) = "$pad$(ket[label]) $(ketstr(label))"
    labelrepr(bra::Bra, label, pad) = "$pad$(bra[label]) $(brastr(label))"

    Base.summary(s::AbstractState) = "$(typeof(s)) with $(length(s)) state(s)"
    Base.show(io::IO, s::AbstractState) = dirac_show(io, s)
    Base.showcompact(io::IO, s::AbstractState) = dirac_showcompact(io, s)
    Base.repr(s::AbstractState) = dirac_repr(s)

####################
# Helper Functions #
####################
    function tensor_state(pairs)
        #pairs structure is: ((label1, value1), (label2, value2)....,)
        return (vcat(map(first, pairs)...), prod(second, pairs))
    end

export Ket,
    Bra,
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
    coeffs