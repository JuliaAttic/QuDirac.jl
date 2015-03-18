###########
# Ket/Bra #
###########
    DEFAULT_INNER() = Orthonormal

    typealias StateDict Dict{Vector,Number}

    type Ket{P<:AbstractInner} <: AbstractState{P}
        dict::StateDict
        Ket() = new(StateDict())
        Ket(dict::Dict) = new(dict)
        Ket(label...) = new(single_dict(StateDict(), collect(label), 1))
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
    Base.get(s::AbstractState, label::Array, default) = haskey(s, label) ? s[label] : default

    Base.delete!(s::AbstractState, label::Array) = (delete!(dict(s), label); return s)

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!(f::Function, ket::Ket) = (filter!(f, dict(ket)); return ket)
    Base.filter!(f::Function, bra::Bra) = (filter!((k,v)->f(k,v'), bra.ket); return bra)

    Base.filter{P}(f::Function, ket::Ket{P}) = Ket{P}(filter(f, dict(ket)))
    Base.filter(f::Function, bra::Bra) = Bra(filter((k,v)->f(k,v'), bra.ket))

    Base.map{P}(f::Function, ket::Ket{P}) = Ket{P}(mapkv(f, dict(ket)))
    Base.map(f::Function, bra::Bra) = mapkv!((k,v)->f(k,v'), similar(bra), bra.ket)

    mapcoeffs{P}(f::Function, ket::Ket{P}) = Ket{P}(mapvals(f, dict(ket)))
    mapcoeffs(f::Function, bra::Bra) = mapvals!(v->f(v'), similar(bra), bra.ket)
    maplabels(f::Function, s::AbstractState) = typeof(s)(mapkeys(f, dict(s)))

##########################
# Mathematical Functions #
##########################
    function inner{A,B}(bra::Bra{A}, ket::Ket{B})
        result = 0
        for (b,c) in dict(bra)
            for (k,v) in dict(ket)
                result += c'*v*inner_eval(A,B,b,k) 
            end
        end
        return result
    end

    function inner{A<:Orthonormal,B<:Orthonormal}(bra::Bra{A}, ket::Ket{B})
        result = 0
        if length(bra) < length(ket)
            for label in keys(dict(bra))
                if haskey(ket, label)
                    result += bra[label]*ket[label]*inner_eval(A,B,label,label)
                end
            end
        else
            for label in keys(dict(ket))
                if haskey(bra, label)
                    result += ket[label]*bra[label]*inner_eval(A,B,label,label)
                end
            end
        end
        return result
    end

    Base.scale!(c::Number, s::AbstractState) = (castvals!(*, c, dict(s)); return s)
    Base.scale!(s::AbstractState, c::Number) = (castvals!(*, dict(s), c); return s)

    Base.scale(c::Number, s::AbstractState) = typeof(s)(castvals(*, c, dict(s)))
    Base.scale(s::AbstractState, c::Number) = typeof(s)(castvals(*, dict(s), c))

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

    filternz!(s::AbstractState) = filter!((k, v) -> v != 0, s)
    filternz(s::AbstractState) = filter((k, v) -> v != 0, s)

    purity(ket::Ket) = trace((ket*ket')^2)
    purity(bra::Bra) = purity(bra')

######################
# Printing Functions #
######################
    labelstr(label) = strip(repr(label)[2:end-1], ',')
    ketstr(label) = "| $(labelstr(label)) $rang"
    brastr(label) = "$lang $(labelstr(label)) |"
    statestr{K<:Ket}(label, ::Type{K}) = ketstr(label)
    statestr{B<:Bra}(label, ::Type{B}) = brastr(label)

    function Base.show(io::IO, s::AbstractState)
        print(io, "$(typeof(s)) with $(length(s)) state(s):")
        pad = "  "
        maxlen = 30
        i = 1
        for label in keys(dict(s))
            if i <= maxlen
                println(io)
                print(io, "$pad$(s[label]) $(statestr(label, typeof(s)))")
                i = i + 1
            else  
                println(io)
                print(io, "$pad$vdots")
                break
            end
        end
    end

####################
# Helper Functions #
####################
    function tensor_state(pairs)
        #pairs structure is: ((label1, value1), (label2, value2)....,)
        return (vcat(map(first, pairs)...), prod(second, pairs))
    end

export Ket,
    Bra,
    maplabels,
    mapcoeffs,
    xsubspace,
    switch,
    permute,
    filternz!,
    filternz,
    purity