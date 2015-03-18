###########
# Ket/Bra #
###########
    typealias StateDict Dict{Tuple,Number}

    type Ket{S} <: AbstractState{S}
        dict::StateDict
        Ket() = new(StateDict())
        Ket(dict::Dict) = new(dict)
        Ket(label...) = new(single_dict(StateDict(), label, 1))
    end

    Ket(items...) = Ket{Orthonormal}(items...)

    type Bra{S} <: AbstractState{S}
        ket::Ket{S}
        Bra(items...) = new(Ket{S}(items...))
        Bra(ket::Ket{S}) = new(ket)
    end

    Bra{S}(ket::Ket{S}) = Bra{S}(ket)
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
    Base.(:(==)){S}(a::Ket{S}, b::Ket{S}) = dict(a) == dict(b)
    Base.(:(==)){S}(a::Bra{S}, b::Bra{S}) = a.ket == b.ket
    Base.hash(s::AbstractState) = hash(dict(s), hash(typeof(s)))

    Base.length(s::AbstractState) = length(dict(s))

    Base.getindex(k::Ket, label::Tuple) = dict(k)[label]
    Base.getindex(b::Bra, label::Tuple) = b.ket[label]'
    Base.getindex(s::AbstractState, i...) = s[i]

    Base.setindex!(k::Ket, c, label::Tuple) = setindex!(dict(k), c, label)
    Base.setindex!(b::Bra, c, label::Tuple) = setindex!(b.ket, c', label)
    Base.setindex!(s::AbstractState, c, i...) = setindex!(s, c, i)

    Base.haskey(s::AbstractState, label::Tuple) = haskey(dict(s), label)
    Base.get(s::AbstractState, label::Tuple, default) = haskey(s, label) ? s[label] : default

    Base.delete!(s::AbstractState, label::Tuple) = (delete!(dict(s), label); return s)

##################################################
# Function-passing functions (filter, map, etc.) #
##################################################
    Base.filter!(f::Function, ket::Ket) = (filter!(f, dict(ket)); return ket)
    Base.filter!(f::Function, bra::Bra) = (filter!((k,v)->f(k,v'), bra.ket); return bra)

    Base.filter{S}(f::Function, ket::Ket{S}) = Ket{S}(filter(f, dict(ket)))
    Base.filter(f::Function, bra::Bra) = Bra(filter((k,v)->f(k,v'), bra.ket))

    Base.map{S}(f::Function, ket::Ket{S}) = Ket{S}(mapkv(f, dict(ket)))
    Base.map(f::Function, bra::Bra) = mapkv!((k,v)->f(k,v'), similar(bra), bra.ket)

    mapcoeffs{S}(f::Function, ket::Ket{S}) = Ket{S}(mapvals(f, dict(ket)))
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

    # This method is more optimized than the general inner() for
    # orthonormal states, but unfortunately sidesteps the error
    # condition for states of differing factor length...
    # function inner{A<:Orthogonal,B<:Orthogonal}(bra::Bra{A}, ket::Ket{B})
    #     result = 0
    #     if length(bra) < length(ket)
    #         for (label, c) in bra.ket
    #             if haskey(ket, label)
    #                 result += c'*ket[label]*inner_eval(A,B,label,label)
    #             end
    #         end
    #     else
    #         for (label, v) in ket
    #             if haskey(bra, label)
    #                 result += v*bra[label]*inner_eval(A,B,label,label)
    #             end
    #         end
    #     end
    #     return result
    # end

    Base.scale!(c::Number, s::AbstractState) = (castvals!(*, c, dict(s)); return s)
    Base.scale!(s::AbstractState, c::Number) = (castvals!(*, dict(s), c); return s)

    Base.scale(c::Number, s::AbstractState) = typeof(s)(castvals(*, c, dict(s)))
    Base.scale(s::AbstractState, c::Number) = typeof(s)(castvals(*, dict(s), c))

    Base.(:+){S}(a::Ket{S}, b::Ket{S}) = Ket{S}(mergef(+, dict(a), dict(b)))
    Base.(:-){S}(a::Ket{S}, b::Ket{S}) = a + (-b)
    Base.(:-){S}(ket::Ket{S}) = mapcoeffs(-, ket)

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

    QuBase.tensor{S}(kets::Ket{S}...) = Ket{S}(mergecart!(tensor_state, StateDict(), map(dict, kets)))
    QuBase.tensor{S}(bras::Bra{S}...) = Bra(tensor(map(ctranspose, bras)...))

    QuBase.normalize(s::AbstractState) = (1/norm(s))*s
    QuBase.normalize!(s::AbstractState) = scale!(1/norm(s), s)

    xsubspace(s::AbstractState, x) = filter((k,v)->sum(k)==x, s)
    switch(s::AbstractState, i, j) = maplabels(label->switch(label,i,j), s)
    permute(s::AbstractState, p::Array{Int}) = maplabels(label->permute(label,p), s)

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
        return (join_tup(map(first, pairs)), prod(second, pairs))
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