##############
# DiracState #
##############
    typealias StateCoeffs Dict{Tuple,Number}

    type Ket{S} <: AbstractState{S}
        coeffs::StateCoeffs
        Ket() = new(StateCoeffs())
        Ket(coeffs::Dict) = new(coeffs)
        Ket(label...) = new(single_dict(StateCoeffs(), label, 1))
    end

    Ket(items...) = Ket{Orthonormal}(items...)

    type Bra{S} <: AbstractState{S}
        ket::Ket{S}
        Bra() = new(Ket())
        Bra(ket::Ket{S}) = new(ket)
        Bra(items...) = new(Ket{S}(items...))
    end

    Bra{S}(ket::Ket{S}) = Bra{S}(ket)
    Bra(items...) = Bra{Orthonormal}(items...)

    coeffs(k::Ket) = k.coeffs
    coeffs(b::Bra) = coeffs(b.ket)

################
# Constructors #
################
    Base.copy(s::AbstractState) = typeof(s)(copy(coeffs(s)))
    Base.similar(s::AbstractState) = typeof(s)(StateCoeffs())

#######################
# Dict-Like Functions #
#######################
    Base.(:(==)){S}(a::Ket{S}, b::Ket{S}) = coeffs(a) == coeffs(b)
    Base.(:(==)){S}(a::Bra{S}, b::Bra{S}) = coeffs(a) == coeffs(b)
    Base.hash(s::AbstractState) = hash(coeffs(s), hash(typeof(s)))

    Base.length(s::AbstractState) = length(coeffs(s))

    Base.getindex(k::Ket, label::Tuple) = coeffs(k)[label]
    Base.getindex(b::Bra, label::Tuple) = coeffs(b)[label]'
    Base.getindex(s::AbstractState, i...) = s[i]

    Base.setindex!(k::Ket, c, label::Tuple) = setindex!(coeffs(k), c, label)
    Base.setindex!(b::Bra, c, label::Tuple) = setindex!(coeffs(b), c', label)
    Base.setindex!(s::AbstractState, c, i...) = setindex!(s, c, i)

    Base.keys(k::Ket) = keys(coeffs(k))
    Base.values(k::Ket) = values(coeffs(k))

    Base.start(k::Ket) = start(coeffs(k))
    Base.next(k::Ket, state) = next(coeffs(k), state)
    Base.done(k::Ket, state) = done(coeffs(k), state)
    
    Base.haskey(s::AbstractState, label) = haskey(coeffs(s), label)
    Base.get(s::AbstractState, label, default) = haskey(s, label) ? s[label] : default

    Base.filter!(f::Function, s::AbstractState) = (filter!(f, coeffs(s)); return s)
    Base.filter(f::Function, s::AbstractState) = typeof(s)(filter(f, coeffs(s)))
    Base.map(f::Function, s::AbstractState) = typeof(s)(mapkv(f, coeffs(s)))
    Base.delete!(s::AbstractState, label) = (delete!(coeffs(s), label); return s)

##########################
# Mathematical Functions #
##########################
    function inner{A,B}(a::Ket{A}, b::Ket{B})
        result = 0
        for (bralabel,c) in a
            for (ketlabel,v) in b
                result += c'*v*inner_eval(A,B,bralabel,ketlabel) 
            end
        end
        return result
    end

    function inner{A<:Orthogonal,B<:Orthogonal}(a::Ket{A}, b::Ket{B})
        result = 0
        if length(a) < length(b)
            for (label, c) in a
                if haskey(b, label)
                    result += c'*b[label]*inner_eval(A,B,label,label)
                end
            end
        else
            for (label, v) in b
                if haskey(a, label)
                    result += v*a[label]'*inner_eval(A,B,label,label)
                end
            end
        end
        return result
    end

    Base.scale!(c::Number, s::AbstractState) = (castvals!(*, c, coeffs(s)); return s)
    Base.scale!(s::AbstractState, c::Number) = (castvals!(*, coeffs(s), c); return s)

    Base.scale(c::Number, s::AbstractState) = typeof(s)(castvals(*, c, coeffs(s)))
    Base.scale(s::AbstractState, c::Number) = typeof(s)(castvals(*, coeffs(s), c))

    +{S}(a::Ket{S}, b::Ket{S}) = Ket{S}(mergef(+, coeffs(a), coeffs(b)))
    -{S}(a::Ket{S}, b::Ket{S}) = a + (-b)
    -{S}(ket::Ket{S}) = Ket{S}(mapvals(-, coeffs(ket)))

    +(a::Bra, b::Bra) = Bra(a.ket+b.ket)
    -(a::Bra, b::Bra) = Bra(a.ket-b.ket)
    -(bra::Bra) = Bra(-bra.ket)

    *(bra::Bra, ket::Ket) = inner(bra.ket,ket)
    *(ket::Ket, bra::Bra) = outer(bra.ket,ket)
    *(a::Ket, b::Ket) = tensor(a,b)
    *(a::Bra, b::Bra) = tensor(a,b)

    *(c::Number, s::AbstractState) = scale(c, s)
    *(s::AbstractState, c::Number) = scale(s, c)
    /(s::AbstractState, c::Number) = scale(s, 1/c)

    Base.ctranspose(k::Ket) = Bra(k)
    Base.ctranspose(b::Bra) = b.ket
    Base.norm(s::AbstractState) = sqrt(sum(v->v^2, values(coeffs(s))))

    QuBase.tensor{S}(kets::Ket{S}...) = Ket{S}(mergecart!(tensor_state, StateCoeffs(), kets))
    QuBase.tensor{S}(bras::Bra{S}...) = Bra(tensor(map(ctranspose, bras)...))

    QuBase.normalize(s::AbstractState) = (1/norm(s))*s
    QuBase.normalize!(s::AbstractState) = scale!(1/norm(s), s)

    xsubspace(s::AbstractState, x) = filter((k,v)->sum(k)==x, s)
    matchlabel_at(s::AbstractState, x, y) = filter((k,v)-> x==k[y], s)
    matchlabel_in(s::AbstractState, x) = filter((k,v)-> x in k, s)

    switch(s::AbstractState, i, j) = mapkeys(k->switch(k,i,j), coeffs(s))
    permute(s::AbstractState, p::Array{Int}) = mapkeys(k->permute(k,p), coeffs(s))

######################
# Printing Functions #
######################
    labelstr(label) = strip(repr(label)[2:end-1], ',')
    ketstr(label) = "| $(labelstr(label)) $rang"
    brastr(label) = "$lang $(labelstr(label)) |"
    statestr(::Ket, label) = ketstr(label)
    statestr(::Bra, label) = brastr(label)

    function Base.show(io::IO, s::AbstractState)
        print(io, "$(typeof(s)) with $(length(s)) state(s):")
        pad = "  "
        maxlen = 30
        i = 1
        for (label,v) in coeffs(s)
            if i <= maxlen
                println(io)
                v = typeof(s)<:Bra ? v' : v
                print(io, "$pad$v $(statestr(s,label))")
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

    mapkv(f::Function, s::AbstractState) = typeof(s)(mapkv(f, coeffs(s)))
    mapvals(f::Function, s::AbstractState) = typeof(s)(mapvals(f, coeffs(s)))
    mapkeys(f::Function, s::AbstractState) = typeof(s)(mapkeys(f, coeffs(s)))

export Ket,
    Bra,
    normalize,
    xsubspace,
    matchlabel_at,
    matchlabel_in,
    switch,
    permute