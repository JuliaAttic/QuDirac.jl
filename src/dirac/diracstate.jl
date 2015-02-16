##############
# DiracState #
##############
    type DiracState{D,S} <: AbstractState{D,S}
        coeffs::ObjectIdDict
        DiracState() = new(ObjectIdDict())
        DiracState(coeffs) = new(coeffs)
    end

    typealias DiracKet{S} DiracState{Ket,S}
    typealias DiracBra{S} DiracState{Bra,S}

################
# Constructors #
################
    single_state{D,S}(::Type{D}, ::Type{S}, label::Tuple) = DiracState{D,S}(singlet_dict(label, 1))

    ket{S}(::Type{S}, label...) = single_state(Ket,S,label)
    ket(label...) = single_state(Ket,Orthonormal,label)
    ket{S}(::Type{S}) = DiracKet{S}()
    ket() = DiracKet{Orthonormal}()

    bra{S}(::Type{S}, label...) = single_state(Bra,S,label)
    bra(label...) = single_state(Bra,Orthonormal,label)
    bra{S}(::Type{S}) = DiracBra{S}()
    bra() = DiracBra{Orthonormal}()

    copy_type{D,S}(::DiracState{D,S}, new_coeffs) = DiracState{D,S}(new_coeffs)

    Base.copy(ds::DiracState) = copy_type(ds, copy(ds.coeffs))
    Base.similar(ds::DiracState) = copy_type(ds, ObjectIdDict())

#######################
# Dict-Like Functions #
#######################
    Base.length(ds::DiracState) = length(ds.coeffs)

    Base.getindex(ds::DiracState, label::Tuple) = ds.coeffs[label]
    Base.getindex(ds::DiracState, i...) = ds[i]

    Base.setindex!(ds::DiracState, c, label::Tuple) = setindex!(ds.coeffs, c, label)
    Base.setindex!(ds::DiracState, c, i...) = setindex!(ds, c, i)

    Base.start(ds::DiracState) = start(ds.coeffs)
    Base.done(ds::DiracState, state) = done(ds.coeffs, state)
    Base.next(ds::DiracState, state) = next(ds.coeffs, state)
    Base.endof(ds::DiracState) = endof(ds.coeffs)
    Base.last(ds::DiracState) = last(ds.coeffs)
    Base.first(ds::DiracState) = first(ds.coeffs)
    Base.collect(ds::DiracState) = collect(ds.coeffs)

    Base.get(ds::DiracState, label, default) = get(ds.coeffs, label, default)
    Base.haskey(ds::DiracState, label) = haskey(ds.coeffs, label)
    Base.keys(ds::DiracState) = keys(ds.coeffs)
    Base.values(ds::DiracState) = values(ds.coeffs)
    Base.filter!(f::Function, ds::DiracState) = (filter!(f, ds.coeffs); return ds)
    Base.filter(f::Function, ds::DiracState) = copy_type(ds, filter(f, ds.coeffs))
    Base.map(f::Function, ds::DiracState) = mapkv(f, ds)
    Base.delete!(ds::DiracState, label) = delete!(ds.coeffs, label)
    
    getstate{D,S}(ds::DiracState{D,S}, label::Tuple) = ds[label] * DiracState{D,S}(singlet_dict(label, 1))
    getstate(ds::DiracState, i...) = getstate(ds, i)

##########################
# Mathematical Functions #
##########################
    function inner{A,B}(db::DiracBra{A}, dk::DiracKet{B})
        result = 0
        for (b,c) in db
            for (k,v) in dk
                result += c*v*inner_eval(A,B,b,k) 
            end
        end
        return result
    end

    function inner{O<:Orthogonal}(db::DiracBra{O}, dk::DiracKet{O})
        if length(db) > length(dk)
            return ortho_inner(dk, db)
        else
            return ortho_inner(db, dk)
        end
    end

    +{D,S}(a::DiracState{D,S}, b::DiracState{D,S}) = DiracState{D,S}(mergef(+, a.coeffs, b.coeffs))
    -{D,S}(a::DiracState{D,S}, b::DiracState{D,S}) = a + (-b)
    -(ds::DiracState) = copy_type(ds, mapvals(-, ds.coeffs))
    *(a::DiracBra, b::DiracKet) = inner(a,b)
    *{D,S}(a::DiracState{D,S}, b::DiracState{D,S}) = tensor(a,b)
    *(c, ds::DiracState) = copy_type(ds, castvals(*, c, ds.coeffs))
    *(ds::DiracState, c) = copy_type(ds, castvals(*, ds.coeffs, c))
    /(ds::DiracState, c) = copy_type(ds, castvals(/, ds.coeffs, c))

    Base.conj(ds::DiracState) = mapvals(conj, ds)
    Base.ctranspose{D,S}(ds::DiracState{D,S}) = DiracState{D',S}(conj(ds.coeffs))
    Base.norm(ds::DiracState) = sqrt(sum(v->v^2, values(ds)))

    QuBase.tensor{D,S}(states::DiracState{D,S}...) = DiracState{D,S}(mergecart!(tensor_state, ObjectIdDict(), states))

    normalize(ds::DiracState) = (1/norm(ds))*ds

    xsubspace(ds::DiracState, x) = filter((k,v)->sum(k)==x, ds)
    matchlabel_at(ds::DiracState, x, y) = filter((k,v)-> x==k[y], ds)
    matchlabel_in(ds::DiracState, x) = filter((k,v)-> x in k, ds)

    switch(ds::DiracState, i, j) = mapkeys(k->switch(k,i,j), ds)
    permute(ds::DiracState, p::Array{Int}) = mapkeys(k->permute(k,p), ds)

######################
# Printing Functions #
######################
    labelstr(label) = strip(repr(label)[2:end-1], ',')
    statestr(label, ::Type{Ket}) = "| $(labelstr(label)) $rang"
    statestr(label, ::Type{Bra}) = "$lang $(labelstr(label)) |"

    function Base.show{D}(io::IO, ds::DiracState{D})
        print(io, "$(summary(ds)) with $(length(ds)) state(s):")
        pad = "  "
        maxlen = 30
        i = 1
        for (k,v) in ds
            if i <= maxlen
                println(io)
                print(io, "$pad$v $(statestr(k,D))")
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
    function ortho_inner(a, b)
        result = 0
        for (label, c) in a
            if haskey(b, label)
                result += c*b[label]
            end
        end
        return result
    end

    function tensor_state(pairs)
        #pairs structure is: ((label1, value1), (label2, value2)....,)
        return (join_tup(map(first, pairs)), prod(second, pairs))
    end

    mapkv(f::Function, ds::DiracState) = copy_type(ds, mapkv(f, ds.coeffs))
    mapvals(f::Function, ds::DiracState) = copy_type(ds, mapvals(f, ds.coeffs))
    mapkeys(f::Function, ds::DiracState) = copy_type(ds, mapkeys(f, ds.coeffs))

export DiracState,
    DiracKet,
    DiracBra,
    ket,
    bra,
    getstate,
    normalize,
    xsubspace,
    matchlabel_at,
    matchlabel_in,
    switch,
    permute