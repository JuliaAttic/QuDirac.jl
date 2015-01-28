import Base: length,
    show

type DiracSet{T}
    dict::OrderedDict{Uint64,T}
    DiracSet(dict) = new(dict)
    DiracSet() = new(OrderedDict(Uint64, T))
end

function loaddict!(dict, arr)
    for i in arr
        dict[labelhash(i)] = i
    end
    return dict
end

DiracSet{T}(arr::AbstractArray{T}) = DiracSet{T}(loaddict!(OrderedDict(Uint64, T), arr))
DiracSet{T}(arr::T...) = DiracSet{T}(loaddict!(OrderedDict(Uint64, T), arr))

function add{D,S,A,B}(a::DiracSet{DiracState{D,S,A}}, b::DiracSet{DiracState{D,S,B}})
    result = DiracSet{DiracState{D,S,promote_type(A,B)}}()
    for av in values(a)
        for bv in values(b)

        end
    end
end

length(ds::DiracSet) = length(ds.dict)

function show(io::IO, ds::DiracSet)
    print(io, "$(summary(ds)):")
    pad = "  "
    maxlen = 30
    i = 1
    for v in values(ds.dict)
        if i <= maxlen
            println(io)
            print(io, "$pad$v")
            i = i + 1
        else  
            println(io)
            print(io, "$pad$vdots")
            break
        end
    end
end


function collide!(f::Function, d::Associative, others::Associative...)
    for other in others
        for (k,v) in other
            if haskey(d, k)
                d[k] = f(d[k], v)
            else   
                d[k] = v
            end
        end
    end
    return d
end

function collide(f::Function, d::Associative, others::Associative...)
    K, V = eltype(d)
    for other in others
        (Ko, Vo) = eltype(other)
        K = promote_type(K, Ko)
        V = promote_type(V, Vo)
    end
    collide!(f, Dict{K,V}(), d, others...)
end


export DiracSet,
    add