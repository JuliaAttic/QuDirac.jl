################
# InnerProduct #
################
immutable InnerProduct{P,B<:AbstractBra,K<:AbstractKet}
    br::B
    kt::K
    InnerProduct(br::AbstractBra{P}, kt::AbstractKet{P}) = new(br, kt)
end

function InnerProduct{P}(br::AbstractBra{P}, kt::AbstractKet{P})
    return InnerProduct{P,typeof(br),typeof(kt)}(br, kt)
end

Base.copy(i::InnerProduct) = i
Base.conj(i::InnerProduct) = InnerProduct(i.k', i.b')
Base.ctranspose(i::InnerProduct) = conj(i)
