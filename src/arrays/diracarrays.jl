import Base: getindex,
    size,
    in,
    summary

# Dirac arrays are subtypes of AbstractQuArrays which generate
# ScaledStates and ScaledOperators as elements. These 
# quantum elements are formed by multiplying together 
# the underlying QuArray's basis states with the associated
# coefficients. 
#
# For example, the nth element of DiracVector{Ket, S} is 
# the nth coefficient times a DiracState{Ket, S} whose label
# is the nth label in the basis. 
#
# Likewise, the (ith, jth) element of DiracMatrix{Ket, S} is 
# the (ith, jth) coefficient times a DiracOperator{S} whose 
# ket label is the ith label in the row basis, and whose 
# bra label is the jth label in the column basis.
#
# This kind of structure means we can do two cool things:
# 
# 1) If operations are being done within a single basis, we can 
# bypass the performance cost of using states/basis and 
# just operate on the coefficient arrays. If we're doing 
# mixed basis operations, however, we can perform operations
# utilizing the bras and kets themselves. 
#
# 2) Storing a basis of labels allows for label-based methods 
# of analysis. An easy example is arbitrary selection/extraction 
# of subspaces using methods like `filter`. (not yet implemented
# here, but examples of which can be found in the old QuDirac repo).

abstract DiracArray{B, T<:AbstractDirac, N} <: AbstractQuArray{B, T, N}


###############
# DiracVector #
###############
    checksize(::Type{Ket}, qa) = size(qa, 2) == 1 
    checksize(::Type{Bra}, qa) = size(qa, 1) == 1 

    type DiracVector{D, 
                     S<:AbstractStructure, 
                     T, 
                     B<:AbstractLabelBasis, 
                     N, 
                     A} <: DiracArray{(B,), ScaledState{D, S, T}, N}
        quarr::QuVector{B, T, N, A}
        function DiracVector{L<:AbstractLabelBasis{S}}(quarr::QuVector{L, T, N, A})
            if checksize(D, quarr)
                new(quarr)
            else 
                error("Coefficient array does not conform to input bases")
            end
        end
    end


    function DiracVector{L<:AbstractLabelBasis,T,N,A}(quarr::QuVector{L,T,N,A}, D::DataType=Ket)
        return DiracVector{D, structure(L), T, L, N, A}(quarr)
    end

    function DiracVector{T,S}(
                        coeffs::AbstractArray{T}, 
                        basis::AbstractLabelBasis{S}, 
                        D::DataType=Ket)
        return DiracVector(QuArray(coeffs, basis), D)    
    end

    function DiracVector{K<:AbstractKet}(arr::AbstractArray{K})
        return DiracVector(map(coeff, arr), LabelBasis(map(state, arr)), Ket)
    end

    function DiracVector{B<:AbstractBra}(arr::AbstractArray{B})
        return DiracVector(map(coeff, arr), LabelBasis(map(state, arr)), Bra)
    end

    typealias KetVector{S<:AbstractStructure, T, B<:AbstractLabelBasis} DiracVector{Ket, S, T, B}
    typealias BraVector{S<:AbstractStructure, T, B<:AbstractLabelBasis} DiracVector{Bra, S, T, B}

    getbasis(dv::DiracVector) = getbasis(dv.quarr, 1)
    getcoeffs(dv::DiracVector) = getcoeffs(dv.quarr)

    ########################
    # Array-like Functions #
    ########################
    size(dv::DiracVector, i...) = size(dv.quarr, i...)

    in(s::AbstractState, dv::DiracVector) = in(label(s), getbasis(dv))
    in(c, dv::DiracVector) = in(c, dv.quarr)

    getpos(dv::DiracVector, s) = getpos(getbasis(dv), s)

    getcoeff(dv::DiracVector, s::AbstractState) = dv.quarr[getpos(dv, s)]
    getcoeff(dv::DiracVector, s::StateLabel) = dv.quarr[getpos(dv, s)]
    getcoeff(dv::DiracVector, s::Tuple) = dv.quarr[getpos(dv, s)]
    getcoeff(dv::DiracVector, i) = dv.quarr[i]

    getstate{D,S<:AbstractStructure}(dv::DiracVector{D,S}, i) = DiracState{D,S}(getbasis(dv)[i])

    getindex(dv::DiracVector, arr::AbstractArray) = DiracVector([dv[i] for i in arr])
    getindex(dv::DiracVector, i::Real) = getcoeff(dv, i) * getstate(dv, i)
    getindex(dv::DiracVector, i) = getcoeff(dv, i) * getstate(dv, i)
    getindex(dv::DiracVector, s::AbstractState) = getcoeff(dv, s) * s
    getindex(dv::DiracVector, label::StateLabel) = getcoeff(dv, label) * s
    getindex(dv::DiracVector, label::Tuple) = getcoeff(dv, label) * s
    
    getindex(dv::KetVector, i, j) = j==1 ? dv[i] : throw(BoundsError())
    getindex(dv::BraVector, i, j) = i==1 ? dv[j] : throw(BoundsError())

    ######################
    # Printing Functions #
    ######################
    summary{S<:AbstractStructure,T}(dv::KetVector{S,T}) = "KetVector{$S} with $(length(dv)) $T entries"
    summary{S<:AbstractStructure,T}(dv::BraVector{S,T}) = "BraVector{$S} with $(length(dv)) $T entries"

###############
# DiracMatrix #
###############
    type DiracMatrix{S<:AbstractStructure, 
                     T, 
                     R<:AbstractLabelBasis, 
                     C<:AbstractLabelBasis,
                     N,
                     A} <: DiracArray{(R,C), ScaledOperator{S, T}, N}
        quarr::QuMatrix{R, C, T, N, A}
        function DiracMatrix{R<:AbstractLabelBasis{S}, C<:AbstractLabelBasis{S}}(arr::QuMatrix{R, C, T, N, A})
            return new(quarr)
        end
    end

############################
# Convenience Constructors #
############################
    one_at_ind!(arr, i) = setindex!(arr, one(eltype(arr)), i)
    single_coeff(i, lens...) = one_at_ind!(zeros(Complex128, lens), i)
    diraccoeffs(i, len, ::Type{Ket}) = single_coeff(i, len)
    diraccoeffs(i, len, ::Type{Bra}) = single_coeff(i, 1, len)

    diracvec(coeffs::AbstractArray, D=Ket, S=AbstractStructure) = DiracVector(coeffs, FockBasis{S}(length(coeffs)), D)
    diracvec(i::Int, b::AbstractLabelBasis, D=Ket) = DiracVector(diraccoeffs(i, length(b), D), b, D)
    diracvec(tup::(Int...), b::AbstractLabelBasis, D=Ket) = DiracVector(diraccoeffs(getpos(b, tup), length(b), D), b, D)

    ketvec(s::(Int...), lens::Int...) = diracvec(s, FockBasis(lens), Ket)
    ketvec(s::(Int...)) = diracvec(s, FockBasis(map(x->x+1, s)), Ket)
    ketvec(s::Int, lens::Int...=s) = diracvec(s, FockBasis(lens), Ket)

    bravec(s::(Int...), lens::Int...) = diracvec(s, FockBasis(lens), Bra)
    bravec(s::(Int...)) = diracvec(s, FockBasis(map(x->x+1, s)), Bra)
    bravec(s::Int, lens::Int...=s) = diracvec(s, FockBasis(lens), Bra)

export DiracArray,
    DiracVector,
    DiracMatrix
    ketvec,
    bravec