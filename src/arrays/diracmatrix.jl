###############
# DiracMatrix #
###############
    type DiracMatrix{S<:AbstractStructure, 
                     T, 
                     B<:AbstractLabelBasis, 
                     C<:AbstractLabelBasis,
                     A} <: DiracArray{B, ScaledOperator{S, T}, 2}
        quarr::QuMatrix{B, T, A}
        function DiracMatrix{B<:AbstractLabelBasis{S}}(arr::QuMatrix{B, T, A})
            return new(quarr)
        end
    end

export DiracMatrix