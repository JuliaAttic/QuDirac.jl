using QuDirac

# Hadamard operator `H` flips the coin; defined as
#
# F | 0 ⟩ = 1/√2 * ( | 0 ⟩ + | 1 ⟩ )
# F | 1 ⟩ = 1/√2 * ( | 0 ⟩ - | 1 ⟩ )
#
const H = func_op(label -> d" 1/√2 * ( | 0 > + (-1)^label[1] * | 1 > ) ", [StateLabel(0), StateLabel(1)]);

# If an operation on a state only alters basis labels, not coefficients,
# it can be simpler and more efficient to use `maplabels` than to construct 
# a whole new operator.
# 
# In this case, we make a function `shift_ket`
# that maps over the basis labels of the input Ket such that:
# 
# shift_state(| 0, j ⟩) -> | 0, j - 1 ⟩
# shift_state(| 1, j ⟩) -> | 0, j + 1 ⟩
#
function shift_map(label::StateLabel)
    if label[1] == 0
        return StateLabel(label[1], label[2] - 1)
    else
        return StateLabel(label[1], label[2] + 1)
    end
end

shift_ket(kt::Ket) = maplabels(shift_map, kt) 

function walk_nsteps(steps)
    steps += 1

    # Allocate space for the results
    results = Array(Ket{KroneckerDelta, 2, Float64}, steps);

    # Initial state
    results[1] = d" 1.0 * | 0,0 > "

    # Calculate the each step by referring to the previous step
    for i=2:steps
        results[i] = shift_ket(act_on(H, results[i-1], 1))
    end

    return results
end
