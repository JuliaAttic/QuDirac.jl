using QuDirac

# This file is an example of a discrete-time quantum random walk 
# implemented using QuDirac. 

# First, we define the Hadamard operator `H`, defined as
#
# H | 0 ⟩ = 1/√2 * (| 0 ⟩ + | 1 ⟩)
# H | 1 ⟩ = 1/√2 * (| 0 ⟩ - | 1 ⟩)
#
@rep_op " H | n > = 1/√2 * (| 0 > + (-1)^n * | 1 >)" 0:1

# If an operation on a state only alters basis labels, not coefficients,
# it can be simpler and more efficient to use `maplabels` than to construct 
# a whole new operator.
# 
# In this case, we make a function `shift_ket`
# that maps over the basis labels of the input Ket such that:
# 
# shift_ket(| 0, j ⟩) -> | 0, j - 1 ⟩
# shift_ket(| 1, j ⟩) -> | 1, j + 1 ⟩
#
function shift_map(label::StateLabel)
    if label[1] == 0
        return StateLabel(0, label[2] - 1)
    else
        return StateLabel(1, label[2] + 1)
    end
end

shift_ket(kt::Ket) = maplabels(shift_map, kt) 

function walk_nsteps(steps)
    steps += 1

    # Allocate space for the results
    results = cell(steps);

    # Initial state
    results[1] = d" 1.0 * | 0,0 > "

    # Calculate the each step by referring to the previous step
    for i=2:steps
        results[i] = shift_ket(act_on(H, results[i-1], 1))
    end

    return results
end

println(""" 
This example provides the `walk_nsteps` function, which takes in the 
number of steps to walk and returns an array containing the result of
each step. The initial state is | 0,0 ⟩.""")
