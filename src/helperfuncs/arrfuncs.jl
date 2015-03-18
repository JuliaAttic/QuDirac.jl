############################
# Combinatorical Functions #
############################
second(t) = t[2]
except(arr, i) = deleteat!(copy(arr), i)
function switch!(arr, i, j)
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp
    return arr
end
switch(arr, i, j) = switch!(copy(arr), i, j)
permute(arr, perm) = permute!(copy(arr), perm)

# Placing here for now - may be useful for
# subscripting label factors in the future
digit_subscript(i::Int) = Base.REPLCompletions.latex_completions("\\_$i",3)[2][1][1]
subscript(i::Int) = reduce(*, reverse(map(digit_subscript, digits(i))))