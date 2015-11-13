using DiracNotation
using Base.Test

##############
# Test State #
##############
T = Float64
L = ASCIIString

c1 = rand()
c2 = rand()
c3 = rand(Complex{Float64})
lbl1 = "x"
lbl2 = "y"
lbl3 = UTF8String("z")
t1 = c1*term(lbl1)
t2 = c2*term(lbl2)
t3 = c3*term(lbl3)
s = t1 + t2

######################
# Property Functions #
######################
@test label(t1) == lbl1
@test coeff(t1) == c1
@test DiracNotation.data(s) == Dict(lbl1=>c1, lbl2=>c2)

@test LabelTerm{labeltype(t1), coefftype(t1)} == eltype(s) == LabelTerm{L,T}

########################
# Promotion/Conversion #
########################
LA, LB = Char, ASCIIString
LC = promote_type(LA, LB)

TA, TB = Complex{Int}, T
TC = promote_type(TA, TB)

@test promote_type(LabelTerm{LA,TA}, LabelTerm{LB,TB}) == LabelTerm{LC,TC}
@test promote_type(LabelTerm{LA,TA}, LabelSum{LB,TB}) == LabelSum{LC,TC}
@test promote_type(LabelSum{LA,TA}, LabelSum{LB,TB}) == LabelSum{LC,TC}

@test convert(LabelTerm{UTF8String, Complex{T}}, t1) == LabelTerm(UTF8String(lbl1), Complex{T}(c1))
@test convert(LabelSum{UTF8String, Complex{T}}, t1) == LabelSum(convert(LabelTerm{UTF8String, Complex{T}}, t1))
@test convert(LabelSum{UTF8String, Complex{T}}, s) == convert(LabelTerm{UTF8String, Complex{T}}, t1) + convert(LabelTerm{UTF8String, Complex{T}}, t2)

############################
# Hashing/Equality/Copying #
############################
t1_copy = copy(t1)
s_copy = copy(s)

@test t1_copy == t1
@test s_copy == s

t1_complex = convert(LabelTerm{UTF8String, Complex{T}}, t1)
s_complex = convert(LabelSum{UTF8String, Complex{T}}, s)

@test hash(t1_complex) == hash(t1)
@test hash(s_complex) == hash(s)

#########################
# Associative Functions #
#########################
@test length(s) == 2

@test s[lbl1] == c1
@test s[lbl2] == c2

@test (s_copy[lbl1] = c2; s_copy[lbl1] == c2)

@test haslabel(s, lbl1)
@test haslabel(s, lbl2)
@test !haslabel(s, lbl3)

default = rand()
@test get(s, lbl1, default) == c1
@test get(s, lbl2, default) == c2
@test get(s, lbl3, default) == default

@test !isempty(s)
@test isempty(empty!(s_copy))

s_copy = copy(s)

@test delete!(s_copy, lbl1) == t2

@test merge(s, t3) == s + t3
@test merge(s, s*s) == s + s*s

#############
# Iteration #
#############
@test s == LabelSum([i for i in s]...)

@test collect(s) == map(LabelTerm, collect(DiracNotation.data(s)))
@test collect(labels(s)) == collect(keys(DiracNotation.data(s)))
@test collect(coeffs(s)) == collect(values(DiracNotation.data(s)))

#####################
# Mapping/Filtering #
#####################
s_copy = copy(s)
filter!(t -> label(t) == lbl1, s_copy)
@test s_copy == t1

@test filter(t -> label(t) == lbl1, s) == t1

s_copy = copy(s)
s_copy[lbl1] = 0
filternz!(s_copy)
@test s_copy == t2

s_copy[lbl2] = 0
@test filternz(s_copy) == LabelSum()

sqr(i) = i*i
@test sqr(t1) + sqr(t2) == map(sqr, s)

##############
# Arithmetic #
##############

# Multiplication #
#----------------#
@test t1*t2 == c1*c2*term(lbl1*lbl2)
@test t2*t1 == c2*c1*term(lbl2*lbl1)
@test t1*s == t1*t1 + t1*t2
@test s*t1 == t1*t1 + t2*t1
@test s*s == t1*t1 + t1*t2 + t2*t1 + t2*t2
@test -t1 == -c1*term(lbl1)
@test -s == -t1 - t2
@test c3*s == c3*t1 + c3*t2

# Division #
#----------#
@test t1/c3 == inv(c3)*t1
@test s/c3 == inv(c3)*s

# Addition #
#----------#
@test t1 + t2 + t3 == LabelSum(t1, t2, t3)
@test t1 + t3 + t3 == t1 + 2*t3
@test s + t1 == 2*t1 + t2
@test s + s == 2*s
@test s + t3 == t1 + t2 + t3

# Subtraction #
#-------------#
@test t1 - t2 - t3 == LabelSum(t1, -t2, -t3)
@test filternz!(t1 + t3 - t3) == t1
@test filternz!(s - t1) == t2
@test filternz!(s - s) == LabelSum()
@test filternz!(s - t3) == filternz!(t1 + t2 - t3)
