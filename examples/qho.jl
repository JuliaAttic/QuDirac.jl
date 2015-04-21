warn("This example requires the Plotly package. See here for more: https://plot.ly/julia/getting-started/")

using QuDirac

####################
# Ladder Operators #
####################
# These functions can be applied to Kets like "lower * k"
@def_op " lower | n > = √n * | n - 1 > "
@def_op " raise | n > = √(n + 1) * | n + 1 > "

#################################
# Hermite Polynomial Evaluation #
#################################
# This is a very naive implementation, 
# likely prone to numerical error for high n. 
# It should suffice for this example, however.
# Implements the expression found here: 
# http://en.wikipedia.org/wiki/Hermite_polynomials#Recursion_relation
function hermite(n, x)
    if n == 0
        return 1
    elseif n == 1
        return 2 * x
    else
        return 2x*hermite(n-1, x) - 2*(n-1)*hermite(n-2, x)
    end
end 

##################
# QHOInner Setup #
##################
# Custom inner product type
# the QHO problem
immutable QHOInner <: AbstractInner end

qho_inner(x::Float64, k::Int) = e^((-x^2)/2) * hermite(k, x) * 1/√(2^k * factorial(k) * √π) # using natural units
qho_inner(n::Int, m::Int) = n==m ? 1.0 : 0.0
qho_inner(pair) = qho_inner(pair[1], pair[2])

QuDirac.inner_rule(::QHOInner, b, k) = mapreduce(qho_inner, *, zip(b, k))
QuDirac.inner_rettype(::QHOInner) = Float64

default_inner(QHOInner())

# With the above, we've defined this behavior for basis states:
#
# < i::Int | * | j::Int > = ∫ ψᵢ'ψⱼ dx =  δᵢⱼ
#
# < x::Float64 | * | n::Int > = ψᵢ(x)

####################
# Make Some Plots! #
####################
# Here we're going to make some plots 
# with Plotly. You'll have to have an 
# account to do this, but it's free and
# totally awesome.

# Given an iterable of x and y points, generate a
# distribution for the state by taking the inner product
gen_z{P}(kt::Ket{P,2}, xpoints, ypoints, s) = inner_rettype(QHOInner())[s*d" < x, y | * kt " for x in xpoints, y in ypoints]

# Generate the distribution above, and package it for a Plotly surface plot
function gen_plot_data{P}(kt::Ket{P,2}, xpoints, ypoints, s=1.0)
    return [
      [
        "z" => gen_z(kt, xpoints, ypoints, s), 
        "x" => xpoints, 
        "y" => ypoints, 
        "type" => "surface"
      ]
    ]
end

# some default stuff
len = 50
max = pi
xpoints = linspace(-max, max, len)
ypoints = copy(xpoints)

info("Loading the Plotly package, this could take a little while if it has to sign in...")

using Plotly

# Generate the distribution, sending it to Plotly.
# Return the response URL, which you can then go to
# to see and interact with your plot.
function plot_wave2D{P}(kt::Ket{P,2}, xpoints, ypoints, s=1.0)
    response = Plotly.plot(gen_plot_data(kt, xpoints, ypoints, s))
    return response["url"]
end

info("Finished loading Plotly. Make sure you're signed in before trying to plot!")

println("""

To generate a plot of a wave function for a 2-factor Ket, 
just call plot_wave2D(kt, xpoints, ypoints). This function
will build your plot and return the URL you should go to to
view the result. Here are some examples: 

# Basis state:
julia> plot_wave2D(d\" | 1, 1 > \", xpoints, ypoints)

# Random superposition of the first 4 basis states:
julia> randkt = normalize!(sum(i -> rand() * ket(i), 0:3))^2
       plot_wave2D(randkt, xpoints, ypoints)""")
