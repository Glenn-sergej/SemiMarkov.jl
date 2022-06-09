using SpecialFunctions
using Unitful
using AdditionalUnits
using LinearAlgebra
using Plots

abstract type AbstractDistribution{N,R} end
abstract type AbstractStochasticProcess end
abstract type AbstractSemiMarkovProcess <: AbstractStochasticProcess end

include("distributions.jl")
include("set_semimarkov.jl")
include("set_A.jl")
include("set_U.jl")
include("set_p.jl")

"""
State transtion diagram of the example of application: reliability of downhole optical monitoring systems from: 
"Mathematical formulation and numerical treatment based on transition frequency densities and quadrature methods for non-homogeneous semi-Markov processes"
- Moura and Droguett
"""

StateTrans = [1 1; 1 2; 2 1; 2 3; 3 1; 3 4]

Settings = Dict("T" => 9000, "dT" => 2, "states" => 4, "trans" => size(StateTrans,1))
Init = [1.0,0.0,0.0,0.0]

t = 0:Settings["dT"]:Settings["T"];

#set empty array for constant probability
constp = ones(size(t))
p1 = constp*0.6.+(-0.00004*t)
p2 = constp*0.4.+(0.00004*t)

StateTransProb = [p1  p2 0.82*constp 0.18*constp 0.62*constp 0.38*constp]
dst = ["Exponential" 1/1e-4 0; "Weibull" 150 1.36; "Exponential" 1/0.05 0; "Lognormal" 2.5 0.25; "Exponential" 1/0.05 0; "Lognormal" 4.0 0.4]

semimarkov = set_semimarkov(Settings, StateTransProb, StateTrans, t, dst);

C = semimarkov[1];
dC = semimarkov[2];
F = semimarkov[3];
replace_nan!(dC);

A = set_A(Settings, StateTrans, Init, dC, t)

U = set_U(Settings, t, dC, StateTrans)

H = U\A

h = zeros(length(t),Settings["states"])

for ns in 1:Settings["states"]
    h[:,ns] = H[ns:Settings["states"]:end]
end

# Determine phi
phi = zeros(Settings["states"],length(t))

for ns in 1:Settings["states"]
    for rt in 1:length(t)
        phi[ns, rt] = set_p(Settings, rt, h, F, Init, ns)
    end
end


plot(t,phi[1,:])
