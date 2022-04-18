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
State transtion diagram of the Semi-Analytical example from: 
"Mathematical formulation and numerical treatment based on transition frequency densities and quadrature methods for non-homogeneous semi-Markov processes"
- Moura and Droguett
"""

StateTrans = [1 2; 2 3]

Settings = Dict("T" => 4500, "dT" => 3, "states" => 3, "trans" => size(StateTrans,1))
Init = [1.0,0.0,0.0]

t = 0:Settings["dT"]:Settings["T"];

#set empty array for constant probability
constp = ones(size(t))

StateTransProb = [1*constp 1*constp]
dst = ["Exponential" 1/1e-3 0; "Weibull" 250 1.5]

semimarkov = set_semimarkov(Settings, StateTransProb, StateTrans, t, dst);

C = semimarkov[1];
dC = semimarkov[2];
F = semimarkov[3];

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

phi_r = phi[1,:].+phi[2,:];

plot(t,phi_r)




