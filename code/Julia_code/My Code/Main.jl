using QuantumOptics
using PyPlot

using Random; Random.seed!(0)
η = 0.9 # Pumping strength
κ = 1 # Decay rate

Ncutoff = 20 # Maximum photon number
T = [0:0.1:10;];