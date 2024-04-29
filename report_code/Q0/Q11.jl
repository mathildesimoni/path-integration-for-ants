using Revise
a = push!(LOAD_PATH, pwd()*"/src")
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons
using Random, Distributions

# Network parameters
N = 300 # number of poisson neurons
x_i = range(start = 0, stop = 2*pi, length = N + 1)[1:N] # equally spaced neurons over the range [0, 2pi)
I_ext = 0.0 # no external input
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution

# simulation parameters
J_values = [2.0, 3.0, 5.0, 10.0]