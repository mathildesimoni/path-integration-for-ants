using Revise
a = push!(LOAD_PATH, pwd()*"/src")
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, BumpAttractor
using Random, Distributions

# Network parameters
N = 300 # number of poisson neurons
x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution
R = 1 # resistance in Mohm
tau = 10.0 # characteristic time in ms
alpha = 2.0 # parameter for the transfer function in 1/mV
beta = 0.5 # parameter for the transfer function in mV
ro = 1 # parameter for the mean firing rate function in 1/ms
J = 5

# simulation parameters
delta_t = 0.1 # timestep for the simulation in ms. MUST BE <= 1
T = 1000 # simulation length in ms
n = Int64(T/delta_t)

# plot parameters
nb_ticks_x = 5

