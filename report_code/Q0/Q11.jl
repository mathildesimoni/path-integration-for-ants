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
I_ext = false # no external input

# simulation parameters
delta_t = 0.1 # timestep for the simulation in ms. MUST BE <= 1
T = 1000 # simulation length in ms
n = Int64(T/delta_t)

# plot parameters
nb_ticks_x = 5

# tryout with a specific J
J = 5
spikes = simulate_network(h_init, x_i, N, delta_t, n, R, tau, I_ext, J, alpha, beta, ro)
p1 = heatmap(transpose(spikes), 
             title="Network Activity", 
             xlabel=L"t"*" (ms)", 
             ylabel= "Neuron", 
             c = reverse(cgrad(:grayC)),
             colorbar=false, 
             right_margin = 5Plots.mm, 
             left_margin = 2Plots.mm, 
             yticks = (range(start = 0, stop = N , length =5), 
             [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), 
             xticks = (Int.(0:n/nb_ticks_x:n), 
             Int.(0:T/nb_ticks_x:T)))

# ranging J
# J_values = [2.0, 3.0, 5.0, 10.0]
# for J in J_values
#     println(J)
#     spikes = simulate_network(h_init, x_i, N, delta_t, n, R, tau, I_ext, J, alpha, beta, ro)
#     heatmap(transpose(spikes), title="Network Activity", xlabel=L"t"*" (ms)", ylabel= "Neuron", c = :grayC, colorbar=false, right_margin = 5Plots.mm, left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), xticks = (Int.(0:n/nb_ticks_x:n), Int.(0:T/nb_ticks_x:T)))
# end

