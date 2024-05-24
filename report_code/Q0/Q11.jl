using Revise
a = push!(LOAD_PATH, pwd()*"/src")
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, SingleBumpAttractor, BumpAttractorUtils
using Random, Distributions


T = 1000.0 # simulation length in ms
delta_t = 0.1 # timestep for the simulation in ms. MUST BE <= 1
n = Int64(T/delta_t)
N = 300

sp = SimulationParameters(
    # Network parameters
    T = T, # simulation length in ms
    I_ext_bool = false, # no external input
    n = n,
    delta_t = delta_t
)

np = NetworkParameters(
    N = N, # number of poisson neurons
    tau = 10.0, # characteristic time in ms
    R = 1, # resistance in Mohm
    alpha = 2.0, # parameter for the transfer function in 1/mV
    beta = 0.5, # parameter for the transfer function in mV
    ro = 1, # parameter for the mean firing rate function in 1/ms
    J = 5
)    

h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution
x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)

# simulation parameters

# plot parameters
nb_ticks_x = 5

spikes = SingleBumpAttractor.simulate_network(h_init, x_i, 0, sp, np)
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

