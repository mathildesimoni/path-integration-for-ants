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
J = 5

# simulation parameters
delta_t = 0.1 # timestep for the simulation in ms. MUST BE <= 1
T = 1000 # simulation length in ms
T = 2000 # to show the drift that happens for larger tmescales
n = Int64(T/delta_t)

# plot parameters
nb_ticks_x = 5

# simulate the network activity
spikes = simulate_network(h_init, x_i, N, delta_t, n, R, tau, I_ext, J, alpha, beta, ro)
heatmap(transpose(spikes), title="Network Activity", xlabel=L"t"*" (ms)", ylabel= "Neuron Location", c = :grayC, colorbar=false, right_margin = 3Plots.mm, left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), xticks = (Int.(0:n/nb_ticks_x:n), Int.(0:T/nb_ticks_x:T)))

# find the location of the bump at every timestep
bump_location = locate_bump.(eachrow(spikes), Ref(x_i))
heatmap(transpose(spikes), title="Network Activity", xlabel=L"t"*" (ms)", ylabel= "neuron location", c = :grayC, colorbar=false, right_margin = 3Plots.mm, left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), xticks = (Int.(0:n/nb_ticks_x:n), Int.(0:T/nb_ticks_x:T)))
scatter!(0:1:n, bump_location * (N/(2*pi)), mc=:red, ms=2, ma=0.2, label = "center of the bump")

# average location of the bump over small bins of time (to have a clearer plot)
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/delta_t)
bump_location_bins = transpose(reshape(bump_location[1:n], bin_length, Int((n)/bin_length)))
# NOT THIS (their reshape function is weird): bump_location_bins = reshape(bump_location[1:n], Int((n)/bin_length), bin_length)
avg_bump_location = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins)) # need to use a circular mean method again
heatmap(transpose(spikes), title="Network Activity", xlabel=L"t"*" (ms)", ylabel= "Neuron Location", c = :grayC, colorbar=false, right_margin = 3Plots.mm, left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), xticks = (Int.(0:n/nb_ticks_x:n), Int.(0:T/nb_ticks_x:T)))
plot!(0:bin_length:n-1, avg_bump_location * (N/(2*pi)), label = "center of the bump")
savefig("data/Q12.pdf")
