using Revise
a = push!(LOAD_PATH, pwd()*"/src")
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, BumpAttractor
using Random, Distributions

# Network parameters
N = 300 # number of poisson neurons
x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
R = 1
tau = 10.0
alpha = 2.0 
beta = 0.5 
ro = 1 
J = 3

# simulation parameters
delta_t = 0.1 # timestep for the simulation in ms. MUST BE <= 1
T = 1000 # simulation length in ms
n = Int64(T/delta_t)

# plot parameters
nb_ticks_x = 5
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/delta_t)

# new parameters
theta = deg2rad(10) # 10 degres
I_ext_bool = false # no external input
h_init_L = rand(Uniform(0,1), N) 
h_init_R = rand(Uniform(0,1), N)

# simulate the network activity
spikes_L, spikes_R = simulate_networks(h_init_L, h_init_R, x_i, N, delta_t, n, R, tau, I_ext_bool, J, alpha, beta, ro, theta)
spikes = (spikes_L + spikes_R)./2
heatmap(transpose(spikes), title="Network Activity", xlabel=L"t"*" (ms)", ylabel= "Neuron Location", c = reverse(cgrad(:grayC)), colorbar=false, right_margin = 3Plots.mm, left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), xticks = (Int.(0:n/nb_ticks_x:n), Int.(0:T/nb_ticks_x:T)))
bump_location_L = locate_bump.(eachrow(spikes_L), Ref(x_i))
bump_location_bins_L = transpose(reshape(bump_location_L[1:n], bin_length, Int((n)/bin_length)))
avg_bump_location_L = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins_L))
plot!(0:bin_length:n-1, avg_bump_location_L * (N/(2*pi)), label = "center of the bump")

# heatmap!(transpose(spikes_R), title="Network Activity", xlabel=L"t"*" (ms)", ylabel= "Neuron Location", c = reverse(cgrad(:grayC)), colorbar=false, right_margin = 3Plots.mm, left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), xticks = (Int.(0:n/nb_ticks_x:n), Int.(0:T/nb_ticks_x:T)))
bump_location_R = locate_bump.(eachrow(spikes_R), Ref(x_i))
bump_location_bins_R = transpose(reshape(bump_location_R[1:n], bin_length, Int((n)/bin_length)))
avg_bump_location_R = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins_R))
plot!(0:bin_length:n-1, avg_bump_location_R * (N/(2*pi)), label = "center of the bump")

last_mean = ((avg_bump_location_R + avg_bump_location_L)/2)[end]
print(last_mean)

savefig("data/Q22.pdf")