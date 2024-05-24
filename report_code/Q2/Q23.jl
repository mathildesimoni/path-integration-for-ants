using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Random, Distributions
using CoupledBumpAttractors, Neurons, BumpAttractorUtils
using Q2

sp.delay = 500
N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t

x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)

# plot parameters
nb_ticks_x = 5
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/delta_t)

# new parameters
theta = deg2rad(10) # 10 degres
I_ext_bool = true # no external input
h_init_L = rand(Uniform(0,1), N) 
h_init_R = rand(Uniform(0,1), N)

spikes_L, spikes_R = CoupledBumpAttractors.simulate_network(h_init_L, h_init_R, Q2.I_ext, Q2.I_ext, x_i, theta, sp, np)
spikes = (spikes_L + spikes_R)./2
# heatmap(transpose(spikes), title="Network Activity", xlabel=L"t"*" (ms)", ylabel= "Neuron Location", c = reverse(cgrad(:grayC)), colorbar=false, right_margin = 3Plots.mm, left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), xticks = (Int.(0:n/nb_ticks_x:n), Int.(0:T/nb_ticks_x:T)))

bump_location_L = locate_bump.(eachrow(spikes_L), Ref(x_i))
bump_location_bins_L = transpose(reshape(bump_location_L[1:n], bin_length, Int((n)/bin_length)))
avg_bump_location_L = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins_L))

bump_location_R = locate_bump.(eachrow(spikes_R), Ref(x_i))
bump_location_bins_R = transpose(reshape(bump_location_R[1:n], bin_length, Int((n)/bin_length)))

avg_bump_location_R = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins_R))

avg_bump_location = (avg_bump_location_L + avg_bump_location_R) / 2


heatmap(transpose(spikes_R), 
    title="Network Activity", 
    xlabel=L"t"*" (ms)", 
    ylabel= "Neuron Location", 
    c = reverse(cgrad(:grayC)), 
    colorbar=false, 
    right_margin = 3Plots.mm, 
    left_margin = 2Plots.mm, 
    yticks = (range(start = 0, stop = N , length =5), 
    [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), 
    xticks = (Int.(0:n/nb_ticks_x:n), 
    Int.(0:T/nb_ticks_x:T)))

plot!(0:bin_length:n-1, avg_bump_location_L * (N/(2*pi)), label = "center of the bump")
plot!(0:bin_length:n-1, avg_bump_location_R * (N/(2*pi)), label = "center of the bump")

# savefig("data/Q23.pdf")