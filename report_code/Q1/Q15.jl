using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, SingleBumpAttractor
using Random, Distributions
using Q1

N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t

x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution

# plot parameters
nb_ticks_x = 5

# new parameter for the connectivity function
phi = 0.05

# no external input
I_ext(x::Real, t::Real) = 0.0

# simulate the network activity
spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, phi, sp, np)
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
heatmap(transpose(spikes), 
        title="Network Activity", 
        xlabel=L"t"*" (ms)", 
        ylabel= "Neuron Location", 
        c = :grayC, 
        colorbar=false, 
        right_margin = 3Plots.mm, 
        left_margin = 2Plots.mm, 
        yticks = (range(start = 0, stop = N , length =5), 
        [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), 
        xticks = (Int.(0:n/nb_ticks_x:n), 
        Int.(0:T/nb_ticks_x:T)))
plot!(0:bin_length:n-1, avg_bump_location * (N/(2*pi)), label = "center of the bump")
# savefig("data/Q15.pdf")