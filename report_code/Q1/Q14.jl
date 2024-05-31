using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, SingleBumpAttractor
using Random, Distributions
using Q1
using BumpAttractorUtils

N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t

x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution

# plot parameters
nb_ticks_x = 5

# plot I_ext for different times
t_values = [350, 500, 650]
p1 = plot()
plot!(xlabel=L"x_i", ylabel=L"I_{ext}")
for t in t_values
    I_ext_t = Q1.I_ext.(x_i, t)
    plot!(x_i, I_ext_t, label="t=$t", lw=2)
end
plot(p1)
savefig("data/Q14_I_ext.pdf")

# simulate the network activity
spikes = SingleBumpAttractor.simulate_network(h_init, x_i, Q1.I_ext, 0.0, sp, np)

# find the location of the bump at every timestep
bump_location = locate_bump.(eachrow(spikes), Ref(x_i))
# average location of the bump over small bins of time (to have a clearer plot)
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/delta_t)
bump_location_bins = transpose(reshape(bump_location[1:n], bin_length, Int((n)/bin_length)))
avg_bump_location = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins)) # need to use a circular mean method again
heatmap(transpose(spikes), title="Network Activity", xlabel=L"t"*" (ms)", ylabel= "Neuron Location", c = :grayC, colorbar=false, right_margin = 3Plots.mm, left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), xticks = (Int.(0:n/nb_ticks_x:n), Int.(0:T/nb_ticks_x:T)))
plot!(0:bin_length:n-1, avg_bump_location * (N/(2*pi)), label = "center of the bump")
# savefig("data/Q14.pdf")