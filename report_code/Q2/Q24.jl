using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Random, Distributions
using CoupledBumpAttractors, Neurons, BumpAttractorUtils
using Q2

sp.delay = 500
sp.rate_neurons = true

N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t
x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)

# plot parameters
nb_ticks_x = 5
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/delta_t)

# coupled bump attractor parameters
theta = deg2rad(10) # 10 degres
h_init_L = rand(Uniform(0,1), N) 
h_init_R = rand(Uniform(0,1), N)

# external input
Io_values = range(-1.5, 1.5, 101)

last_means = zeros(size(Io_values)[1])
for (i, Io) in enumerate(Io_values)
    I_ext_L(x,t) = Q2.I_ext(x,t,-Io)
    I_ext_R(x,t) = Q2.I_ext(x,t,Io)
    spikes_L, spikes_R = CoupledBumpAttractors.simulate_network(h_init_L, h_init_R, I_ext_L, I_ext_R, x_i, theta, sp, np)
    avg_bump_location_L = Utils.spikes_to_average_bump_location(spikes_L, x_i, bin_size, sp)
    avg_bump_location_R = Utils.spikes_to_average_bump_location(spikes_R, x_i, bin_size, sp)
    avg_bump_location = (avg_bump_location_L + avg_bump_location_R) / 2
    last_means[i] = avg_bump_location[end]
end

p = plot()
p = Utils.plot_segments(collect(Io_values), last_means, tol = 1, xlabel = L"I_o", ylabel = L"\theta_{mean}")
yticks!(range(start = minimum(last_means), stop = maximum(last_means) , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"])
vline!([-0.35,0.35], legend = false, color = :Purple)
display(p)
savefig("./data/Q24.pdf")

