# will take a little moment to install the packages on the first run
using Pkg
Pkg.add("Plots")
Pkg.add("LaTeXStrings") 

using Plots, LaTeXStrings

# include the functions from the other files
include("poisson_neurons.jl")

# ------------
# Question 0.1
# ------------

# define the parameters
h_array = range(-5, 5, length=1000)
alpha_values = [0.5 1.0 2.0 100.0]
beta_values = [-2.0 0.0 2.0]

# ranging alpha
beta = 0.0
p1 =plot()
# plot!(title = "Activation \$g\$ vs. Potential \$h\$", xlabel=L"x", ylabel=L"g(x)")
plot!(xlabel=L"h", ylabel=L"g")
for alpha in alpha_values
    # the dot . here is used to broadcast the function g over the array h_arr (g is defined for scalars)
    result = g.(h_array, alpha, beta)
    plot!(h_array, result, label=L"\alpha=%$(alpha)", lw=2)
end

# ranging beta
p2 =plot()
alph = 1.0
plot!(xlabel=L"h", ylabel=L"g")
for beta in beta_values
    result = g.(h_array, alph, beta)
    plot!(h_array, result, label=L"\beta=%$(beta)", lw=2)
end

# plot(p1, p2, layout=(2,1), size=(800, 800))
plot(p1)
plot(p2)

# save
# savefig("data/sigmoid_alpha.png")
# savefig("data/sigmoid_beta.png")

# plot(h, g_values, label=labels, lw=2)


# ------------
# Question 0.2
# ------------

# define the parameters
Io = 2.0 # parameter for the input function in nA
omega = 1.0 # parameter for the input function in 1/ms
alpha = 2.0 # parameter for the transfer function in 1/mV
beta = 0.5 # parameter for the transfer function in mV
tau = 10.0 # characteristic time in ms
R = 1 # resistance in Mohm
ro = 1 # parameter for the mean firing rate function in 1/ms

N = 100 # number of poisson neurons
h_init = 0.0 # initial potential in mV
T = 100 # simulation length in ms
delta_t = 0.1 # timestep for the simulation in ms. MUST BE <= 1
n = Int64(T/delta_t)
t = range(0, step = delta_t, length = n+1) 

# SIMULATION FOR 1 NEURON ========================

I_t = I.(t, Io, omega) # input to the neurons
plot(t, I_t, label="Input Current", lw=2, xlabel=L"t" * "(ms)", ylabel=L"I"*" (nA)")

h_arr_for_loop = h_for_loop(h_init, delta_t, n, R, I_t, tau) # evolution of neuron potential
plot(t, h_arr_for_loop, label="Neuron Potential", lw=2, xlabel=L"t" * " (ms)", ylabel=L"h" * " (mV)")

h_arr = h(h_init, delta_t, n, R, I_t, tau) # evolution of neuron potential
plot(t, h_arr, label="Neuron Potential", lw=2, xlabel=L"t" * " (ms)", ylabel=L"h" * " (mV)")

if h_arr_for_loop != h_arr
    throw(AssertionError("Problem in the computation of potential function!"))
end

r_arr = (r.(h_arr, alpha, beta, ro) * delta_t)[1:n]
plot(t[1:n], r_arr, label="Spike probability", lw=2, xlabel=L"t"*" (ms)")

# sample from the distribution
rand_nb = rand(Float64, n)
spikes = rand_nb .<= r_arr
plot(t[1:n], spikes, label=L"Spikes", lw=2, xlabel=L"t" * " (ms)")

# calculate total number of spikes
nb_spikes = sum(spikes)
println("Total number of spikes: ", nb_spikes)

# calculate average number of spikes per ms
bin_size = 1/delta_t
nb_spikes_per_ms = sum(reshape(spikes, Int64(n/bin_size), Int64(bin_size)), dims = 2)
avg_spikes = sum(nb_spikes_per_ms) / length(nb_spikes_per_ms)
println("Average number of spikes per ms: ", avg_spikes)

# simulation for N neurons
I_t = I.(t, Io, omega) # input to the neurons
plot(t, I_t, label=L"Input \: Current", lw=2, xlabel=L"t \: (ms)", ylabel=L"I \: (nA)")
h_init_N = h_init * ones(N)

h_arr = h(h_init, delta_t, n, R, I_t, tau) # evolution of neuron potential
plot(t, h_arr, label=L"Neuron \: Potential", lw=2, xlabel=L"t \: (ms)", ylabel=L"h \: (mV)")

r_arr = (r.(h_arr, alpha, beta, ro) * delta_t)[1:n]
plot(t[1:n], r_arr, label=L"Spike \: probability", lw=2, xlabel=L"t \: (ms)")

# sample from the distribution
rand_nb = rand(Float64, (n, N))

spikes = rand_nb .<= r_arr

# calculate average number of spikes per ms
bin_length = 1
bin_size = bin_length/delta_t
nb_spikes_per_ms = sum(reshape(spikes, (Int64(n/bin_size), Int64(bin_size), N)), dims = 2) ./ bin_length
avg_population_spikes_per_ms = sum(nb_spikes_per_ms, dims=3) / N
avg_population_spikes_per_ms = avg_population_spikes_per_ms[:,1,:]
t = range(0, stop = T, length = Int64(n/bin_size))
plot(t, avg_population_spikes_per_ms, label="Population Rate", lw=2, xlabel=L"t"*" (ms)")

# Compare with the theoretical instantaneous rate
t = range(0, stop = T, length = n+1)
plot(t, r.(t, alpha, beta, ro, Io, omega), label="Theoretical Rate", lw=2, xlabel=L"t"* " (ms)")
t = range(0, stop = T, length = Int64(n/bin_size))
plot!(t, avg_population_spikes_per_ms, label="Population Rate", lw=2, xlabel=L"t"*" (ms)")