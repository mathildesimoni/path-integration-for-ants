using Revise
a = push!(LOAD_PATH, pwd()*"/src")
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons

# SIMULATION FOR 1 NEURON =====================================================

# define the parameters
Io = 2.0 # parameter for the input function in nA
omega = 1.0 # parameter for the input function in 1/ms
alpha = 2.0 # parameter for the transfer function in 1/mV
beta = 0.5 # parameter for the transfer function in mV
tau = 0.01 # characteristic time in ms
# tau = 10.0 # characteristic time in ms
R = 1 # resistance in Mohm
ro = 1 # parameter for the mean firing rate function in 1/ms

N = 1000000 # number of poisson neurons
h_init = 0.0 # initial potential in mV
T = 10 # simulation length in ms
delta_t = 0.01 # timestep for the simulation in ms. MUST BE <= 1
n = Int64(T/delta_t)
t = range(0, step = delta_t, length = n+1) 

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


