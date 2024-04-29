using Revise
a = push!(LOAD_PATH, pwd()*"/src")
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons

# SIMULATION FOR N NEURONS ========================================================

# define the parameters
R = 1 # resistance in Mohm
ro = 1 # parameter for the mean firing rate function in 1/ms
alpha = 2.0 # parameter for the transfer function in 1/mV
beta = 0.5 # parameter for the transfer function in mV
delta_t = 0.1 # timestep for the simulation in ms. MUST BE <= 1
tau = 10.0 # characteristic time in ms

N = 100 # number of poisson neurons
T = 1000 # simulation length in ms
Io = 2.0 # parameter for the input function in nA
omega = 1.0 # parameter for the input function in 1/ms
h_init = 0.0 # initial potential in mV

bin_length = 1 # in ms
bin_size = bin_length/delta_t # number of timesteps in a bin
n = Int64(T/delta_t)

spikes = simulate_spikes(N, h_init, delta_t, T, R, Io, omega, alpha, beta, tau, ro)
avg_population_spikes_per_ms = compute_rate(spikes, delta_t, bin_length, N, n)

# Increase N
new_N = 1000
new_spikes = simulate_spikes(new_N, h_init, delta_t, T, R, Io, omega, alpha, beta, tau, ro)
new_avg_population_spikes_per_ms = compute_rate(new_spikes, delta_t, bin_length, new_N, n)

bin_t = range(0, T, Int64(T/bin_length))
plot(bin_t, avg_population_spikes_per_ms, title="Population firing rate", lw=1, xlabel=L"t"*" (ms)", label=L"N = %$(N)")
plot!(bin_t, new_avg_population_spikes_per_ms, lw=1, xlabel=L"t"*" (ms)", label=L"N = %$(new_N)")
theoretical_t = range(0, T, n+1)
plot!(theoretical_t, theoretical_r.(theoretical_t, alpha, beta, ro, Io, omega), lw=1, xlabel=L"t"* " (ms)", label="Theoretical", alpha=0.2)
savefig("data/Q02_N$(N)_vs_N$(new_N)_rate_vs_theoretical_T$(T).pdf")

# Zoom in 
new_T = 50
n = Int64(new_T/delta_t)
spikes = simulate_spikes(N, h_init, delta_t, new_T, R, Io, omega, alpha, beta, tau, ro)
avg_population_spikes_per_ms = compute_rate(spikes, delta_t, bin_length, N, n)
new_spikes = simulate_spikes(new_N, h_init, delta_t, new_T, R, Io, omega, alpha, beta, tau, ro)
new_avg_population_spikes_per_ms = compute_rate(new_spikes, delta_t, bin_length, new_N, n)
theoretical_t = range(0, stop = new_T, length = n+1)
bin_t = range(0, new_T, Int64(new_T/bin_length))
plot(bin_t, avg_population_spikes_per_ms, title="Population firing rate", lw=2, xlabel=L"t"*" (ms)", label=L"N = %$(N)")
plot!(bin_t, new_avg_population_spikes_per_ms, lw=2, xlabel=L"t"*" (ms)", label=L"N = %$(new_N)")
theoretical_t = range(0, new_T, n+1)
plot!(theoretical_t, theoretical_r.(theoretical_t, alpha, beta, ro, Io, omega), lw=2, xlabel=L"t"* " (ms)", label="Theoretical", alpha=0.8)
savefig("data/Q02_sim_vs_theoretical_rate_N$(N)_T$(new_T).png")

# bin length
plot(bin_t, new_avg_population_spikes_per_ms, title="Population firing rate", lw=2, xlabel=L"t"*" (ms)", label="bin length"*L"= %$bin_length"*" ms")
new_bin_length = delta_t
avg_population_spikes_per_ms = compute_rate(spikes, delta_t, new_bin_length, N, n)
new_avg_population_spikes_per_ms = compute_rate(new_spikes, delta_t, new_bin_length, new_N, n)
theoretical_t = range(0, stop = new_T, length = n+1)
bin_t = range(0, new_T, Int64(ceil(new_T/new_bin_length)))
theoretical_t = range(0, new_T, n+1)
plot!(bin_t, new_avg_population_spikes_per_ms, title="Population firing rate", lw=1, xlabel=L"t"*" (ms)", label="bin length"*L"= \Delta t")
plot!(theoretical_t, theoretical_r.(theoretical_t, alpha, beta, ro, Io, omega), lw=2, xlabel=L"t"* " (ms)", label="Theoretical", alpha=0.8)
savefig("data/Q02_simulated_vs_theoretical_adjusted_bin_length.pdf")

# initial rate
avg_population_spikes_per_ms = new_avg_population_spikes_per_ms
new_ro = 0.4
new_spikes = simulate_spikes(N, h_init, delta_t, new_T, R, Io, omega, alpha, beta, tau, new_ro)
new_avg_population_spikes_per_ms = compute_rate(spikes, delta_t, new_bin_length, N, n)
plot(bin_t, avg_population_spikes_per_ms, title="Population firing rate", lw=2, xlabel=L"t"*" (ms)", label=L"r_0= %$ro")
plot!(bin_t, new_avg_population_spikes_per_ms, title="Population firing rate", lw=2, xlabel=L"t"*" (ms)", label=L"r_0= %$new_ro")
plot!(theoretical_t, theoretical_r.(theoretical_t, alpha, beta, new_ro, Io, omega), lw=2, xlabel=L"t"* " (ms)", label="Theoretical", alpha=0.8)

# decrease tau
new_T = 50
tau1 = 1
tau2 = 0.1
new_bin_length = delta_t
n = Int64(ceil(new_T/delta_t))
spikes = simulate_spikes(new_N, h_init, delta_t, new_T, R, Io, omega, alpha, beta, tau, ro)
avg_population_spikes_per_ms = compute_rate(spikes, delta_t, new_bin_length, new_N, n)
spikes_tau1 = simulate_spikes(new_N, h_init, delta_t, new_T, R, Io, omega, alpha, beta, tau1, ro)
avg_population_spikes_per_ms_tau1 = compute_rate(spikes_tau1, delta_t, new_bin_length, new_N, n)
spikes_tau2 = simulate_spikes(new_N, h_init, delta_t, new_T, R, Io, omega, alpha, beta, tau2, ro)
avg_population_spikes_per_ms_tau2 = compute_rate(spikes_tau2, delta_t, new_bin_length, new_N, n)
theoretical_t = range(0, stop = new_T, length = n+1)
bin_t = range(0, new_T, Int64(new_T/new_bin_length))
theoretical_t = range(0, new_T, n+1)
plot(bin_t, avg_population_spikes_per_ms, title="Population firing rate", lw=1, xlabel=L"t"*" (ms)", label=L"\tau = %$(tau)")
plot!(bin_t, avg_population_spikes_per_ms_tau1, lw=1, xlabel=L"t"*" (ms)", label=L"\tau = %$(tau1)")
plot!(bin_t, avg_population_spikes_per_ms_tau2, lw=1, xlabel=L"t"*" (ms)", label=L"\tau = %$(tau2)")
plot!(theoretical_t, theoretical_r.(theoretical_t, alpha, beta, ro, Io, omega), lw=2, xlabel=L"t"* " (ms)", label="Theoretical", alpha=0.8)
savefig("data/Q02_rate_vs_varying_tau.pdf")
