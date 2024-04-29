module Neurons
    
    using LinearAlgebra
    
    export sigmoid, g, h_for_loop, h, r, I, theoretical_r, simulate_spikes, compute_rate
    
    function sigmoid(z::Real)
        return 1.0 / (1.0 + exp(-z))
    end

    function g(x::Real, alpha::Real=1.0, beta::Real=0.0)
        return sigmoid(2.0 * alpha * (x - beta))
    end

    function h_for_loop(h_init::Real, delta_t::Real, n::Real, R::Real, I_t::Array, tau::Real)
        # I_t is the input already computed for every timestep (because otherwise we would have to pass the arguments the function I needs to this function too)
        h_t = zeros(Float64, n+1)
        h_t[1] = h_init
        for i in range(1, length=n)
            h_t[i+1] = (1 - delta_t/tau) * h_t[i] + (delta_t/tau) * R * I_t[i]
        end
        return h_t
    end

    function h(h_init::Real, delta_t::Real, n::Real, R::Real, I_t::Array, tau::Real)
        lower_diag = (delta_t/tau - 1) * ones(n)
        diag = ones(n+1)
        A = Bidiagonal(diag, lower_diag, :L)
        b = vcat([h_init], (delta_t * R / tau) * I_t[1:n])
        return A\b
    end

    function r(h::Real, alpha::Real=1.0, beta::Real=0.0, ro::Real=1.0)
        return ro * g(h, alpha, beta)
    end

    # slow oscillating input for the poisson neurons
    function I(t::Real, Io::Real, omega::Real)
        return Io * sin(omega * t)
    end

    # Theoretical instantaneous firing rate
    function theoretical_r(t::Real, alpha::Real=1.0, beta::Real=0.0, ro::Real=1.0, Io::Real=1.0, omega::Real=10.0)
        return ro * g(I(t, Io, omega), alpha, beta)
    end

    # Simulate the spikes of N poisson neurons
    function simulate_spikes(
                N::Int=100, 
                h_init::Real=0, 
                delta_t::Real=0.1, 
                T::Real=1000, 
                R::Real=1, 
                Io::Real=2.0, 
                omega::Real=10, 
                alpha::Real=2, 
                beta::Real=0.5, 
                tau::Real=10, 
                ro::Real=1)
        n = Int64(round(T/delta_t))
        t = range(0, step = delta_t, length = n+1) 
        I_t = I.(t, Io, omega) # input to the neurons
        h_arr = h(h_init, delta_t, n, R, I_t, tau) # evolution of neuron potential
        r_arr = r.(h_arr, alpha, beta, ro) # instantaneous firing rate 
        P_arr = (r_arr * delta_t)[1:n] # spike probability (exclude last timestep)
        rand_nb = rand(Float64, (n, N))
        spikes = rand_nb .<= P_arr
        return spikes
    end

    # Compute the average population firing rate
    function compute_rate(
                spikes::BitMatrix, 
                delta_t::Real=0.1, 
                bin_length::Real=1,
                N::Int=100,
                n::Int=1000)
        bin_size = Int64(ceil(bin_length/delta_t)) # number of timesteps in a bin
        nb_spikes_per_ms = sum(reshape(spikes, bin_size, Int64(ceil(n/bin_size)), N), dims = 1) ./ bin_length
        avg_population_spikes_per_ms = sum(nb_spikes_per_ms, dims=3) / N
        avg_population_spikes_per_ms = avg_population_spikes_per_ms[1,:,:] # remove the second dimension
        return avg_population_spikes_per_ms
    end

end