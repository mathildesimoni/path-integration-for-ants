module SingleBumpAttractor


    using Neurons
    using Random, Distributions

    using BumpAttractorUtils
    
    # first collective variable
    function m_cos(N::Int, x_i::Array, S_i::Array)
        return 1/N * sum(cos.(x_i) .* S_i)
    end

    # second collective variable
    function m_sin(N::Int, x_i::Array, S_i::Array)
        return 1/N * sum(sin.(x_i) .* S_i)
    end

    # Input function
    function I(J::Real, x::Real, m_cos::Real, m_sin::Real, I_ext::Real, phi::Real = 0.0)
        return J * (cos(x - phi) * m_cos + sin(x - phi) * m_sin) + I_ext
    end

    function simulate_network(h_init::Array,
                    x_i::Array,
                    I_ext::Function,
                    phi::Real,
                    sp::SimulationParameters,
                    np::NetworkParameters)
        
        n = sp.n
        delta_t = sp.delta_t

        N = np.N
        R = np.R
        tau = np.tau
        alpha = np.alpha
        beta = np.beta
        ro = np.ro
        J = np.J

        # initialization
        h_t = zeros(Float64, (n+1, N)) # TODO: no need to store the full matrix if we don't return it
        h_t[1, :] = h_init

        S_i = zeros(Int64, (n+1, N))
        S_i[1, :] = zeros(Int64, N) # no spikes initially

        for i in range(1, length=n)
            t = i * delta_t # current time

            # calculate collective variables (only once per timestep)
            m_cos_t = m_cos(N, x_i, S_i[i, :])
            m_sin_t = m_sin(N, x_i, S_i[i, :])

            # calulate input to each neuron
            I_ext_t = I_ext.(x_i, t)
            I_t = I.(J, x_i, m_cos_t, m_sin_t, I_ext_t, phi) # TODO: check it broadcasts the x_i AND the I_ext_t

            # update potential at next timestep
            h_t[i+1, :] = (1 - delta_t/tau) * h_t[i, :] + (delta_t/tau) * R * I_t

            # update spike array by using the poisson neuron model
            r_arr = r.(h_t[i+1, :], alpha, beta, ro) # instantaneous firing rate 
            P_arr = (r_arr * delta_t) # spike probability
            rand_nb = rand(Float64, N)
            S_i[i+1, :] = Int.(rand_nb .<= P_arr)/delta_t
        end

        # replace the non zero values in S_i by 1
        S_i[S_i.>0] .= 1

        return S_i
    end

end