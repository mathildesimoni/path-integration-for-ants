module CoupledBumpAttractors

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

    # input function for the left/right bump attractor network
    # need to pass -1 * theta for the right bump attractor
    function I(J::Real, x::Real, m_cos_L::Real, m_sin_L::Real, m_cos_R::Real, m_sin_R::Real, I_ext::Real, theta::Real = 0.0)
        return J * (cos(x + theta) * (m_cos_L + m_cos_R) + sin(x + theta) * (m_sin_L + m_sin_R)) + I_ext
    end

    function I_ext_init(x::Real, mean)
        return pdf(Normal(mean, pi/8), x)
    end

    function simulate_network(
        h_init_L::Array,
        h_init_R::Array,
        I_ext_L::Base.Callable,
        I_ext_R::Base.Callable,
        x_i::Array,
        theta::Real,
        sp::SimulationParameters,
        np::NetworkParameters,
    )
        
        # unpack parameters
        n = sp.n
        delta_t = sp.delta_t
        delay = sp.delay
        
        N = np.N
        R = np.R
        tau = np.tau
        alpha = np.alpha
        beta = np.beta
        ro = np.ro
        J = np.J

        # initialization
        h_t_L = zeros(Float64, (n+1+delay, N)) # TODO: no need to store the full matrix if we don't return it
        h_t_R = zeros(Float64, (n+1+delay, N)) 
        h_t_L[1, :] = h_init_L
        h_t_R[1, :] = h_init_R

        S_i_L = zeros(Float64, (n+1+delay, N))
        S_i_R = zeros(Float64, (n+1+delay, N))
        S_i_L[1, :] = zeros(Float64, N) # no spikes initially
        S_i_R[1, :] = zeros(Float64, N) # no spikes initially

        for i in range(1, length=n + delay)
            t = (i-delay) * delta_t # current time

            # calculate collective variables (only once per timestep)
            m_cos_t_L = m_cos(N, x_i, S_i_L[i, :])
            m_sin_t_L = m_sin(N, x_i, S_i_L[i, :])
            m_cos_t_R = m_cos(N, x_i, S_i_R[i, :])
            m_sin_t_R = m_sin(N, x_i, S_i_R[i, :])

            # calulate input to each neuron
            
            if i < delay # initial input
                I_ext_t_L = I_ext_init.(x_i, pi)
                I_ext_t_R = I_ext_init.(x_i, pi)
                current_theta = 0.0
            elseif !sp.I_ext_bool # no external input
                I_ext_t_L = 0.0
                I_ext_t_R = 0.0
                current_theta = theta
            else
                I_ext_t_L = I_ext_L.(x_i, t)
                I_ext_t_R = I_ext_R.(x_i, t)

                current_theta = theta
            end
            I_t_L = I.(J, x_i, m_cos_t_L, m_sin_t_L, m_cos_t_R, m_sin_t_R, I_ext_t_L, current_theta) # TODO: check it broadcasts the x_i AND the I_ext_t
            I_t_R = I.(J, x_i, m_cos_t_L, m_sin_t_L, m_cos_t_R, m_sin_t_R, I_ext_t_R, -current_theta) 

            # update potential at next timestep
            h_t_L[i+1, :] = (1 - delta_t/tau) * h_t_L[i, :] + (delta_t/tau) * R * I_t_L
            h_t_R[i+1, :] = (1 - delta_t/tau) * h_t_R[i, :] + (delta_t/tau) * R * I_t_R

            # update spike array by using the poisson neuron model
            r_arr_L = r.(h_t_L[i+1, :], alpha, beta, ro) # instantaneous firing rate 
            r_arr_R = r.(h_t_R[i+1, :], alpha, beta, ro)

            P_arr_L = (r_arr_L * delta_t) # spike probability
            P_arr_R = (r_arr_R * delta_t) 
            rand_nb_L = rand(Float64, N)
            rand_nb_R = rand(Float64, N)
            if !sp.rate_neurons
                S_i_L[i+1, :] = Int.(rand_nb_L .<= P_arr_L)/delta_t
                S_i_R[i+1, :] = Int.(rand_nb_R .<= P_arr_R)/delta_t
            else
                S_i_L[i+1, :] = r_arr_L
                S_i_R[i+1, :] = r_arr_R
            end
        end

        if !sp.rate_neurons
            # replace the non zero values in S_i by 1
            S_i_L[S_i_L.>0] .= 1
            S_i_R[S_i_R.>0] .= 1
        end

        return S_i_L[delay+1:n+1+delay,:], S_i_R[delay+1:n+1+delay,:]
    end

end