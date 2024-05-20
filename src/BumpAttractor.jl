module BumpAttractor

    export simulate_network, locate_bump, locate_bump_avg, I_ext, simulate_networks, I_ext_2
    using Neurons
    using Random, Distributions
    
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

    # input function for the left/right bump attractor network
    # need to pass -1 * theta for the right bump attractor
    function I(J::Real, x::Real, m_cos_L::Real, m_sin_L::Real, m_cos_R::Real, m_sin_R::Real, I_ext::Real, theta::Real = 0.0)
        return J * (cos(x + theta) * (m_cos_L + m_cos_R) + sin(x + theta) * (m_sin_L + m_sin_R)) + I_ext
    end

    function I_ext(x::Real, t::Real)
        # using cases, k is 1 if 300 < t < 600, 2 if 600 < t < 900, 3 if 900 < t < 1200
        k = 300 <= t < 400 ? 1 :
            600 <= t < 700 ? 2 :
            return 0.0
        std = pi/8
        mean = k == 1 ? 2*pi/3 : 4*pi/3 
        return pdf(Normal(mean, std), x)
    end

    function I_ext_2(x::Real, t::Real, Io::Real)
        return 300 <= t < 600 ?  Io : 0.0
    end

    function I_ext_init(x::Real, mean)
        return pdf(Normal(mean, pi/8), x)
    end

    function simulate_network(
                    h_init::Array,
                    x_i::Array,
                    N::Int,
                    delta_t::Real,
                    n::Real,
                    R::Real,
                    tau::Real,
                    I_ext_bool::Bool, 
                    J::Real,
                    alpha::Real, 
                    beta::Real, 
                    ro::Real,
                    phi::Real = 0.0)
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
            if !I_ext_bool # no external input
                I_t = I.(J, x_i, m_cos_t, m_sin_t, 0.0, phi)
            else
                I_ext_t = I_ext.(x_i, t)
                I_t = I.(J, x_i, m_cos_t, m_sin_t, I_ext_t, phi) # TODO: check it broadcasts the x_i AND the I_ext_t
            end

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

    function locate_bump(S_i::AbstractArray, x_i::Array)
        return mod(weighted_circular_mean(S_i, x_i), 2*pi) # within [0, 2pi] range
    end

    # the only difference with the previous function is the signature (need it to be different for broadcasting purposes)
    function locate_bump_avg(S_i::Array, x_i::AbstractArray)
        return mod(weighted_circular_mean(S_i, x_i), 2*pi) # within [0, 2pi] range
    end

    function weighted_circular_mean(w, x)
        # implement the circular mean technique
        # https://en.wikipedia.org/wiki/Circular_mean#Using_complex_arithmetic
        cos_weighted_sum = sum(w .* cos.(x))
        sin_weighted_sum = sum(w .* sin.(x))
        return atan(sin_weighted_sum, cos_weighted_sum) 
    end

    function simulate_networks(
        h_init_L::Array,
        h_init_R::Array,
        x_i::Array,
        N::Int,
        delta_t::Real,
        n::Real,
        R::Real,
        tau::Real,
        I_ext_bool::Bool, 
        J::Real,
        alpha::Real, 
        beta::Real, 
        ro::Real,
        theta::Real = 0.0,
        I_ext_L::Base.Callable = I_ext,
        I_ext_R::Base.Callable = I_ext,
        rate_neurons::Bool = False)
        delay = 500
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
            elseif !I_ext_bool # no external input
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
            if !rate_neurons
                S_i_L[i+1, :] = Int.(rand_nb_L .<= P_arr_L)/delta_t
                S_i_R[i+1, :] = Int.(rand_nb_R .<= P_arr_R)/delta_t
            else
                S_i_L[i+1, :] = r_arr_L
                S_i_R[i+1, :] = r_arr_R
            end
        end

        if !rate_neurons
            # replace the non zero values in S_i by 1
            S_i_L[S_i_L.>0] .= 1
            S_i_R[S_i_R.>0] .= 1
        end

        return S_i_L[delay:n+1+delay,:], S_i_R[delay:n+1+delay,:]
    end

end