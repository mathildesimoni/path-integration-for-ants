module BumpAttractor

    export simulate_network, locate_bump, locate_bump_avg, I_ext
    using Neurons
    
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

    function I_ext(x::Real, t)
        if t < 300 || t >= 700
            return 0.0
        elseif t < 400
            u_k = (2*pi)/3
        elseif t < 600
            return 0.0
        else
            u_k = (4*pi)/3
        end
        sigma = pi/8
        return 1/(sqrt(2*pi)*sigma) * exp((-(x - u_k)^2 / (2*sigma^2)))
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

end