module Q2
    using BumpAttractorUtils
    using Distributions

    export np, sp

    # simulation parameters
    T = 1000 # simulation length in ms
    delta_t = 0.1 # timestep for the simulation in ms. MUST BE <= 1
    n = Int64(T/delta_t)

    sp = SimulationParameters(
        # Network parameters
        T = T,
        n = n,
        delta_t = delta_t,
        I_ext_bool = true,
    )

    # Network parameters
    np = NetworkParameters(
        N = 300, # number of poisson neurons
        R = 1,
        tau = 10.0,
        alpha = 2.0, 
        beta = 0.5,
        ro = 1,
        J = 3
    )

    function I_ext(x::Real, t::Real)
        # using cases, k is 1 if 300 < t < 600, 2 if 600 < t < 900, 3 if 900 < t < 1200
        k = 300 <= t < 400 ? 1 :
            600 <= t < 700 ? 2 :
            return 0.0
        std = pi/8
        mean = k == 1 ? 2*pi/3 : 4*pi/3 
        return pdf(Normal(mean, std), x)
    end
    
end