module Q1
    
    using BumpAttractorUtils
    using Distributions, Random

    export np, sp

    T = 1000 # simulation length in ms
    delta_t = 0.1 # timestep for the simulation in ms. MUST BE <= 1

    # simulation parameters
    sp = SimulationParameters(
        T=T,
        delta_t=delta_t,
        n = Int64(T/delta_t),
        I_ext_bool = false # no external input
    )

    # Network parameters
    np = NetworkParameters(
        N = 300, # number of poisson neurons
        tau = 10.0, # characteristic time in ms
        R = 1, # resistance in Mohm
        alpha = 2.0, # parameter for the transfer function in 1/mV
        beta = 0.5, # parameter for the transfer function in mV
        ro = 1, # parameter for the mean firing rate function in 1/ms
        J = 5 
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