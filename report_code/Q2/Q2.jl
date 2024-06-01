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
        delta_t = delta_t
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

    function I_ext(x::Real, t::Real, Io::Real)
        return 300 <= t < 600 ?  Io : 0.0
    end
end