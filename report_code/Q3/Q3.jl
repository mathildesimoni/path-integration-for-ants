module Q3
    using BumpAttractorUtils, CoupledBumpAttractors
    using Distributions

    export np, sp

    # simulation parameters
    T = 1000 # simulation length in ms
    delta_t = 0.1 # timestep for the simulation in ms. MUST BE <= 1
    n = Int64(T/delta_t)
    speed = 0.001  # m/ms
    volatility = 0.1
    # volatility = 0.01
    theta_coupled = deg2rad(10) # 10 degres 
    Io = 1
    
    global sp = SimulationParameters()
    
    # Network parameters
    global np = NetworkParameters()
    
    function init_params()
        global sp = SimulationParameters(
            # Network parameters
            T = T,
            n = n,
            delta_t = delta_t
        )
            
        # Network parameters
        global np = NetworkParameters(
            N = 300, # number of poisson neurons
            R = 1,
            tau = 10.0,
            alpha = 2.0, 
            beta = 0.5,
            ro = 1,
            J = 5
        )
    end     
    
    init_params()

    t_to_idx(t) = Int(floor(t/sp.delta_t)) + 1
    
    I_ext_head(x::Real, t::Real, Io::Real, theta::Function) = Io*cos(x-theta(t))
    
    function get_I_head_cos(S_i::Array, t::Real)
        x_i = collect(range(start = 0, stop = 2*pi, length = np.N + 1)[1:np.N])
        I_head_val = (S_i/sp.delta_t) * cos.(x_i)
        return I_head_val[t_to_idx(t)]
    end
    
    function get_I_head_sin(S_i::Array, t::Real)
        x_i = collect(range(start = 0, stop = 2*pi, length = np.N + 1)[1:np.N])
        I_head_val = (S_i/sp.delta_t) * sin.(x_i)
        return I_head_val[t_to_idx(t)]
    end
    
    function integrate_position(J_head::Real, I_head::Function)
        h_init_L = rand(Uniform(0,1), np.N) 
        h_init_R = rand(Uniform(0,1), np.N)
        
        I_ext_L(x, t) = -J_head/np.N * I_head(t)
        I_ext_R(x, t) = J_head/np.N * I_head(t)
        
        x_i = collect(range(start = 0, stop = 2*pi, length = np.N + 1)[1:np.N])

        return CoupledBumpAttractors.simulate_network(h_init_L, h_init_R, I_ext_L, I_ext_R, x_i, theta_coupled, sp, np)

    end 

end