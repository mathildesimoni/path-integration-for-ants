using Revise
a = push!(LOAD_PATH, pwd()*"/src")
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using RandomTrajectory

T = sp.T
delta_t = sp.delta_t
speed = Q3.speed
volatility = Q3.volatility
angles, pos = RandomTrajectory.random_trajectory(speed, T, delta_t, volatility)
plot(pos[:,1], pos[:,2], label=false)




