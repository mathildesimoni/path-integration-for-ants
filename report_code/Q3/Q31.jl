using Revise
a = push!(LOAD_PATH, pwd()*"/src")
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, BumpAttractor
using Random, Distributions
using RandomTrajectory

T = 1000
delta_t = 0.1 
speed = 0.001  # m/ms
volatility = 0.05
angles, pos = RandomTrajectory.random_trajectory(speed, T, delta_t, volatility)
plot(pos[:,1], pos[:,2], label=false)

