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
# plot(pos[:,1], pos[:,2], label=false)
p = RandomTrajectory.plot_trajectory(pos[:,1], pos[:,2], angles)
display(p)

savefig("./data/Q31.pdf")

# using GMT

# p = plot()
# plot(rand(5,2))
# psimage!("report_code/Q3/icons/ant.png", D="g0.2/0.5+jCM+w2c", fmt=:png, show=1)

# # Start a new GMT session
# GMT.gmt("begin")

# # Create your plot
# GMT.plot(pos[:,1], pos[:,2])

# first_x = pos[:,1][1] * 4.8
# first_y = pos[:,2][1] * 4.8
# last_x = pos[:,1][end] * 4.8
# last_y = pos[:,2][end] * 4.8
# angle = pi
# # Overlay the image using psimage
# # GMT.psimage("report_code/Q3/icons/ant.png", D="x1/1+w2c", O=true, R = pi/2)
# GMT.psimage("report_code/Q3/icons/house.png", D="x$first_x/$first_y+w2c", O=true, R=angle)
# GMT.psimage("report_code/Q3/icons/ant.png", D="x$last_x/$last_y+w2c", O=true, R=angle)

# # Finalize the plot
# GMT.gmt("end show")



