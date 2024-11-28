include("Fluid2d.jl")
using .Fluid2d

NI = 18
NJ = 18
NB = 4
gamma = 1.4
dir = "D:\\physics\\pressure_drag\\julia\\KH_grid400\\"
dir_o = dir * "data\\"
f_settings = dir * "settings.dat"
f_coord = dir * "coordinate.dat"

start = Dates.now()

# create fluid object
fluid = IdealGas(NI,NJ,NB,gamma)

# input coordinate
input(f_coord, fluid.coord)
# set metrices/Jacobian/dx for CFL
calc_metrices_dx(fluid.coord)
# input initial condition
input(dir_o * "b0000000.dat", fluid.basic)
t = 0.0
# output settings
print(fluid.coord.x)