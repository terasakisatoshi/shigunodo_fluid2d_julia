using Test

using Fluid2d
using Dates

@testset "Fluid2d" begin
    NI = 18
    NJ = 18
    NB = 4
    gamma = 1.4
    dir = joinpath(pkgdir(Fluid2d), "data")
    dir_o = dir
    f_settings = joinpath(dir, "settings.dat")
    f_coord = joinpath(dir, "coordinate.dat")

    start = Dates.now()

    # create fluid object
    fluid = IdealGas(NI,NJ,NB,gamma)

    # input coordinate
    input(f_coord, fluid.coord)
    # set metrices/Jacobian/dx for CFL
    calc_metrices_dx(fluid.coord)
    # input initial condition
    input(joinpath(dir_o, "b0000000.dat"), fluid.basic)
    t = 0.0
    # output settings
    @test size(fluid.coord.x) == (NI, NJ)
end

