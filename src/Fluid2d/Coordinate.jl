module Coordinate

import ..GenStructCoord

# calc metrices/Jacobian/dx for CFL
function calc_metrices_dx(coord::GenStructCoord)

  NI = coord.NI
  NJ = coord.NJ
  NB = coord.NB

  # metrices devided by Jacobian
  # evaluated at cell-points (i,j)
  # returning (NI-2*NB+2,NJ-2*NB+2) arrays
  # using 2nd order central difference
  coord.ixs .= 0.5 .* (coord.y[NB:NI-NB+1, NB+1:NJ-NB+2] .- coord.y[NB:NI-NB+1, NB-1:NJ-NB])
  coord.iys .= (-0.5) .* (coord.x[NB:NI-NB+1, NB+1:NJ-NB+2] .- coord.x[NB:NI-NB+1, NB-1:NJ-NB])
  coord.jxs .= (-0.5) .* (coord.y[NB+1:NI-NB+2, NB:NJ-NB+1] .- coord.y[NB-1:NI-NB, NB:NJ-NB+1])
  coord.jys .= 0.5 .* (coord.x[NB+1:NI-NB+2, NB:NJ-NB+1] .- coord.x[NB-1:NI-NB, NB:NJ-NB+1])

  # inverse of Jacobian
  # evaluated at cell-nodes (i+0.5,j+0.5)
  # returning (NI-1,NJ-1) array
  coord.s .= ((coord.x[2:NI, 2:NJ] .- coord.x[1:NI-1, 1:NJ-1]) .* (coord.y[1:NI-1, 2:NJ] .- coord.y[2:NI, 1:NJ-1]) .- (coord.y[2:NI, 2:NJ] .- coord.y[1:NI-1, 1:NJ-1]) .* (coord.x[1:NI-1, 2:NJ] .- coord.x[2:NI, 1:NJ-1])) .* 0.5

  # dx for CFL condition
  # minimum interval distance among adjacent 4 cells
  for j = 1:NJ-2*NB
    for i = 1:NI-2*NB
      dx1 = coord.x[NB+i+1, NB+j] - coord.x[NB+i, NB+j]
      dy1 = coord.y[NB+i+1, NB+j] - coord.y[NB+i, NB+j]
      dx = dx1 * dx1 + dy1 * dy1
      dx1 = coord.x[NB+i, NB+j+1] - coord.x[NB+i, NB+j]
      dy1 = coord.y[NB+i, NB+j+1] - coord.y[NB+i, NB+j]
      dd  = dx1 * dx1 + dy1 * dy1
      dx = min(dx, dd)
      dx1 = coord.x[NB+i-1, NB+j] - coord.x[NB+i, NB+j]
      dy1 = coord.y[NB+i-1, NB+j] - coord.y[NB+i, NB+j]
      dd  = dx1 * dx1 + dy1 * dy1
      dx = min(dx, dd)
      dx1 = coord.x[NB+i, NB+j-1] - coord.x[NB+i, NB+j]
      dy1 = coord.y[NB+i, NB+j-1] - coord.y[NB+i, NB+j]
      dd  = dx1 * dx1 + dy1 * dy1
      dx = min(dx, dd)
      coord.dx[i,j] = sqrt(dx)
    end
  end

  return
end # function metrices/Jacobian/dx for CFL

# export calc_metrices_dx

end # module Coordinate