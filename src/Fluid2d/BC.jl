module BC
import ..BasicVarHD

# reflecting selected boundary condition
function reflect_bc(bc_type, basic, eos)
  if bc_type == "periodical_in_i"
    bc_periodical_in_i(basic)
  else
    throw(DomainError(bc_type, "boundary condition failed to be identified"))
  end
  return
end

# periodical in i-direction
function bc_periodical_in_i(basic::BasicVarHD)
  NI = basic.NI
  NJ = basic.NJ
  NB = basic.NB
  # periodical in i-direction
  # for i = 0
  basic.rho[1:NB,:] .= basic.rho[NI-2*NB+1:NI-NB,:]
  basic.u[1:NB,:] .= basic.u[NI-2*NB+1:NI-NB,:]
  basic.v[1:NB,:] .= basic.v[NI-2*NB+1:NI-NB,:]
  basic.e[1:NB,:] .= basic.e[NI-2*NB+1:NI-NB,:]
  # for i = NI
  basic.rho[NI-NB+1:NI,:] .= basic.rho[NB+1:2*NB,:]
  basic.u[NI-NB+1:NI,:] .= basic.u[NB+1:2*NB,:]
  basic.v[NI-NB+1:NI,:] .= basic.v[NB+1:2*NB,:]
  basic.e[NI-NB+1:NI,:] .= basic.e[NB+1:2*NB,:]
  # no updates for j-boundaries since Dirichlet
  return
end

# export reflect_bc

end # module BC