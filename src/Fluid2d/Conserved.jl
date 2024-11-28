module Conserved
# construction of each terms in conservative equation
using ..Eos: calc_p 

# conservative var from basic var
function calc_conserv(rho, u, v, e, s)
  return [rho*s, rho*u*s, rho*v*s, e*s]
end

# convective flux from basic var
function calc_flux_conv(rho,u,v,e,ixs,iys,eos)
  bigu = ixs*u + iys*v
  p = calc_p(rho,u,v,e,eos)
  return [rho*bigu, rho*u*bigu + ixs*p, rho*v*bigu + iys*p, (e+p)*bigu]
end

# basic var from conservative var
function calc_basic(vec_q,s)
  rho = vec_q[1]/s
  u = vec_q[2]/rho/s
  v = vec_q[3]/rho/s
  e = vec_q[4]/s
  return rho,u,v,e
end

# export calc_conserv, calc_flux_conv, calc_basic

end # module Conserved