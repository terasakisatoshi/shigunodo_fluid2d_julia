module FluxScheme

import ..BasicVarHD
using ..Fnd: muscl_minmod, mp5
using ..Eos: calc_p, calc_e, calc_eigen
using ..Conserved: calc_conserv, calc_flux_conv

# reconstruction schemes for BasicVarHD

# reconstruct by basic variables with MUSCL-minmod
function reconst_by_basic_muscl(rho,u,v,e)
  rho_l, rho_r = muscl_minmod(rho[1],rho[2],rho[3],rho[4])
  u_l, u_r = muscl_minmod(u[1],u[2],u[3],u[4])
  v_l, v_r = muscl_minmod(v[1],v[2],v[3],v[4])
  e_l, e_r = muscl_minmod(e[1],e[2],e[3],e[4])
  return rho_l, u_l, v_l, e_l, rho_r, u_r, v_r, e_r
end

# reconstruct by basic variables with MP5
@inline function reconst_by_basic_mp5(rho,u,v,e)
  # use @inbounds macro to avoid boundary check
  rho1 = @inbounds rho[1]
  rho2 = @inbounds rho[2]
  rho3 = @inbounds rho[3]
  rho4 = @inbounds rho[4]
  rho5 = @inbounds rho[5]
  rho6 = @inbounds rho[6]

  u1 = @inbounds u[1]
  u2 = @inbounds u[2]
  u3 = @inbounds u[3]
  u4 = @inbounds u[4]
  u5 = @inbounds u[5]
  u6 = @inbounds u[6]

  v1 = @inbounds v[1]
  v2 = @inbounds v[2]
  v3 = @inbounds v[3]
  v4 = @inbounds v[4]
  v5 = @inbounds v[5]
  v6 = @inbounds v[6]

  e1 = @inbounds e[1]
  e2 = @inbounds e[2]
  e3 = @inbounds e[3]
  e4 = @inbounds e[4]
  e5 = @inbounds e[5]
  e6 = @inbounds e[6]

  rho_l, rho_r = mp5(rho1,rho2,rho3,rho4,rho5,rho6)
  u_l, u_r = mp5(u1,u2,u3,u4,u5,u6)
  v_l, v_r = mp5(v1,v2,v3,v4,v5,v6)
  e_l, e_r = mp5(e1,e2,e3,e4,e5,e6)

  #rho_l, rho_r = mp5(rho[1],rho[2],rho[3],rho[4],rho[5],rho[6])
  #u_l, u_r = mp5(u[1],u[2],u[3],u[4],u[5],u[6])
  #v_l, v_r = mp5(v[1],v[2],v[3],v[4],v[5],v[6])
  #e_l, e_r = mp5(e[1],e[2],e[3],e[4],e[5],e[6])
  return rho_l, u_l, v_l, e_l, rho_r, u_r, v_r, e_r
end

# calc numerical flux with selected scheme
function calc_num_flux(flux_scheme,rho_l,u_l,v_l,e_l,rho_r,u_r,v_r,e_r,ixs,iys,s,eos)
  if flux_scheme == "Roe_FDS"
    vec_fc = roe_fds(rho_l, u_l, v_l, e_l, rho_r, u_r, v_r, e_r, ixs, iys, s, eos)
  else
    throw(DomainError(flux_scheme, "flux scheme failed to be identified"))
  end
end

# Roe average
function roe_average(rho_l,u_l,v_l,e_l,rho_r,u_r,v_r,e_r,eos)
  p_l = calc_p(rho_l,u_l,v_l,e_l,eos)
  p_r = calc_p(rho_r,u_r,v_r,e_r,eos)
  h_l = (e_l + p_l)/rho_l
  h_r = (e_r + p_r)/rho_r
  sqr_l = sqrt(rho_l)
  sqr_r = sqrt(rho_r)

  rho_a = sqr_l * sqr_r
  u_a = (sqr_l*u_l + sqr_r*u_r)/(sqr_l + sqr_r)
  v_a = (sqr_l*v_l + sqr_r*v_r)/(sqr_l + sqr_r)
  h_a = (sqr_l*h_l + sqr_r*h_r)/(sqr_l + sqr_r)
  e_a = calc_e(rho_a,u_a,v_a,h_a,eos)
  return rho_a,u_a,v_a,e_a
end

# classical Roe-type FDS
@inline function roe_fds(rho_l, u_l, v_l, e_l, rho_r, u_r, v_r, e_r, ixs, iys, s, eos)
  ix = ixs/s
  iy = iys/s

  rho_a,u_a,v_a,e_a = roe_average(rho_l,u_l,v_l,e_l,rho_r,u_r,v_r,e_r,eos)

  vec_q_l = calc_conserv(rho_l,u_l,v_l,e_l,s)
  vec_q_r = calc_conserv(rho_r,u_r,v_r,e_r,s)
  vec_f_l = calc_flux_conv(rho_l,u_l,v_l,e_l,ixs,iys,eos)
  vec_f_r = calc_flux_conv(rho_r,u_r,v_r,e_r,ixs,iys,eos)
  dia_lam, mat_r, mat_rinv = calc_eigen(rho_a,u_a,v_a,e_a,ix,iy,eos)

  vec_fc = 0.5 .* (vec_f_r .+ vec_f_l  .- mat_r * abs.(dia_lam) * mat_rinv * (vec_q_r .- vec_q_l))

  return vec_fc
end

# export reconst_by_basic_muscl, reconst_by_basic_mp5, calc_num_flux

end # module FluxScheme