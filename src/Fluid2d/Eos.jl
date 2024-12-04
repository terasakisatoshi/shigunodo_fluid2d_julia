module Eos

using LinearAlgebra
using StaticArrays: @SMatrix

import ..IdealEoS
# Equation of state

# EoS for ideal gas

@inline function calc_p(rho,u,v,e,eos::IdealEoS)
  return (eos.gam-1.0)*(e - 0.5*rho*(u*u+v*v))
end
@inline function calc_temp(rho,u,v,e,eos::IdealEoS)
  return e/rho - 0.5*(u*u+v*v)
end
@inline function calc_cs(rho,u,v,e,eos::IdealEoS)
  return sqrt(eos.gam*(eos.gam-1.0)*(e/rho - 0.5*(u*u+v*v)))
end
@inline function calc_e(rho,u,v,h,eos::IdealEoS)
  return rho*(h + (eos.gam-1.0)*0.5*(u*u+v*v))/eos.gam
end
@inline function calc_rho_e(p,temp,u,v,eos::IdealEoS)
  rho = p/temp/(eos.gam-1.0)
  e = rho*(temp + 0.5*(u*u+v*v))
  return rho, e
end
@inline function calc_e_wp(rho,u,v,p,eos::IdealEoS)
  return p/(eos.gam-1.0) + 0.5*rho*(u*u+v*v)
end
# calc eigen values/vectors of flux Jacobian
@inline function calc_eigen(rho,u,v,e,ix,iy,eos::IdealEoS)
  cs = calc_cs(rho,u,v,e,eos)
  p = calc_p(rho,u,v,e,eos)
  h = (e+p)/rho
  #sqr = sqrt(ix*ix + iy*iy)
  sqr = hypot(ix,iy)
  inv_sqr = inv(sqr)
  ixb = ix * inv_sqr
  iyb = iy * inv_sqr
  bigu = ix*u + iy*v
  bigub = bigu * inv_sqr
  inv_cs² = inv(cs*cs)
  b1 = 0.5*(u*u+v*v)*(eos.gam-1.0) * inv_cs²
  b2 = (eos.gam-1.0) * inv_cs²
  # diagonal matrix of eigen values
  # dia_lam = Diagonal([bigu-cs*sqr, bigu, bigu+cs*sqr, bigu])

  dia_lam = @SMatrix [
    bigu-cs*sqr 0.0 0.0 0.0;
    0.0 bigu 0.0 0.0;
    0.0 0.0 bigu+cs*sqr 0.0;
    0.0 0.0 0.0 bigu
  ]
  # right eigen matrix
  mat_r = @SMatrix [1.0          1.0              1.0             0.0;
          (u - ixb*cs)  u                (u + ixb*cs)    (-iyb);
          (v - iyb*cs)  v                (v + iyb*cs)     ixb;
          (h - cs*bigub)  0.5*(u*u+v*v)  (h + cs*bigub)  (-(iyb*u - ixb*v))]
  # inverse of righr eigen matrix
  mat_rinv = @SMatrix [0.5*(b1 + bigub/cs)  (-0.5*(ixb/cs + b2*u))  (-0.5*(iyb/cs + b2*v))  0.5*b2;
              (1.0 - b1)           b2*u                    b2*v                    (-b2);
              0.5*(b1 - bigub/cs)  0.5*(ixb/cs - b2*u)     0.5*(iyb/cs - b2*v)     0.5*b2;
              (iyb*u - ixb*v)      (-iyb)                  ixb                     0.0]
  return dia_lam, mat_r, mat_rinv
end

# export calc_p, calc_cs, calc_e, calc_temp, calc_e_wp, calc_rho_e, calc_eigen

end # module Eos