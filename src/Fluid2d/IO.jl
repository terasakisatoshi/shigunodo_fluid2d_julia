module IO

import ..BasicVarHD
import ..EulerEq
import ..IdealEoS
import ..GenStructCoord
using Printf

function input(f_name, basic::BasicVarHD)

  NI = basic.NI
  NJ = basic.NJ

  f = open(f_name,"r")
  body = readlines(f)
  close(f)

  k = 1
  for j = 1:NJ
    for i = 1:NI
      basic.rho[i,j] = parse(Float64, body[(k-1)*NI*NJ+(j-1)*NI+i])
    end
  end

  k = 2
  for j = 1:NJ
    for i = 1:NI
      basic.u[i,j] = parse(Float64, body[(k-1)*NI*NJ+(j-1)*NI+i])
    end
  end

  k = 3
  for j = 1:NJ
    for i = 1:NI
      basic.v[i,j] = parse(Float64, body[(k-1)*NI*NJ+(j-1)*NI+i])
    end
  end

  k = 4
  for j = 1:NJ
    for i = 1:NI
      basic.e[i,j] = parse(Float64, body[(k-1)*NI*NJ+(j-1)*NI+i])
    end
  end

  return
end # function input{BasicVarHD}

function input(f_name, coord::GenStructCoord)

  NI = coord.NI
  NJ = coord.NJ

  f = open(f_name,"r")
  body = readlines(f)
  close(f)

  # print(body[10], coord.x[1,1])

  k = 1
  for j = 1:NJ
    for i = 1:NI
      coord.x[i,j] = parse(Float64, body[(k-1)*NI*NJ+(j-1)*NI+i])
    end
  end

  k = 2
  for j = 1:NJ
    for i = 1:NI
      coord.y[i,j] = parse(Float64, body[(k-1)*NI*NJ+(j-1)*NI+i])
    end
  end

  return
end # function input{GenStructCoord}

function time_to_hms(time::Integer)
  hrs = div(time, 3600)
  secs = time - hrs * 3600
  mins = div(secs, 60)
  secs = secs - mins * 60
  return hrs, mins, secs
end

function output(dir_o, f_settings, tstep, t, iter, basic::BasicVarHD, cpu_time, rest_time)

  NI = basic.NI
  NJ = basic.NJ

  # output basic var to file
  f_name = dir_o * @sprintf("b%07d.dat", tstep)
  f = open(f_name,"w")

  for j = 1:NJ
    for i = 1:NI
      @printf(f,"%.18e\n", basic.rho[i,j])
    end
  end

  for j = 1:NJ
    for i = 1:NI
      @printf(f,"%.18e\n", basic.u[i,j])
    end
  end

  for j = 1:NJ
    for i = 1:NI
      @printf(f,"%.18e\n", basic.v[i,j])
    end
  end

  for j = 1:NJ
    for i = 1:NI
      @printf(f,"%.18e\n", basic.e[i,j])
    end
  end

  close(f)

  # output status to settings-file
  hh, mm, ss = time_to_hms(cpu_time)
  hr, mr, sr = time_to_hms(rest_time)
  line = @sprintf("elapsed: %3d h %02d m %02d s | tstep = %5d | t = %10.4f | iter = %7d | rest: %3d h %02d m %02d s\n", hh, mm, ss, tstep, t, iter, hr, mr, sr)
  f = open(f_settings,"a")
  print(f, line)
  close(f)
  print(line)

  return
end # output

function output(f_settings, t_max, n_out, basic, eos::IdealEoS, eq::EulerEq)
  f = open(f_settings,"w")
  print(f,"Euler equation with ideal EoS")
  @printf(f,"NI = %10d\n", basic.NI)
  @printf(f,"NJ = %10d\n", basic.NJ)
  @printf(f,"NB = %10d\n", basic.NB)
  @printf(f,"gamma = %10.4f\n", eos.gam)
  @printf(f,"Tmax = %10.4f\n", t_max)
  @printf(f,"Nout = %10d\n", n_out)
  print(f, "\n")
  close(f)

  return
end

# export input, output

end # module IO