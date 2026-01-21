       program Labinit
       implicit none
       integer N, Lo, Lp, Lm, j, has_microphy
       integer idum, iy, ir(97)
       double precision nu, kT, D_v, ks
       double precision g, alpha, Tdif, Tdif_e
       double precision H, lambda, Uo, f
       double precision C2, ZC2, chi
       double precision dt, tmax
       double precision bT, Pr, Sc
       double precision Up, fp
       double precision z, u, v, w, T, pv, s
       double precision Tv_b, Tv_t, pv_t, pv_b
       double precision es_t, es_b, T_o, p
       double precision pv_o, pvdif
       double precision small, random, s_limit
       character(len=:),allocatable ::      input_dir

 10    format(3i12)
 20    format(3d15.4)
 30    format(5d16.8)

  
c      multiplicative factor
       small = 2.d-10
c      Kinematic viscocity
       nu = 1.41d-5
c      Heat Diffusivity
       kT = 1.96d-5
c      Mass diffusivity
       D_v = 2.2705d-05
c      gravity
       g = 9.81d0
c      thermal expansion coeff
       alpha = 3.5d-3
c      Temperature difference
       Tdif = 12.d0
c      Reference temperature in celsius
       T_o = 15.d0
c
       s_limit = 30.d0
c      pressure
       p = 1.01d5
c      height
       H = 1.d0
c      velocity difference / shear
       Uo = 0.d0
c      coriolis factor
       f = 0.d0
c      C^2 model parameter
       C2 = 1.5d3
c      non-dim kolmogorov scale
       ZC2 = 1.d5
       chi = 0.d0
c      number of grid cells
       N = 6000
c      1/3rd of smallest eddy size
       Lo = 25
c      1/3rd of highly probable eddy size
       Lp = 100
c      1/3 of largest eddy size
       Lm = N/3
c      time step
       dt = 1.d0/(1.d0*N*N)
c      max simulation time
       tmax = 1.d-3
c      vapor pressure at top plate Tt = (T0-Tdiff) in celsius
       es_t = 6.112d0*dexp(17.67d0*(T_o-Tdif/2.d0)
     +        /(T_o-Tdif/2.d0+243.5d0))*100.d0
c      vapor pressure at bottom plate Tt = (T0+Tdiff) in celsius
       es_b = 6.112d0*dexp(17.67d0*(T_o+Tdif/2.d0)
     +        /(T_o+Tdif/2.d0+243.5d0))*100.d0
c      mixing ratios
       pv_t = 0.622d0*es_t/(p-es_t)
       pv_b = 0.622d0*es_b/(p-es_b)
c      Virtual temperature
c      Mani Convert temperature to Kelvin
       Tv_t = (T_o-Tdif/2.d0+273.15)*(1.d0+pv_t/0.622d0)/(1.d0+pv_t)
       Tv_b = (T_o+Tdif/2.d0+273.15)*(1.d0+pv_b/0.622d0)/(1.d0+pv_b)
c      difference b/w virtual temperature and mixing ratio
       Tdif_e = Tv_b-Tv_t
       pvdif = pv_b-pv_t
c      ??? reference mixing ratio; Shouldn't this be calculated from T_o
       pv_o  = (pv_t+pv_b)/2
       bT = (8.d0*g*alpha*Tdif_e*C2*H*H*H)/(27.d0*nu*nu)
c      Prandtl number
       Pr = nu/kT
c      Schmidt number
       Sc = nu/D_v
c      non-dimensional velocity difference
       Up = dsqrt(C2)*Uo*H/nu
c      non-dimensional coriolis factor
       fp = f*H*H/nu
c      has droplet microphysics       
       has_microphy = 0   
       
       input_dir  = "input/default/" 
       call system("mkdir -p "// input_dir)
       call chdir(input_dir)
       write(*,*) "Writing file LabExppar.dat with init params to"
       call system("pwd")


       open(100,file="LabExppar.dat", status="unknown")
       write(100,10) N
       write(100,20) bT, Pr, Sc
       write(100,20) Up, fp
       write(100,20) dt, tmax
       write(100,20) ZC2, chi
       write(100,10) Lo, Lp, Lm
       write(100,10) -1234
       write(100,20) p, T_o, Tdif, s_limit, pv_o, pvdif
       !used in microphysics for vol calc
       write(100,20) H
       write(100,10) has_microphy

       close(100)
       idum = -5555
       open(100, file="U.dat", status="unknown")
       open(200, file="V.dat", status="unknown")
       open(300, file="W.dat", status="unknown")
       open(400, file="T.dat", status="unknown")
       open(500, file="pv.dat", status="unknown")
       open(600, file="s.dat", status="unknown")

       do j=1, N
        z = (1.d0*j)/(1.d0*N)
        if (j .eq. N) small = 0.d0
c       assumes shear in x direction only
        u = small*(random(idum,iy,ir) - 0.5d0) + (z*Up)
        v = small*(random(idum,iy,ir) - 0.5d0)
        w = small*(random(idum,iy,ir) - 0.5d0)
        T = (T_o-Tdif/2.d0) + Tdif  *z*H
        pv = pv_b           + pvdif *z*H
        s = 0.d0
        write(100, 30) z, u
        write(200, 30) z, v
        write(300, 30) z, w
        write(400, 30) z, z
        write(500, 30) z, z
        write(600, 30) z, s
       enddo
       close(100)
       close(200)
       close(300)
       close(400)
       close(500)
       close(600)

       stop
       end

c produces uniformly distributed random values between 0 and 1
       function random(idum,iy,ir)
       implicit none
       integer idum, iy, m, ia, ic, j
       double precision random, rm
       integer ir(97)
       m = 714025
       ia = 1366
       ic = 150889
       rm = 1.4005112d-6
       if (idum .lt. 0) then
        idum = mod(ic-idum,m)
        do j=1,97
         idum = mod(ia*idum+ic,m)
         ir(j) = idum
        enddo
        idum = mod(ia*idum+ic,m)
        iy = idum
       endif
       j = 1 + (97*iy)/m
       if (j .gt. 97 .or. j .lt. 1) then
        j = mod(j,97)
        if (j .ge. 0) j=-j
        j = j + 1
       endif
       iy = ir(j)
       idum = mod(ia*idum+ic,m)
       ir(j) = idum
       random = iy*rm
       return
       end
