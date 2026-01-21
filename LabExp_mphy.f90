program odt_exp
    use microphysics, only :init_microphysics, add_prtcl_to_dom, move_prtcl_by_triplet_map,&
                            prtcl_growth, prtcl_settling, inject_stime,inject_etime, compute_microphy_stats,&
                            final_step_microphy,sat_mixing_ratio, tscale_nu 
    
    implicit none
    integer M, L, idum, iy, Length
    integer N, Lo, Lp, Lm, j, has_microphy
    integer NL(10000), ir(97), num_args, realz_int, fhdl_mp_log, debug
    character(len=32), dimension(:), allocatable :: args
    character(len=:), allocatable :: in_dir,cmd_str,exp_name, out_dir, mphy_nml, realization
    character(len=256) :: parent_dir
    double precision Nt, Np, Na, Nd, Co, Cm, uK, vK, wK, PE, mp_iter
    double precision pmax, p, prob, pa, Sc
    double precision Tdif, T_o, p_atm, s_limit, s_m, pv_o, pvdif, H
    double precision dt, time, td, tmax, te, ts, random
    double precision u(10000), ua(10000), ur(10000), eu(10000)
    double precision v(10000), va(10000), vr(10000), ev(10000)
    double precision w(10000), wa(10000), wr(10000), ew(10000)
    double precision T(10000), Ta(10000), Tr(10000), eT(10000)
    double precision pv(10000), pva(10000), pvr(10000),&
                       epv(10000)
    double precision s(10000), sa(10000), sr(10000),&
                       esat(10000)
    double precision PL(10000), a(14), pdf(5,1000),&
                       pdf_pv(5,1000), pdf_s(5,1000)
    double precision T_pv_covar(10000), index_arr(10000)
    !bulk mean
    double precision T_bmean(1000),T_bvar(1000),pv_bmean(1000),pv_bvar(1000),&
                     ss_bmean(1000),ss_bvar(1000), time_bevol(1000), deltat_bevol(1000),&
                     ss_bmin(1000), ss_bmax(1000)
    !domain mean
    double precision T_dmean(1000),T_dvar(1000),pv_dmean(1000),pv_dvar(1000),&
                     ss_dmean(1000),ss_dvar(1000), time_devol(1000), deltat_devol(1000),&
                     ss_dmin(1000), ss_dmax(1000)
    integer          bulk_sidx, bulk_eidx, bulk_len,wclock_tstart, wclock_tend, err_mphy    
    double precision cpu_tstart, cpu_tend, bulk_sfrac, bulk_efrac


    !clock 
    call cpu_time(cpu_tstart)
    call get_unix_time(wclock_tstart)

    !debug_mode
    debug     = 0
    !microphy iteration
    mp_iter   = 0
    !microphysics error
    err_mphy  = 0

    !Command line arguments
    num_args = command_argument_count()
    allocate(args(num_args))
    do j = 1, num_args
     call get_command_argument(j,args(j))
    end do

    !get in_dir, experiment name and realization value
    if (num_args   < 4)then
        write(*,*) "Incorrect number of arguments"
        write(*,*) " call as <program executable>  <in_dir> <namelist_file> <exp_name> <realization>"
        return
    else
        in_dir         = "input/"//trim(args(1)) // "/"
        mphy_nml       = trim(args(2))
        exp_name       = trim(args(3))
        realization    = trim(args(4))
    endif
    read(realization,*) realz_int
    write(*,'(8A)') " input:", in_dir," mphy nml:",mphy_nml,&
          " exp:",exp_name," realization:", realization

   call getcwd(parent_dir)
   write(*,*) "parent_dir: ", parent_dir

   !change to input directory
   cmd_str             = in_dir
   write(*,*) "changing to input dir: "
   call chdir(cmd_str)
   call system("pwd")

   10    format(5g16.6)
   20    format(g20.1,3g12.3,2i8)
    pmax = 0.1d0
    Nd   = 0.d0

    !read model params and initialization values
    call readpar(idum,N,Lo,Lp,Lm,dt,td,tmax,a,Sc,p_atm,&
                    T_o,Tdif,s_limit,pv_o,pvdif,H, has_microphy)
    idum   = idum + realz_int;
    write(*,"(A,I8)") "Random seed Idum:", idum

    !define bulk region
    bulk_sfrac = 0.2 
    bulk_efrac = 0.8
    bulk_sidx  = int(bulk_sfrac * N) +1
    bulk_eidx  = int(bulk_efrac * N)
    bulk_len   = bulk_eidx - bulk_sidx + 1
    write(*,*) ""
    write(*,*) "Bulk sidx:", bulk_sidx, "bulk eidx:", bulk_eidx, "bulk len:", bulk_len

    call init(N,u,v,w,T,pv,s)
    call zeroparam(NL,Nt,Np,Na,pdf,pdf_pv,pdf_s,time,te,ts,s_m)
    call zerovars(N,ua,ur,eu,va,vr,ev,wa,wr,ew,Ta,&
                    pva,Tr,pvr,eT,epv,sa,sr,esat,T_pv_covar)

    !change to parent directory
    call chdir(trim(parent_dir))
    write(*,*) "now at parent dir:"
    call system("pwd");
    
    !different output directory for different realization
    out_dir     = "output/"//exp_name//"/"//realization

    !create output dir
    cmd_str     = "mkdir -p "//out_dir
    call system(cmd_str)
    
    !cp config file to output
    call system("cp "//trim(in_dir)//"/LabExppar.dat "//trim(out_dir))
    !change to output directory
    
    

    !integrate microphysics
    !subroutine init_microphysics(idir,nml_fname, N, tmax, td, T_o, Tdif, press,fhdl_log,debug)
    if (has_microphy == 1)then            
      write(*,*) "************ ODT simulation with microphysics code"
      fhdl_mp_log = 1001
      open(fhdl_mp_log, file=out_dir//"/microphysics.log")
      call init_microphysics(in_dir,mphy_nml,out_dir,exp_name, H, N, tmax, td, T_o, Tdif, p_atm,fhdl_mp_log,realization,debug)
      flush(fhdl_mp_log)
    end if 

    !change to output directory
    !Since ODT output files are scatterred across various function.
    !It is better to change the current working directory to output directory
    call chdir(out_dir)
    write(*,*) "now at output dir:"
    call system("pwd");

    !compute eddy probability function parameters
    call LenProb(Lo,Lm,PL,a(14),Co,Cm)
    open(101, file="EddyDiagram.dat", status="unknown")
    open(102, file="Warnings.dat",    status="unknown")
    open(201, file="evol_dom_scalar.txt")
    write(201,'(A,A15,10A16)') '#', "time_ndim", "time_dim", "deltat_dim","T_mean","T_var",&
    "qv_mean", "qv_var", "ss_mean", "ss_var",     "ss_min" , "ss_max"  
    write(201,'(A,A15,10A16)') '#', "ndim", "secs", "secs","ndim","ndim",&
    "ndim", "ndim", "%", "%",  "%", "%"  

    open(202, file="evol_bulk_scalar.txt")
    write(202,'(A,A15,10A16)') '#', "time_ndim", "time_dim", "deltat_dim","T_mean","T_var",&
    "qv_mean", "qv_var", "ss_mean", "ss_var",     "ss_min" , "ss_max"  
    write(202,'(A,A15,10A16)') '#', "ndim", "secs", "secs","ndim","ndim",&
    "ndim", "ndim", "%", "%",  "%", "%"      
    write(6,*) "tmax,t", tmax, dt

    !Mani:Use random.dat to find random number generator's probability distribution
    !Mani: random produces a uniform distribution 
    !open(1001,file="out_dir//random.dat");
    !do j = 1, 10000
    !    write(1001,*) random(idum,iy,ir)
    !end do
    !close(1001);

    !test module scope
    !call test_scope()
    
    !begin ODT simulation
    do while (time .le. tmax)
        !write(*,'(A, 4E16.4)') "time, dt, td, time-te :", time, dt, td, time-te
        time = time + dt
        Nt   = Nt + 1.d0
       !write(*,*) "Nt:", Nt, " time:", time
       !track time and time step
        !long simuilations: 2000000   
        if (mod(Nt,2000000.0) == 0) then
          write(*,'(A24, F12.0, E12.4)') "Nt, time:",  Nt, time
          write(*,'(A24, 5F12.0)') "Nt, Np, Nd, Na, mp_iter:",  Nt, Np, Nd, Na, mp_iter
          write(*,'(A24, 2E12.4)') "tmax, dt:",        tmax, dt
          write(*,*)
        endif
        
        ! if (mod(Nt,3000000.0) == 0) then
        !   exit
        ! endif  

        !call diffusion  and compute stats
        if ((time-te) .ge. td) then
            !write(*,'(A36,F12.1)') "**** Diffusion Nd:" , Nd+1
            !write(*,'(A36,2E12.4)') "time, dt:", time, time-te
            !flush(6)

            call vis(N,u,v,w,time-te,a(9),a(10))
            call difT(N,T,pv,time-te,a(4),Sc)
            Nd = Nd + 1.d0
            !compute supersat before calling microphysics
            call compsup(N,T,pv,s,Tdif,T_o,p_atm)

            if (has_microphy == 1 .and. time >= inject_stime)then 
              !inject_stime and inject_etime are expressed as non-dim simulation time
              !call test_scope() 
              mp_iter  = mp_iter + 1.0d0
              call microphy_processes(N,T,pv,s, p_atm,time,te,mp_iter,fhdl_mp_log,Nd+1.0,0, err_mphy)              
              call compute_microphy_stats(T(1:N), pv(1:N), s(1:N), N, time, time-te, 0)
            end if

            !compute stats
            call statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,&
                            time-te)
            call statsdens(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr,eT,&
                            epv,esat,time-te)
            call statscovar(N,time-te,T,pv,T_pv_covar);
            call statsevol(201,0,N,T(1:N),pv(1:N),s(1:N),time,time-te,Nd,&
                      time_devol,deltat_devol,T_dmean,T_dvar,pv_dmean,pv_dvar,ss_dmean,ss_dvar,ss_dmin,ss_dmax)
            call statsevol(202,0,bulk_len,T(bulk_sidx:bulk_eidx),pv(bulk_sidx:bulk_eidx),s(bulk_sidx:bulk_eidx),time,time-te,Nd,& 
                      time_bevol,deltat_bevol,T_bmean,T_bvar, pv_bmean,pv_bvar,ss_bmean,ss_bvar,ss_bmin,ss_bmax)

            call comppdf(N,T,pdf,time-te)
            call comppdf(N,pv,pdf_pv,time-te)
            call comppdfsup(N,s,pdf_s,time-te,s_limit,T_o,Tdif,s_m)
            
            te = time

            if (err_mphy == -1) exit

        endif

        !get random eddy size, L
        L = 3*Length(PL,a(14),Co,Cm,idum,iy,ir)
        !get random position, M such that M+L doesn't go beyond domain boundary
        !hence multiply by (N-L)
        M = 1 + int(random(idum,iy,ir)*(N-L))
        !probability for eddy of given size and is it energetically possible (based buoyancy and kolmogorov scales)
        p = dt*prob(N,M,L,u,v,w,T,pv,a,uK,vK,wK,PE,T_o,Tdif,pv_o,pvdif)
        !if real i.e gt 0
        if (p .gt. 0.d0) then
            !if probability gt pmax then change the dt by same factor
            if (p .gt. pmax) then
                write(102,20) "Time step warning: ", time, dt, p, M, L
                dt = dt*pmax/p
                p = pmax
            endif
            pa = pa + p
            Np = Np + 1.d0
        endif

        !eddy is not called at every time step; so some iterations are wasted
        !w/o call to eddy or diffusion
        if (random(idum,iy,ir) .lt. p) then
            !write(*,'(A36,F12.1)') "**** Eddy mapping Na:" , Na+1
            !write(*,'(A36,2E12.4)') "time, dt: ", time, time-te
            !flush(6)
            !write(6,*) "random value < p", p
            !call energy(N,M,M+L,a(2),u,v,w,T)
            call eddy(N,M,L,u,v,w,T,pv,uK,vK,wK,PE,a(12))
            !call energy(N,M,M+L,a(2),u,v,w,T)

            !integrate microphysics into ODT : droplets moved by eddies
            if (has_microphy == 1 .and. time >= inject_stime )then 
              !move droplets by triplet map
              call reset_index_arr(index_arr, N)
              call triplet(N,M,L,index_arr)
              call move_prtcl_by_triplet_map(index_arr)
            end if

            !call diffusion
            call vis(N,u,v,w,time-te,a(9),a(10))
            call difT(N,T,pv,time-te,a(4),Sc)
            Nd = Nd + 1.d0
            !compute supersat before calling microphysics
            call compsup(N,T,pv,s,Tdif,T_o,p_atm)

            !compute stats
            if (has_microphy == 1 .and. time >= inject_stime)then 
              !inject_stime and inject_etime are expressed as non-dim simulation time
              !call test_scope() 
              mp_iter  = mp_iter + 1.0d0
              call microphy_processes(N,T,pv,s, p_atm,time,te,mp_iter,fhdl_mp_log,Nd+1.0,0, err_mphy)               
              call compute_microphy_stats(T(1:N), pv(1:N), s(1:N), N, time, time-te, 0)
            end if

            call statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time-te)
            call statsdens(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr,eT,&
                         epv,esat,time-te)
            call statscovar(N,time-te,T,pv,T_pv_covar);
            call statsevol(201,0,N,T(1:N),pv(1:N),s(1:N),time,time-te,Nd,&
                      time_devol,deltat_devol,T_dmean,T_dvar,pv_dmean,pv_dvar,ss_dmean,ss_dvar,ss_dmin,ss_dmax)
            call statsevol(202,0,bulk_len,T(bulk_sidx:bulk_eidx),pv(bulk_sidx:bulk_eidx),s(bulk_sidx:bulk_eidx),time,time-te,Nd,& 
                      time_bevol,deltat_bevol,T_bmean,T_bvar, pv_bmean,pv_bvar,ss_bmean,ss_bvar,ss_bmin,ss_bmax)
            call comppdf(N,T,pdf,time-te)
            call comppdf(N,pv,pdf_pv,time-te)
            call comppdfsup(N,s,pdf_s,time-te,s_limit,T_o,Tdif,s_m)
            !write time eddy position (midpoint) and eddy half length in non-dimensional form
            write(101,10) time, 1.d0*(M+(L/2))/(1.d0*N),&
                      (1.d0*L)/(2.d0*N)
            Na = Na + 1.d0
            
            NL(L/3) = NL(L/3) + 1
            te = time

            if (err_mphy == -1) exit
        endif

        if (Np > 1.d4) then
          call raisedt(Np,dt,pa)
        endif
        flush(6)

    enddo

    write (*,*) "outside loop"

    if (err_mphy /= -1) then

      !integrate microphysics into ODT: condensation and gravitational settling
      if (has_microphy == 1 .and. time >= inject_stime )then       
        !inject_stime and inject_etime are expressed as non-dim simulation time
        !call test_scope() 
        mp_iter  = mp_iter + 1.0d0
        call microphy_processes(N,T,pv,s, p_atm,time,te,mp_iter,fhdl_mp_log,Nd+1.0,1,err_mphy)
        call compute_microphy_stats(T(1:N), pv(1:N), s(1:N), N, time, time-te, 1)
      end if

      call vis(N,u,v,w,time-te,a(9),a(10))
      call difT(N,T,pv,time-te,a(4),Sc)
      Nd = Nd + 1.d0
      !compute supersat before calling microphysics
      call compsup(N,T,pv,s,Tdif,T_o,p_atm)

      call statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time-te)
      call statsdens(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr,eT,&
                      epv,esat,time-te)
      call statscovar(N,time-te,T,pv,T_pv_covar);
      call statsevol(201,1,N,T(1:N),pv(1:N),s(1:N),time,time-te,Nd,&
      time_devol,deltat_devol,T_dmean,T_dvar,pv_dmean,pv_dvar,ss_dmean,ss_dvar,ss_dmin,ss_dmax)
      call statsevol(202,1,bulk_len,T(bulk_sidx:bulk_eidx),pv(bulk_sidx:bulk_eidx),s(bulk_sidx:bulk_eidx),time,time-te,Nd,& 
      time_bevol,deltat_bevol,T_bmean,T_bvar, pv_bmean,pv_bvar,ss_bmean,ss_bvar,ss_bmin,ss_bmax)
      call comppdf(N,T,pdf,time-te)
      call comppdf(N,pv,pdf_pv,time-te)
      call comppdfsup(N,s,pdf_s,time-te,s_limit,T_o,Tdif,s_m)

      !write data to usual output dir
      call outputstats(N,Na,Nt,Lo,Lm,NL,PL,pdf,pdf_pv,pdf_s,&
                        time,dt,s_limit,s_m)
      call saveveldata(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time)
      call savedensdata(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr,&
                        eT,epv,esat,time)
      call savecovar(N,time,T_pv_covar,Ta,Tr,pva,pvr);
      if (has_microphy == 1 ) then 
        call final_step_microphy()
      endif  
      !call test_scope()
      write (*,*) "Save data completed"

    else
      !Write model state to dump directory
      !current dir changed to dump directory in microphysics module
      call saveveldata(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time)
      call savedensdata(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr,&
                        eT,epv,esat,time)
      write (*,'(42x,A)') "**** Error during microphysics call"
      write (*,'(42x,A)') "ODT model state dumped"

    endif
    
    close(101)
    close(102)
    close(201)
    close(202)

    !clock 
    write(*,*)""
    write(*,*)"************************ Execution time ****"
    call cpu_time(cpu_tend)
    call get_unix_time(wclock_tend)
    write(*,'(A,F8.2,A)') "       CPU time: ", (cpu_tend - cpu_tstart)/60.0d0 , " minutes" 
    write(*,'(A,F8.2,A)') "Wall clock time: ", (wclock_tend - wclock_tstart)/60.0d0 , " minutes" 

    stop
end

subroutine get_unix_time(tvar)
  !I am creating this function since I there is a variable defined as time
  !so I couldn't invoke this function
  integer   tvar

  tvar = time()

end subroutine

subroutine test_scope()
  use microphysics, only :    inject_stime, inject_etime,  m0_aero,n_prtcl_curr,n_prtcl_alltime
  implicit none

  write(*,*) "test module: inject_stime     : ", inject_stime 
  write(*,*) "test module: inject_etime     : ", inject_etime
  write(*,*) "test module: m0_aero          : ", m0_aero
  write(*,*) "test module: n_prtcl_curr     : ", n_prtcl_curr
  write(*,*) "test module: n_prtcl_alltime  : ", n_prtcl_alltime
  flush(6)

end subroutine  

subroutine reset_index_arr(index_arr, N)
  implicit none
  double precision , intent(inout) :: index_arr(10000)
  integer,intent(in)      ::  N

  !local variable
  integer                 :: j

  do j = 1, N
    index_arr(j) = j * 1.0d0
  end do

end subroutine

subroutine microphy_processes(N,T,pv,s, p_atm,time,te,mp_iter,fhdl_mp_log, Nd, flush_arr, err_stat)
  use microphysics, only : inject_stime, inject_etime,inject_n_tdiff, inject_onetime, is_injected_onetime,&
                           add_prtcl_to_dom, prtcl_growth, prtcl_settling,compute_microphy_stats,&
                           debug_gbl, tscale_nu , n_prtcl_curr
                           
  implicit none                         

  integer,intent(in)            :: N,fhdl_mp_log, flush_arr
  double precision,intent(in)   :: time, te, p_atm, Nd, mp_iter
  double precision,intent(inout):: T(N),pv(N), s(N)
  integer, intent(inout)        :: err_stat
  integer                       :: Nd_int,iter_debug, mp_iter_int

  !if time step is zero then return
  if ((time-te) < 1e-16) return

  iter_debug   = 1000 !test: change to 1000 for long simulations
  Nd_int       = nint(Nd)
  mp_iter_int  = nint(mp_iter)
  debug_gbl    = 0

  if( (mod(mp_iter_int,iter_debug) == 0) ) then
    debug_gbl  = 1
    write(fhdl_mp_log,'(A,A,F12.0,A)') repeat("=",36)," Calling microphysics: Iteration:", mp_iter, repeat("=",36)
    write(fhdl_mp_log,'(A,E16.4,A,E16.4)') "non-dim Time:", time, " delta_time:", (time-te)
    write(fhdl_mp_log,'(A,F16.6,A,F16.6)') "    dim Time:", time * tscale_nu , " delta_time:", (time-te)*tscale_nu
    write(fhdl_mp_log,'(A,I16)')           "          Nd:", Nd_int  
  end if

  if (inject_onetime == 1) then
    !write(*,'(A,E16.4,A)') "time_ndim:", time,  " === Adding droplets ONE TIME to the domain ..."
    if ( is_injected_onetime == 0) then
      call add_prtcl_to_dom(1.0d0, time)
    end if  

  else 
    if ( time >= inject_stime .and. time <= inject_etime .and. mod(Nd_int,inject_n_tdiff) == 0 ) then
    !write(*,'(A,E16.4,A)') "time_ndim:", time,  "=== Adding droplets to the domain ..."
      call add_prtcl_to_dom(time-te, time)
    end if

  end if
  
  !call microphysical process only if particles are present
  if (n_prtcl_curr > 0) then
    call prtcl_settling(T,pv,p_atm,time-te, N)
    call prtcl_growth(T,pv,s,p_atm,te,time, N, err_stat)
    if (err_stat == -1)then
      return
    endif  
  endif

  !compute stats independent of presence or absence of particles
  !call compute_microphy_stats(T(1:N), pv(1:N), s(1:N), N, time, time-te, flush_arr)

  if( (mod(mp_iter_int,iter_debug) == 0) ) then
    flush(fhdl_mp_log)
  end if  

end subroutine

function prob(N,M,L,u,v,w,T,pv,a,uK,vK,wK,PE,T_o,Tdif,pv_o,pvdif)
    implicit none
    integer N, M, L, j
    double precision prob, u(10000), v(10000), w(10000)
    double precision T(10000), pv(10000), a(14), uK, vK,&
                       wK, PE, KE
    double precision Tv(10000)
    double precision TK, p, x, psiK
    double precision Tdif, T_o, pv_o, pvdif
    double precision T_b,pv_b,T_t,pv_t,T_cell,pv_cell
    double precision Tv_b,Tv_t,Tvdif

    x = (1.d0*L)/(1.d0*N)
    T_b     = T_o  + Tdif/2.0d0
    T_t     = T_o  - Tdif/2.0d0
    pv_b    = pv_o + pvdif/2.0d0
    pv_t    = pv_o - pvdif/2.0d0
    Tv_b    = (T_b+273.15) * (1.0d0 + pv_b/0.622d0)/(1.0d0+pv_b)
    Tv_t    = (T_t+273.15) * (1.0d0 + pv_t/0.622d0)/(1.0d0+pv_t)
    Tvdif   = (Tv_b - Tv_t) 

    do j = 1, N
      !Tv(j)=( ( ((T(j)-0.5d0)*Tdif+T_o) * &
      !          (1.d0+((pv(j)-0.5d0)*pvdif+pv_o)/0.622d0)/ &
      !          (1.d0+((pv(j)-0.5d0)*pvdif+pv_o)) &
      !        ) - T_o + Tdif*0.5d0 &
      !     )/Tdif
      T_cell  = T_b  - T(j) * Tdif 
      pv_cell = pv_b - pv(j) * pvdif
      !compute virtual temperature
      Tv(j)   = (T_cell+ 273.15)*(1.0d0 + pv_cell/0.622d0)/(1.0d0+pv_cell)
      !convert to value between 0 and 1
      Tv(j)   = (Tv_b - Tv(j))/Tvdif

    enddo

    uK = psiK(N,M,L,u)
    vK = psiK(N,M,L,v)
    wK = psiK(N,M,L,w)
    TK = -psiK(N,M,L,Tv)
    PE = a(2)*TK*x
    KE = ((1.d0 - a(12))*wK*wK) + (a(12)*((uK*uK) + &
          (vK*vK))/2.d0)
    p  = ((KE + PE)*x*x) - a(11)

    if (p .gt. 0.d0) then

    prob = a(1)*(1.d0-x)*dsqrt(p)*dexp(3.d0*a(14)/(x*N))/(x*x)
    else
    prob = 0.d0
    endif
   return
end

function psiK(N,M,L,psi)
    implicit none
    integer N, M, L, j
    double precision psiK, sum, z
    double precision psi(10000)

    sum = psi(M)*L/2.d0
    do j=1, L-1
        z = L -( 2*j)
        sum = sum + (psi(j+M)*z)
    enddo
    sum = sum - (psi(M+L)*L/2.d0)
    psiK = 4.d0*sum/(9.d0*L*L)
  return
end

subroutine eddy(N,M,L,u,v,w,T,pv,uK,vK,wK,PE,c)
    implicit none
    integer N, M, L
    double precision u(10000), v(10000), w(10000)
    double precision T(10000), pv(10000), uK, vK, wK, PE, c
    double precision qu, qv, qw, cu, cv, cw

    qu = dsqrt((c*((vK*vK) + (wK*wK))/2.d0) + ((1.d0-c)*uK*uK))
    if (uK .gt. 0.d0) then
        cu = 6.75d0*(qu - uK)/(1.d0*L)
    else
        cu = -6.75d0*(qu + uK)/(1.d0*L)
    endif
    qv = dsqrt((c*((uK*uK) + (wK*wK))/2.d0) + ((1.d0-c)*vK*vK))
    if (vK .gt. 0.d0) then
        cv = 6.75d0*(qv - vK)/(1.d0*L)
    else
        cv = -6.75d0*(qv + vK)/(1.d0*L)
    endif
    qw = dsqrt((c*((uK*uK) + (vK*vK))/2.d0) + ((1.d0-c)*wK*wK)+&
                PE)
    if (wK .gt. 0.d0) then
        cw = 6.75d0*(qw - wK)/(1.d0*L)
    else
        cw = -6.75d0*(qw + wK)/(1.d0*L)
    endif
    call triplet(N,M,L,u)
    call triplet(N,M,L,v)
    call triplet(N,M,L,w)
    call triplet(N,M,L,T)
    call triplet(N,M,L,pv)
    call addK(N,M,L,u,cu)
    call addK(N,M,L,v,cv)
    call addK(N,M,L,w,cw)
    return
end

subroutine triplet(N,M,L,psi)
    implicit none
    integer N, M, L, Lo, j, k
    double precision psi(10000), x(10000)
    Lo = L/3
    do j = 1, Lo
        k = M + 3*(j-1)
        x(j) = psi(k)
    enddo
    do j=1, Lo
        k = M + L + 1 - (3*j)
        x(j+Lo) = psi(k)
    enddo
    do j = 1, Lo
        k = M + (3*j) - 1
        x(j+Lo+Lo) = psi(k)
    enddo
    do j=1, L
        k = M+j-1
        psi(k) = x(j)
    enddo
    return
end

subroutine addK(N,M,L,u,c)
    implicit none
    integer N, M, L, Lo, j, j1, j2, j3
    double precision c, y1, y2, y3
    double precision u(10000)
    Lo = L/3
    do j = 1, Lo
        y1 = -(2.d0*j)
        y2 = (4.d0*(j+Lo)) - (2.d0*L)
        y3 = (2.d0*L) - (2.d0*(j+Lo+Lo))
        j1 = M + j
        j2 = M + j + Lo
        j3 = M + j + Lo + Lo
        u(j1) = u(j1) + (c*y1)
        u(j2) = u(j2) + (c*y2)
        u(j3) = u(j3) + (c*y3)
    enddo
    return
end

subroutine difT(N,T,pv,dt,Pr,Sc2)
    implicit none
    integer N, j
    double precision T(10000), pv(10000), x(10000),&
                        xpv(10000)
    double precision l(10000), lpv(10000), d(10000),&
                        dpv(10000)
    double precision r(10000), rpv(10000)
    double precision dt, Pr, Sc2, De, Dv
    De = (dt*N*N)/(2.d0*Pr)
    Dv = (dt*N*N)/(2.d0*Sc2)
    l(1) = 0.d0
    d(1) = 1.d0 + (2.d0*De)
    r(1) = -De

    lpv(1) = 0.d0
    dpv(1) = 1.d0 + (2.d0*Dv)
    rpv(1) = -Dv

    do j = 2, N-1
        l(j) = -De
        d(j) = 1.d0 + (2.d0*De)
        r(j) = -De
        lpv(j) = -Dv
        dpv(j) = 1.d0 + (2.d0*Dv)
        rpv(j) = -Dv
    enddo
    d(N) = 1.d0
    l(N) = 0.d0
    r(N) = 0.d0
    dpv(N) = 1.d0
    lpv(N) = 0.d0
    rpv(N) = 0.d0

    x(1) = ((1.d0 - (2.d0*De))*T(1)) + (De*T(2))
    xpv(1) = ((1.d0 - (2.d0*Dv))*pv(1)) + (Dv*pv(2))

    do j = 2, N-1
        x(j) = ((1.d0-(2.d0*De))*T(j)) + (De*(T(j+1) + T(j-1)))
        xpv(j) = ((1.d0-(2.d0*Dv))*pv(j)) + (Dv*(pv(j+1)+&
                 pv(j-1)))
    enddo
    x(N) = T(N)
    xpv(N) = pv(N)
    call tridiagonal(N,l,d,r,x,T)
    call tridiagonal(N,lpv,dpv,rpv,xpv,pv)
    return
end

subroutine vis(N,u,v,w,dt,Uo,f)
    implicit none
    integer N, j
    double precision u(10000), v(10000), w(10000),&
                      dt, Uo, f, De
    double precision l(10000), d(10000), r(10000)
    double precision xu(10000), xv(10000), xw(10000)
    double precision duf(10000), dvf(10000), cft, sft
    De = (dt*N*N)/(2.d0)
    l(1) = 0.d0
    d(1) = 1.d0 + (2.d0*De)
    r(1) = -De
    do j = 2, N-1
        l(j) = -De
        d(j) = 1.d0 + (2.d0*De)
        r(j) = -De
    enddo
    d(N) = 1.d0
    l(N) = 0.d0
    r(N) = 0.d0
    cft = dcos(f*dt) - 1.d0
    sft = dsin(f*dt)
    do j = 1, N-1
        duf(j) = ((u(j) - Uo)*cft) + (v(j)*sft)
        dvf(j) = (v(j)*cft) - ((u(j) - Uo)*sft)
    enddo
    xu(1) = ((1.d0-(2.d0*De))*u(1)) + (De*u(2))
    xv(1) = ((1.d0-(2.d0*De))*v(1)) + (De*v(2))
    xw(1) = ((1.d0-(2.d0*De))*w(1)) + (De*w(2))
    do j = 2, N-1
        xu(j) = ((1.d0-(2.d0*De))*u(j)) + (De*(u(j+1)+u(j-1)))
        xv(j) = ((1.d0-(2.d0*De))*v(j)) + (De*(v(j+1)+v(j-1)))
        xw(j) = ((1.d0-(2.d0*De))*w(j)) + (De*(w(j+1)+w(j-1)))
    enddo
    xu(N) = u(N)
    xv(N) = v(N)
    xw(N) = w(N)
    call tridiagonal(N,l,d,r,xu,u)
    call tridiagonal(N,l,d,r,xv,v)
    call tridiagonal(N,l,d,r,xw,w)
    do j = 1, N-1
        u(j) = u(j) + duf(j)
        v(j) = v(j) + dvf(j)
    enddo
    return
end

subroutine tridiagonal(N,l,d,r,x,y)
    implicit none
    integer N, j
    double precision l(10000), d(10000), r(10000)
    double precision x(10000), y(10000), b, g(10000)
    b = d(1)
    y(1) = x(1)/b
    do j = 2, N
        g(j) = r(j-1)/b
        b = d(j) - (l(j)*g(j))
        if (b .eq. 0.d0) then 
          write(*,*) 'Tridiagonal:  failure'
          stop
        end if
        y(j) = (x(j) - (l(j)*y(j-1)))/b
    enddo
    do j = N-1, 1, -1
        y(j) = y(j) - (g(j+1)*y(j+1))
    enddo
    return
end

subroutine raisedt(Np,dt,p)
    implicit none
    double precision Np, dt, p, pmin
    pmin = 1.d-3
    p = p/Np
    if (p .lt. (pmin/2.d0)) then
        dt = dt*2.d0;
    else
        dt = dt*pmin/p
    endif
    p = 0.d0
    Np = 0.d0
    return
end

function Length(PL,xp,Co,Cm,idum,iy,ir)
    implicit none
    integer Length, idum, iy, n
    integer ir(97)
    double precision xp, Co, Cm, r, random, x
    double precision PL(10000)
    r = random(idum,iy,ir)
    x = -xp/dlog((Co*r)+(Cm*(1.d0-r)))
    n = int(x)-1
    if (r .gt. PL(n)) then
        10     n = n + 1
        if (r .gt. PL(n)) goto 10
    endif
    if (r .lt. PL(n-1)) then
        20     n = n - 1
        if (r .lt. PL(n-1)) goto 20
    endif
    Length = n
    return
end

subroutine LenProb(Lo,Lm,P,xp,Co,Cm)
    implicit none
    integer Lo, Lm, L
    double precision P(10000), xp, Co, Cm, C, z
    Co = dexp(-xp/(1.d0*Lo))
    Cm = dexp(-xp/(1.d0*Lm))
    C = 0.d0
    do L = Lo, Lm
        z = dexp(-xp/(1.d0*L))*(dexp(xp/(L*(L+1.d0)))-1.d0)
        C = C + z
    enddo
    C = 1.d0/C
    do L = 1, Lo-1
        P(L) = 0.d0
    enddo
    do L = Lo, Lm
        z = dexp(-xp/(1.d0*L))*(dexp(xp/(L*(L+1.d0)))-1.d0)
        P(L) = P(L-1) + (C*z)
    enddo
    do L = Lm+1, 10000
        P(L) = 0.d0
    enddo
    return
end

subroutine comppdf(N,T,pdf,dt)
    implicit none
    integer N, j, k(5), i, w
    double precision dt, T(10000), pdf(5,1000), X
    ! removed '-' in X
    X = 1.d3
    w = 5
    k(1) = N/2
    k(2) = N/4
    k(3) = N/8
    k(4) = 5*w
    k(5) = 2*w
    do j = -w, w
        i = int(X*(T(j+k(1))))
        pdf(1,i) = pdf(1,i) + dt
        i = int(X*(T(j+k(2))))
        pdf(2,i) = pdf(2,i) + dt
        i = int(X*(T(j+k(3))))
        pdf(3,i) = pdf(3,i) + dt
        i = int(X*(T(j+k(4))))
        pdf(4,i) = pdf(4,i) + dt
        i = int(X*(T(j+k(5))))
        pdf(5,i) = pdf(5,i) + dt

    enddo
    return
end

subroutine comppdfsup(N,s,pdf,dt,s_limit,T_o,Tdif,s_m)
    implicit none
    integer N, j, k(5), i, w
    double precision dt, s(10000), pdf(5,1000), X
    double precision s_limit, T_o, Tdif, es_b, es_t
    double precision es_m, s_m

    es_b = 6.112*dexp(17.67*(T_o+Tdif/2)/(T_o+Tdif/2+243.5))
    es_t = 6.112*dexp(17.67*(T_o-Tdif/2)/(T_o-Tdif/2+243.5))
    es_m = 6.112*dexp(17.67*T_o/(T_o+243.5))
    s_m = ((es_b+es_t)/(2*es_m)-1)*100

    X = 1.d3
    w = 5
    k(1) = N/2
    k(2) = N/4
    k(3) = N/8
    k(4) = 5*w
    k(5) = 2*w

    !grid cells +/- 5 levels above and below kth cell
    do j = -w, w
        !Mani: I don't know why the PDF is shifted by 500 position
        !I think this allows negative supersaturation values i.e sub-saturation
       i = int(X/2+X*(s(j+k(1))-s_m)/s_limit+1)
       pdf(1,i) = pdf(1,i) + dt

       i = int(X/2+X*(s(j+k(2))-s_m)/s_limit+1)
       pdf(2,i) = pdf(2,i) + dt

       i = int(X/2+X*(s(j+k(3))-s_m)/s_limit+1)
       pdf(3,i) = pdf(3,i) + dt

       i = int(X/2+X*(s(j+k(4))-s_m)/s_limit+1)
       pdf(4,i) = pdf(4,i) + dt

       i = int(X/2+X*(s(j+k(5))-s_m)/s_limit+1)
       pdf(5,i) = pdf(5,i) + dt

    enddo
    return
end

subroutine statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,dt)
    implicit none
    integer N, j
    double precision u(10000), ua(10000), ur(10000),&
                       eu(10000)
    double precision v(10000), va(10000), vr(10000),&
                       ev(10000)
    double precision w(10000), wa(10000), wr(10000),&
                       ew(10000)
    double precision dt, dudz, dvdz, dwdz
    do j = 1, N
        ua(j) = ua(j)+(u(j)*dt)
        va(j) = va(j)+(v(j)*dt)
        wa(j) = wa(j)+(w(j)*dt)
        ur(j) = ur(j)+(u(j)*u(j)*dt)
        vr(j) = vr(j)+(v(j)*v(j)*dt)
        wr(j) = wr(j)+(w(j)*w(j)*dt)
    enddo
    ! dz = 1/N
    dudz = (u(2) - u(1))*N
    dvdz = (v(2) - v(1))*N
    dwdz = (w(2) - w(1))*N
    eu(1) = eu(1)+(dudz*dudz*dt)
    ev(1) = ev(1)+(dvdz*dvdz*dt)
    ew(1) = ew(1)+(dwdz*dwdz*dt)
    do j = 2, N-1
        dudz = (u(j+1) - u(j-1))*N/2.d0
        dvdz = (v(j+1) - v(j-1))*N/2.d0
        dwdz = (w(j+1) - w(j-1))*N/2.d0
        eu(j) = eu(j)+(dudz*dudz*dt)
        ev(j) = ev(j)+(dvdz*dvdz*dt)
        ew(j) = ew(j)+(dwdz*dwdz*dt)
    enddo
    dudz = (u(N) - u(N-1))*N
    dvdz = (v(N) - v(N-1))*N
    dwdz = (w(N) - w(N-1))*N
    eu(N) = eu(N)+(dudz*dudz*dt)
    ev(N) = ev(N)+(dvdz*dvdz*dt)
    ew(N) = ew(N)+(dwdz*dwdz*dt)
    return
end

subroutine statsdens(N,T,pv,s,Ta,pva,sa,Tr,pvr,&
                       sr,eT,epv,esat,dt)
    implicit none
    integer N, j
    double precision T(10000), Ta(10000), Tr(10000),&
                       eT(10000)
    double precision pv(10000), pva(10000), pvr(10000),&
                       epv(10000)
    double precision s(10000), sa(10000), sr(10000),&
                       esat(10000)


    double precision dt, delT, delpv, dels
    do j = 1, N
        Ta(j) = Ta(j)+(T(j)*dt)
        Tr(j) = Tr(j)+(T(j)*T(j)*dt)
        pva(j) = pva(j)+(pv(j)*dt)
        pvr(j) = pvr(j)+(pv(j)*pv(j)*dt)
        sa(j) = sa(j)+(s(j)*dt)
        sr(j) = sr(j)+(s(j)*s(j)*dt)
    enddo
    delT = T(2) - T(1)
    delpv = pv(2) - pv(1)
    dels  = s(2) - s(1)
    eT(1) = eT(1)+(delT*delT*N*N*dt)
    epv(1) = epv(1)+(delpv*delpv*N*N*dt)
    esat(1) = esat(1)+(dels*dels*N*N*dt)
    do j = 2, N-1
        delT = T(j+1) - T(j-1)
        eT(j) = eT(j)+(delT*delT*N*N*dt/4.d0)
        delpv = pv(j+1) - pv(j-1)
        epv(j) = epv(j)+(delpv*delpv*N*N*dt/4.d0)
        dels = s(j+1) - s(j-1)
        esat(j) = esat(j)+(dels*dels*N*N*dt/4.d0)

    enddo
    delT = T(N) - T(N-1)
    eT(N) = eT(N)+(delT*delT*N*N*dt)
    delpv = pv(N) - pv(N-1)
    epv(N) = epv(N)+(delpv*delpv*N*N*dt)
    dels = s(N) - s(N-1)
    esat(N) = esat(N)+(dels*dels*N*N*dt)

    return
end

subroutine statscovar(N,dt,T,pv,T_pv_covar)
    implicit none
    integer N,j
    double precision T(10000), pv(10000),T_pv_covar(10000)
    double precision dt

    do j = 1, N
        T_pv_covar(j)  = T_pv_covar(j) + T(j)*pv(j)*dt
    enddo

    return
end

subroutine statsevol(fhdl,flush_arr,alen,T,pv,s,time,deltat,Nd,time_evol,deltat_evol,&
                     T_smean,T_svar,pv_smean,pv_svar,ss_smean,ss_svar,ss_min,ss_max)
    use microphysics, only: tscale_nu
    implicit none
    
    !evolution array len < 10,000
    integer,parameter :: elen = 1000
    integer           :: idx , j, ite_start, alen, fhdl, ite, flush_arr
    double precision  :: T(alen),pv(alen),s(alen), Nd, time,deltat
    double precision  :: T_smean(elen),T_svar(elen),pv_smean(elen),pv_svar(elen),ss_smean(elen),ss_svar(elen)
    double precision  :: time_evol(elen), deltat_evol(elen),ss_min(elen),ss_max(elen)
     
    !set
    ite            = int(Nd)
    idx            = mod(ite,elen)
    if (idx == 0) idx = elen

    !calc stats
    T_smean(idx)   = sum(T)/alen
    T_svar(idx)    = sum( (T - T_smean(idx))**2.d0 )/alen
    pv_smean(idx)  = sum(pv)/alen
    pv_svar(idx)   = sum( (pv - pv_smean(idx))**2.d0 )/alen

    do j = 1,alen
      !to avoid negative zeros affecting stats 
      if (s(j) > -1.0e-8 .and. s(j) < 1.0e-8  ) s(j) = 0.0d0
    end do

    ss_smean(idx)  = sum(s)/alen
    ss_svar(idx)   = sum( (s - ss_smean)**2.d0 )/alen
    ss_min(idx)    = minval(s)
    ss_max(idx)    = maxval(s)
    time_evol(idx) = time
    deltat_evol(idx)   = deltat

    !every 10,000the iteration dump the values into a file
    !if file flush is 1
    if (idx == elen .or. flush_arr == 1) then
      
      ite_start     = ite - idx 
      do j = 1, idx
        write(fhdl, '(11E16.8)') time_evol(j), time_evol(j)*tscale_nu, deltat_evol(j)*tscale_nu, T_smean(j), T_svar(j),&
                     pv_smean(j), pv_svar(j), ss_smean(j), ss_svar(j),ss_min(j),ss_max(j)
      end do
      flush(fhdl)
      
    end if  

    return
end

subroutine saveveldata(N,u,ua,ur,eu,v,va,vr,ev,w,&
                               wa,wr,ew,tf)
    implicit none
    integer N, j, k
    double precision u(10000), ua(10000), ur(10000),&
                    eu(10000)
    double precision v(10000), va(10000), vr(10000),&
                    ev(10000)
    double precision w(10000), wa(10000), wr(10000),&
                    ew(10000)
    double precision tf, z, small, dis, Re, etot

    10     format(5g16.6)
    20     format(g30.1,g15.3)
    small = 1.d-10
    dis = 0.d0
    do j = 1, N
        ua(j) = ua(j)/tf
        va(j) = va(j)/tf
        wa(j) = wa(j)/tf
        ur(j) = (ur(j)/tf)-(ua(j)*ua(j))
        vr(j) = (vr(j)/tf)-(va(j)*va(j))
        wr(j) = (wr(j)/tf)-(wa(j)*wa(j))
        eu(j) = eu(j)/tf
        ev(j) = ev(j)/tf
        ew(j) = ew(j)/tf
        dis = dis + eu(j) + ev(j) + ew(j)
    enddo
    dis = dis/(1.d0*N)
    Re = 0.d0
    k = 3*N/8
    do j = 1, N/4
        Re = Re + ur(j+k) + vr(j+k) + wr(j+k)
    enddo
    Re = dsqrt(4.d0*Re/(1.d0*N))
    write(6,20) "Core Reynolds Number, ReC: ", Re
    write(6,20) "Total KE dissipation:      ", dis
    write(6,20) "Wall stress, u:            ", ua(1)*N
    write(6,20) "                           ", (ua(N)-ua(N-1))*N
    write(6,20) "Wall stress, v:            ", va(1)*N
    write(6,20) "                           ", (va(N)-va(N-1))*N
    write(6,20) "Wall stress, w:            ", wa(1)*N
    write(6,20) "                           ", (wa(N)-wa(N-1))*N
    write(6,20)

    open(100, file="U.dat", status="unknown")
    open(200, file="V.dat", status="unknown")
    open(300, file="W.dat", status="unknown")
    open(400, file="eu.dat", status="unknown")

    do j = 1, N
        z = (1.d0*j)/(1.d0*N)
        etot = eu(j) + ev(j) + ew(j)
        if (abs(ua(j)) .lt. small) ua(j) = 0.d0
        if (abs(va(j)) .lt. small) va(j) = 0.d0
        if (abs(wa(j)) .lt. small) wa(j) = 0.d0
        if (ur(j) .lt. 0.d0) ur(j) = 0.d0
        if (vr(j) .lt. 0.d0) vr(j) = 0.d0
        if (wr(j) .lt. 0.d0) wr(j) = 0.d0
        write(100,10) z, u(j), ua(j), dsqrt(ur(j))
        write(200,10) z, v(j), va(j), dsqrt(vr(j))
        write(300,10) z, w(j), wa(j), dsqrt(wr(j))
        write(400,10) z, eu(j), ev(j), ew(j), etot
    enddo

    close(100)
    close(200)
    close(300)
    close(400)
    return
end

subroutine savedensdata(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr,&
                                eT,epv,esat,tf)
    implicit none
    integer N, j, k
    double precision T(10000), Ta(10000), Tr(10000),&
                    eT(10000)
    double precision pv(10000), pva(10000), pvr(10000),&
                    epv(10000)
    double precision s(10000), sa(10000), sr(10000),&
                    esat(10000)

    double precision tf, z, small, dis, dis_pv, x, x_pv

    10     format(5g16.8)
    20     format(g30.1,g15.4)

    small = 1.d-10
    dis = 0.d0
    dis_pv = 0.d0
    do j = 1, N
        Ta(j) = Ta(j)/tf
        Tr(j) = (Tr(j)/tf)-(Ta(j)*Ta(j))
        eT(j) = eT(j)/tf
        dis = dis + eT(j)
        pva(j) = pva(j)/tf
        pvr(j) = (pvr(j)/tf)-(pva(j)*pva(j))
        epv(j) = epv(j)/tf
        dis_pv = dis_pv + epv(j)
        sa(j) = sa(j)/tf
        sr(j) = (sr(j)/tf)-(sa(j)*sa(j))
        esat(j) = esat(j)/tf

    enddo

    dis = dis/(1.d0*N)
    dis_pv = dis_pv/(1.d0*N)
    x = 0.d0
    x_pv = 0.d0
    k = 3*(N/8)
    do j = 1, N/4
        x = x + Tr(j+k)
        x_pv = x_pv + pvr(j+k)
    enddo
    x = dsqrt(4.d0*x/(1.d0*N))
    x_pv = dsqrt(4.d0*x_pv/(1.d0*N))

    open(100, file="T.dat", status="unknown")
    open(200, file="pv.dat", status="unknown")
    open(300, file="Nu.dat", status="unknown")
    open(400, file="s.dat", status="unknown")

    write(300,20) "Core Temperature and pv Fluct.:   ", x, x_pv
    write(300,20) "Total Temperature and pv dissip.: ", dis, dis_pv
    write(300,20) "Nusselt Number:            ", Ta(1)*N
    write(300,20) "                           ", (Ta(N)-Ta(N-1))*N
    write(300,20) "Sherwood Number:            ", pva(1)*N
    write(300,20) "                         ", (pva(N)-pva(N-1))*N

    do j = 1, N
        z = (1.d0*DBLE(j))/(1.d0*N)
        if (abs(Ta(j)) .lt. small) Ta(j) = 0.d0
        if (Tr(j) .lt. 0.d0) Tr(j) = 0.d0
        if (abs(pva(j)) .lt. small) pva(j) = 0.d0
        if (pvr(j) .lt. 0.d0) pvr(j) = 0.d0
        write(100,10) z, T(j), Ta(j), dsqrt(Tr(j))
        write(200,10) z, pv(j), pva(j), dsqrt(pvr(j))
        if (abs(sa(j)) .lt. small) sa(j) = 0.d0
        if (sr(j) .lt. 0.d0) sr(j) = 0.d0
        write(400,10) z, s(j), sa(j), dsqrt(sr(j))

    enddo


    close(100)
    close(200)
    close(300)
    close(400)

    return
end

subroutine savecovar(N,tf,T_pv_covar,Ta,Tr,pva,pvr)
    integer N,j
    double precision tf,Z,numer,denom
    double precision T_pv_covar(10000)
    double precision Ta(10000), Tr(10000)
    double precision  pva(10000),pvr(10000)

    do j = 1, N
        numer        = (T_pv_covar(j)/tf)-(Ta(j)*pva(j))
        denom        =  dsqrt(Tr(j))*dsqrt(pvr(j))
        if (abs(denom) < 1e-9) then
          T_pv_covar(j) = 0
        else
          T_pv_covar(j)= numer/denom;
        end if  
    end do

    open(100, file="T_qv_covar.dat")
    do j = 1, N
        z = (1.d0*DBLE(j))/(1.d0*N);
        write(100,'(2g16.6)') z,T_pv_covar(j)
    enddo
    close(100)

    return
end

subroutine outputstats(N,Na,Nt,Lo,Lm,NL,PL,pdf,&
                       pdf_pv,pdf_s,tf,dt,s_limit,s_m)

    implicit none
    integer N, NL(10000), Lo, Lm, j
    double precision Na, Nt, tf, dt, x, xo, x_sup
    double precision PL(10000), pdf(5,1000), pdf_pv(5,1000)
    double precision pdf_s(5,1000), s_limit, s_m


    10    format(6g16.8)
    20    format(g30.1,g15.4)

    open(100, file="PDF_T.dat", status="unknown")
    open(200, file="PDF_pv.dat", status="unknown")
    open(300, file="PDF_s.dat", status="unknown")
    open(500, file="PL.dat", status="unknown")
    open(600, file="out_ODT.dat", status="unknown")

    do j = Lo, Lm
        if (NL(j) .gt. 0) then
            x = (1.d0*NL(j))/Na
            xo = PL(j) - PL(j-1)
            write(500,10) dlog((3.d0*j)/(1.d0*N)), dlog(x), dlog(xo)
        endif
    enddo


    write(600,20) "Final Time Step:           ", dt
    write(600,20) "Average Eddy Time Step:    ", tf/Nt
    write(600,20) "Time Between Eddies:       ", tf/Na
    write(600,20) "Total Number of Eddies:    ", Na
    write(600,20) "Eddie Acceptance Rate:     ", Na/Nt
    write(600,20)


    do j = 1 , 1000
        x = j/1.d3
        !removed '-' above
        pdf(1,j) = pdf(1,j)/tf
        pdf(2,j) = pdf(2,j)/tf
        pdf(3,j) = pdf(3,j)/tf
        pdf(4,j) = pdf(4,j)/tf
        pdf(5,j) = pdf(5,j)/tf
        write(100,10) x, pdf(1,j), pdf(2,j), pdf(3,j), pdf(4,j),&
                    pdf(5,j)
        pdf_pv(1,j) = pdf_pv(1,j)/tf
        pdf_pv(2,j) = pdf_pv(2,j)/tf
        pdf_pv(3,j) = pdf_pv(3,j)/tf
        pdf_pv(4,j) = pdf_pv(4,j)/tf
        pdf_pv(5,j) = pdf_pv(5,j)/tf
        write(200,10) x, pdf_pv(1,j), pdf_pv(2,j), pdf_pv(3,j),&
                    pdf_pv(4,j), pdf_pv(5,j)

        pdf_s(1,j) = pdf_s(1,j)/tf
        pdf_s(2,j) = pdf_s(2,j)/tf
        pdf_s(3,j) = pdf_s(3,j)/tf
        pdf_s(4,j) = pdf_s(4,j)/tf
        pdf_s(5,j) = pdf_s(5,j)/tf

        x_sup = s_m-s_limit/2.d0+DBLE(j)*s_limit/1.d3
        write(300,10) x_sup, pdf_s(1,j), pdf_s(2,j), pdf_s(3,j),&
                       pdf_s(4,j), pdf_s(5,j)

    enddo
    close(100)
    close(200)
    close(300)
    close(500)
    close(600)

    return
end

subroutine init(N,u,v,w,T,pv,s)
    implicit none
    integer N, j
    double precision u(10000), v(10000), w(10000)
    double precision T(10000), pv(10000), s(10000), z
    10    format(5g16.8)

    open(100, file="U.dat", status="old")
    open(200, file="V.dat", status="old")
    open(300, file="W.dat", status="old")
    open(400, file="T.dat", status="old")
    open(500, file="pv.dat", status="old")
    open(600, file="s.dat", status="old")



    do j = 1, N
        read(100,10) z, u(j)
        read(200,10) z, v(j)
        read(300,10) z, w(j)
        read(400,10) z, T(j)
        read(500,10) z, pv(j)
        read(600,10) z, s(j)
    enddo

    close(100)
    close(200)
    close(300)
    close(400)
    close(500)
    close(600)

    return
end

subroutine zeroparam(NL,Nt,Np,Na,pdf,pdf_pv,pdf_s,time,te,&
                        ts,s_m)

    implicit none
    integer NL(10000), j

    double precision Nt, Np, Na, pdf(5,1000), pdf_pv(5,1000),&
                         pdf_s(5,1000), time, te, ts, s_m
    Nt = 0.d0
    Np = 0.d0
    Na = 0.d0
    time = 0.d0
    te = 0.d0
    ts = 0.d0
    s_m = 0.d0
    do j=1, 1000
        pdf(1,j) = 0.d0
        pdf(2,j) = 0.d0
        pdf(3,j) = 0.d0
        pdf(4,j) = 0.d0
        pdf(5,j) = 0.d0
        pdf_pv(1,j) = 0.d0
        pdf_pv(2,j) = 0.d0
        pdf_pv(3,j) = 0.d0
        pdf_pv(4,j) = 0.d0
        pdf_pv(5,j) = 0.d0
        pdf_s(1,j) = 0.d0
        pdf_s(2,j) = 0.d0
        pdf_s(3,j) = 0.d0
        pdf_s(4,j) = 0.d0
        pdf_s(5,j) = 0.d0
    enddo
    do j=1, 10000
        NL(j) = 0
    enddo
    return
end

subroutine zerovars(N,ua,ur,eu,va,vr,ev,wa,wr,ew,&
                   Ta,pva,Tr,pvr,eT,epv,sa,sr,esat,T_pv_covar)
    implicit none
    integer N, j, n_k
    double precision ua(10000), ur(10000), eu(10000)
    double precision va(10000), vr(10000), ev(10000)
    double precision wa(10000), wr(10000), ew(10000)
    double precision Ta(10000), Tr(10000), eT(10000)
    double precision pva(10000), pvr(10000), epv(10000)
    double precision sa(10000), sr(10000), esat(10000)
    double precision T_pv_covar(10000)

    do j=1, N
        ua(j) = 0.d0
        ur(j) = 0.d0
        eu(j) = 0.d0
        va(j) = 0.d0
        vr(j) = 0.d0
        ev(j) = 0.d0
        wa(j) = 0.d0
        wr(j) = 0.d0
        ew(j) = 0.d0
        Ta(j) = 0.d0
        Tr(j) = 0.d0
        eT(j) = 0.d0
        pva(j) = 0.d0
        pvr(j) = 0.d0
        epv(j) = 0.d0
        sa(j) = 0.d0
        sr(j) = 0.d0
        esat(j) = 0.d0
        T_pv_covar(j) = 0.d0
    enddo
    return
end

subroutine readpar(idum,N,Lo,Lp,Lm,dt,td,tmax,a,Sc,p_atm,T_o,Tdif,&
                           s_limit,pv_o,pvdif,H, has_microphy)
    implicit none
    integer idum, N, Lo, Lp, Lm, has_microphy
    double precision dt, tmax, td, a(14), Sc
    double precision p_atm, T_o, Tdif, s_limit, pv_o, pvdif, H
    10    format(3i12)
    20    format(3d15.4)
    30    format(g30.1,g15.4)
    40    format(g30.1,i15)
    

    open(100,file="LabExppar.dat", status="old")

    read(100,10) N
    read(100,20) a(2), a(4), Sc
    read(100,20) a(9), a(10)
    read(100,20) dt, tmax
    read(100,20) a(11), a(12)
    read(100,10) Lo, Lp, Lm
    read(100,10) idum
    read(100,20) p_atm, T_o, Tdif,s_limit, pv_o, pvdif
    read(100,20) H
    read(100,10) has_microphy

    close(100)
    a(14) = 2.d0*Lp
    a(1) = dexp(-a(14)/(1.d0*Lm))-dexp(-a(14)/(1.d0*Lo))
    a(1) = a(1)*N/(3.d0*a(14))
    td = 1.d1/(1.d0*N*N)
    if (a(10)*td .gt. 0.1d0) td = 0.1d0/a(10)
    write(6,30) "Dimensionless Parameters:  "
    write(6,30) "   Buoyancy, Temperature:  ", a(2)
    write(6,30) "   Prandtl Number and Sc:  ", a(4), Sc
    write(6,30) "   Wind Speed:             ", a(9)
    write(6,30) "   Coriolis Parameter:     ", a(10)
    write(6,30)
    write(6,30) "Model Parameters:          "
    write(6,30) "   ZC2 =                   ", a(11)
    write(6,30) "   KE mixing ratio =       ", a(12)
    write(6,30) "   a1 =                    ", a(1)
    write(6,30) "   Smallest eddy:      ", (3.d0*Lo)/(1.d0*N)
    write(6,30) "   Most likely eddy:   ", (3.d0*Lp)/(1.d0*N)
    write(6,30) "   Largest eddy:       ", (3.d0*Lm)/(1.d0*N)
    write(6,30)
    write(6,40) "Grid Points:               ", N
    write(6,30) "Total Simulation Time:     ", tmax
    write(6,30) "Diffusive Time Step:       ", td
    write(6,30) "Initial Eddy Time Step:    ", dt
    return
end

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

subroutine energy(N,M1,M2,bT,u,v,w,T)
    implicit none
    integer N, M1,M2, j
    double precision bT, Eu, Ev, Ew, ET, EK
    double precision u(10000), v(10000), w(10000), T(10000)
    10    format(i12,3g16.6)
    Eu = 0.d0
    Ev = 0.d0
    Ew = 0.d0
    ET = 0.d0
    do j = M1, M2
        Eu = Eu + (u(j)*u(j))
        Ev = Ev + (v(j)*v(j))
        Ew = Ew + (w(j)*w(j))
        ET = ET + (T(j)*T(j))
    enddo
    Eu = Eu/(2.d0*N)
    Ev = Ev/(2.d0*N)
    Ew = Ew/(2.d0*N)
    EK = Eu + Ev + Ew
    ET = -ET*27.d0*bT/(8.d0*N*N)
    write(6, 10) M2-M1, EK+ET, EK, ET
    return

end

subroutine compsup(N,T,pv,s,Tdif,T_o,p_atm)
    use microphysics, only: sat_mixing_ratio
    implicit none
    
    integer N,j, xx
    double precision T(10000), pv(10000), s(10000)
    double precision Tdif, T_o, p_atm, T_t, T_b
    double precision pvT(2)
    double precision es_b, es_t, pv_t, pv_b, pv_sat
    
    10    format(g16.8)

    pvT(1)= 0.d0
    pvT(2)= 0.d0
    es_b= 0.d0
    es_t= 0.d0
    pv_t= 0.d0
    pv_b= 0.d0
    pv_sat= 0.d0

    T_b  = T_o + Tdif/2.0d0
    T_t  = T_o - Tdif/2.0d0
    es_b = 6.112d0*dexp(17.67d0*T_b/(T_b+243.5))*100.d0
    es_t = 6.112d0*dexp(17.67d0*T_t/(T_t+243.5))*100.d0
    pv_b = 0.622d0*es_b/(p_atm-es_b)
    pv_t = 0.622d0*es_t/(p_atm-es_t)

    !write(*,'(A, 1F12.6)') 'Sat mix ratio bottom:', pv_b
    !write(*,'(A, 1F12.6)') 'Sat mix ratio top   :', pv_t
    !write(*,*)""
    ! write(*,'(A)') "Supersaturation1:"      
    ! write(*,'(A)') "Temperature:"      
    ! write(*,'(10F12.6)')  T(1:10)
    ! write(*,'(A)') "water vapor:"
    ! write(*,'(10F12.6)')  pv(1:10)
    ! write(*,'(A)') "SS :"


    do j = 1, N
        !Mani: since T and pv have values between 0 and 1;
        !multiplying by Tdiff will scale them back; then add to Ttop, which is 0
        !T in Celsius, pv_sat in Pascals (hence the factor 100.d0)
        !pv_sat = 6.112d0*dexp(17.67d0*(T(j)*Tdif+T_o-Tdif/2.d0)/&
        !         (T(j)*Tdif+T_o-Tdif/2+243.5d0))*100.d0
        pv_sat = 6.112d0*dexp(17.67d0*(T_b - T(j)*Tdif)/&
                 ((T_b - T(j)*Tdif)+243.5d0))*100.d0
        pv_sat = 0.622d0*pv_sat/(p_atm-pv_sat)
        
        !any scalar can be obtained from scaling back for values between 0 and 1
        !super saturation in percentage 
        !s(j) = ((pv(j)*(pv_b-pv_t)+pv_t)/pv_sat-1.d0)*100.d0
         s(j) = ((pv_b - pv(j)*(pv_b-pv_t))/pv_sat-1.d0)*100.d0

        ! !debug 
        ! if (j > 2000 .and. j < 2002) then
        !   write(*,*) ""
        !   ! write(*,'(A,F12.6,A,F12.6)') "T_nondim :", T(j), " qv_nondim:" , pv(j)
        !   ! write(*,'(A,F12.6,A,F12.6)') "T_b :", T_b, " T_t:" , T_t
        !   ! write(*,'(A,F12.6,A,F12.6)') "T_b :", pv_b, " T_t:" , pv_t
        !   ! write(*,'(A,F12.6,A,F12.6)') "T_dim :", T_b - T(j)*Tdif, " qv_dim:" , pv_b - pv(j)*(pv_b - pv_t)
        !   !bot 1 and top 0
        !   !pv_sat = 6.112d0*dexp(17.67d0*(T(j)*Tdif+T_o-Tdif/2.d0)/&
        !   !         (T(j)*Tdif+T_o-Tdif/2+243.5d0))*100.d0          
        !   !pv_sat = 0.622d0*pv_sat/(p_atm-pv_sat)
        !   !s(j) = ((pv(j)*(pv_b-pv_t)+pv_t)/pv_sat-1.d0)*100.d0
        !   !bot 0 and top 1
        !   ! write(*,'(A,F12.6,A,F12.6,A,F12.6)') "SS1 :", s(j), " qv1:", (pv(j)*(pv_b-pv_t)+pv_t), ' qvs1 :', pv_sat
        !   !pv_sat = sat_mixing_ratio(T_o+Tdif/2 - T(j)*Tdif,p_atm)
        !   ! s(j) = ((pv_b - pv(j)*(pv_b-pv_t))/pv_sat-1.d0)*100.d0
        !   write(*,'(A,F12.6,A,F12.6,A,F12.6)') "SS2 :", s(j), " qv2:", (pv_b-pv(j)*(pv_b-pv_t)), ' qvs2 :', pv_sat
        ! end if

    enddo



  return

end
