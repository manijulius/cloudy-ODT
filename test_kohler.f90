!Mani Rajagopal
!Nov 22nd 2023 
!subroutines and functions are from the test program and microphysics implemented earlier
!Hence somelines/parts of program might seem unneccessary

!This program helps find the Kohler curve (x, and y values) - equlibrium droplet radius for various supersaturation 
!Kohler curve is obtained using two different methods for comparison
!1)Using the equations from Roger and Yau textbook (chap 6) Eqs. 6.7 and 6.8  
!2)Using droplet growth model used in EMPM and the new version ODT w/ microphysics
!The droplet growth model makes a droplet grow/shrink based on the environmental conditions, 
!which also leads to changing the the environment, particularly temperature and humidity
!Therefore, I wrote a driver to iteratively reset the environment allowing for the droplet to grow 
!to the equilibrium for constant environment condition

!compile as 
!gfortran -o test_kohler.exe  const.f90 array.f90 ew.f90   fcnkb0.f90 fcnkb.f90  rkqs.f90 rkck.f90 odeint.f90 test_kohler.f90
!ensure that you have above program files

!execute 
!./test_kohler.exe

!output
!test_kohler.log -> log file
!kohler_RY2.txt  -> col1 is saturation ratio in percentage, col2 is equilibrium radius in meters
!kohlerdgm2.txt  -> col1 is saturation ratio in percentage, col2 is equilibrium radius in meters

program test_kohler
  implicit none
  double precision :: mass_aero, Temp_env, rho_aero, inv_vol_gcell, molar_mass_aero
  integer          :: fhdl_log, aero_type, nions_aero

  !write(7,"(/,A)") "======================= Unit test kohler curve ============"
  fhdl_log            = 7
  open(7,file="test_kohler.log")

  mass_aero           = 2.212E-18
  Temp_env            = 288.0d0

  aero_type           = 1 !1 -> Nacl, 2->Ammonium sulphate, 3-> ammonium bisulphate
  rho_aero            = 2.163d3
  nions_aero          = 2
  molar_mass_aero     = 58.4428e-3

  inv_vol_gcell       = 1.0d0/6000*1e-6 !1 cm3 is divided into 6000 grid cells
  

  call kohler_curve_from_dgm(mass_aero,Temp_env)
  call kohler_curve_from_RY(mass_aero,Temp_env)

  close(7);
  stop
  contains 

  function sat_mixing_ratio(T_in_C, press)
    implicit none
      double precision :: sat_mixing_ratio , T_in_C, es, press
      double precision, EXTERNAL :: ew
  
      es = 6.112d0 * dexp(17.67d0*(T_in_C)/(T_in_C+243.5)) * 100.d0
      !es  = ew(T_in_C + 273.15) * 100
      sat_mixing_ratio = 0.622d0 * es/(press - es)
  
      return
  end
 
  subroutine set_odeint_dgm_params(gcell_sidx, aerosol_type, m0_aerosol, invvol_gcell)
    !Subroutines odeint, fcnkb, fcnkb0 for particle growth is used as is from the EMPM model
    !these subroutines uses parameters from the array.f90 module
    ! to avoid changing code. I am setting those variables from here  
    use array, only :  ndmax, solute_mass, dmaxa, grid_scale, err_dgm
    implicit none

    integer, intent(in)    :: gcell_sidx, aerosol_type
    double precision, intent(in) ::  m0_aerosol, invvol_gcell

    !starting index of array
    ndmax                 = gcell_sidx
    !aerosol type
    dmaxa                 = aerosol_type 
    !solute mass in kg
    solute_mass           = m0_aerosol
    !grid scale is actually inverse of grid volume
    grid_scale            = invvol_gcell
    
    !reset the error flag
    err_dgm = 0

  end subroutine

  subroutine indiv_prtcl_growth(T_env, qv_env, p_env, ss_env, ql_prtcl, t1, t2, r_prtcl, m0_aerosol, aerosol_type,&
      invvol_gcell, has_cs_effect, is_iterative, debug)
    !This subroutine is a driver for the droplet growth model using in EMPM and the new version of ODT  
    use const, only                   : pi43, rho_w
    use array, only                   : err_dgm
    implicit none

    double precision, intent(inout)   :: T_env, qv_env, p_env, r_prtcl, ss_env, ql_prtcl
    double precision, intent(in)      :: t1, t2, m0_aerosol, invvol_gcell 
    integer, intent(in)               :: has_cs_effect, debug, aerosol_type
    !external subroutines
    EXTERNAL fcnkb,rkqs,fcnkb0, odeint

    !local variables 
    integer, parameter  :: n_yvar = 8
    double precision    :: y_arr(n_yvar), y_before(n_yvar), dt, eps_err, dt_min, n_ok, n_bad, qvs_env
    double precision    :: t1_intg, t2_intg, tstep_intg
    integer             :: debug_lcl, iter_count, is_iterative

    debug_lcl         = debug
    tstep_intg        = 1e-2

    !integration time - conv from non-dim to dim
    dt                = (t2 - t1)
    !acceptable margin of error 
    eps_err           = 1.0e-6

    !minimum time step between t1 and t2 in RK method
    !RK 5th order integration between t1 and t2 will use timesteps smaller than  (t2-t1) to get accurate integral
    dt_min            = 0.0d0

    y_arr(1)          = r_prtcl
    y_arr(2)          = qv_env
    y_arr(3)          = T_env
    y_arr(4)          = ss_env
    y_arr(5)          = p_env
    y_arr(6)          = 0 !vertical velocity
    y_arr(7)          = 0 !height
    y_arr(8)          = ql_prtcl !(pi43 * rho_w * r_prtcl**3.0  - m_aero )* inv_vol_dom
    y_before          = y_arr

    !I use EMPM model's droplet growth model encapsulated in subprogram called odeint
    !some of the subprogram parameters are set in a file called array.f90 particularly aero mass
    !if each particle has different aerosol mass then reset the value before calling odeint
    call set_odeint_dgm_params(1, aerosol_type, m0_aerosol, invvol_gcell)

    !error flag for inegration function
    err_dgm   = 0

    !iteratively integrate a droplet  
    if ( is_iterative == 1) then
      t1_intg       = t1
      t2_intg       = t1_intg + tstep_intg
      iter_count    = 0
      do while(t2_intg < t2)
        iter_count  = iter_count + 1
        if (has_cs_effect == 1) then
          call odeint(y_arr,n_yvar,t1_intg,t2_intg,eps_err,dt,dt_min,n_ok,n_bad,fcnkb,rkqs)
        else
          call odeint(y_arr,n_yvar,t1_intg,t2_intg,eps_err,dt,dt_min,n_ok,n_bad,fcnkb0,rkqs)
        end if
          t1_intg     = t2_intg
          t2_intg     = t2_intg + tstep_intg
        if (isnan(y_arr(1)) .or. err_dgm == -1) exit
      enddo

      if ( t2 > t1_intg .and. t2 < t2_intg .and. err_dgm == 0 .and. isnan(y_arr(1)) .eqv. .false.)then
        iter_count  = iter_count + 1
        if (has_cs_effect == 1) then
          call odeint(y_arr,n_yvar,t1_intg,t2,eps_err,dt,dt_min,n_ok,n_bad,fcnkb,rkqs)
        else
          call odeint(y_arr,n_yvar,t1_intg,t2,eps_err,dt,dt_min,n_ok,n_bad,fcnkb0,rkqs)
        end if
      endif  

    else
      iter_count  = 1
      if (has_cs_effect == 1) then
        call odeint(y_arr,n_yvar,t1,t2,eps_err,dt,dt_min,n_ok,n_bad,fcnkb,rkqs)
      else
        call odeint(y_arr,n_yvar,t1,t2,eps_err,dt,dt_min,n_ok,n_bad,fcnkb0,rkqs)
      end if
    endif    


    !control debug_lcl t through parent subroutine prtcl_growth()
    if (debug_lcl == 1 .or. isnan(y_arr(1)) .or. err_dgm == -1) then
      write(fhdl_log,*) ""
      if (isnan(y_arr(1))) write(fhdl_log,'(A8,A48,A)')"pgrow:"," Error:", "droplet growth - Nan value"
      if (err_dgm == -1)   write(fhdl_log,'(A8,A48,A)')"pgrow:"," Error:"," droplet growth - error flag"
      write(fhdl_log,'(A8,A48,F16.6)') "pgrow:", "time intvl:", t2-t1
      write(fhdl_log,'(A8,A48,I16)') "pgrow:", " droplet grew iteratively: ", iter_count
      write(fhdl_log,'(A8,A24,8E16.6)') "pgrow:","Before ode integral:", y_before
      write(fhdl_log,'(A8,A24,8E16.6)') "pgrow:","After ode integral :", y_arr
      flush(fhdl_log)
    end if    

    r_prtcl    = y_arr(1)
    qv_env     = y_arr(2)
    T_env      = y_arr(3)
    !ss and ql from integration is slightly different from calc ss  abd ql
    !y_arr(4) is also used as error flag
    ss_env     = y_arr(4)
    ql_prtcl   =  y_arr(8)


  end subroutine

  subroutine prtcl_growth_const_env(T_env, ss_env, p_env, m_aero, radius, time_start,time_end)
    use const, only :  pi43, rho_w
    implicit none
    !grow particle at constant temperature and supersaturation
    double precision,intent(in)     :: T_env,  p_env, m_aero, ss_env
    double precision,intent(inout)  :: radius
    !local variables
    double precision                :: qv_env, ql_env, qv_sat, t1, t2, time_start, time_end
    integer                         :: has_cs_effect, debug_lcl
    double precision                :: T_temp, qv_temp, p_temp, ss_temp, ql_temp,r_aero, tstep

    
    debug_lcl           = 0
    has_cs_effect       = 1
    tstep               = 0.01

    qv_sat              = sat_mixing_ratio( T_env - 273.15, p_env)
    qv_env              = (1.0d0 + ss_env/100) * qv_sat
    
    r_aero              = m_aero/(rho_aero*pi43)
    r_aero              = r_aero**(1/3.0d0)
    ql_env              = pi43 * (radius-r_aero)**3 * rho_w * inv_vol_gcell
    
    
    t1                  = time_start
    if ( debug_lcl == 1) then
      write(7,"(A,F12.2,A,F12.2)")  "time start:", time_start,   " end:", time_end
      write(7,"(A,E12.6)")   "Before prtcl growth r      = :", radius
      write(7,"(A,F12.6,A)") "Before prtcl growth ss     = :", ss_env, " %" 
      write(7,"(A,F12.6)")   "Before prtcl growth T      = :", T_env
      write(7,"(A,E12.6)")   "Before prtcl growth qv     = :", qv_env
    end if   

    do while(t1 < time_end)
      t2                = t1 + tstep
      !reset the environment to constant value
      T_temp            = T_env
      qv_temp           = qv_env
      p_temp            = p_env
      ss_temp           = (ss_env)/100.0d0

      ql_temp           = pi43 * (radius-r_aero)**3 * rho_w * inv_vol_gcell
      call indiv_prtcl_growth(T_temp, qv_temp, p_temp, ss_temp, ql_temp, t1, t2, radius, m_aero, aero_type, inv_vol_gcell, &
           has_cs_effect,0,0)

      t1                = t2
    end do

    if ( debug_lcl == 1) then
      write(7,"(A,E12.6)")   "After prtcl growth r      = :", radius
      write(7,"(A,F12.6,A)") "After prtcl growth ss     = :", ss_temp*100, " %" 
      write(7,"(A,F12.6)")   "After prtcl growth T      = :", T_temp
      write(7,"(A,E12.6)")   "After prtcl growth qv     = :", qv_temp
    endif
      
  end subroutine
 
  subroutine kohler_curve_from_dgm(m_aero,T_env)
    use const, only : pi43, rho_w 
    implicit none
    double precision :: m_aero, s_max, p_env, s_min, s_env, r_aero , T_env
    double precision :: time_start, time_end, radius1, radius2, s_step, ss_env
    integer          :: is_equ, debug_lcl, fhdl_kcurve

    debug_lcl     = 0
    fhdl_kcurve   = 601
    open(fhdl_kcurve,file="kohler_dgm2.txt")

    r_aero        = m_aero/rho_aero
    r_aero        = (r_aero/pi43)**(1.d0/3.d0)
    p_env         = 0.101000E+06
    
    s_min        = 40.0d0 !%
    s_max        = 101    !%
    s_step       = 0.01!%
    !s_min        = 99.0d0 !%
    !s_max        = s_min + 0.01d0    !%
    
    s_env        = s_min 
    do while (s_env < s_max)
      radius1     = r_aero
      time_start  = 0
      !Allow 10 seconds for droplet to reach equilibrium
      time_end    = time_start + 10
      ss_env      = s_env - 100
      call prtcl_growth_const_env(T_env, ss_env, p_env, m_aero, radius1, time_start,time_end)

      if ( debug_lcl == 1 )then
        write(7,'(A,F12.6,A)') "==== ss env :", s_env, " %"
        write(7,"(A,E12.6)")   "prtcl radius r after 10 secs     = :", radius1
      end if 

      radius2     = r_aero
      time_start  = 0
      !Allow 11 seconds for droplet to reach equilibrium
      time_end    = time_start + 11
      call prtcl_growth_const_env(T_env, ss_env, p_env, m_aero, radius2, time_start,time_end)
      if ( debug_lcl == 1 )then
        write(7,"(A,E12.6)")   "prtcl radius r after 11 secs     = :", radius2
        
      end if  
      
      write(7,"(A,E12.6)")   "diff radius1-radius2               = :", abs(radius1-radius2)*1e6

      !If the particle is in stable equilibrium then radius1 and radius2 will be same (very close)
      if( abs((radius1-radius2)*1e6) < 1e-2)then
        is_equ    = 1
        write(7,'(A,F12.6,A)') "**** Alert: droplet growth is in  stable equilibrium at saturation :", s_env , " %" 
        write(fhdl_kcurve,'(2E16.8)') s_env, radius1
        
      else
        is_equ    = 0
        write(7,'(A,F12.6,A)') "**** Alert: droplet growth is in unstable equilibrium at saturation :", s_env , " %" 
        exit

      end if
      s_env      = s_env + s_step
    end do  

    write(7,'(A,F12.6,A)') "Loop ended at saturation : ", s_env , " %"
    
    close(fhdl_kcurve)
  end subroutine  

  subroutine kohler_curve_from_RY(m_aero,T_env)
   ! Critical radius and supersaturation (for diagnostic output)
   ! Use Rogers & Yau Eqs. 6.7 and 6.8 
   
    implicit none
   ! T_env: temperature (K) and  m_aero:   solute mass (kg)
    double precision,intent(in)       :: T_env, m_aero

    !local variables 
    double precision,parameter        :: m_per_cm = 0.01
    double precision                  :: a_RY, b_RY
    double precision                  :: r_min, r_max, r_prtcl, r_step, s_equ
    integer                           :: fhdl_kcurve

    fhdl_kcurve    = 602
    !T_RY: temperature (K)
    a_RY      = 3.3e-5/T_env ! cm
    !a_RY      = 3.2038e-5/T_env ! cm
    a_RY      = a_RY * m_per_cm
  
    b_RY      = 4.3 * m_aero * nions_aero/ molar_mass_aero ! cm^3
    b_RY      = b_RY * m_per_cm**3.0d0
    
    open(fhdl_kcurve,file="kohler_RY2.txt")

    r_min     = 0.07e-6
    r_max     = 10.0e-6
    r_step    = 0.001e-6

    r_prtcl   = r_min
    do while (r_prtcl < r_max)
      !write(7,'(A,E12.6)') "r_prtcl :" , r_prtcl
      s_equ     = 1+ a_RY/r_prtcl - b_RY/r_prtcl**3.0d0
      write(fhdl_kcurve,'(2E16.8)') s_equ*100.0d0, r_prtcl
      r_prtcl = r_prtcl + r_step
    end do  

    write(7,'(A)') "Kohler curve from RY completed"
    close(fhdl_kcurve)
    
  end subroutine

end program