program test_microphysics
  use microphysics, only : init_microphysics, add_prtcl_to_dom, move_prtcl_by_triplet_map,prtcl_growth, prtcl_settling,&
                           inject_n_tdiff, inject_stime, inject_etime, &
                            indiv_prtcl_term_vel, indiv_prtcl_growth, sat_mixing_ratio, m0_aero, &
                           inv_vol_gcell, create_bedges_xpos, create_bedges_r_prtcl, bin_data, verify_prtcl_count,&
                           print_nprtcl_in_gcell, find_r_crit, find_r_crit_dgm, rho_aero, nions_aero, molar_mass_aero, &
                           vol_gcell
  implicit none
  character(len=:), allocatable :: namelist_fname, idir,odir, realz_str, cmd_str, exp_name
  integer, parameter :: N = 6000
  integer            ::  j, n_bedges_xpos, n_bedges_rprtcl
  double precision   :: T_o, Tdif, press, index_arr(N), ref_index(N), nu, H, tmax, td, time1, time_diff, test_data(N), rseed1
  double precision,allocatable :: bedges_xpos(:), bedges_rprtcl(:), freq_xpos(:), freq_rprtcl(:)
  double precision   :: T_env,m_aero,r_crit,ss_crit,r_crit_dgm, ss_crit_dgm

  !external function
  double precision    ::  rand1

  idir           = "input/steady_state/";
  namelist_fname = "nml_mphy_inj10_vol0.01";
  exp_name       = "test_microphy"
  realz_str      = "001" 
  odir           = "./output/test/"//realz_str
  nu             = 1.41E-5
  H              = 1.0 
  tmax           = 1.0d-3
  td             = 1.0d1/(1.0d0 * N * N)
  T_o            = 15
  Tdif           = 12
  press          = 1013.10 * 100
  time_diff      = tmax/1e5
  time1          = 0.0d0

  cmd_str        = "mkdir -p "//odir
  call system(cmd_str)
  open(7,file="test_microphysics.log") 

  !call init_microphysics(idir,nml_fname,odir, H, N, tmax, td, T_o, Tdif, press, fhdl_mp_log,realz_str, debug_mode)
  call init_microphysics(idir,namelist_fname,odir, exp_name,  H, N, tmax, td, T_o, Tdif, press,7,realz_str,1)
  ! !adding ccn/droplets to the domain
  ! write(7,"(/,A)") "====================== Unit test for adding droplets ============"
  ! write(7,"(/,A)") "************ Adding droplets 1st time ************"
  ! call add_prtcl_to_dom(time_diff,time1)
  ! call verify_prtcl_count()
  ! call print_nprtcl_in_gcell()

  ! write(7,"(/,A)") "************ Adding droplets 2nd time ************"
  ! call add_prtcl_to_dom(time_diff, time1 + time_diff)
  ! call verify_prtcl_count()
  ! call print_nprtcl_in_gcell()

  ! move droplets by eddy mapping; set following values in namelist to make sure all gcells have a droplet
  ! match the gcell pos difference and prtcl pos difference
  ! max_prtcl_dom        = 60000, max_prtcl_gcell      = 100, max_prtcl_alltime    = 120000,
  ! n_prtcl_uvol_utime  = 4.0E011,   add_drop_vol_frac  = 1.0
  

  ! write(7,"(/,A)") "====================== Unit test for moving droplets ============"
  ! write(7,"(/,A)") "************ Eddy size - 12 grid cells ************"
  ! do j = 1, N
  !   index_arr(j) = j * 1.0d0
  !   ref_index(j) = j * 1.0d0
  ! end do
  ! write(7,*) "Before triplet map - sum of differences (old pos - new):" , sum(abs(index_arr-ref_index))
  ! write(7,*) "old index position:"
  ! write(7,'(12F6.1)') index_arr(1:12)
  ! call triplet(N,1,12,index_arr);
  ! write(7,*) "new index position:"
  ! write(7,'(12F6.1)') index_arr(1:12)
  ! write(7,*) "After triplet map - sum of differences:" , sum(abs(index_arr-ref_index))
  ! call move_prtcl_by_triplet_map(index_arr)

  ! write(7,"(/,A)") "************ Eddy size - 2000 grid cells ************"
  ! do j = 1, N
  !   index_arr(j) = j * 1.0d0
  !   ref_index(j) = j * 1.0d0
  ! end do
  ! write(7,*) "Before triplet map - sum of differences(old pos - new):" , sum(abs(index_arr-ref_index))
  ! call triplet(N,1,900,index_arr);
  ! write(7,*) "After triplet map - sum of differences (old pos - new):" , sum(abs(index_arr-ref_index))
  ! call move_prtcl_by_triplet_map(index_arr)

  !write(7,"(/,A)") "====================== Unit test r_prtcl vs terminal velocity ============"
  !call test_rprtcl_vs_termvel()
  
  
  write(7,"(/,A)") "======================= Unit test for single particle growth ============"
  call test_prtcl_growth()
  
  !write(7,"(/,A)") "======================= Unit test iterative particle growthh ============"
  !call test_prtcl_growth_iterate()

  !write(7,"(/,A)") "======================= Unit test bin edges and bins ==================="
  ! bedges_xpos     = create_bedges_xpos()
  ! n_bedges_xpos   = size(bedges_xpos,dim=1)
  ! write(7,=) "n_bedges_xpos:" , n_bedges_xpos
  ! write(7,'(A,200F12.2)') "bedges_xpos:" , bedges_xpos

  ! bedges_rprtcl    = create_bedges_r_prtcl()
  ! n_bedges_rprtcl  = size(bedges_rprtcl,dim=1)
  ! write(7,*) "n_bedges_rprtcl:" , n_bedges_rprtcl
  ! write(7,'(A,200E12.2)') "bedges_rprtcl:" , bedges_rprtcl*1.0d6

  ! rseed1          = time() + 10000000.0d0
  ! do j = 1, N
  !   test_data(j)  = rand1(rseed1)
  ! end do  
  ! freq_xpos       = bin_data(bedges_xpos, n_bedges_xpos, test_data, N)
  ! open(10, file="output/test_microphy/001/PDF_xpos_prtcl_test.txt")
  ! write(10, '(A)') "#This file has the frequency of droplet position as fraction of domain length"
  ! write(10, '(A,I6,A)') "The position data is grouped into ", n_bedges_xpos-1, " bins."
  ! write(10, '(A)')   "First line of data is bin center"
  ! write(10, '(A)')   "Second line of data is bin width"
  ! write(10,'(200F12.4)')  (bedges_xpos(2:n_bedges_xpos) + bedges_xpos(1:n_bedges_xpos-1)) * 0.5d0 
  ! write(10,'(200F12.4)')  (bedges_xpos(2:n_bedges_xpos) - bedges_xpos(1:n_bedges_xpos-1)) 
  ! write(10,'(200F12.4)')  freq_xpos
  ! close(10)
  ! write(7,*) "Total frequency of prtcl pos bins : ", sum(freq_xpos)
  
  ! !create test data that is 1 to 100 micron diameter
  ! test_data       = 1.0d-6 + test_data *  39.d0 * 1.0d-6 
  ! freq_rprtcl      = bin_data(bedges_rprtcl, n_bedges_rprtcl, test_data, N)
  ! open(10, file="output/test_microphy/001/PDF_r_prtcl_test.txt")
  ! write(10, '(A)')      "#This file has the frequency of droplet radius"
  ! write(10, '(A,I6,A)') "#The radius data is grouped into ", n_bedges_rprtcl-1, " bins."
  ! write(10, '(A)')      "#First line of data is bin center"
  ! write(10, '(A)')      "#Second line of data is bin width"
  ! write(10,'(200E12.4)')  (bedges_rprtcl(2:n_bedges_rprtcl) + bedges_rprtcl(1:n_bedges_rprtcl-1)) * 0.5d0 
  ! write(10,'(200E12.4)')  (bedges_rprtcl(2:n_bedges_rprtcl) - bedges_rprtcl(1:n_bedges_rprtcl-1)) 
  ! write(10,'(200F12.2)')  freq_rprtcl
  ! close(10)
  ! write(7,*) "Total frequency of prtcl radius bins: ", sum(freq_rprtcl)

  ! write(7,"(/,A)") "======================= Unit test r_crit for given aerosol mass ============"
  ! !call find_aero_r_crit(T_env,m_aero,r_crit,ss_crit)
  !T_env  = 288
  !case study value
  !m_aero  = 2.212E-18
  ! r_crit  = 0.9232E-6
  ! ss_crit  = 0.8274E-1
  ! r_crit_dgm     =  0.66256512618537078     !use r_prtcl, rho_prtcl
  ! s_crit_dgm     =  0.11193110845725038     !use r_prtcl, rho_prtcl

  ! !from steve
  ! ! m_aero      = 2.4847e-18
  ! ! r_crit      = 0.9784E-6
  ! ! ss_crit     = 0.7807E-1
  
  ! call find_r_crit(T_env,m_aero,r_crit,ss_crit)
  ! write(7,*) 'mass of the solute =', m_aero
  ! write(7,*) 'critical radius (um)         =', r_crit  * 1.0E6
  ! write(7,*) 'critical supersaturation (%) =', ss_crit * 100. 

  ! call find_r_crit_dgm(T_env, m_aero, r_crit, r_crit_dgm, ss_crit_dgm)
  ! write(7,*) 'mass of the solute =', m_aero
  ! write(7,*) 'critical radius from dgm (um)         =', r_crit_dgm  * 1.0E6
  ! write(7,*) 'critical supersaturation from dgm (%) =', ss_crit_dgm * 100. 

  !write(7,"(/,A)") "======================= Unit test kohler curve ============"
  ! m_aero        = 2.212E-18
  ! T_env         = 288.0d0
  ! call kohler_curve_from_dgm(m_aero,T_env)
  ! call kohler_curve_from_RY(m_aero,T_env)

  close(7);
  stop
  contains 

  subroutine triplet(N,M,L,psi)
    ! this subprogram is a subroutine in ODT code and I just copied it here w/o changes
    implicit none
    integer N, M, L, Lo, j, k
    double precision psi(N), x(N)
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

  subroutine test_rprtcl_vs_termvel()
    implicit none
    
    !Subprogram to test the terminal velocity function
    !for droplet radius from 1 to 40 microns, compute the terminal velocity

    double precision   ::  T_env, qv_env, p_env, radius, term_vel, vel_factor
    integer            :: j

    write(7,*) "test_rprtcl_vs_termvel"
    open(200,file = "prtcl_vs_termvel.txt")


    p_env         = 1013.10 * 100
    T_env         = 300 
    qv_env        = sat_mixing_ratio(T_env-273.15,p_env)
    vel_factor    = 1.0d0
    do j = 1, 40
      radius = j * 1.0d-6 !micron 
      !call indiv_prtcl_term_vel(T_env, qv_env, p_env, r_prtcl,vel_factor, term_vel)
      call indiv_prtcl_term_vel(T_env, qv_env, p_env,radius, vel_factor, term_vel)
      write(200,"(2F12.4)") radius*1.0d6, term_vel
    end do
    
    close(200)

  end subroutine

  subroutine test_prtcl_growth()
    !subprogram to find the droplet equilibrium radius
    use const, only : pi43, rho_w
    double precision   ::  T_env, qv_env, p_env, radius, time1, time2, qv_sat, delta_qv, r1, qv1, r2, qv2, delta_ql, eps_err
    double precision   :: ss_env, ss1, ss2, delta_ss, T1, T2, ss2_calc, ql_env, ql1, ql2, ql2_calc, r0_aero, ss_env_calc
    integer            :: has_cs_effect

    write(7,"(/,A)") "************************************ test indvidual droplet  growth"
    
    !dgm model works good;
    ! time1         = 0.0
    ! time2         = 0.1d0
    ! radius        = 1.0d-6
    ! has_cs_effect = 1
    ! p_env         = 1013.10 * 100
    ! T_env         = 300 
    ! qv_sat        = sat_mixing_ratio(T_env-273.15,p_env)
    ! write(7,"(A,F12.1,A,E12.4)") "Sat mixing ration at T", T_env, " is :" , qv_sat * 1.0d3
    ! qv_env        = 1.01 * qv_sat
    ! ss_env        = (qv_env/qv_sat - 1 )
    ! ql_env        = pi43 * rho_w  * inv_vol_gcell * radius**3;

    !values dgm model works but produces large temperature and humidity perturbations for large time steps
    !a sign that its failing

    !yarr values (r, qv, T, ss, p, w, z, ql)
    !0.1000E-05  0.1142E-01  0.2882E+03  0.6982E-01  0.1010E+06  0.0000E+00  0.0000E+00  0.2513E-04
    !
    !print temperature value in fcbkb.f90 for all calls to it helped me see what is going on
    !delta_time of 2   seconds produces the numerical error since temperature swings from small to large values
    !amplitudes around initial value finally reaches negative value throwing error
    
    !delta time of 1   seconds doesn't produce numerical error but temperature swings 
    !outside Tbot and Ttop but reaches a stable value between Tbot and Top
    
    !delta time of 0.5 seconds  doesn't produce numerical error but small temperature swings(0.5 K) within Tbot and Top

    !delta time of 0.1 seconds doesn't produce numerical error but smaller temperature swings(0.02K) within Tbot and Top

    !delta time of 0.05 seconds doesn't produce numerical error but smaller temperature swings(0.01K) within Tbot and Top

    !delta time of 0.01 seconds doesn't produce numerical error but smaller temperature swings(0.001K) within Tbot and Top
    
    !delta time of 0.001 seconds doesn't produce numerical error but smaller temperature swings(0.0001K) within Tbot and Top

    !a timestep < 0.1 seems suitable for dgm model  
    ! time1         = 0.0
    ! time2         = 0.001
    ! has_cs_effect = 1
    ! radius        = 0.1000E-05
    ! p_env         = 0.1010E+06
    ! T_env         = 0.2882E+03
    ! qv_env        = 0.1142E-01
    ! ss_env        = 0.6982E-01
    ! ql_env        = pi43 * rho_w  * inv_vol_gcell * radius**3;
    
    !testing some error situation in ODT
    ! t1, t2      1.8669      3.7339 
    !yarr values (r, qv, T, ss, p, w, z, ql)
    ! 0.20000000E-04  0.11542476E-01  0.28832524E+03  0.70088446E-01  0.10100000E+06  0.00000000E+00  0.00000000E+00  0.20106193E+00
    !
    !error when injecting dry aerosol
    ! t1, t2:          6.8126          6.8256
    !!yarr values (r, qv, T, ss, p, w, z, ql)
    ! 0.625000E-07    0.113124E-01    0.288009E+03    0.696259E-01    0.101000E+06    0.000000E+00    0.000000E+00   -0.118935E-11
    
    !error when injecting dry aerosol
    !I could not reproduce error here
    ! t1, t2:          2186.6031          2186.6126
    !!yarr values (r, qv, T, ss, p, w, z, ql)
    !0.630000E-07    0.112612E-01    0.288963E+03    0.390479E-03    0.101000E+06    0.000000E+00    0.000000E+00    0.104739E-11

    !error when injecting dry aerosol
    !I could not reproduce error
    !t1, t2:          3426.3858          3426.3880
    !!yarr values (r, qv, T, ss, p, w, z, ql)
    !0.630000E-07    0.105923E-01    0.288006E+03    0.172437E-02    0.101000E+06    0.000000E+00    0.000000E+00    0.104739E-11

    !error changing the area frac - Nan value
    !!yarr values (r, qv, T, ss, p, w, z, ql)
    !0.630000E-07    0.109306E-01    0.288474E+03    0.254828E-02    0.101000E+06    0.000000E+00    0.000000E+00    0.523608E-16
    !t1, t2: 38.7115, 38.7302

    !0.630000E-07    0.108602E-01    0.287421E+03    0.672762E-01    0.101000E+06    0.000000E+00    0.000000E+00    0.523608E-16
    !t1, t2: 16.7109 , 16.7295
    !The run failed for vol_cell 0.167d-11  with inv_vol_gcell = 0.6d12
    !the domain volume is 0.01 * 0.01 * 1.0 * 0.0001 (area_frac) = 1d-8 m3;   vol_gcell = 1d-8/6000 = 0.1667d-11
    !the droplet growth model failed if the gcel volume is smaller than  0.5d-11

    !don't forget to change mphy namelist file name at the top before testing
    eps_err       = 1e-4
    m_aero        = 2.212E-18

    time1         = 16.7070
    time2         = 16.7256
    radius        = 0.630000E-07
    T_env         = 0.288106E+03
    p_env         = 0.101000E+06
    qv_env        = 0.113889E-01
    ql_env        = 0.523608E-16
    ss_env        = 0.700425E-01 

    qv_sat        = sat_mixing_ratio(T_env-273.15,p_env)
    write(7,"(A,F12.1,A,E12.4)") "Sat mixing ratio at T", T_env, " is :" , qv_sat * 1.0d3
    ss_env_calc   = (qv_env/qv_sat-1)
    write(7,"(A,F12.1,A,E12.4,A)") "SS calc:", T_env, " is :" , ss_env_calc * 1.0d2 , " %"
    !ql_env       = pi43*radius**3.0 * rho_w
    ! r0_aero     = m0_aero/rho_aero
    ! r0_aero     = (r0_aero/pi43)**(1.d0/3.d0) 
    has_cs_effect = 1

    write(7,"(A,E12.4)") "Solute  mass         :", m0_aero
    write(7,*) ""

    r1            = radius;
    T1            = T_env; 
    qv1           = qv_env
    ss1           = ss_env 
    ql1           = ql_env
    write(7,"(A,E12.6)")   "Before prtcl growth r      = :", r1
    write(7,"(A,F12.6,A)") "Before prtcl growth ss     = :", ss1 *1e2, " %" 
    write(7,"(A,F12.6)")   "Before prtcl growth T      = :", T1
    write(7,"(A,E12.6)")   "Before prtcl growth qv     = :", qv1
    
    write(7,"(A,E12.6)")   "Before prtcl growth ql     = :", ql1
    write(7,"(A,E12.6,A)") "inv_vol_gcell              = :", inv_vol_gcell
    write(7,"(A,E12.6,A)") "vol_gcell                  = :", vol_gcell

    call indiv_prtcl_growth(T_env, qv_env, p_env, ss_env, ql_env, time1, time2, radius, m0_aero, 1,inv_vol_gcell,&
            has_cs_effect,eps_err,1,1)
    r2      = radius; 
    T2      = T_env; 
    qv2     = qv_env
    ss2     = ss_env 
    ql2     = ql_env
    qv_sat  = sat_mixing_ratio(T2-273.15,p_env)
    ss2_calc = (qv2/qv_sat -1 )
    ql2_calc = pi43 * rho_w  * inv_vol_gcell * r2**3 

    write(7,*) ""
    write(7,"(A,E12.6)")   "After prtcl growth r                  = :", r2
    write(7,"(A,F12.6,A)") "After prtcl growth ss_from_qv_n_T     = :", ss2_calc*1e2 , " %" 
    write(7,"(A,F12.6,A)") "After prtcl growth ss_from_integration= :", ss2*1e2 , " %" 
    write(7,"(A,F12.6)")   "After prtcl growth T                  = :", T2
    write(7,"(A,E12.6)")   "After prtcl growth qv                 = :", qv2
    
    write(7,"(A,E12.6,A)") "After prtcl growth ql_from_integration= :", ql2
    write(7,"(A,E12.6,A)") "After prtcl growth ql_from_r          = :", ql2_calc

    delta_qv  = qv2 - qv1
    !divide by cell volume to covert it to kg/m^3 similar units as delta_qv (kg/m^3 or Kg/kg)
    delta_ql  = pi43 * rho_w  * inv_vol_gcell * (r2**3  -  r1**3)
    !delta_ql  = ql2 - ql1
    delta_ss  = ss2  - ss1 
    write(7,*) "delta qv:" , delta_qv
    write(7,*) "delta ql:" , delta_ql , "%"
    write(7,*) "delta ss:" , delta_ss , "%"

  end subroutine

  subroutine prtcl_growth_const_env(T_env, ss_env, p_env, m_aero, radius, time_start,time_end)
    use const, only :  pi43, rho_w
    implicit none

    double precision,intent(in)     :: T_env,  p_env, m_aero, ss_env
    double precision,intent(inout)  :: radius
    !local variables
    double precision                :: qv_env, ql_env, qv_sat, t1, t2, time_start, time_end, eps_err
    integer                         :: has_cs_effect, debug_lcl
    double precision                :: T_temp, qv_temp, p_temp, ss_temp, ql_temp,r_aero

    
    debug_lcl           = 0
    eps_err             = 1.0d-4    
    qv_sat              = sat_mixing_ratio( T_env - 273.15, p_env)
    qv_env              = (1.0d0 + ss_env/100) * qv_sat
    
    r_aero              = m_aero/(rho_aero*pi43)
    r_aero              = r_aero**(1/3.0d0)
    
    ql_env              = pi43 * radius**3 * rho_w * inv_vol_gcell
    has_cs_effect       = 1
    
    t1                  = time_start
    if ( debug_lcl == 1) then
      write(7,"(A,F12.2,A,F12.2)")  "time start:", time_start,   " end:", time_end
      write(7,"(A,E12.6)")   "Before prtcl growth r      = :", radius
      write(7,"(A,F12.6,A)") "Before prtcl growth ss     = :", ss_env, " %" 
      write(7,"(A,F12.6)")   "Before prtcl growth T      = :", T_env
      write(7,"(A,E12.6)")   "Before prtcl growth qv     = :", qv_env
    end if   

    do while(t1 < time_end)
      t2                = t1 + 0.01
      !call indiv_prtcl_growth(T_env, qv_env, p_env, ss_env, ql_prtcl, t1, t2, r_prtcl, m0_aerosol, aerosol_type,&
      ! invvol_gcell, has_cs_effect, debug)
      T_temp            = T_env
      qv_temp           = qv_env
      p_temp            = p_env
      ss_temp           = (ss_env)/100.0d0
      ql_temp           = pi43 * (radius)**3 * rho_w * inv_vol_gcell
      call indiv_prtcl_growth(T_temp, qv_temp, p_temp, ss_temp, ql_temp, t1, t2, radius, m_aero, 1, inv_vol_gcell, &
           has_cs_effect,eps_err,1,0)

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
    double precision :: m_aero, s_max, p_env, s_min, qv_env, qv_sat, s_env, r_aero , T_env
    double precision :: time_start, time_end, radius1, radius2, s_step, ss_env
    integer          :: is_equ, debug_lcl, fhdl_kcurve

    debug_lcl     = 0
    fhdl_kcurve   = 601
    open(fhdl_kcurve,file="kohler_dgm.txt")

    r_aero        = m_aero/rho_aero
    r_aero        = (r_aero/pi43)**(1.d0/3.d0)
    ! r_crit      = 0.92315830E-06  
    ! ss_crit     = 8.274733E-4
    ! r_hdlt      = 0.4E-6
    !ss_hdplt     = 0.035E-2
    p_env         = 0.101000E+06
    
    s_min        = 99.0d0 !%
    s_max        = 101    !%
    !s_min        = 99.0d0 !%
    !s_max        = s_min + 0.01d0    !%
    
    s_env        = s_min 
    s_step       = 0.001!%
        
    do while (s_env < s_max)
      radius1     = r_aero
      time_start  = 0
      time_end    = time_start + 10
      ss_env      = s_env - 100
      call prtcl_growth_const_env(T_env, ss_env, p_env, m_aero, radius1, time_start,time_end)

      if ( debug_lcl == 1 )then
        write(7,'(A,F12.6,A)') "==== ss env :", s_env, " %"
        write(7,"(A,E12.6)")   "prtcl radius r after 10 secs     = :", radius1
      end if 

      radius2     = r_aero
      time_start  = 0
      time_end    = time_start + 11
      call prtcl_growth_const_env(T_env, ss_env, p_env, m_aero, radius2, time_start,time_end)
      if ( debug_lcl == 1 )then
        write(7,"(A,E12.6)")   "prtcl radius r after 11 secs     = :", radius2
        
      end if  
      
      write(7,"(A,E12.6)")   "diff radius1-radius2               = :", abs(radius1-radius2)*1e6

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
    
    open(fhdl_kcurve,file="kohler_RY.txt")

    r_min     = 0.1e-6
    r_max     = 10.0e-6
    r_prtcl   = r_min
    r_step    = 0.01e-6
    do while (r_prtcl < r_max)
      write(7,'(A,E12.6)') "r_prtcl :" , r_prtcl
      s_equ     = 1+ a_RY/r_prtcl - b_RY/r_prtcl**3.0d0
      write(fhdl_kcurve,'(2E16.8)') s_equ*100.0d0, r_prtcl
      r_prtcl = r_prtcl + r_step
    end do  

    write(7,'(A)') "Kohler curve from RY completed"
    close(fhdl_kcurve)
  end subroutine

end program