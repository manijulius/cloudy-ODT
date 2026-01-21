!Mani Rajagopal
!Aug 2023 - Jan 2024
!This module adds the following microphysical processes to the ODT model
!Move particle by eddies, condensational growth, and gravitational settling

module microphysics
  implicit none
  !Make all variables and functions private
  !Make only variables and functions functions used in ODT as public
  !any variable added should be allocated and initialized right away
  private 
  save

  !Buffer constants - store stats in a buffer to reduce number of write to file to improve speed
  integer, parameter:: n_buff_evol = 1000, n_buff_pdf =1000, n_buff_prof=1000 
  !Five particle types:  dry aero, haze dplt, cloud dplt, intru obs, all
  integer, parameter:: n_ptype = 5, iaero=1, ihdplt=2, icdplt=3, i_instru=4, iprtcl=5


  !domain variables;
  !len_dom = height H from ODT
  double precision              :: area_frac, width_dom
  double precision              :: len_dom, vol_dom, vol_gcell, inv_vol_gcell, len_gcell, inv_vol_dom

  character(len=:),allocatable  :: input_dir, namelist_file, output_dir, exp_name, realz
  character(len=256)            :: parent_dir
   
  !random seed
  double precision              :: rand_seed1, rand_seed2, rand_seed3

  !domain boundary conditions
  double precision              :: T_bot, T_mid, T_top, qvs_bot, qvs_top, T_diff, qvs_diff
  !time realted variables
  double precision              :: tscale_nu, sim_time, sim_time_ndim, tintvl_diff, tintvl_diff_ndim

  !debug flag to print statement in each processes/function
  integer                       :: debug_gbl

  !aerosol effects
  !Three aero types 1) Nacl, 2) Ammonium Sulphate, 3) Ammonium bisulphate
  !These number codes should match fcnkb.f90
  integer                       :: has_curv_sol_eff, aero_type

  !aero related properties
  double precision              :: molar_mass_aero
  double precision              :: rho_aero, r_min_instru
  integer                       :: nions_aero

  !different methods to inject particles in to domain
  !1)monodisperse, 2) Gamma distribution 3) from a files with aero mas, particle size and count
  integer                       :: add_prtcl_mthd
  
  !monodisperse particle injection parameters
  double precision              :: m0_aero, r0_prtcl

  !Gamma distribution particle injection parameters
  double precision              :: lwc0_prtcl_sum
  integer                       :: nbins_prtcl_dis

  !User defined particle distribution with nbins_prtcl_dis 
  character(len=64)             :: dis_prtcl_file
  
  !injection rate
  integer                       :: is_injected_onetime,  inject_onetime, is_inj_blob_rand, inject_n_tdiff
  double precision              :: n_prtcl_uvol_utime,  n_prtcl_uvol_onetime
  integer                       :: max_prtcl_dom, max_prtcl_gcell
  double precision              :: inject_vol,  inject_at, inject_stime, inject_etime
  

  !injection partcile size distribution
  double precision, allocatable :: dis_r_prtcl(:), dis_r_crit(:),  dis_m_aero(:), dis_r_aero(:)
  double precision, allocatable :: dis_n_uvol_utime(:), dis_n_bvol_utime(:)
  double precision, allocatable :: dis_n_residue(:), dis_n_alltime(:)
  
  !book keeping of injection
  integer                       :: n_inj_alltime
  double precision              :: inject_time(n_buff_evol)
  integer                       :: inject_prtcl_count(n_buff_evol)

  !properties of individual particles currently in the domain
  !Throughout the simulation new particles are added to the domain through injection
  !while  some of the particles are removed through settling
  !properties of particles currenlty active in the domain  are stored in the following variables 
  integer                       :: n_prtcl_curr, n_inj_curr, n_fallout_curr, n_actvd_curr, n_deactvd_curr
  double precision, allocatable :: r_prtcl_curr(:), m_aero_curr(:), x_prtcl_curr(:)
  double precision, allocatable :: r_crit_curr(:), ql_prtcl_curr(:), r_aero_curr(:)
  integer, allocatable          :: gcell_prtcl_curr(:), unqid_prtcl_curr(:)

  !Store activation or deactivation of a particle 
  !(1) -> activation (-1) -> deactivation  (0) -> no change in particle type
  !sum of values in this array will give net activation/deactivation
  integer, allocatable          :: is_actvd_curr(:)

  !properties of particles that fallout at bottom boundaries
  double precision, allocatable :: r_prtcl_fallout(:)
  double precision, allocatable :: r_crit_fallout(:), ql_prtcl_fallout(:), r_aero_fallout(:)

  !total count
  integer                       :: n_prtcl_alltime, prtcl_counter, n_fallout_alltime

  !variable related to grid cells and particle moved by eddies - eddy mapping
  integer                       :: n_gcell, gcell_1_idx, gcell_n_idx
  integer, allocatable          :: gcell_prtcl_count(:),  gcell_prtcl_currid_list(:,:)
  integer, allocatable          :: did_edymov_prtcl(:)

  !terminal velocity factor for particle settling
  !default value 1.0. It can be tuned to a different value to adjust for tubrulence
  double precision              :: term_vel_factor
  

  !recording pdf, prof, and particle history
  integer                       :: do_rec_prtcl_hist,    do_rec_pdf_n_prof
  double precision              :: rec_stime_pdf_n_prof, rec_etime_pdf_n_prof
  double precision              :: rec_stime_prtcl_hist, rec_etime_prtcl_hist
  integer                       :: rec_intvl_prtcl_hist  
  integer                       :: write_intvl_prtcl_hist

  !buffer for particle history
  integer                       :: max_prtcl_hist, bidx_rec_hist
  !rec particle history - their properties and env (T, qv, ss) at each diffusion time step
  integer, allocatable          :: uniqid_prtcl_hist(:)
  double precision, allocatable :: r_prtcl_hist(:), r_aero_hist(:), r_crit_hist(:), x_prtcl_hist(:)
  double precision, allocatable :: T_prtcl_hist(:), qv_prtcl_hist(:), ss_prtcl_hist(:)
  double precision, allocatable :: time_prtcl_hist(:), deltat_prtcl_hist(:)
  double precision              :: time_next_phist_rec, time_next_phist_write
  

  !data arrays
  !pdf bins for data
  integer                       :: n_bedges_r_prtcl, n_bins_r_prtcl, n_bedges_pprof, n_bins_pprof
  integer                       :: n_bedges_T, n_bins_T, n_bedges_qv, n_bins_qv, n_bedges_ss, n_bins_ss
  double precision, allocatable :: bedges_T(:), bedges_qv(:), bedges_ss(:)
  double precision, allocatable :: bedges_r_prtcl(:), bedges_pprof(:)

  double precision              :: pdf_stime, pdf_etime, pdf_deltat_sum, prof_stime, prof_etime
  double precision,allocatable  :: pdf_dom_T(:),pdf_bulk_T(:), pdf_dom_qv(:), pdf_bulk_qv(:)
  double precision,allocatable  :: pdf_dom_ss(:), pdf_bulk_ss(:)
  double precision,allocatable  :: pdf_r_prtcl(:), pdf_r_aero(:), pdf_r_hdplt(:), pdf_r_cdplt(:), pdf_r_instru(:)

  !buffer array for microphysics evolution data
  double precision              :: rmin_ptype_evol(n_ptype,n_buff_evol), rmax_ptype_evol(n_ptype,n_buff_evol)
  double precision              :: rmean_ptype_evol(n_ptype,n_buff_evol), ql_ptype_evol(n_ptype,n_buff_evol)
  double precision              :: n_ptype_evol(n_ptype,n_buff_evol), r2mean_ptype_evol(n_ptype,n_buff_evol)

  double precision              :: time_prtcl_evol(n_buff_evol),deltat_prtcl_evol(n_buff_evol)
  double precision              :: n_prtcl_alltime_evol(n_buff_evol)
  double precision              :: n_inj_evol(n_buff_evol),  n_fallout_alltime_evol(n_buff_evol)
  double precision              :: n_fall_ptype_evol(n_ptype,n_buff_evol), rmean_fall_ptype_evol(n_ptype,n_buff_evol)
  double precision              :: r2mean_fall_ptype_evol(n_ptype,n_buff_evol), ql_fall_ptype_evol(n_ptype,n_buff_evol)
  double precision              :: n_actvd_evol(n_buff_evol), n_deactvd_evol(n_buff_evol)

  !profiles of E(x) , and E(x^2) over time : 
  !store sum and sum of squares; then compute mean and variance later
  double precision,allocatable :: prof_Tsum(:),prof_T2sum(:), prof_qvsum(:), prof_qv2sum(:)
  double precision,allocatable :: prof_sssum(:), prof_ss2sum(:)
  
  !profile of prtcl properties
  double precision,allocatable :: prof_prtcl_nsum(:), prof_prtcl_rsum(:), prof_prtcl_r2sum(:), prof_prtcl_r3sum(:)
  double precision,allocatable :: prof_cdplt_nsum(:), prof_cdplt_rsum(:), prof_cdplt_r2sum(:), prof_cdplt_r3sum(:)
  double precision,allocatable :: prof_hdplt_nsum(:), prof_hdplt_rsum(:), prof_hdplt_r2sum(:), prof_hdplt_r3sum(:)
  double precision,allocatable :: prof_inj_nsum(:), prof_actvd_nsum(:), prof_deactvd_nsum(:)

  !buffer indices to buffer data and then write it to a file
  integer                      :: bidx_inj, bidx_prtcl_evol, bidx_pdf, bidx_prof
  
  !file handles
  integer                      :: fhdl_pdf_bT, fhdl_pdf_bqv, fhdl_pdf_bss
  integer                      :: fhdl_pdf_r_prtcl, fhdl_pdf_r_aero, fhdl_pdf_r_hdplt, fhdl_pdf_r_cdplt, fhdl_pdf_r_instru
  integer                      :: fhdl_pdf_dT, fhdl_pdf_dqv, fhdl_pdf_dss, fhdl_pdf_time
  integer                      :: fhdl_prof_Tsum, fhdl_prof_T2sum, fhdl_prof_qvsum, fhdl_prof_qv2sum
  integer                      :: fhdl_prof_sssum, fhdl_prof_ss2sum, fhdl_prof_time
  integer                      :: fhdl_profprtcl_nsum, fhdl_profprtcl_rsum, fhdl_profprtcl_r2sum, fhdl_profprtcl_r3sum
  integer                      :: fhdl_profcdplt_nsum, fhdl_profcdplt_rsum, fhdl_profcdplt_r2sum, fhdl_profcdplt_r3sum
  integer                      :: fhdl_profhdplt_nsum, fhdl_profhdplt_rsum, fhdl_profhdplt_r2sum, fhdl_profhdplt_r3sum
  integer                      :: fhdl_profinj_nsum, fhdl_profactvd_nsum, fhdl_profdeactvd_nsum 
  integer                      :: fhdl_inject
  integer                      :: fhdl_evol_aero, fhdl_evol_hdplt, fhdl_evol_cdplt, fhdl_evol_instru, fhdl_evol_prtcl  
  integer                      :: fhdl_evol_mphy, fhdl_evol_fallout
  integer                      :: fhdl_log

  !flag to denote if header is written
  integer                      :: is_fhdr_evol_prtcl, is_write_fhdr_pdf, is_write_fhdr_prof
  
  namelist /domain_params/  area_frac, width_dom
  namelist /microphy_params/max_prtcl_dom, max_prtcl_gcell, &
                            add_prtcl_mthd, nbins_prtcl_dis, inject_stime, inject_etime,  inject_onetime,&
                            n_prtcl_uvol_onetime, n_prtcl_uvol_utime,&
                            inject_n_tdiff, is_inj_blob_rand , inject_at,  inject_vol,  & 
                            m0_aero, r0_prtcl, r_min_instru, lwc0_prtcl_sum, dis_prtcl_file,& 
                            aero_type, has_curv_sol_eff,term_vel_factor,& 
                            do_rec_pdf_n_prof, rec_stime_pdf_n_prof, rec_etime_pdf_n_prof,&
                            max_prtcl_hist, do_rec_prtcl_hist, rec_stime_prtcl_hist, rec_etime_prtcl_hist, rec_intvl_prtcl_hist, &
                            write_intvl_prtcl_hist
                             
  !variables and functions made public for ODT  
  public   init_microphysics, add_prtcl_to_dom, move_prtcl_by_triplet_map, prtcl_growth, prtcl_settling, &
           inject_n_tdiff, inject_stime, inject_etime, rec_stime_prtcl_hist, rec_etime_prtcl_hist, &
           compute_microphy_stats, final_step_microphy, n_prtcl_curr, n_prtcl_alltime,&
           is_injected_onetime, inject_onetime, debug_gbl, tscale_nu

  !variables and functions made public for unit testing         
  public   indiv_prtcl_term_vel, indiv_prtcl_growth, sat_mixing_ratio, m0_aero, inv_vol_gcell,&
           create_bedges_r_prtcl, create_bedges_xpos, bin_data, verify_prtcl_count, print_nprtcl_in_gcell, find_r_crit,&
           find_r_crit_dgm, rho_aero, nions_aero, molar_mass_aero, vol_gcell

  contains 

  !initialize microphysics based on ODT parameters, and boundary conditions
  !allocate arrays, initialize variables to zero, and open files
  subroutine init_microphysics(idir,nml_fname,odir,ename,H, N, tmax, td, T_o, Tdif, press, fhdl_mp_log,realz_str, debug_mode)
    use const, only               : dyn_vis                    
    implicit none

    integer, intent(in)           :: N
    integer                       :: fhdl_mp_log 
    character(len=:), allocatable, intent(in) :: idir, nml_fname, odir, realz_str, ename
    double precision, intent(in)  :: T_o, Tdif, press, tmax, td, H

    !local variable
    double precision              :: es_bot, es_top, rshift
    integer                       :: debug_mode
    character(len=:), allocatable :: cmd_str, rshift_str, prtcl_pfile
    integer                       :: iostat, is_sim_contd,fhdl_prtcl_prop

    !debug_mode
    debug_gbl       = debug_mode
    is_sim_contd    = 0

    !namelist dir and filename
    call getcwd(parent_dir)
    input_dir       = idir
    namelist_file   = nml_fname
    output_dir      = odir
    exp_name        = ename
    prtcl_pfile     = input_dir//"prtcl_prop_curr.txt"
    
    !When there is an error particle state is dumped before
    !future change: maybe the domain state can also be dumped
    !if prtcl_prop_curr.txt is present then particle state is initialized from this file
    !The simulation can continue from that state 
    fhdl_prtcl_prop = 901
    open(fhdl_prtcl_prop,file=prtcl_pfile,iostat=iostat,status="old")
    if (iostat == 0 )  is_sim_contd  = 1

    if (fhdl_mp_log > 0) then
      fhdl_log = fhdl_mp_log
    else
      fhdl_log    = 902
      open(fhdl_log, file=output_dir//"/microphysics.log")
    endif

    call read_namelist();
    call print_namelist();
     
    !domain length
    len_dom        = H
    !domain volume
    vol_dom        = area_frac * width_dom**2 * len_dom;  
    
    write(fhdl_log,'(A8,A48, F16.4)') "init:", "len_dom    :", len_dom
    write(fhdl_log,'(A8,A48, F16.4)') "init:", "width_dom  :", width_dom
    write(fhdl_log,'(A8,A48, F16.4)') "init:", "area_frac  :", area_frac
    write(fhdl_log,'(A8,A48, E16.4)') "init:", "vol_dom    :", vol_dom
    
    !grid cell info
    n_gcell       = N
    len_gcell     = (1.0d0*len_dom)/n_gcell
    vol_gcell     = len_gcell * area_frac * width_dom**2
    inv_vol_gcell = 1.0d0/vol_gcell
    inv_vol_dom   = 1.0d0/vol_dom
    gcell_1_idx   = 1
    gcell_n_idx   = n_gcell      
    
    !temperature and water vapor
    T_diff        = Tdif
    T_mid         = T_o 
    T_bot         = T_o + Tdif/2.0d0
    T_top         = T_o - Tdif/2.0d0
    es_bot        = 6.112d0 * dexp(17.67d0*(T_bot)/(T_bot+243.5)) * 100.d0
    es_top        = 6.112d0 * dexp(17.67d0*(T_top)/(T_top+243.5)) * 100.d0
    qvs_bot       = 0.622d0 * es_bot/(press - es_bot)
    qvs_top       = 0.622d0 * es_top/(press - es_top)
    qvs_diff      = qvs_bot - qvs_top
    write(fhdl_log,'(A8, A48, F16.4)') "init:", "pressure (pa)  :", press
    write(fhdl_log,'(A8, A48, F16.4)') "init:", "T_bot (K)      :", T_bot + 273.15
    write(fhdl_log,'(A8, A48, F16.4)') "init:", "T_top (K)      :", T_top + 273.15
    write(fhdl_log,'(A8, A48, F16.4)') "init:", "qvs_bot (g/kg) :", qvs_bot * 1e3
    write(fhdl_log,'(A8, A48, F16.4)') "init:", "qvs_top (g/kg) :", qvs_top * 1e3

    !velocity diffusion time scale     
    tscale_nu     = len_dom**2/dyn_vis
    write(fhdl_log,*)""
    write(fhdl_log,'(A8, A48, F16.4)') "init:", "tscale_nu:", tscale_nu

    !time steps
    sim_time_ndim = tmax
    sim_time      = tmax * tscale_nu 
    write(fhdl_log,'(A8,A48,E16.4)') "init:", "simulation non-dim time:", sim_time_ndim
    write(fhdl_log,'(A8,A48,E16.4)') "init:", "simulation time        :", sim_time

    tintvl_diff_ndim = td   
    tintvl_diff      = td   * tscale_nu

    write(fhdl_log,'(A8,A48,I16)') "init:", "Particle records array size:", max_prtcl_hist
    
    !random seed
    !rshift            =  1e10 
    !use realization number padded with 0, repeated four times
    ! if realization number = 001, then rshift = 001001001001
    realz             = realz_str
    rshift_str        = repeat(realz_str,4)
    read(rshift_str,*) rshift
    write(fhdl_log,*) ""
    write(fhdl_log,'(A8,A48,F16.0)') "init:", "random seed shifted by  :", rshift
    rand_seed1        =  rshift + 1e10 
    rand_seed2        =  rshift + 2e10 
    rand_seed3        =  rshift + 3e10 

    call init_variables();
    
    call open_files();
    !Open files before creating bins, since bin information are written to files
    call create_pdf_bins()
    call aero_prop();
    
    call allocate_variables();
    !call allocate variables before creating  prtcl dist Since distribution variables need to be allocated
    call create_prtcl_dist();
    
    write(fhdl_log, '(A8,A48,E16.4)')"init:", " Min instrument radius :" , r_min_instru
  
    !copy config files from input to output dir
    cmd_str =  "cp "//input_dir//"/"//namelist_file//" "//output_dir//"/"
    call system(cmd_str)
    cmd_str =  "cp "//input_dir//"/LabExppar.dat"//" "//output_dir//"/"
    call system(cmd_str)

    if( is_sim_contd == 1)then
      call set_curr_prtcl_prop(fhdl_prtcl_prop)
      !inject_stime  = 0
    endif  

    flush(fhdl_log)
  end subroutine

  !read the namelist for microphysical processes
  subroutine read_namelist()
    implicit none

    character(len=128) :: ifile, msg
    integer     :: ios
    
    if (len_trim(namelist_file) == 0)then
      write(fhdl_log,'(A8,A)') "init:", "Error: namelist file unknown"
      stop
    end if  
      
    ifile     = input_dir // "/"// namelist_file 
    write(fhdl_log,'(A8,A48,A)') "init:","namelist file:" , ifile
    open(201, file=ifile, form="FORMATTED", iostat=ios, iomsg=msg);
    if (ios < 0 ) then
      write(fhdl_log,'(A8,A)') "init:","Error: while opening namelist file"
      write(fhdl_log,'(A8,A,A)') "init:","Error msg:", msg
      flush(fhdl_log)
      stop   
    end if  
    
    read(201,NML=domain_params);
    read(201,NML=microphy_params);

    
    close(201)

  end subroutine
  
  !print the namelist for microphysical processes
  subroutine print_namelist()

    write(fhdl_log,'(A8,A)') "init:", "------------ namelist values ----"    
    if (len_trim(namelist_file) == 0)then
      write(fhdl_log,'(A8,A)') "init:", "Error: namelist filename variable is empty"
      flush(fhdl_log)
      stop

    else
      write(fhdl_log,  nml=domain_params   )
      write(fhdl_log,  nml=microphy_params )

    end if
    flush(fhdl_log)

  end subroutine

  !initialize scalar variables
  subroutine init_variables()
    !Bookkeeping - for particle currently in the domain
    is_injected_onetime     = 0
    n_inj_alltime           = 0
    prtcl_counter           = 0
    n_prtcl_curr            = 0
    n_prtcl_alltime         = 0
    n_fallout_curr          = 0
    n_fallout_alltime       = 0
    n_inj_curr              = 0
    n_actvd_curr            = 0
    n_deactvd_curr          = 0
    !n_prtcl_alltime = n_inj_alltime  = prtcl_counter

    !arrays for book keeping of particle injection
    inject_time             = 0.0d0
    inject_prtcl_count      = 0
    
    !set buffer index to 0
    bidx_inj                = 0 
    bidx_prtcl_evol         = 0
    bidx_pdf                = 0

    !init buffer array for evol of particle prop 
    time_prtcl_evol         = 0.0d0 
    deltat_prtcl_evol       = 0.0d0
    n_inj_evol              = 0
    n_prtcl_alltime_evol    = 0
    n_fallout_alltime_evol  = 0
    n_actvd_evol            = 0
    n_deactvd_evol          = 0
    n_fall_ptype_evol       = 0
    rmean_fall_ptype_evol   = 0.0d0
    r2mean_fall_ptype_evol  = 0.0d0
    ql_fall_ptype_evol      = 0.0d0

    n_ptype_evol            = 0
    rmin_ptype_evol         = 0.0d0
    rmax_ptype_evol         = 0.0d0 
    rmean_ptype_evol        = 0.0d0
    r2mean_ptype_evol       = 0.0d0
    ql_ptype_evol           = 0.0d0 
    
    !aerosols 
    molar_mass_aero         = 0.0d0
    nions_aero              = 0

    !particle history
    bidx_rec_hist           = 0
    time_next_phist_rec     = rec_stime_prtcl_hist
    time_next_phist_write   = rec_stime_prtcl_hist +  write_intvl_prtcl_hist

  end subroutine 

  !allocate arrays and initialize to 0
  subroutine allocate_variables()
    implicit none

    !current particles to grid cell mapping
    !There are two indices for tracking particles 
    !one index for current particles in the domain and other index for all particles injected so far
    !some would have fallenout
    ! colum "0": tell number of particles  in a grid cell, n
    ! column "1 to n" has the current id of the particle
    ! number of rows = n_gcell, which is the number of grid cells 
    allocate(gcell_prtcl_currid_list(0:max_prtcl_gcell,n_gcell));
    gcell_prtcl_currid_list         = 0;

    !This variable has count of particles in a grid cell
    !same as col "0" of above variable
    allocate(gcell_prtcl_count(n_gcell));
    gcell_prtcl_count              = 0;

    !Following variables hold the injected particle size distribution
    !which could be dry aerosol, haze particle, or cloud particle
    allocate(dis_r_prtcl(nbins_prtcl_dis))
    allocate(dis_m_aero(nbins_prtcl_dis))
    allocate(dis_r_aero(nbins_prtcl_dis))
    allocate(dis_r_crit(nbins_prtcl_dis))
    allocate(dis_n_bvol_utime(nbins_prtcl_dis))
    allocate(dis_n_uvol_utime(nbins_prtcl_dis))
    allocate(dis_n_residue(nbins_prtcl_dis));
    allocate(dis_n_alltime(nbins_prtcl_dis));
    
    !These don't take lot of space 
    !so I am not checking for memory allocation errors
    dis_r_prtcl           = 0.0d0
    dis_m_aero            = 0.0d0
    dis_r_aero            = 0.0d0
    dis_r_crit            = 0.0d0
    dis_n_bvol_utime      = 0.0d0
    dis_n_uvol_utime      = 0.0d0
    dis_n_residue         = 0.0d0
    dis_n_alltime         = 0.0d0

    !allocate pdf buffer variables
    allocate(pdf_dom_T(n_bins_T))
    allocate(pdf_bulk_T(n_bins_T))
    allocate(pdf_dom_qv(n_bins_qv))
    allocate(pdf_bulk_qv(n_bins_qv))
    allocate(pdf_dom_ss(n_bins_ss))
    allocate(pdf_bulk_ss(n_bins_ss))

    allocate(pdf_r_prtcl(n_bins_r_prtcl))
    allocate(pdf_r_aero(n_bins_r_prtcl))
    allocate(pdf_r_hdplt(n_bins_r_prtcl))
    allocate(pdf_r_cdplt(n_bins_r_prtcl))
    allocate(pdf_r_instru(n_bins_r_prtcl))

    !These don't take lot of space 
    !So I am not checking for memory allocation errors
    pdf_dom_T             = 0.0d0
    pdf_bulk_T            = 0.0d0
    pdf_dom_qv            = 0.0d0
    pdf_bulk_qv           = 0.0d0
    pdf_dom_ss            = 0.0d0
    pdf_bulk_ss           = 0.0d0

    pdf_r_prtcl           = 0.0d0
    pdf_r_aero            = 0.0d0
    pdf_r_hdplt           = 0.0d0
    pdf_r_cdplt           = 0.0d0
    pdf_r_instru          = 0.0d0

    !allocate profile buffer variables
    allocate(prof_Tsum(n_gcell))
    allocate(prof_T2sum(n_gcell))
    allocate(prof_qvsum(n_gcell))
    allocate(prof_qv2sum(n_gcell))
    allocate(prof_sssum(n_gcell))
    allocate(prof_ss2sum(n_gcell))

    allocate(prof_prtcl_nsum(n_bins_pprof))
    allocate(prof_prtcl_rsum(n_bins_pprof))
    allocate(prof_prtcl_r2sum(n_bins_pprof))
    allocate(prof_prtcl_r3sum(n_bins_pprof))
    allocate(prof_cdplt_nsum(n_bins_pprof))
    allocate(prof_cdplt_rsum(n_bins_pprof))
    allocate(prof_cdplt_r2sum(n_bins_pprof))
    allocate(prof_cdplt_r3sum(n_bins_pprof))
    allocate(prof_hdplt_nsum(n_bins_pprof))
    allocate(prof_hdplt_rsum(n_bins_pprof))
    allocate(prof_hdplt_r2sum(n_bins_pprof))
    allocate(prof_hdplt_r3sum(n_bins_pprof))
    allocate(prof_inj_nsum(n_bins_pprof))
    allocate(prof_actvd_nsum(n_bins_pprof))
    allocate(prof_deactvd_nsum(n_bins_pprof))

    !These don't take lot of space 
    !so I am not checking for memory allocation errors
    prof_Tsum             = 0.0d0
    prof_T2sum            = 0.0d0
    prof_qvsum            = 0.0d0
    prof_qv2sum           = 0.0d0
    prof_sssum            = 0.0d0
    prof_ss2sum           = 0.0d0
    prof_prtcl_nsum       = 0.0d0
    prof_prtcl_rsum       = 0.0d0
    prof_prtcl_r2sum      = 0.0d0
    prof_prtcl_r3sum      = 0.0d0
    prof_cdplt_nsum       = 0.0d0
    prof_cdplt_rsum       = 0.0d0
    prof_cdplt_r2sum      = 0.0d0
    prof_cdplt_r3sum      = 0.0d0
    prof_hdplt_nsum       = 0.0d0
    prof_hdplt_rsum       = 0.0d0
    prof_hdplt_r2sum      = 0.0d0
    prof_hdplt_r3sum      = 0.0d0
    prof_inj_nsum         = 0.0d0
    prof_actvd_nsum       = 0.0d0
    prof_deactvd_nsum     = 0.0d0

    !=== particles are from index 1 to n_prtcl_curr
    !radius of each particle
    allocate( r_prtcl_curr(max_prtcl_dom) );
    !mass of solute in each particle
    allocate( m_aero_curr(max_prtcl_dom) );
    !radius of dry aerosol
    allocate( r_aero_curr(max_prtcl_dom) );
    !critical radius for cloud droplet activation
    allocate( r_crit_curr(max_prtcl_dom) );
    !grid cell location of each particle
    allocate( gcell_prtcl_curr(max_prtcl_dom) );
    !position of particle (in non-dimensional space 0 to 1)
    allocate( x_prtcl_curr(max_prtcl_dom) );
    !particle liquid water
    allocate(ql_prtcl_curr(max_prtcl_dom)); 
    !particle id 
    allocate(unqid_prtcl_curr(max_prtcl_dom));
    !Flag to verify particle movement by eddy
    allocate(did_edymov_prtcl(max_prtcl_dom));
    !actvation/deactivation or same type; 
    !1 -> activation -1 -> deactivation  0-> no change in particle type
    allocate(is_actvd_curr(max_prtcl_dom));

    !These arrays take lot of space so I am checking if they are successsfully allocated
    if (.not. (allocated(r_prtcl_curr) .and. allocated(m_aero_curr)   .and.  &
               allocated(r_aero_curr)  .and. allocated(r_crit_curr)   .and.  & 
               allocated(gcell_prtcl_curr).and. &
               allocated(x_prtcl_curr) .and. allocated(ql_prtcl_curr) .and.  & 
               allocated(unqid_prtcl_curr) .and. allocated(did_edymov_prtcl) .and. &
               allocated(is_actvd_curr)) )  then
       write(fhdl_log,'(A8,A)') "init:", "Error during memory allocation for active particles"         
       flush(fhdl_log)         
       stop
    end if

    r_prtcl_curr       = 0.0d0
    m_aero_curr        = 0.0d0
    r_aero_curr        = 0.0d0
    r_crit_curr        = 0.0d0
    gcell_prtcl_curr   = 0
    x_prtcl_curr       = 0.0d0
    ql_prtcl_curr      = 0.0d0
    unqid_prtcl_curr   = 0 
    did_edymov_prtcl   = 0
    is_actvd_curr      = 0

    !properties of particles that fallout
    allocate( r_prtcl_fallout(max_prtcl_gcell) );
    !radius of dry aerosol
    allocate( r_aero_fallout(max_prtcl_gcell) );
    !critical radius for cloud droplet activation
    allocate( r_crit_fallout(max_prtcl_gcell) );
    !particle liquid water
    allocate( ql_prtcl_fallout(max_prtcl_gcell)); 

    if (.not. (allocated(r_prtcl_fallout) .and. allocated(ql_prtcl_fallout) .and. &
               allocated(r_aero_fallout)  .and. allocated(r_crit_fallout) ) )  then
      write(fhdl_log,'(A8,A)') "init:", "Error during memory allocation for fallout particles"         
      flush(fhdl_log)         
      stop
    end if

    r_prtcl_fallout       = 0.0d0
    r_aero_fallout        = 0.0d0
    r_crit_fallout        = 0.0d0
    ql_prtcl_fallout      = 0.0d0

    !Throughout the simulation particles are removed settling
    !So the particle position, radius, supersaturation, 
    !is recorded for all particles that ever existed during the simulation
    !The data is sampled every second and stored in a following buffer variables.
    !written every 30 or 60 seconds as configured in the namelisr


    allocate(r_prtcl_hist(max_prtcl_hist))
    allocate(x_prtcl_hist(max_prtcl_hist))
    allocate(r_aero_hist(max_prtcl_hist))
    allocate(r_crit_hist(max_prtcl_hist))
    allocate(T_prtcl_hist(max_prtcl_hist))
    allocate(qv_prtcl_hist(max_prtcl_hist))
    allocate(ss_prtcl_hist(max_prtcl_hist))
    allocate(time_prtcl_hist(max_prtcl_hist))
    allocate(deltat_prtcl_hist(max_prtcl_hist))    
    allocate(uniqid_prtcl_hist(max_prtcl_hist))

    !These arrays take lot of space so I am checking if they successsfully allocated
    if (.not. (allocated(r_prtcl_hist) .and. allocated(x_prtcl_hist) .and. &
               allocated(r_aero_hist)  .and. allocated(r_crit_hist) .and. &
               allocated(T_prtcl_hist) .and. allocated(qv_prtcl_hist) .and. & 
               allocated(ss_prtcl_hist).and. allocated(time_prtcl_hist)  .and. &
               allocated(deltat_prtcl_hist) .and. allocated(uniqid_prtcl_hist) ) )  then
       write(fhdl_log,'(A8,A)') "init:", "Error during memory allocation for particle history"         
       flush(fhdl_log)         
       stop
    end if


    r_prtcl_hist         = 0.0d0 
    x_prtcl_hist         = 0.0d0 
    r_aero_hist          = 0.0d0 
    r_crit_hist          = 0.0d0 
    T_prtcl_hist         = 0.0d0 
    qv_prtcl_hist        = 0.0d0 
    ss_prtcl_hist        = 0.0d0 
    uniqid_prtcl_hist    = 0
    time_prtcl_hist      = 0.0d0
    deltat_prtcl_hist    = 0.0d0

  end subroutine
  
  !open all files at one place to keep track of it easily
  subroutine open_files()
    implicit none

    !for each file handle opened
    !please close it in the function below
      

    fhdl_evol_mphy  = 803
    fhdl_inject     = 804  



    fhdl_pdf_dT     = 811
    fhdl_pdf_bT     = 812
    fhdl_pdf_dqv    = 813
    fhdl_pdf_bqv    = 814
    fhdl_pdf_dss    = 815
    fhdl_pdf_bss    = 816
    fhdl_pdf_time   = 817

    fhdl_prof_Tsum   = 818
    fhdl_prof_T2sum  = 819
    fhdl_prof_qvsum  = 820
    fhdl_prof_qv2sum = 821
    fhdl_prof_sssum  = 822
    fhdl_prof_ss2sum = 823
    fhdl_prof_time  = 824
    

    fhdl_pdf_r_prtcl = 825
    fhdl_pdf_r_aero  = 826
    fhdl_pdf_r_hdplt = 827
    fhdl_pdf_r_cdplt = 828
    fhdl_pdf_r_instru= 829

    fhdl_evol_aero   = 830
    fhdl_evol_hdplt  = 831
    fhdl_evol_cdplt  = 832
    fhdl_evol_instru = 833
    fhdl_evol_prtcl  = 834
    fhdl_evol_fallout= 835

    fhdl_profprtcl_nsum   = 836
    fhdl_profprtcl_rsum   = 837
    fhdl_profprtcl_r2sum  = 838
    fhdl_profprtcl_r3sum  = 839
    fhdl_profcdplt_nsum   = 840
    fhdl_profcdplt_rsum   = 841
    fhdl_profcdplt_r2sum  = 842
    fhdl_profcdplt_r3sum  = 843
    fhdl_profhdplt_nsum   = 844
    fhdl_profhdplt_rsum   = 845
    fhdl_profhdplt_r2sum  = 846
    fhdl_profhdplt_r3sum  = 847
    fhdl_profinj_nsum     = 848
    fhdl_profactvd_nsum   = 849
    fhdl_profdeactvd_nsum = 850


    open(fhdl_evol_mphy,  file =output_dir//"/evol_mphy_prop.txt")
    open(fhdl_inject,      file =output_dir//"/injection_prtcl.txt")
    
    open(fhdl_pdf_dT,    file =output_dir//"/pdf_dom_T.txt")
    open(fhdl_pdf_bT,    file =output_dir//"/pdf_bulk_T.txt")
    open(fhdl_pdf_dqv,   file =output_dir//"/pdf_dom_qv.txt")
    open(fhdl_pdf_bqv,   file =output_dir//"/pdf_bulk_qv.txt")
    open(fhdl_pdf_dss,   file =output_dir//"/pdf_dom_ss.txt")
    open(fhdl_pdf_bss,   file =output_dir//"/pdf_bulk_ss.txt")
    open(fhdl_pdf_time,  file =output_dir//"/pdf_time.txt")


    open(fhdl_prof_Tsum,  file =output_dir//"/prof_Tsum.txt")
    open(fhdl_prof_T2sum, file =output_dir//"/prof_T2sum.txt")
    open(fhdl_prof_qvsum, file =output_dir//"/prof_qvsum.txt")
    open(fhdl_prof_qv2sum,file =output_dir//"/prof_qv2sum.txt")
    open(fhdl_prof_sssum, file =output_dir//"/prof_sssum.txt")
    open(fhdl_prof_ss2sum,file =output_dir//"/prof_ss2sum.txt")
    open(fhdl_prof_time, file =output_dir//"/prof_time.txt")

    open(fhdl_pdf_r_prtcl, file =output_dir//"/pdf_radius_prtcl.txt")
    open(fhdl_pdf_r_aero,  file =output_dir//"/pdf_radius_aero.txt")
    open(fhdl_pdf_r_hdplt, file =output_dir//"/pdf_radius_hdplt.txt")
    open(fhdl_pdf_r_cdplt, file =output_dir//"/pdf_radius_cdplt.txt")
    open(fhdl_pdf_r_instru,file =output_dir//"/pdf_radius_instru.txt")

    open(fhdl_evol_aero,  file =output_dir//"/evol_aero_prop.txt")
    open(fhdl_evol_hdplt,  file =output_dir//"/evol_hdplt_prop.txt")
    open(fhdl_evol_cdplt,  file =output_dir//"/evol_cdplt_prop.txt")
    open(fhdl_evol_instru,  file =output_dir//"/evol_instru_prop.txt")
    open(fhdl_evol_prtcl,  file =output_dir//"/evol_prtcl_prop.txt")
    open(fhdl_evol_fallout,  file =output_dir//"/evol_fallout_prop.txt")

    open(fhdl_profprtcl_nsum, file =output_dir//"/prof_prtcl_nsum.txt")
    open(fhdl_profprtcl_rsum, file =output_dir//"/prof_prtcl_rsum.txt")
    open(fhdl_profprtcl_r2sum, file =output_dir//"/prof_prtcl_r2sum.txt")
    open(fhdl_profprtcl_r3sum, file =output_dir//"/prof_prtcl_r3sum.txt")
    open(fhdl_profcdplt_nsum, file =output_dir//"/prof_cdplt_nsum.txt")
    open(fhdl_profcdplt_rsum, file =output_dir//"/prof_cdplt_rsum.txt")
    open(fhdl_profcdplt_r2sum, file =output_dir//"/prof_cdplt_r2sum.txt")
    open(fhdl_profcdplt_r3sum, file =output_dir//"/prof_cdplt_r3sum.txt")
    open(fhdl_profhdplt_nsum, file =output_dir//"/prof_hdplt_nsum.txt")
    open(fhdl_profhdplt_rsum, file =output_dir//"/prof_hdplt_rsum.txt")
    open(fhdl_profhdplt_r2sum, file =output_dir//"/prof_hdplt_r2sum.txt")
    open(fhdl_profhdplt_r3sum, file =output_dir//"/prof_hdplt_r3sum.txt")
    open(fhdl_profinj_nsum, file =output_dir//"/prof_inj_nsum.txt")
    open(fhdl_profactvd_nsum, file =output_dir//"/prof_actvd_nsum.txt")
    open(fhdl_profdeactvd_nsum, file =output_dir//"/prof_deactvd_nsum.txt")

    !flag for writing file headers
    is_fhdr_evol_prtcl   = 1
    is_write_fhdr_pdf    = 1
    is_write_fhdr_prof   = 1

  end subroutine  
  
  !Close all files at one place to keep track of it easily
  subroutine close_files()



    close(fhdl_evol_mphy)
    close(fhdl_inject)

    close(fhdl_pdf_dT)
    close(fhdl_pdf_bT)
    close(fhdl_pdf_dqv)
    close(fhdl_pdf_bqv)
    close(fhdl_pdf_dss)
    close(fhdl_pdf_bss)
    close(fhdl_pdf_time)

    close(fhdl_prof_time)
    close(fhdl_prof_Tsum)
    close(fhdl_prof_T2sum)
    close(fhdl_prof_qvsum)
    close(fhdl_prof_qv2sum)
    close(fhdl_prof_sssum)
    close(fhdl_prof_ss2sum)
    
    close(fhdl_pdf_r_prtcl)
    close(fhdl_pdf_r_aero)
    close(fhdl_pdf_r_hdplt)
    close(fhdl_pdf_r_cdplt)
    close(fhdl_pdf_r_instru)
    
    close(fhdl_evol_aero)
    close(fhdl_evol_hdplt)
    close(fhdl_evol_cdplt)
    close(fhdl_evol_instru)
    close(fhdl_evol_prtcl)
    close(fhdl_evol_fallout)
       
    close(fhdl_profprtcl_nsum)
    close(fhdl_profprtcl_rsum)
    close(fhdl_profprtcl_r2sum)
    close(fhdl_profprtcl_r3sum)
    close(fhdl_profcdplt_nsum)
    close(fhdl_profcdplt_rsum)
    close(fhdl_profcdplt_r2sum)
    close(fhdl_profcdplt_r3sum)
    close(fhdl_profhdplt_nsum)
    close(fhdl_profhdplt_rsum)
    close(fhdl_profhdplt_r2sum)
    close(fhdl_profhdplt_r3sum)
    close(fhdl_profinj_nsum)
    close(fhdl_profactvd_nsum)
    close(fhdl_profdeactvd_nsum)

  end subroutine

  !create particle size distribution based on input parameters from namelist
  !uses subprograms - set_monodisp_prtcl_dist, set_gamma_prtcl_dist, read_prtcl_dist_file
  subroutine create_prtcl_dist()
    !This subprogram creates particle size distribution using one of three methods
    !1) Add monodisperse sized particles  
    !   input: mthd = 1
    !   input: n_prtcl_uvol_utime -> total number conc per unit vol    
    !   input: r0_prtcle  -> radius is same for all particles
    !   input: m0_aero  -> solute mass is same for all particles  

    !2) Add particles based on gamma particle distribution with bin size of 0.5 micrometer
    !   input: mthd = 2
    !   input: n_prtcl_add_uvol -> total number conc per unit vol        - dist parameter   
    !   input: lwc0_prtcl_sum -> total liquid water content of all particles - dist parameter
    !   bin size = 0.5 micrometer, number of bins  = nbins_dis
    
    !3) Add particles based on size distribution from file 
    !   input: mthd = 3  
    !   input: dis_file - filename in the input_dir that has number per unit volume, solute mass for each size bine
    !   bin size = 0.5 micrometer, number of bins  = nbins_dis
    
    use const, only : pi43
    implicit none
    !local variables
    double precision      :: vol_blob, n_prtcl_uv_ut, r_crit, ss_crit
    integer               :: j, bin_count
    character(len=:), allocatable :: ifile

    vol_blob              = vol_dom * inject_vol 
    write(fhdl_log,*) ''
    write(fhdl_log,'(A8,A48, E16.4)') "init:","domain vol in (m^3)                :", vol_dom
    write(fhdl_log,'(A8,A48, E16.5)') "init:",'injection blob as domain vol frac  :', inject_vol
    write(fhdl_log,'(A8,A48, E16.5)') "init:",'injection blob vol   (m^3)         :', vol_blob
    write(fhdl_log,'(A8,A48, E16.5)') "init:",'inj blob vol center at (non-dim)   :', inject_at
    dis_r_prtcl           = 0.0d0
    dis_m_aero            = 0.0d0
    dis_n_bvol_utime      = 0.0d0
    dis_n_uvol_utime      = 0.0d0

    if (inject_onetime == 1)then
      n_prtcl_uv_ut        = n_prtcl_uvol_onetime
    else
      n_prtcl_uv_ut        = n_prtcl_uvol_utime
    end if

    if (add_prtcl_mthd     == 1) then
      call set_monodisp_prtcl_dist(n_prtcl_uv_ut,vol_blob,r0_prtcl,m0_aero,dis_r_prtcl,dis_m_aero,&
                                  dis_n_bvol_utime, dis_n_uvol_utime)
      !Change number of bins for the monodisperse distribution as 1
      nbins_prtcl_dis      = 1
    else if(add_prtcl_mthd == 2) then  
      call set_gamma_prtcl_dist(n_prtcl_uv_ut,vol_blob,lwc0_prtcl_sum,m0_aero,nbins_prtcl_dis,&
                                dis_r_prtcl,dis_m_aero,dis_n_bvol_utime, dis_n_uvol_utime)
    
    else if(add_prtcl_mthd == 3) then
      ifile   = input_dir // trim(dis_prtcl_file)
      call read_prtcl_dist_file(ifile, n_prtcl_uv_ut, vol_blob, nbins_prtcl_dis,dis_r_prtcl,dis_m_aero,&
                               dis_n_bvol_utime, dis_n_uvol_utime)
    
    else
      write(fhdl_log,'(A8,A)') "init:", "Error: Incorrect particle size distribution method " 
      stop

    endif
    
    write(fhdl_log,*) ''
    write(fhdl_log,'(A8,A)') "init:",'---- particle distribution in the injection volume ----'
    write(fhdl_log,'(A8,A)') "init:",'aerosol_radius    crit_radius   particle_radius    aerosol_mass     particle_concentration '  
    write(fhdl_log,'(A8,A)') "init:",'  (m)                 (m)           (m)               (kg)            (# in inj vol/sec)' 
    do j = 1, nbins_prtcl_dis
      bin_count           = nint(dis_n_bvol_utime(j))
      dis_r_aero(j)       = dis_m_aero(j)/rho_aero
      dis_r_aero(j)       = (dis_r_aero(j)/pi43)**(1.d0/3.d0)   
      r_crit              = 0.0d0
      ss_crit             = 0.0d0
      call find_r_crit(T_mid+273.15,dis_m_aero(j),r_crit,ss_crit);
      dis_r_crit(j)       =  r_crit
      !dis_r_crit(j)       =  r_min_instru
      
      if ( dis_n_bvol_utime(j) > 0.0d0) then
        write(fhdl_log,'(A8,5(e12.5,4x))') "init:", dis_r_aero(j), dis_r_crit(j), dis_r_prtcl(j), dis_m_aero(j), dis_n_bvol_utime(j)    
      end if

      if ( dis_r_prtcl(j) < dis_r_aero(j) )then
        write(fhdl_log, '(A8,A)') "init:", "Error: particle radius smaller than dry aerosol radius"
        flush(fhdl_log)
        stop
      endif

    end do
    write(fhdl_log,*) ''

    flush(fhdl_log)
  end subroutine

  !create monodisperse particle size distribution
  subroutine set_monodisp_prtcl_dist(nprtcl_uv_ut,bvol,r0_add, m0_add,dis_r,dis_m,dis_n_bv_ut,dis_n_uv_ut)
    !set first bin with monodisperse particle properties
    !==== input
    !nprtcl_uv_ut -> n prtcles in unit vol and unit time - from namelist
    !bvol         -> injection blob vol,
    !r0_add       -> injected particles' radius
    !m0_add       -> aerosol mass in the prtcl, rest is liquid water
    !=== output
    !dis_r,dis_m,dis_n_bv_ut,dis_n_uv_ut  have same array size i.e n bins
    ! monodisperse prtcl size dis  have only one radius bin
    !dis_r        -> prtcl radius (bin center - average of left and right radius bin edges)
    !dis_m        -> aerosol mass
    !dis_n_bv_ut  -> inject rate expressed as n prtcl in each radius bin in injection blob vol in unit time
    !dis_n_uv_ut  -> inject rate expressed as n prtcl in each radius bin per unit vol in unit time

    implicit none
    
    !input
    double precision, intent(in)  :: r0_add, m0_add, bvol, nprtcl_uv_ut
    !output
    double precision, allocatable, intent(inout) ::dis_r(:),dis_m(:), dis_n_bv_ut(:), dis_n_uv_ut(:)

    !local variable
    double precision              :: n_bv_ut
    
    write(fhdl_log,*) ''
    write(fhdl_log,'(A8,A48)') "init:", '---- prtcl_dist_type = 1 monodisperse dis'
    write(fhdl_log,'(A8,A48, E16.4)') "init:", 'radius       r0     :', r0_add
    write(fhdl_log,'(A8,A48, E16.4)') "init:", 'aerosol mass m0     :', m0_add
    !n particles injected = (n_prtcl_uvol_utime * inj_vol * time_step)
    !n particles inj unit time = (n_prtcl_uvol_utime * inj_vol )
    !inj_vol - blob_vol
    n_bv_ut               = nprtcl_uv_ut * bvol
    write(fhdl_log,'(A8,A48, E16.4)') "init:", 'n_prtcl in unit vol unit time  :', nprtcl_uv_ut
    write(fhdl_log,'(A8,A48, F16.4)') "init:", 'n_prtcl in blob vol  unit time :', n_bv_ut
    
    dis_r(1)              = r0_add
    dis_m(1)              = m0_add
    dis_n_bv_ut(1)        = n_bv_ut
    dis_n_uv_ut(1)        = nprtcl_uv_ut

    flush(fhdl_log)
  end subroutine

  !create gamma particle size distribution
  subroutine set_gamma_prtcl_dist(nprtcl_uv_ut,bvol,lwc_add,m0_add,nbins_dis,dis_r,dis_m,dis_n_bv_ut,dis_n_uv_ut)
    !Uses subroutine size_dis" from EMPM to find the size distribution
    !Use gamma distribution to determine size distribution
    !bin size = 0.5 micrometer, number of bins  = nbins_dis
    !==== input
    !nprtcl_uv_ut -> n prtcles in unit vol and unit time - from namelist
    !bvol         -> injection blob vol,
    !lwc_add      -> total liquid water content of all particles - gamma distribution parameter
    !m0_add       -> aerosol mass in the prtcl - this imply that all particles in the dist have same aero mass
    !nbins_dis    -> n bins in the distribution from the namelist
    !=== output
    !dis_r,dis_m,dis_n_bv_ut,dis_n_uv_ut  have same array size i.e n bins
    !dis_r        -> prtcl radius (bin center - average of left and right radius bin edges)
    !dis_m        -> aerosol mass
    !dis_n_bv_ut  -> inject rate expressed as n prtcl in each radius bin in injection blob vol in unit time
    !dis_n_uv_ut  -> inject rate expressed as n prtcl in each radius bin per unit vol in unit time
    implicit none

    !input
    integer, intent(in)           :: nbins_dis
    double precision, intent(in)  :: lwc_add, bvol, nprtcl_uv_ut, m0_add
    !output
    double precision, allocatable, intent(inout) ::dis_r(:),dis_m(:), dis_n_bv_ut(:), dis_n_uv_ut(:)
    
    !local variable 
    integer                       :: j

    !external subroutine 
    EXTERNAL                      size_dis

    !calculate particle size distribution using gamma function based on total number concentration (n_prtcl_uvol_utime) and 
    !lwc0_prtcl_sum given in the namelist file
    !call size_dis(N_i,  lwc0_prtcl_sum,radius_drop,con_drop,     mass_solute,vol_dom, m_i,       ne)
    write(fhdl_log,*) ''
    write(fhdl_log,'(A8,A)') "init:", '---- prtcl_dist_type = 2 from gamma dis'
    write(fhdl_log,'(A8,A48, E16.4)') "init:", ' distribution parameters: total particles:', nprtcl_uv_ut
    write(fhdl_log,'(A8,A48, E16.4)') "init:", ' distribution parameters: total lwc     :', lwc_add
    call size_dis(nprtcl_uv_ut,lwc_add,dis_r,dis_n_uv_ut,dis_m,1.0d0,m0_add,nbins_dis)
    write(fhdl_log,'(A8,A48,I16)') "init:", "nbins_dis :", nbins_dis

    do j = 1, nbins_dis
      dis_n_bv_ut(j)    = dis_n_uv_ut(j) * bvol
    end do

    flush(fhdl_log)
  end subroutine

  !read size distribution from the file
  subroutine read_prtcl_dist_file(dis_file,n_prtcl_uv_ut, bvol, nbins_dis, dis_r,dis_m,dis_n_bv_ut,dis_n_uv_ut)
    !read the particle size distribution from the file and stored in the array
    !==== input
    !dis_file     -> file with particle size distribution 
    !nprtcl_uv_ut -> n prtcles in unit vol and unit time - from namelist
    !bvol         -> injection blob vol,
    !nbins_dis     -> number of bins from the namelist
    !=== output
    !dis_r,dis_m,dis_n_bv_ut,dis_n_uv_ut  have same array size i.e n bins
    !dis_r        -> prtcl radius (bin center - average of left and right radius bin edges)
    !dis_m        -> aerosol mass
    !dis_n_bv_ut  -> inject rate expressed as n prtcl in each radius bin in injection blob vol in unit time
    !dis_n_uv_ut  -> inject rate expressed as n prtcl in each radius bin per unit vol in unit time
    implicit none

    !input
    character(len=:), allocatable, intent(in)    ::  dis_file
    double precision, intent(in)                 :: bvol, n_prtcl_uv_ut
    integer, intent(in)                          :: nbins_dis
    !output
    double precision, allocatable, intent(inout) ::dis_r(:),dis_m(:), dis_n_bv_ut(:), dis_n_uv_ut(:)

    ! -- read particle data from file (the output of Dr. Austin's DGM code)
    
    !local variables
    integer                       :: j
    double precision              :: sum_prtcl_dis

    write(fhdl_log,*) ''
    write(fhdl_log,'(A8,A48)') "init:",'---- prtcl_dist_type = 3, dis from a file'
    write(fhdl_log,'(A8,A48,A)') "init:",'reading particle data from file:', dis_file
    
    !read particle size distribution from the file
    open (201,  FILE=dis_file)
    do j = 1, nbins_dis
      read(201,*) dis_r(j)
    end do
    sum_prtcl_dis      = 0.0d0
    do j = 1, nbins_dis
      read(201,*) dis_m(j),  dis_n_uv_ut(j)
      sum_prtcl_dis    = sum_prtcl_dis + dis_n_uv_ut(j)
    end do
    close(201)

    !update particle size distribution from file
    !call accu(vol_dom,           dis,       con_drop,     ndrop,            ne,           ne1)
    !call accu(bvol,dis_nprtcl,dis_n_uv_ut,n_prtcl_add,nbins_dis,nbins_dis)
    do j = 1, nbins_dis
      if ( dis_n_uv_ut(j) > 0.0d0 ) then
        !convert total freq of distribution equal to injection rate
        dis_n_uv_ut(j)  = dis_n_uv_ut(j)/sum_prtcl_dis * n_prtcl_uv_ut
        dis_n_bv_ut(j)  = dis_n_uv_ut(j) * bvol
      end if
    end do

    flush(fhdl_log) 
    
  end subroutine
  
  !set aerosol properties based on aerosol type chosen in the namelist
  subroutine aero_prop()
    implicit none

    write(fhdl_log,*) ""
    if (aero_type .eq. 1) then
      !NaCl
      molar_mass_aero = 58.4428e-3 ! kg per mole 
      rho_aero        = 2.163d3 !kg/m3  
      nions_aero      = 2

    else if (aero_type .eq. 2) then
      ! Ammonium sulphate (NH4)2SO4
      molar_mass_aero = 132.1395e-3  ! kg per mole
      rho_aero        = 1.770d3 !kg/ m3 
      nions_aero      = 3

    else if (aero_type .eq. 3) then
      ! Ammonium bisulfate (NH4)HSO4
      molar_mass_aero = 115.11e-3 ! kg per mole
      rho_aero        = 1.780d3   !kg/ m3 
      nions_aero      = 2
      
    else
      write(fhdl_log,'(A8,A)') "init:", 'Error: aero value must be either 1 or 2'
      flush(fhdl_log)
      stop 

    end if

    flush(fhdl_log)
  end subroutine

  !find critical radius and supersaturation for a given dry aerosol mass
  !based on kohler theory in Roger and Yau book
  !critical radius is typically used to classify if a particle is a cloud droplet or a haze droplet
  !if r_prtcl < r_crit then it is haze droplet else cloud droplet
  subroutine find_r_crit(T_env,m_aero,r_crit,ss_crit)
   
   ! Use Rogers & Yau Eqs. 6.7 and 6.8 

   !input 
   !T_env   -> temperature
   !m_aero  -> mass of aerosol
   !output 
   !r_crit  -> critical radius  
   !ss_crit -> critical supersaturation

   implicit none
   ! T_env: temperature (K) and  m_aero:   solute mass (kg)
    double precision,intent(in)       :: T_env, m_aero
    double precision,intent(inout)    :: ss_crit, r_crit

    !local variables 
    double precision,parameter        :: m_per_cm = 0.01
    double precision                  :: a_RY, b_RY
  
    !T_RY: temperature (K)
    a_RY      = 3.3e-5/T_env ! cm
    !a_RY      = 3.2038e-5/T_env ! cm
    a_RY      = a_RY * m_per_cm
  
    b_RY      = 4.3 * m_aero * nions_aero/ molar_mass_aero ! cm^3
    b_RY      = b_RY * m_per_cm**3.0d0
  
    ss_crit   = sqrt(4.0d0 * a_RY**3.0d0 / (27.0d0 * b_RY) )
    r_crit    = sqrt(3.0d0 * b_RY / a_RY) 
  
  end subroutine  

  !find critical radius and supersaturation based on droplet growth model used in this program
  subroutine find_r_crit_dgm(T_env,m_aero,r_prtcl,r_crit,ss_crit)
   !use a and b constants from droplet growth model in fcnkb.f90
   use const, only: Rv, pi43, Mw, rho_w 
   implicit none

   ! T_env: temperature (K) and  m_aero:   solute mass (kg)
   double precision,intent(in)       :: T_env, m_aero,r_prtcl
   double precision,intent(inout)    :: ss_crit, r_crit

   !local variables 
   double precision,parameter        :: m_per_cm = 0.01
   double precision                  :: a, b, sigma, Ms, rho_prtcl, c7am, c7nacl, c7
   sigma       = 7.392730e-2  
   Ms          = molar_mass_aero
   c7am        = 0.4363021
   c7nacl      = 0.5381062
   if (aero_type == 1) then
    c7  = c7nacl
   else
    c7  = c7am
   endif

   rho_prtcl   = (r_prtcl**3 * pi43 * rho_w + m_aero*c7)/(r_prtcl**3 * pi43)
   !a          = 2*sigma/(Rv*rho_prtcl*T_env)
   !b          = (Mw/Ms)*m_aero/(pi43 * r_prtcl**3 * rho_prtcl - m_aero) * r_prtcl**3
   a           = 2 * sigma / (Rv * rho_w * T_env)
   b           = (Mw/Ms)* m_aero *nions_aero  /(pi43 * rho_w)
   ss_crit     = sqrt(4.0d0 * a**3.0d0 / (27.0d0 * b) )
   r_crit      = sqrt(3.0d0 * b / a) 
  
  end subroutine  

  !Add particle to dom either one time or through continuous injection
  !use other subpprograms
  !get particle size based chosen size distbution, assign uniq id, 
  !add or compute particle properties, position particles randomly in the domain 
  subroutine add_prtcl_to_dom(deltat_ndim, time_ndim)

    implicit none
    double precision, intent(in)    :: deltat_ndim, time_ndim
    !local variables
    integer                         :: n_prtcl_add, debug_lcl
    double precision                :: tintvl 

    debug_lcl                       = 0
    n_prtcl_add                     = 0
    n_inj_curr                     = 0

    !call update_prtcl_gcell_from_xpos()
    !call print_nprtcl_in_gcell()
    if (debug_lcl == 1 .or. debug_gbl == 1) then
      write(fhdl_log,*) ''
      write(fhdl_log,'(A8,42x,A)') "inj:", "************ Adding particles to the domain ****"
    end if

    !tscale_nu is computed during the call to init_microphysics
    if (inject_onetime == 1) then
      tintvl        = 1
      write(fhdl_log,'(A8,42x,A)') "inj:","************ Onetime injection"
    else 
      tintvl        = deltat_ndim * tscale_nu
      if (debug_lcl == 1 .or. debug_gbl == 1) then
        write(fhdl_log,'(A8,42x,A)') "inj:","************ Continuous injection"
      endif  
    end if

    if (debug_lcl == 1 .or. debug_gbl == 1) then
      write(fhdl_log,'(A8,A48,F16.4,A)') "inj:", "time interval  (dim)     :", tintvl,   " secs"
    endif  

    !==== Start of executable unit
    !==== Treat this as a group of operations , all global variables (various prtcl and gcell properties)
    !==== will tally after all statements in this unit are executed. Breaking this unit will create inconsistent particle properties
    !==== Any verification across variables should be done after this unit else they will fail since they are not updated yet
    call add_prtcl_prop_from_dis(n_prtcl_curr, nbins_prtcl_dis, dis_r_prtcl, dis_m_aero, dis_r_aero,dis_n_bvol_utime, &
                                    max_prtcl_dom, r_prtcl_curr, m_aero_curr, r_aero_curr,n_prtcl_add, tintvl, dis_n_residue);
    
    if (debug_lcl == 1 .or. debug_gbl == 1) then                                
      write(fhdl_log,'(A8,A48,I16)') "inj:", "Number of particles added:", n_prtcl_add
    endif

    !Check if adding particles exceeds any array thresholds                                   
    if ( ( n_prtcl_add + n_prtcl_curr) > max_prtcl_dom)then
      write(fhdl_log, '(A8,A)') "inj:", "Error: # particles added exceeds # of inst max particles  (array alloc)" 
      write(fhdl_log, '(A8,A,I16)') "inj:", "Error: max particles allowed in the domain                :", max_prtcl_dom
      write(fhdl_log, '(A8,A,I16)') "inj:", "Error: number of active particles currently in the domain :", n_prtcl_curr
      write(fhdl_log, '(A8,A,I16)') "inj:", "Error: number of  particles  added                        :", n_prtcl_add
      flush(fhdl_log);
      stop
    end if

    if(n_prtcl_add > 0) then
      call assign_unqid_for_prtcl(n_prtcl_curr, n_prtcl_add, prtcl_counter, unqid_prtcl_curr);
      call calc_prtcl_ql(n_prtcl_curr+n_prtcl_add, r_prtcl_curr,r_aero_curr,inv_vol_dom, ql_prtcl_curr);
      call assign_prtcl_position(is_inj_blob_rand, inject_at, inject_vol, n_prtcl_curr, n_prtcl_add,&
                                rand_seed1, rand_seed2, x_prtcl_curr)
      
      !This variable keeps track of all particles added since the start of simulation
      n_prtcl_alltime                = prtcl_counter;

      
      !this variable keeps track of number of particles currently in the domain
      !some of the particles add to the domain maybe removed through settling
      n_prtcl_curr                   = n_prtcl_curr + n_prtcl_add;
      
      ! do j  = 1, n_prtcl_curr
      !   write(fhdl_log, '(A48,E12.4)') " prtcl radius:" , r_prtcl_curr(j)
      !   write(fhdl_log, '(A48,E12.4)') " prtcl ql    :" , ql_prtcl_curr(j)
      ! end do

      !save the injection details
      bidx_inj                       = bidx_inj+1
      inject_prtcl_count(bidx_inj)   = n_prtcl_add
      inject_time(bidx_inj)          = time_ndim
      n_inj_alltime                  = n_inj_alltime + n_prtcl_add
      n_inj_curr                     = n_prtcl_add

      call update_prtcl_gcell_from_xpos()
      call print_nprtcl_in_gcell()                                               
      !This creates 2D array ; list of particle ids at each grid cell
      !plus checks if array limit set by max_prtcl_gcell is exceeded or not  
      call update_gcell_prtcl_map() 
    endif 
    !==== end of executable unit 

    if (inject_onetime == 1) then
      is_injected_onetime = 1
    endif  

    if (debug_lcl == 1 .or. debug_gbl == 1) then
      write(fhdl_log,'(A8,A48,I16)') "inj:","Number of particles currently in the domain   :" , n_prtcl_curr
      write(fhdl_log,'(A8,A48,I16)') "inj:","Number of particles all time                  :" , n_prtcl_alltime
      flush(fhdl_log)
    end if 

  end subroutine 

  !inject particle based size distribution configured in the namelist
  subroutine add_prtcl_prop_from_dis(n_prtcl_currt, nbins_dis, dis_r, dis_m, dis_raero, dis_n_bv_ut, max_prtcl, r_prtcl, m_aero, &
                                    r_aero,n_prtcl_add, tintvl, dis_n_res)
    !Add individual particle and its properties to the domain 
    !step 1) find the number of particles to add from particle distribution in the given tintvl
    !step 2) add properties to current particles list

    !input
    implicit none                                    
    double precision, allocatable, intent(in)   :: dis_r(:),dis_m(:), dis_n_bv_ut(:), dis_raero(:)
    integer, intent(in)                         :: n_prtcl_currt, nbins_dis, max_prtcl
    double precision, intent(in)                :: tintvl    
    
    !output
    double precision, allocatable, intent(inout):: r_prtcl(:), m_aero(:), dis_n_res(:),r_aero(:)  
    integer, intent(inout)                      :: n_prtcl_add

    !local variable
    integer             :: j, k, idx,  bin_count_int, debug_lcl
    double precision    :: bin_count_float, n_prtcl_add_float

    debug_lcl              = 0
    n_prtcl_add            = 0  
    n_prtcl_add_float      = 0.0

    do  j = 1, nbins_dis
      bin_count_float      = dis_n_bv_ut(j) * tintvl + dis_n_res(j)
      bin_count_int        = int(bin_count_float)
      n_prtcl_add_float    = n_prtcl_add_float + bin_count_float
      n_prtcl_add          = n_prtcl_add + bin_count_int
    end do
    
    if (debug_lcl == 1 .or. debug_gbl == 1) then
      write(fhdl_log,'(A8,A48,I16)')   "inj:","Particles currently in the domain - n_prtcl  :", n_prtcl_currt
      write(fhdl_log,'(A8,A48,E16.4)') "inj:","Adding particles to domain - n_prtcl_add(float):", n_prtcl_add_float
      write(fhdl_log,'(A8,A48,I16)')   "inj:","Adding particles to domain - n_prtcl_add(int)  :", n_prtcl_add
      
    endif

    if ( (n_prtcl_currt + n_prtcl_add) > max_prtcl ) then
      write(fhdl_log,'(A8,A)') "inj:", "Error: added particles exceeds domain max  (array limit/memory allocation)"
      flush(fhdl_log)
      stop
    end if

    idx = n_prtcl_currt
    do j  = 1 , nbins_dis
      !convert number of particles per unit time to number per time interval
      !If a bin count is < 0.5, using int or nint will make freq zero. This results in no particles of this size 
      !added to domain. By mainitaining residue from past addition, soon they will add up and when freq > 1, 
      !that particular particle size will be added to the domain       
      bin_count_float     =  dis_n_bv_ut(j) * tintvl + dis_n_res(j)
      bin_count_int       =  int(bin_count_float)
      dis_n_res(j)        =  bin_count_float - (1.0d0 * bin_count_int)         
      dis_n_alltime(j)    =  dis_n_alltime(j) + bin_count_int

      if (bin_count_int > 0) then
        do k = 1, bin_count_int
          idx             =  idx + 1
          r_prtcl(idx)    =  dis_r(j)
          m_aero(idx)     =  dis_m(j)
          r_aero(idx)     =  dis_raero(j)
          if (debug_lcl > 0)then
            write(fhdl_log,'(A8,A48,E16.4)') "inj:", "Adding prtcl radius:", r_prtcl(idx)
          endif    
        end do
      end if
    end do

    if (debug_lcl == 1 .or. debug_gbl == 1) then
      flush(fhdl_log)
    end if  
  end subroutine

  !assign unique id to each particle
  subroutine assign_unqid_for_prtcl(n_prtcl_currt,n_prtcl_add,id_counter, id_prtcl)

      !input
      integer, intent(in)                 :: n_prtcl_currt,n_prtcl_add
      !output
      integer, allocatable, intent(inout) :: id_prtcl(:)
      integer, intent(inout)              :: id_counter
       
      !local variables
      integer           :: j, idx

      !id_counter keeps count of particles injected so far
      !n_prtcl_currt keeps count of particles currently in the domain
      idx   = n_prtcl_currt
      
      do j = 1, n_prtcl_add
        idx = idx + 1
        id_counter = id_counter + 1
        id_prtcl(idx) = id_counter  
      end do

  end subroutine

  !compute ql for each particle
  subroutine calc_prtcl_ql(n_prtcl,r_prtcl,r_aero, inv_vol, ql_prtcl)
    use const, only: pi43, rho_w

    !input
    integer, intent(in)                :: n_prtcl
    double precision, intent(in)       :: inv_vol
    double precision, allocatable, intent(in) :: r_prtcl(:),r_aero(:)
    !output
    double precision, allocatable, intent(inout) ::  ql_prtcl(:)

    !local variables
    integer           :: j
    do j  = 1, n_prtcl
      !write(fhdl_log,'(A,E12.2)') "radius:" , r_prtcl(idx)
      ql_prtcl(j) = pi43 * (r_prtcl(j)-r_aero(j))**3.0 * rho_w * inv_vol
    end do
    
  end subroutine

  !assign particle position 
  !There are two steps involved in assigning position
  !1) the center of injected blob  - this can be at a random location or at a chosen position(bpos)
  !2) particle is positioned randomly within the injected blob volume 
  !blob volume (blen) is configured in the namelist
  subroutine assign_prtcl_position( is_rnd_bpos, bpos, blen, n_prtcl_currt, n_prtcl_add, rseed1, rseed2, x_prtcl)
    !input
    double precision, intent(in)  :: bpos, blen, rseed1, rseed2 
    integer, intent(in)           :: is_rnd_bpos, n_prtcl_add, n_prtcl_currt
    
    !output
    double precision, allocatable, intent(inout)  :: x_prtcl(:)
     
    !local  variables
    double precision              :: bcenter, bstart, bend
    integer                       :: idx, j, debug_lcl
    !external function
    double precision              :: rand1
    debug_lcl             = 0

    if (is_rnd_bpos == 0) then
      bcenter             = bpos
    else
      bcenter             = (blen * 0.5d0) + rand1(rseed1) *  (1 - blen)
    end if  

    bstart                = bcenter - (blen * 0.5d0); 
    bend                  = bcenter + (blen * 0.5d0);

    !print rand1(rseed2) randseed2
    if (bstart < 0 .or. bend > 1) then
      write(fhdl_log,'(A8,A)') "inj:", "Error: particle blob outside the domain boundary"
      write(fhdl_log,'(A8,A)') "inj:", "Error: blob_center_pos:",  bpos
      write(fhdl_log,'(A8,A)') "inj:", "Error: blob_size      :",  blen
      flush(fhdl_log)
      stop
    end if

    idx     = n_prtcl_currt
    do j    = 1, n_prtcl_add 
      idx   = idx + 1  
      !use a different random seed for particle position
      x_prtcl(idx)  = bstart + (rand1(rseed2) * blen )
      !write(fhdl_log,'(A48,F12.8)') "particle position  :", x_prtcl(idx)
              !verification
      if(x_prtcl(idx) < 0 .or. x_prtcl(idx) > len_dom) then
        write(fhdl_log,'(A8,A)') "inj:","Error: Assigning pos to new particle but pos outside domain"
        write(fhdl_log,'(A8,A)') "inj:","Error: particle curr idx    :", idx
        write(fhdl_log,'(A8,A)') "inj:","Error: particle new position:", x_prtcl(idx)
        flush(fhdl_log)
        stop
      end if

    end do

    !call the program to calc gcell from parent subroutine
    !call update_prtcl_gcell_from_xpos()
    if (n_prtcl_add > 0)then
      if ( debug_lcl == 1 .or. debug_gbl == 1 .or. inject_onetime == 1) then
        write(fhdl_log,'(A8,42x,A)') "inj:","************ Assigning particle position: ****"  
        write(fhdl_log,'(A8,A48,F16.2)') "inj:","bstart         :", bstart
        write(fhdl_log,'(A8,A48,F16.2)') "inj:","bend           :", bend
        write(fhdl_log,'(A8,A48,E16.8)') "inj:","Min particle pos:", minval(x_prtcl(n_prtcl_currt+1:n_prtcl_currt+n_prtcl_add))
        write(fhdl_log,'(A8,A48,E16.8)') "inj:","Max particle pos:", maxval(x_prtcl(n_prtcl_currt+1:n_prtcl_currt+n_prtcl_add))
        flush(fhdl_log)
      endif
    end if

  end subroutine
  
  !update the grid cell index for each particle based on its current position
  !current pos is modified when a particle is first injected, moved by eddy, or through gravitional settling
  subroutine update_prtcl_gcell_from_xpos()  
    implicit none
    integer   :: j, gcell_idx
    
    do j = 1, n_prtcl_curr
      gcell_idx  = int(x_prtcl_curr(j)/len_gcell) + 1
      if (gcell_idx < 1 .or. gcell_idx > n_gcell) then
        write(fhdl_log,'(A)')       "Error: gcell idx outside 1 and N"
        write(fhdl_log,'(A,F12.8)') "Error: particle position x       :", x_prtcl_curr(j)
        write(fhdl_log,'(A,I12)')   "Error: particle gcell    idx     :", gcell_idx
        flush(fhdl_log)
        stop 
      end if
      gcell_prtcl_curr(j) = gcell_idx
    end do

  end subroutine

  !This subprogram is used for debugging  
  !called after injection, eddy mapping, or settling of aerosols
  !call calc_prtcl_gcell_from_pos to update gcell_prtcl_count variable before calling this function
  subroutine print_nprtcl_in_gcell()

    implicit none
    integer           :: j, gcell_idx, ngt4, nis4, nis3, nis2, nis1, max_prtcl
    integer           :: debug_lcl
    
    debug_lcl         = 0

    ngt4              = 0 
    nis4              = 0
    nis3              = 0
    nis2              = 0
    nis1              = 0

    gcell_prtcl_count  = 0
    do j = 1, n_prtcl_curr
      gcell_idx                   = gcell_prtcl_curr(j)
      gcell_prtcl_count(gcell_idx) = gcell_prtcl_count(gcell_idx) + 1;
    end do

    max_prtcl   = 0 
    do j = 1, n_gcell

      if ( gcell_prtcl_count(j) > max_prtcl) max_prtcl = gcell_prtcl_count(j)

      if (gcell_prtcl_count(j)       > 4) then
          ngt4 = ngt4 + 1
      else if (gcell_prtcl_count(j) == 4) then
          nis4 = nis4 + 1
      else if (gcell_prtcl_count(j) == 3) then
          nis3 = nis3 + 1
      else if (gcell_prtcl_count(j) == 2) then
          nis2 = nis2 + 1
      else if (gcell_prtcl_count(j) == 1) then
          nis1 = nis1 + 1
      end if      
    end do

    if (debug_lcl == 1 .or. debug_gbl == 1) then
      write(fhdl_log,*) 'Maximum particles in a grid cell:',   max_prtcl
      write(fhdl_log,*) ngt4+nis4+nis3+nis2+nis1, ' cells have particles'
      write(fhdl_log,*) ngt4, ' cells have >4 particles'
      write(fhdl_log,*) nis4, ' cells have 4 particles'
      write(fhdl_log,*) nis3, ' cells have 3 particles'
      write(fhdl_log,*) nis2, ' cells have 2 particles'
      write(fhdl_log,*) nis1, ' cells have 1 particles'
      flush(fhdl_log)
    endif

  end subroutine

 !Mapping of particles in each grid cell  
 !when an eddy moves a grid cell, all particles in that grid cell is moved with it
 !This mapping makes moving of particle quick and efficient
 subroutine update_gcell_prtcl_map()

    implicit none
    integer                :: j, gcell_idx, prtcl_count
    
    !Find particles' grid cell position and updates the variable "gcell_prtcl_curr"
    call update_prtcl_gcell_from_xpos()

    !clear the exisiting map
    gcell_prtcl_currid_list = 0

    do j = 1, n_prtcl_curr
      gcell_idx                           = gcell_prtcl_curr(j)
      prtcl_count                         = gcell_prtcl_currid_list(0,gcell_idx)
      prtcl_count                         = prtcl_count + 1
      gcell_prtcl_currid_list(0,gcell_idx) = prtcl_count;
      if (prtcl_count > max_prtcl_gcell) then
        write(fhdl_log,'(A)')           "Error: Number of particles  in a gcell > max_prtcl_gcell parameter"
        write(fhdl_log,'(A,I12,A,I12)') "Error: Number of particles  at grid cell# ", gcell_idx, " is " , prtcl_count
        write(fhdl_log,'(A,I12)')       "Error:  max_prtcl_gcell parameter = ", max_prtcl_gcell
        flush(fhdl_log);
        stop
      end if
      !The map has current particle index and not the unique id  
      gcell_prtcl_currid_list(prtcl_count,gcell_idx) = j;
    end do
    
    
 end subroutine
 
 !Move particles by an eddy
 !uses subprogram update_gcell_prtcl_map
 subroutine move_prtcl_by_triplet_map(triplet_map_idx)
   implicit none
   !triplet_map_idx(i) has the old location for the grid cell i 
   double precision,intent(in)     :: triplet_map_idx(:)
   integer                         :: n_prtcl_gcell, k, new_idx, old_idx, dpos_gcell, curr_idx_prtcl
   integer                         :: debug_lcl
   double precision                :: dx_gcell , dpos_sum1, dpos_sum2, tpos

   debug_lcl                       = 0

   !write(fhdl_log,*) "------------ Before moving particle by eddy map ------------"
   !call drop_map(      n,       i1,        m1,  ncell,flag,index,   ibound1,ibound2)
   !call drop_map(len_eddy,sidx_eddy,gcell_1_idx,n_gcell,flag,index,gcell_1_idx,n_gcell)
   
   !update grid cell and curr particle mapping
   call update_gcell_prtcl_map()

   did_edymov_prtcl  = 0
   dpos_sum1         = 0.0
   dpos_sum2         = 0.0
   do new_idx        =   1, n_gcell
    old_idx          = int(triplet_map_idx(new_idx)) 
    !change in gcell position due to eddy mapping
    dpos_gcell       = (new_idx - old_idx)
    
    !change in gcell position in meters due to eddy mapping
    dx_gcell         = dpos_gcell * len_gcell

    n_prtcl_gcell    = gcell_prtcl_currid_list(0,old_idx)
    
    tpos          = 0.0
    if (dpos_gcell /= 0 .and.  n_prtcl_gcell > 0) then
      !write(fhdl_log,*) "gcell_idx, n_prtcl_gcell        :", new_idx, n_prtcl_gcell
      !verfication
      dpos_sum1      = dpos_sum1 + abs(dpos_gcell)
      do k = 1, n_prtcl_gcell
        !gcell_prtcl_currid_list has curr particle idx (1:n_prtcl_curr) not uniq particle id
        curr_idx_prtcl                   = gcell_prtcl_currid_list(k,old_idx);
        x_prtcl_curr(curr_idx_prtcl)     = x_prtcl_curr(curr_idx_prtcl) + dx_gcell
        did_edymov_prtcl(curr_idx_prtcl) = 1
        tpos                             = tpos + abs(dx_gcell)
        !verification
        if(x_prtcl_curr(curr_idx_prtcl) < 0 .or. x_prtcl_curr(curr_idx_prtcl) > len_dom) then
          write(fhdl_log,*) ""
          write(fhdl_log,'(A8,A)') "edymp:", "Error: particle position outside domain"
          write(fhdl_log,'(A8,A)') "edymp:", "Error: particle active idx  :", curr_idx_prtcl
          write(fhdl_log,'(A8,A)') "edymp:", "Error: particle new position:", x_prtcl_curr(curr_idx_prtcl)
          flush(fhdl_log)
          stop
        end if

      end do
      !write(fhdl_log,*) "tpos:", tpos, " n_prtcl_gcell:", n_prtcl_gcell, " len_gcell:", len_gcell
      dpos_sum2       = dpos_sum2 + (tpos/n_prtcl_gcell)/len_gcell
    end if
    
   end do

  
   !After particle repositioning by eddy; update gcell idx
   call update_prtcl_gcell_from_xpos()
   call update_gcell_prtcl_map()

   if ( abs(dpos_sum1-dpos_sum2) > 1e-4 )then
    write(fhdl_log,'(A8,A)')"edymp:", "Error: particle movment by eddies dpos_sum1 /= dpos_sum2"
    flush(fhdl_log)
   endif

   if (debug_lcl == 1 .or. debug_gbl == 1 .or. abs(dpos_sum1-dpos_sum2) > 1e-4) then
    write(fhdl_log,'(A)') " "
    write(fhdl_log,'(A8,42x,A)')"edymp:", "************ Droplets moved by eddies ****"
    write(fhdl_log,'(A8,A48,F16.2)') "edymp:","Sum of change in prtcl pos, sum1:" , dpos_sum1
    write(fhdl_log,'(A8,A48,F16.2)') "edymp:","Sum of change in prtcl pos, sum2:" , dpos_sum2
    write(fhdl_log,'(A)') " "

    flush(fhdl_log)
   endif


 end subroutine  

 !Iterate through each particle in the domain 
 !call indiv_prtcl_term_vel
 subroutine prtcl_settling(T_arr,qv_arr, press, deltat_ndim, N)
    use const, only : rho_w, dyn_vis, g , Rd

    !input
    double precision, intent(in) :: T_arr(N), qv_arr(N), press, deltat_ndim
    integer             :: N

    !local variables
    double precision    :: term_vel, tintvl
    integer             :: j,   idx, n_fallout, gidx, debug_lcl
    double precision    :: T_dim(N), qv_dim(N)
    
    !toggle debug mode for this subprogram; on or off represented as 1 or 0
    debug_lcl           = 0 !debug_mode    
    tintvl              = deltat_ndim * tscale_nu
    T_dim               = T_bot   - T_arr  * T_diff + 273.15
    qv_dim              = qvs_bot - qv_arr * qvs_diff
    
    if (n_prtcl_curr < 1) then
      return
    end if
    
    if (debug_lcl == 1 .or. debug_gbl == 1) then
      write(fhdl_log,*)""
      write(fhdl_log,'(A8,42x,A)') "psetl:","************ Particle settling ****"
      write(fhdl_log,'(A8,A48,F16.4,A)') "psetl:","time interval          :", tintvl , " secs"
      write(fhdl_log,'(A8,A48,I16)') "psetl:","before settling: number of particles active :", n_prtcl_curr
    endif

    !update the particles position
    call update_prtcl_gcell_from_xpos()

    do j = 1,n_prtcl_curr
      !to debug all particles, change condition to j>0
      if ( (debug_lcl == 1 .or. debug_gbl == 1) .and. (j == 1) ) then
        write(fhdl_log,*)""
        write(fhdl_log,'(A8,A48,I16)')      "psetl:","particle#...............:", j
        write(fhdl_log,'(A8,I3,A45,I16)')   "psetl:", j, " old pos from eddy movement:", did_edymov_prtcl(j)
        write(fhdl_log,'(A8,I3,A45,F16.4)') "psetl:", j, " old pos  (non-dim)        :", x_prtcl_curr(j) 
        
      end if

      !v(i) = c*(dsize(i)**2.0)
      
      gidx              = gcell_prtcl_curr(j)
      term_vel          = 0.0d0
      call indiv_prtcl_term_vel(T_dim(gidx), qv_dim(gidx), press, r_prtcl_curr(j), term_vel_factor, term_vel)

      !bug fix convert to non-dim
      x_prtcl_curr(j)   = x_prtcl_curr(j) + term_vel * tintvl/len_dom

      if (debug_lcl == 1 .or. debug_gbl == 1 .and. (j == 1) ) then
        write(fhdl_log,*)""
        write(fhdl_log,'(A8,I3,A45,F16.4,A)') "psetl:",j, " fall velocity          :", term_vel , " m/s"
        write(fhdl_log,'(A8,I3,A45,F16.4,A)') "psetl:",j, " displacement           :", term_vel * tintvl,  " m"
        write(fhdl_log,'(A8,I3,A45,F16.4)')   "psetl:",j, " displacement(non-dim)  :", term_vel * tintvl/len_dom
        write(fhdl_log,'(A8,I3,A45,F16.4)')   "psetl:",j, " new position (non-dim) :", x_prtcl_curr(j)
      end if

    end do

    !remove particles outside the domain
    !if particle pos is below domain domain boundary remove it
    idx                   = 0
    n_fallout             = 0
    r_prtcl_fallout       = 0.0d0
    r_aero_fallout        = 0.0d0
    r_crit_fallout        = 0.0d0
    ql_prtcl_fallout      = 0.0d0

    do j = 1,n_prtcl_curr
      if (x_prtcl_curr(j) < 0 ) then
        
        n_fallout = n_fallout + 1;
        !store properties of fallout partciles to compute stats 
        r_prtcl_fallout(n_fallout)    = r_prtcl_curr(j)
        r_aero_fallout(n_fallout)     = r_aero_curr(j)
        r_crit_fallout(n_fallout)     = r_crit_curr(j)
        ql_prtcl_fallout(n_fallout)   = ql_prtcl_curr(j)
                       

      else
        idx                   = j - n_fallout
        !shift particle properties by n_fallout - neat technique
        !This will overwrite the properties of particle that settled at bottom boundary 
        !since idx is always < j
        x_prtcl_curr(idx)      = x_prtcl_curr(j)
        r_prtcl_curr(idx)      = r_prtcl_curr(j)
        m_aero_curr(idx)       = m_aero_curr(j)
        r_aero_curr(idx)       = r_aero_curr(j)
        r_crit_curr(idx)       = r_crit_curr(j)
        ql_prtcl_curr(idx)     = ql_prtcl_curr(j)
        unqid_prtcl_curr(idx)  = unqid_prtcl_curr(j)
        gcell_prtcl_curr(idx)  = gcell_prtcl_curr(j)
        is_actvd_curr(idx)     = is_actvd_curr(j)
      end if
    end do

    !change number of particles in the domain
    if (n_fallout > 0) then
      n_prtcl_curr               = n_prtcl_curr - n_fallout

      !reset the rest of the particle property array to zero
      x_prtcl_curr(n_prtcl_curr+1:max_prtcl_dom)       = 0.0d0
      r_prtcl_curr(n_prtcl_curr+1:max_prtcl_dom)       = 0.0d0
      r_aero_curr(n_prtcl_curr+1:max_prtcl_dom)        = 0.0d0
      r_crit_curr(n_prtcl_curr+1:max_prtcl_dom)        = 0.0d0
      m_aero_curr(n_prtcl_curr+1:max_prtcl_dom)        = 0.0d0
      ql_prtcl_curr(n_prtcl_curr+1:max_prtcl_dom)      = 0.0d0
      unqid_prtcl_curr(n_prtcl_curr+1:max_prtcl_dom)   = 0.0d0
      gcell_prtcl_curr(n_prtcl_curr+1:max_prtcl_dom)   = 0.0d0
      is_actvd_curr(n_prtcl_curr+1:max_prtcl_dom)      = 0

    end if
    
    !recalculate particle gcell idx after setting
    call update_prtcl_gcell_from_xpos()
    !recalculate mapping between gcell-particles and verify if max_prtcl_gcell is not exceeded
    call update_gcell_prtcl_map()

    !write(fhdl_log,*) "------------ After settling ------------"
    !call print_nprtcl_in_gcell()

    n_fallout_alltime   = n_fallout_alltime + n_fallout
    n_fallout_curr      = n_fallout

    if (debug_lcl == 1 .or. debug_gbl == 1) then
      write(fhdl_log,'(A8,A48,I12)') "psetl:","after settling: number of particles in the dom :", n_prtcl_curr
      write(fhdl_log,'(A8,A48,I12)') "psetl:","after settling: number of particles fellout:", n_fallout
      
      write(fhdl_log,'(A8,A48,I12)') "psetl:","number of particles injected all time:", n_prtcl_alltime
      write(fhdl_log,'(A8,A48,I12)') "psetl:","number of particles fellout all time :", n_fallout_alltime
    endif  

    !This happened during a particular simulation but I was not able to reproduce the error
    !because the random number generation was based on time() for random seed
    if (n_prtcl_curr > n_prtcl_alltime ) then
      write(fhdl_log, *) "Error: n_prtcl_curr > n_prtcl_alltime"
      write(fhdl_log, '(A8,A48,I12)') "psetl:","Error: n_prtcl_curr        :", n_prtcl_curr
      write(fhdl_log, '(A8,A48,I12)') "psetl:","Error: n_prtcl_alltime     :", n_prtcl_alltime
      flush(fhdl_log)
      stop
    end if

    if (debug_lcl == 1 .or. debug_gbl == 1) then
      flush(fhdl_log)
    end if
    
 end subroutine

 !A particle's position is changed based on its terminal velocity
 !computed from stokes law for particles with radius  < 40 micrometers
 subroutine  indiv_prtcl_term_vel(T_env, qv_env, p_env, r_prtcl,vel_factor, term_vel)
  !the formula used here is applicable only for particles with radius  < 40 micrometers

  use const, only : rho_w, dyn_vis, Rd, g
  double precision, intent(in)      :: T_env, qv_env, p_env, r_prtcl, vel_factor
  double precision, intent(inout)   :: term_vel
  !local variable
  double precision      :: rho_air, Tv, c

  Tv                    = T_env * ( 1 + 0.622 * qv_env )
  rho_air               = p_env/ (Rd * Tv)
  c                     = (2.0 * g * rho_w) / (9.0 * dyn_vis * rho_air)
  !v(i) = c*(dsize(i)**2.0)
  !-1 for negative Z direction
  term_vel              = -1.0d0 * vel_factor * c * (r_prtcl**2.0)

 end subroutine

 !subprogram for verification and diagnostics
 function prtcl_count_in_a_gcell(gcell_idx)
  implicit none
   
  integer , intent(in) :: gcell_idx

  !local variable
  integer              :: j, prtcl_count_in_a_gcell, count

  count = 0

  do j = 1, n_prtcl_curr
    if (gcell_prtcl_curr(j) == gcell_idx) then
      count = count + 1
    endif
  end do 

  prtcl_count_in_a_gcell = count

 end

 !Iterate through particles in the domain and call droplet growth model
 !implemented in subprogram indiv_prtcl_growth
 subroutine prtcl_growth(T_arr, qv_arr, ss_arr, press, t1,t2, N, err_stat)
  use const, only                   : pi43, rho_w, Rd
  use array, only                   : err_dgm 
  implicit none 

  !input
  double precision, intent(in)                 :: press
  double precision, intent(in)                 :: t1, t2 
  integer, intent(in)                          :: N 
  double precision,  intent(inout)             :: T_arr(N), qv_arr(N), ss_arr(N)
  integer, intent(inout)                       :: err_stat
   
  !external subroutines
  EXTERNAL fcnkb,rkqs,fcnkb0, odeint

  !local variables 
  integer             :: j, gcell_idx, debug_lcl, n_same_prtcl_type, sum_actv_states
  double precision    :: dt, dt_min, r_mean, ql_sum
  double precision    :: T_dim(N), qv_dim(N), time1, time2
  double precision    :: T_env, qv_env, p_env, ss_env, ql_prtcl, r_prtcl , qvs_env, ql_verf
  double precision    :: r_crit, ss_crit, r_prtcl_before, eps_err, inv_mass_gcell, Tv, rho_air, gscale
  integer             :: is_iter  
  
  !turn debug mode for this subprogram; on or off represented as 1 or 0
  debug_lcl           = 0 !debug_mode
  r_mean              = 0.0d0
  ql_sum              = 0.0d0
  err_stat            = 0
  n_actvd_curr        = 0
  n_deactvd_curr      = 0
  n_same_prtcl_type   = 0
  is_iter             = 0
  is_actvd_curr       = 0
  eps_err             = 1e-4

  if (n_prtcl_curr < 1) then 
    return
  end if
  
  r_mean              = sum(r_prtcl_curr(1:n_prtcl_curr))/n_prtcl_curr
  ql_sum              = sum(ql_prtcl_curr(1:n_prtcl_curr))
  

  if (debug_lcl == 1 .or. debug_gbl == 1) then
    write(fhdl_log,*) ""
    write(fhdl_log, "(A8,42x,A)")    "pgrow:", "************ All particles growth ****"
    write(fhdl_log,'(A8,A48,I16)')   "pgrow:", "number of particles :", n_prtcl_curr
    write(fhdl_log,'(A8,A48,E16.4)') "pgrow:", "before prtcl growth: Mean radius:", r_mean
    write(fhdl_log,'(A8,A48,E16.4)') "pgrow:", "before prtcl growth: total particle water :", ql_sum
  endif

  !Convert T and qv from non-dim to dim
  !T_arr has value 0 to 1 starting from bottom to top; hence subtract it
  !similarly for qvs
  T_dim             = T_bot   - T_arr * T_diff + 273.15
  qv_dim            = qvs_bot - qv_arr * qvs_diff

  !integration time - conv from non-dim to dim 
  time1             = t1 * tscale_nu  
  time2             = t2 * tscale_nu 
  dt                = (time2 - time1)


  if (debug_lcl == 1 .or. debug_gbl == 1) then
    write(fhdl_log,*) ""
    
    write(fhdl_log, "(A8,A48,F16.4,A)") "pgrow:","dim dt:", dt , " secs"
    write(fhdl_log, "(A8,A48,F16.4)")   "pgrow:","Time t1:", time1
    write(fhdl_log, "(A8,A48,F16.4)")   "pgrow:","Time t2:", time2
    write(fhdl_log, "(A8,A48,F16.4,A)")   "pgrow:","nondim dt:", t2-t1, " secs"
    write(fhdl_log, "(A8,A48,F16.4,A)") "pgrow:","Time scale:", tscale_nu/60 , " mins"
    
    flush(fhdl_log)
  end if
  

  !minimum time step between t1 and t2 in RK method
  !RK 5th order integration between t1 and t2 will use timesteps smaller than  (t2-t1) to get accurate integral
  dt_min              = 0.0d0

  do j = 1, n_prtcl_curr
    
    r_prtcl           = r_prtcl_curr(j)
    r_prtcl_before    = r_prtcl_curr(j)
    gcell_idx         = gcell_prtcl_curr(j)

    T_env             = T_dim(gcell_idx)
    qv_env            = qv_dim(gcell_idx)
    p_env             = press
    qvs_env           = sat_mixing_ratio(T_env-273.15,p_env)
    !ss calculated in ODT is incorrect so compute directly from qv and qvs
    !ss_env           = (qv_env/qvs_env-1.0d0)
    !conv percentage to fraction
    ss_env            = ss_arr(gcell_idx)/100 
    ql_prtcl          = ql_prtcl_curr(j)
    Tv                = T_env * ( 1 + 0.622 * qv_env )
    rho_air           = p_env/ (Rd * Tv)
    !mass of dry air in gcell to compute mixing ratio
    inv_mass_gcell    = 1.0d0/(rho_air*vol_gcell) 
    !inv_mass_gcell    = 1.0d0/(rho_air*vol_gcell*(1-qv_env)) 
    gscale            = inv_mass_gcell
    !gscale            = inv_vol_gcell

    !since dgm gives error if time steps are close >= 1 sec
    !It is better to call dgm itervatively through small time steps
    !I chose the threshold 0.1 just to be safe 
    is_iter     = 0
    if (dt >= 0.01) then 
      is_iter   = 1
      !substepping helps solve the growth error than using eps
      !eps_err   = 1e-6
    endif  

    if ( (debug_lcl == 1 .or. debug_gbl == 1) .and. j == 1 ) then
      !debug a particular particle (j == 58) to avoid cluttering the log file
      !to debug all particles, change condition to (j > 0)
      write(fhdl_log,*) ""
      write(fhdl_log, '(A8,A48,I16)') "pgrow:","-------- Droplet# ", j
      write(fhdl_log, '(A8,A48,I16)') "pgrow:","gcell idx:", gcell_idx
      write(fhdl_log, '(A8,A48,I16)') "pgrow:","particles count in the gcell:", prtcl_count_in_a_gcell(gcell_prtcl_curr(j))
      flush(fhdl_log)
      !subroutine indiv_prtcl_growth(T_env, qv_env, p_env, ss_env, ql_prtcl, t1, t2, r_prtcl, m0_aerosol, aerosol_type,&
      !gscale, has_cs_effect, debug)
      call indiv_prtcl_growth(T_env, qv_env, p_env, ss_env, ql_prtcl, time1,time2, &
                r_prtcl,m_aero_curr(j),aero_type,gscale, has_curv_sol_eff,eps_err,is_iter,1)
    
    else
      call indiv_prtcl_growth(T_env, qv_env, p_env, ss_env, ql_prtcl, time1,time2, &
                r_prtcl,m_aero_curr(j),aero_type,gscale, has_curv_sol_eff,eps_err,is_iter,0)

    end if
    
    if ( isnan(r_prtcl) .or. err_dgm == -1 ) then
      write(fhdl_log,*) ""
      write(fhdl_log,'(A16,A48,I16)') "pgrow: Error:", "particle curr idx: ", j
      write(fhdl_log, "(A16,A48,F16.4,A)") "pgrow: Error","Time interval:", (time2 - time1), " secs"
      write(fhdl_log, "(A16,A48,F16.4)") "pgrow: Error","Time t1:", time1
      write(fhdl_log, "(A16,A48,F16.4)") "pgrow: Error","Time t2:", time2
      write(fhdl_log,*) ""
      write(fhdl_log,'(A16,A48,E16.8)') "pgrow: Error:", "particle gcell T (k): ", T_dim(gcell_idx)
      write(fhdl_log,'(A16,A48,E16.8)') "pgrow: Error:", "particle gcell qv (kg/kg): ", qv_dim(gcell_idx)
      write(fhdl_log,'(A16,A48,E16.8)') "pgrow: Error:", "particle gcell ss (%): ", ss_arr(gcell_idx)

      call dump_prtcl_state()
      err_stat = -1

      flush(fhdl_log)
      flush(6)
      return
    endif  

    !check if the particle size < dry aerosol size
    !don't update the radius, ql, T_env, qv_env, ss_env; 
    if (r_prtcl >= r_aero_curr(j)) then
      r_prtcl_curr(j)    = r_prtcl
      
      T_dim(gcell_idx)   = T_env
      qv_dim(gcell_idx)  = qv_env
      !override the value from the ode int
      !ss_arr(gcell_idx)  = ss_env * 100
      qvs_env            = sat_mixing_ratio(T_env-273.15,p_env)
      ss_env             = (qv_env/qvs_env - 1.0d0)
      ss_arr(gcell_idx)  = ss_env * 100

     !override ql value from ode int
     !ql_prtcl_curr(j)   = ql_prtcl
     !ql_prtcl is computed from the radius minus the mass of aerosol 
     !ql_prtcl_curr(j)   = (pi43 * r_prtcl**3 * rho_w) * inv_vol_dom

      !r_crit depends on the temperature
      r_crit             = 0.0d0
      ss_crit            = 0.0d0
      call find_r_crit(T_env,m_aero_curr(j),r_crit,ss_crit)
      r_crit_curr(j)    = r_crit  
      !r_crit_curr(j)   = r_min_instru

    else 
        write(fhdl_log, '(A8,A,I8)') "pgrow:", "error: prctl:", j, "growth and its effects ignored"
        write(fhdl_log, '(A8,A)') "pgrow:", "error: Since prtcl radius became smaller than dry aero radius"
        flush(fhdl_log)

    endif  

    !Find number of activated and deactivated particles before and after calling droplet growth model
    if (r_prtcl_before < r_crit .and. r_prtcl >= r_crit ) then
      n_actvd_curr        =  n_actvd_curr  + 1
      is_actvd_curr(j)    = 1
    elseif(r_prtcl_before >= r_crit .and. r_prtcl < r_crit )  then
      n_deactvd_curr      = n_deactvd_curr + 1
      is_actvd_curr(j)    = -1
    else
      n_same_prtcl_type   = n_same_prtcl_type + 1
      is_actvd_curr(j)    = 0
    endif  


  end do

  !calc ql based on radius 
  call calc_prtcl_ql(n_prtcl_curr,r_prtcl_curr,r_aero_curr,inv_vol_dom,ql_prtcl_curr);

  !conv T, and qv from dim to non dim
  T_arr                = (T_bot  - (T_dim - 273.15))/T_diff
  qv_arr               = (qvs_bot - qv_dim)/qvs_diff
  r_mean               = sum(r_prtcl_curr(1:n_prtcl_curr))/n_prtcl_curr
  ql_sum               = sum(ql_prtcl_curr(1:n_prtcl_curr))
  
  sum_actv_states = (n_actvd_curr+n_same_prtcl_type+n_deactvd_curr)
  if (debug_lcl == 1 .or. debug_gbl == 1 .or. sum_actv_states /= n_prtcl_curr) then
    write(fhdl_log,*) ""
    write(fhdl_log,'(A8,A48,E16.4)') "pgrow:","afer prtcl growth: Mean radius:", r_mean
    write(fhdl_log,'(A8,A48,E16.4)') "pgrow:","after prtcl growth: sum particle water:", ql_sum

    ql_verf = sum(pi43*r_prtcl_curr(1:n_prtcl_curr)**3.0d0*rho_w*inv_vol_dom)
    write(fhdl_log,'(A8,A48,E16.4)') "pgrow:","after prtcl growth: sum particle water+m_aero:", ql_verf
    write(fhdl_log,*) ""
    write(fhdl_log,'(A8,A48,I16)') "pgrow:","# prtcl activated:",   n_actvd_curr
    write(fhdl_log,'(A8,A48,I16)') "pgrow:","# prtcl deactivated:", n_deactvd_curr
    write(fhdl_log,'(A8,A48,I16)') "pgrow:","# prtcl same regime:", n_same_prtcl_type
    
    write(fhdl_log,'(A8,A48,I16)') "pgrow:","# prtcl sum above:",     sum_actv_states
    write(fhdl_log,'(A8,A48,I16)') "pgrow:","# prtcl in the domain:", n_prtcl_curr
    flush(fhdl_log)

    if (sum_actv_states /= n_prtcl_curr)then
      write(fhdl_log,'(A8,A48,I16)') "pgrow:", "Error: sum of prtcls by actvd states is not equal to prtcl in the dom"
      flush(fhdl_log)
    endif

  endif

 end subroutine

 !uses droplet growth model (dgm) code from EMPM
 !Following fortran files from EMPM should be compiled with this file
 !array.f90 const.f90 rand1.f90 rkqs.f90 fcnkb.f90 fcnkb0.f90 rkck.f90 ew.f90 odeint.f90 size_dis.f90
 subroutine indiv_prtcl_growth(T_env, qv_env, p_env, ss_env, ql_prtcl, t1, t2, r_prtcl, m0_aerosol, aerosol_type,&
                               gscale, has_cs_effect,eps_err, is_iterative, debug)
  use const, only                   : pi43, rho_w
  use array, only                   : err_dgm
  implicit none
  
  double precision, intent(inout)   :: T_env, qv_env, p_env, r_prtcl, ss_env, ql_prtcl
  double precision, intent(in)      :: t1, t2, m0_aerosol, gscale , eps_err
  integer, intent(in)               :: has_cs_effect, debug, aerosol_type
  !external subroutines
  EXTERNAL fcnkb,rkqs,fcnkb0, odeint

  !local variables 
  integer, parameter  :: n_yvar = 8
  double precision    :: y_arr(n_yvar), y_before(n_yvar), dt, dt_min, n_ok, n_bad, qvs_env
  double precision    :: t1_intg, t2_intg, tstep_intg
  integer             :: debug_lcl, iter_count, is_iterative

  debug_lcl         = debug
  tstep_intg        = 1e-2

  !integration time - conv from non-dim to dim
  dt                = (t2 - t1)


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
  call set_odeint_dgm_params(1, aerosol_type, m0_aerosol, gscale)

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
 
 !Subroutines odeint, fcnkb, fcnkb0 for particle growth is used as it is from the EMPM model
 !these subroutines uses parameters from the array.f90 module. Theses parameters could have beenpassed as arguments but don't know
 !why it is implemented this way. I am just following EMPM implementation and setting those variables here 
 !to avoid changing these files. 
 !This function should be called before calling the odeint (droplet growth function)
 subroutine set_odeint_dgm_params(gcell_sidx, aerosol_type, m0_aerosol, gscale)

  use array, only :  ndmax, solute_mass, dmaxa, grid_scale, err_dgm
  implicit none

  integer, intent(in)    :: gcell_sidx, aerosol_type
  double precision, intent(in) ::  m0_aerosol, gscale

  !starting index of array
  ndmax                 = gcell_sidx
  !aerosol type
  dmaxa                 = aerosol_type 
  !solute mass in kg
  solute_mass           = m0_aerosol
  !grid scale is actually inverse of grid mass/volume
  grid_scale            = gscale
  
  !reset the error flag
  err_dgm = 0

end subroutine
 
!if there is any error the particle state is dumped 
!This allows for program to restart from the same model state for debugging purposes.
!The steady state profiles of scalars from input directory is used for scalar state
!Future change: This subroutine and set_curr_prtcl_prop can be changed to store and restore 
!the scalar state too rather than using initial state
 subroutine dump_prtcl_state()
  implicit none

  
  !local variables
  character(len=:),allocatable      :: dump_dir
  integer                           :: j

  write(*,*) ""
  write(*,'(A8,A48,A)') "pgrow:", "**** Error:","during microphysics call"
  write(*,'(A8,A48,A)') "pgrow:", "**** Error:","particle properties is dumped to a file" 

  !change to parent directory
  call chdir(trim(parent_dir))
  write(*,'(A8,A48,A)') "pgrow:", "now at parent dir:",parent_dir
  !call system("pwd");

  !create dump directory
  dump_dir     = "input/"//exp_name//"/"//realz//"/"
  call system("mkdir -p "//dump_dir)
  call system("cp "//input_dir//"/"//"LabExppar.dat"//" "//dump_dir)
  call system("cp "//input_dir//"/"//namelist_file//" "//dump_dir)
  call chdir(dump_dir)

  write(*,'(A8,A48,A)') "pgrow:", "dir :", dump_dir 
  write(*,'(A8,A48,A)') "pgrow:", "file:", "prtcl_prop_curr.txt"
  write(*,*) ""


  !open scalar profile files
  open(701,file="prtcl_prop_curr.txt")

  do j = 1, n_prtcl_curr
    write(701,'(2I16,5E16.8)') j,  unqid_prtcl_curr(j), x_prtcl_curr(j), r_prtcl_curr(j),&  
          m_aero_curr(j), r_aero_curr(j), ql_prtcl_curr(j)
  enddo  
  !close scalar profile files
  close(701)

 end subroutine 

!restore the particle properties to the model state when error occurred.
 subroutine set_curr_prtcl_prop(fhdl_prtcl_prop)
  !restore the particle prop from the file dump
  !This feature to restore the model state and continuing it to run
  !is not tested rigourosly. So use it with care
  implicit none

  integer         :: fhdl_prtcl_prop, iostat, j, cidx
  
  write(fhdl_log,*) ""
  write(fhdl_log, '(A8,42x,A)') "init:", "************ Reading the particle property file and setting the variables"

  iostat    = 0
  j         = 0
  do 
    j   = j + 1
    read(fhdl_prtcl_prop,'(2I16,5E16.8)',iostat=iostat) cidx,  unqid_prtcl_curr(j), x_prtcl_curr(j), &
          r_prtcl_curr(j), m_aero_curr(j), r_aero_curr(j), ql_prtcl_curr(j)
    if (iostat < 0) then 
      j = j-1
      exit
    endif  
  enddo   

  n_prtcl_curr = j
  write(fhdl_log, '(A8,A48, I16)') "init:", "from file - n_prtcl_curr:", n_prtcl_curr

  if ( n_prtcl_curr > max_prtcl_dom)then
    write(fhdl_log, '(A8,A)') "init:", "Error: # particles added exceeds # of inst max particles  (array alloc)" 
    write(fhdl_log, '(A8,A,I16)') "init:", "Error: max particles allowed in the domain                :", max_prtcl_dom
    write(fhdl_log, '(A8,A,I16)') "init:", "Error: number of active particles currently in the domain :", n_prtcl_curr
    flush(fhdl_log);
    stop
  endif   
  
  unqid_prtcl_curr(n_prtcl_curr+1:max_prtcl_dom) = 0
  x_prtcl_curr(n_prtcl_curr+1:max_prtcl_dom)     = 0.0d0
  r_prtcl_curr(n_prtcl_curr+1:max_prtcl_dom)     = 0.0d0
  m_aero_curr(n_prtcl_curr+1:max_prtcl_dom)      = 0.0d0
  r_aero_curr(n_prtcl_curr+1:max_prtcl_dom)      = 0.0d0
  ql_prtcl_curr(n_prtcl_curr+1:max_prtcl_dom)    = 0.0d0
  is_actvd_curr(n_prtcl_curr+1:max_prtcl_dom)    = 0

  prtcl_counter   = maxval(unqid_prtcl_curr)
  n_prtcl_alltime = prtcl_counter
  n_inj_alltime   = prtcl_counter  
  call update_prtcl_gcell_from_xpos()
  call update_gcell_prtcl_map()
  
  if( inject_onetime == 1 ) then
    is_injected_onetime = 1
  endif  

 end subroutine  

!store lagrangian particle properties and particle environment
 subroutine rec_prtcl_hist(T_arr, qv_arr,ss_arr,N, t_ndim, deltat_ndim, flush_arr)
      implicit none
      integer, intent(in)           :: N
      double precision, intent(in)  :: T_arr(N), qv_arr(N), ss_arr(N), t_ndim, deltat_ndim
      integer                       :: flush_arr
      !local variable
      integer                       :: debug_lcl
      double precision              :: time_dim, deltat_dim

      debug_lcl         = 0  
      time_dim          = t_ndim * tscale_nu
      deltat_dim        = deltat_ndim * tscale_nu

      !The checks need in the following order else it will produce unexpected result

      !if user doesn't want to rec prtcl hist then do nothing
      if (do_rec_prtcl_hist /= 1) return

      !check the start and end time for recording prtcl hist
      if (time_dim < rec_stime_prtcl_hist ) return

      !Keep this as last check since the buffer needs to be flush
      if (time_dim > rec_etime_prtcl_hist)then
        !don't save current prtcl properties to buffer
        call write_prtcl_hist(time_next_phist_rec - rec_intvl_prtcl_hist)
        !Stop further recording or writing to file
        do_rec_prtcl_hist  = 0
        return
      endif 

      !flush_arr is typically used at end of simulation
      !Hence the current data need to stored and then written to the file
      !This if block should be executed even if there is no particle, not a rec invtvl, not a write intervel 

      if (flush_arr == 1) then
        call write_prtcl_hist(time_next_phist_rec)
        !returning since I don't want to do other checks
        return
      endif  

      !if there are no particles then return 
      !this will skip only current time step
      if (n_prtcl_curr < 1) return
 
      !Don't record curr prtcl data if rec intvl is not exceeded
      if( time_dim < time_next_phist_rec   ) return

      !check if time intvl for writing prtcl hist is exceeded 
      if( time_dim > time_next_phist_write )  then
        !Write the data to a file and then store the current prtcl data in to buffer
        call write_prtcl_hist(time_next_phist_write)
        time_next_phist_write = time_next_phist_write + write_intvl_prtcl_hist
      endif 

      !if buffer limit is exceeded 
      if( (bidx_rec_hist + n_prtcl_curr) > max_prtcl_hist ) then
        call write_prtcl_hist(time_next_phist_write-rec_intvl_prtcl_hist)
      end if

      !==== Store the curr prtcl data
      call store_prtcl_hist_to_buff(T_arr, qv_arr, ss_arr,N, t_ndim, deltat_ndim)
      if (debug_lcl == 1) then
        write(fhdl_log,'(A8,A48,F16.4)') "stats:", "phist: time of the record:", (time_next_phist_rec)
        write(fhdl_log,'(A8,A48,I16)')   "stats:", "phist: bidx_rec_hist :", bidx_rec_hist
        flush(fhdl_log) 
      endif
      time_next_phist_rec   = time_next_phist_rec + rec_intvl_prtcl_hist
      
 end subroutine
 
 !To avoid writing to files often, store the particle properties and its environment
 !into a buffer and then write it to a file, when the buffer is full or 
 !write interval (from namelist) has elaspsed
 subroutine store_prtcl_hist_to_buff(T_arr, qv_arr,ss_arr,N, t_ndim, deltat_ndim)
  implicit none
  integer, intent(in)           :: N
  double precision, intent(in)  :: T_arr(N), qv_arr(N), ss_arr(N), t_ndim, deltat_ndim
  !local variable
  integer                       :: j, gid
  double precision              :: T_dim(N), qv_dim(N)

  !==== store the particle properties and its env into buffer
  !Convert T and qv from non-dim to dim
  !T_arr has value 0 to 1 starting from bottom to top; hence subtract it
  !similarly for qvs
  T_dim             = T_bot   - T_arr * T_diff    + 273.15
  qv_dim            = qvs_bot - qv_arr * qvs_diff 

  do j  = 1, n_prtcl_curr
    bidx_rec_hist                      = bidx_rec_hist + 1

    time_prtcl_hist(bidx_rec_hist)     = t_ndim
    deltat_prtcl_hist(bidx_rec_hist)   = deltat_ndim
    uniqid_prtcl_hist(bidx_rec_hist)   = unqid_prtcl_curr(j)
    r_prtcl_hist(bidx_rec_hist)        = r_prtcl_curr(j)
    r_aero_hist(bidx_rec_hist)         = r_aero_curr(j)
    r_crit_hist(bidx_rec_hist)         = r_crit_curr(j)
    x_prtcl_hist(bidx_rec_hist)        = x_prtcl_curr(j)
    
    gid                                = gcell_prtcl_curr(j)
    T_prtcl_hist(bidx_rec_hist)        = T_dim(gid)
    qv_prtcl_hist(bidx_rec_hist)       = qv_dim(gid)
    ss_prtcl_hist(bidx_rec_hist)       = ss_arr(gid)

  end do

 end subroutine 

 !write the particle properties and its environment to a file
 subroutine write_prtcl_hist(time_write)
  implicit none

  double precision,intent(in)     :: time_write

  !local variables 
  integer                         :: fhdl_prtcl_hist, debug_lcl , j
  character(len=64)               :: fname

  debug_lcl                       = 1 
  fhdl_prtcl_hist                 = 951

  call system("mkdir -p hist_prtcl")
  write(fname,'(A,I5.5,A)')  "hist_prtcl_", nint(time_write), ".txt"
  open(fhdl_prtcl_hist, file="hist_prtcl/"//trim(fname))

  if( debug_lcl == 1 )then
    write(fhdl_log,'(A8,A48,A)') "stats:", "hist_prtcl filename:  ", "hist_prtcl/"//trim(fname)
    write(fhdl_log,'(A8,A48,F16.4)') "stats:", "time_write:", time_write
    write(fhdl_log,'(A8,A48,I16)') "stats:", "bidx_rec_hist:", bidx_rec_hist
    flush(fhdl_log)
  endif  

  write(fhdl_prtcl_hist,'(A,A15,10A16)') "#",   "time_ndim", "time_dim",    "delta_time", &
          "uniq_prtcl_id",    "r_prtcl",        "r_aero",     "r_crit",     "x_pos_ndim", &
          "Temperature",      "water_vapor",    "supersat"
  write(fhdl_prtcl_hist,'(A,A15,10A16)') "#",   "ndim",       "secs",    "secs", &
          "num",    "meters",   "meters",       "meters",     "ndim", &
          "K",      "kg/kg",    "%"
          
  do j = 1, bidx_rec_hist
    write(fhdl_prtcl_hist,'(E16.8,2F16.6,I16,7E16.8)') time_prtcl_hist(j), time_prtcl_hist(j)*tscale_nu, &
      deltat_prtcl_hist(j)*tscale_nu, uniqid_prtcl_hist(j),  r_prtcl_hist(j),  r_aero_hist(j), r_crit_hist(j), &
      x_prtcl_hist(j),                      T_prtcl_hist(j),       qv_prtcl_hist(j), ss_prtcl_hist(j)
  end do
  
  bidx_rec_hist       = 0
  time_prtcl_hist     = 0.0d0
  deltat_prtcl_hist   = 0.0d0
  uniqid_prtcl_hist   = 0
  r_prtcl_hist        = 0.0d0
  r_aero_hist         = 0.0d0
  r_crit_hist         = 0.0d0
  T_prtcl_hist        = 0.0d0
  qv_prtcl_hist       = 0.0d0
  ss_prtcl_hist       = 0.0d0

  flush(fhdl_prtcl_hist)
  close(fhdl_prtcl_hist)

 end subroutine 

 !calc saturation mixing ration in kg/kg for a given temperature and pressure
 function sat_mixing_ratio(T_in_C, press)
  implicit none
    double precision :: sat_mixing_ratio , T_in_C, es, press
    double precision, EXTERNAL :: ew

    es = 6.112d0 * dexp(17.67d0*(T_in_C)/(T_in_C+243.5)) * 100.d0
    !es  = ew(T_in_C + 273.15) * 100
    sat_mixing_ratio = 0.622d0 * es/(press - es)

    return
 end
 
 !create bin edges for particle position
 !Currently these bins are not used but created this to study how the particles
 !spread within the domain if the injected at the center of the domain
 function create_bedges_xpos() 
  implicit none

  integer,parameter                               :: n_bedges = 51
  double precision,dimension(n_bedges)            :: xpos_bedges, create_bedges_xpos
  integer                                         :: j
  double precision                                :: bwidth
  
  xpos_bedges       = 0.0d0
  bwidth            = 1.0d0/(n_bedges-1)
  do j  = 1, n_bedges
    xpos_bedges(j)  =  (j-1) * bwidth
  end do

  create_bedges_xpos = xpos_bedges
 end

 !create bin edges for particle size distribution ( based on particle radius) 
 function create_bedges_r_prtcl()
  implicit none
  integer,parameter                               :: n_bedges = 73
  double precision,dimension(n_bedges)            :: rprtcl_bedges, create_bedges_r_prtcl
  integer                                         :: j

  rprtcl_bedges(1:18)     = (/ (j*0.01d0*0.5,  j=2,19) /)
  rprtcl_bedges(19:36)    = (/ (j*0.1d0*0.5,   j=2,19) /)
  rprtcl_bedges(37:54)    = (/ (j*1.0d0*0.5,   j=2,19) /)
  rprtcl_bedges(55:72)    = (/ (j*10.0d0*0.5,  j=2,19) /)
  rprtcl_bedges(73)       = 100.0d0
  !radius binedges in microns
  rprtcl_bedges           = rprtcl_bedges  * 1.0d-6

  !write(fhdl_log,*) "bedges for r_prtcl:"
  !write(fhdl_log,'(500E12.4)') rprtcl_bedges

  create_bedges_r_prtcl = rprtcl_bedges

 end 
 
 !Create bin edges for scalar (T, qv) values from 0 to 1,since they are in the non-dimesional form
 function create_bedges_scalar()
  implicit none
  integer,parameter                               :: n_bedges = 100 + 1
  double precision,dimension(n_bedges)            :: scalar_bedges, create_bedges_scalar
  integer                                         :: j

  scalar_bedges      = 0.0d0

  do j = 1, n_bedges
    scalar_bedges(j) =  (j-1) * 1.0d0 / (n_bedges-1.0d0)
  end do
  
  !write(fhdl_log,*) "bedges for scalar:"
  !write(fhdl_log,'(500E12.4)') scalar_bedges

  create_bedges_scalar = scalar_bedges

 end 

 !create bin edges from supersaturation 
 function create_bedges_supersat()
  implicit none

       
  integer,parameter                     :: n_bedges = 111, n_by_2=55
  double precision,dimension(n_bedges)  :: ss_bedges, create_bedges_supersat
  double precision,dimension(n_by_2)    ::temp
  integer                                         :: j

  ss_bedges      = 0.0d0

  temp(1:18)     = (/ (j*0.01d0*0.5,  j=2,19) /)
  temp(19:36)    = (/ (j*0.1d0*0.5,   j=2,19) /)
  temp(37:54)    = (/ (j*1.0d0*0.5,   j=2,19) /)
  temp(55)       = 10.0d0
  
  !create symmetrical negative and postive ss log bins
  do j = 1, n_by_2
    ss_bedges(j)              = -1.0d0 * temp(n_by_2+1-j)
    ss_bedges(n_by_2+1+j)     =  1.0d0 * temp(j)
  enddo  
  ss_bedges(n_by_2+1)         = 0.0d0      

  
  !write(fhdl_log,*) "bedges for supersaturation:"
  !write(fhdl_log,'(500E12.4)') ss_bedges

  create_bedges_supersat = ss_bedges

 end 
 
 !If the number of grid cells are in thousands
 !creating a profile of particle statistic at each grid cell is computationally expensive
 !Also there may not be many particles to create a reliable statistic
 !for this reason the vertical axis(xpos) is divided into bins with small bins close to boundary layer to resolve variability there
 !and large bins at the core. The particle statistic is computed for particles within the same xpos bin
 function create_bedges_pprof()
  implicit none
  
  integer,parameter                :: max_bedges=1000
  double precision                 :: x_bedges(max_bedges)
  double precision,allocatable     :: create_bedges_pprof(:)
  integer                          :: bidx, k, n_half_bedges, n_bedges

  double precision                 :: diff_depth, bl_depth, bl_dx, bulk_center, bulk_dx , len_gcell_ndim, diff_dx

  len_gcell_ndim     = len_gcell/len_dom
  !diffusion layer is 10% of domain 
  diff_depth         = 0.1  
  !diffusion layer bin width
  diff_dx            = 0.001d0
  !boundary layer is 20% of the domain
  bl_depth           = 0.2d0
  !boundary layer bin width 
  bl_dx              = 0.01d0 
  !domain center
  bulk_center        = 0.5d0
  bulk_dx            = 0.01d0

  
  x_bedges           = 0.0d0
  bidx               = 1
  x_bedges(bidx)     = 0.0d0

  !diff layer  
  do while(x_bedges(bidx) <  diff_depth)
    bidx             = bidx + 1
    x_bedges(bidx)   = x_bedges(bidx-1) + diff_dx
  end do
  write(fhdl_log,'(A,F16.4)') "last bedge:" , x_bedges(bidx)

  !boundary layer
  do while(x_bedges(bidx) <  bl_depth)
    bidx             = bidx + 1
    x_bedges(bidx)   = x_bedges(bidx-1) + bl_dx
  end do
  write(fhdl_log,'(A,F16.4)') "last bedge:" , x_bedges(bidx)

  !bulk layer
  do while( x_bedges(bidx) <  bulk_center)
    bidx             = bidx + 1
    x_bedges(bidx)   = x_bedges(bidx-1) + bulk_dx
  end do
  write(fhdl_log,'(A,F16.4)') "last bedge:" , x_bedges(bidx)

  n_half_bedges      = bidx

  if (  (2*n_half_bedges) > max_bedges) then
    write(fhdl_log,'(A,A64)')  "init:",  "Error creating bedges for particle prof - small array size"  
    write(fhdl_log,'(A,A64,I16)')  "init:",  "Error n_half_bedges:", n_half_bedges
    write(fhdl_log,'(A,A64,I16)')  "init:",  "Error max_bedges   :", max_bedges
    stop 
  endif

  !make bedges symmeterical; but results in 
  !central bin encompossing a twice the volume
  bidx = bidx - 1
  n_half_bedges = n_half_bedges - 1
  write(fhdl_log,'(A,F16.4)') "center bedge:" , x_bedges(bidx)

  !create symmetrical bedges for other half of the domain
  do k = n_half_bedges, 1, -1
    bidx             = bidx + 1
    x_bedges(bidx)   = 1.0d0 - x_bedges(k)
  enddo  
  n_bedges            = bidx
  write(fhdl_log,'(A,F16.4)') "last bedge:" , x_bedges(bidx)
  allocate(create_bedges_pprof(n_bedges))

  create_bedges_pprof = x_bedges(1:n_bedges)

 end 
 
 !master subroutine that calls other subprograms to create the bin edges
 subroutine create_pdf_bins()
  implicit none
  character(len=:),allocatable      :: pname          

  !====bins for particle radius
  bedges_r_prtcl        = create_bedges_r_prtcl()
  n_bedges_r_prtcl      = size(bedges_r_prtcl, dim=1)
  n_bins_r_prtcl        = n_bedges_r_prtcl - 1

  pname                 = "particle"
  call write_rdist_fhdr(fhdl_pdf_r_prtcl,pname)
  pname                 = "aerosol"
  call write_rdist_fhdr(fhdl_pdf_r_aero,pname)
  pname                 = "haze droplet"
  call write_rdist_fhdr(fhdl_pdf_r_hdplt,pname)
  pname                 = "cloud droplet"
  call write_rdist_fhdr(fhdl_pdf_r_cdplt,pname)
  pname                 = "instrument measured particle"
  call write_rdist_fhdr(fhdl_pdf_r_instru,pname)

  !====bins for scalar
  pname                 =  "non-dim temperature" 
  bedges_T              = create_bedges_scalar()
  n_bedges_T            = size(bedges_T, dim=1)    
  n_bins_T              = n_bedges_T-1         
  call write_pdf_fhdr(fhdl_pdf_dT,bedges_T,n_bedges_T, pname)
  call write_pdf_fhdr(fhdl_pdf_bT,bedges_T,n_bedges_T, pname)

  pname                 =  "non-dim water vapor" 
  bedges_qv             = create_bedges_scalar()  
  n_bedges_qv           = size(bedges_qv, dim=1)
  n_bins_qv             = n_bedges_qv-1
  call write_pdf_fhdr(fhdl_pdf_dqv,bedges_qv,n_bedges_qv, pname)
  call write_pdf_fhdr(fhdl_pdf_bqv,bedges_qv,n_bedges_qv, pname)

  pname                 =  "supersaturation" 
  bedges_ss             = create_bedges_supersat()
  n_bedges_ss           = size(bedges_ss, dim=1)
  n_bins_ss             = n_bedges_ss-1
  call write_pdf_fhdr(fhdl_pdf_dss,bedges_ss,n_bedges_ss, pname)
  call write_pdf_fhdr(fhdl_pdf_bss,bedges_ss,n_bedges_ss, pname)  

  !====bins for prtcl vertical prof
  bedges_pprof            = create_bedges_pprof()
  n_bedges_pprof          = size(bedges_pprof, dim=1)
  n_bins_pprof            = n_bedges_pprof - 1
  call write_prof_fhdr1(fhdl_profprtcl_nsum)
  call write_prof_fhdr1(fhdl_profprtcl_rsum)
  call write_prof_fhdr1(fhdl_profprtcl_r2sum)
  call write_prof_fhdr1(fhdl_profprtcl_r3sum)
  call write_prof_fhdr1(fhdl_profcdplt_nsum)
  call write_prof_fhdr1(fhdl_profcdplt_rsum)
  call write_prof_fhdr1(fhdl_profcdplt_r2sum)
  call write_prof_fhdr1(fhdl_profcdplt_r3sum)
  call write_prof_fhdr1(fhdl_profhdplt_nsum)
  call write_prof_fhdr1(fhdl_profhdplt_rsum)
  call write_prof_fhdr1(fhdl_profhdplt_r2sum)
  call write_prof_fhdr1(fhdl_profhdplt_r3sum)
  call write_prof_fhdr2(fhdl_profinj_nsum)
  call write_prof_fhdr2(fhdl_profactvd_nsum)
  call write_prof_fhdr2(fhdl_profdeactvd_nsum)
 end subroutine
 
 !write file header and bin edges for particle size distribution 
 subroutine write_rdist_fhdr(fhdl,ptype_name)
  implicit none
  character(len=:),allocatable      :: ptype_name
  integer                :: fhdl

  write(fhdl, '(A)')      "#This file has the timestep weighted freq dist of "//trim(ptype_name)//" radius"
  write(fhdl, '(A,I6,A)') "#The radius data is grouped into ", n_bins_r_prtcl, " bins."
  write(fhdl, '(A)')   "#First line of data is bin edges hence has an additional column"
  write(fhdl,'(A,A)') "#","Each col corresponds to radius bin and Each row corresponds to time period in pdf_time.txt" 
  write(fhdl,'(500E12.4)')  bedges_r_prtcl
   
 end subroutine

 !write file header and vertical bin edges for time-weighted particle profile 
 subroutine write_prof_fhdr1(fhdl)
  implicit none

  integer                :: fhdl
  write(fhdl, '(A)')      "#This file has profile of a statistic (statistic name given by the filename)"
  write(fhdl, '(A,I6,A,I6,A)') "#The data is grouped into ", n_bins_pprof, " vertical bins.THe statistic is computed on each bin"
  write(fhdl, '(A)')   "#and weighted by timestep. First line of data is non-dim bin edges in vertical direction. "
  write(fhdl,'(A,A)')  "#","Each col corresponds to a vertical bin. Each row corresponds to time period in prof_time.txt " 
  write(fhdl,'(500E12.4)')  bedges_pprof
  

 end subroutine
 
 !write file header and vertical bin edges for particle profile 
 subroutine write_prof_fhdr2(fhdl)
  implicit none

  integer                :: fhdl
  write(fhdl, '(A)')      "#This file has profile of a statistic (statistic name given by the filename)"
  write(fhdl, '(A,I6,A,I6,A)') "#The data is grouped into ", n_bins_pprof, " vertical bins.THe statistic is computed on each bin"
  write(fhdl, '(A)')        "#First line of data is non-dim bin edges in vertical direction. "
  write(fhdl,'(A,A)')  "#","Each col corresponds to a vertical bin. Each row corresponds to time period in prof_time.txt " 
  write(fhdl,'(500E12.4)')  bedges_pprof
  

 end subroutine

 !write file header and bin edges for pdf
 subroutine write_pdf_fhdr(fhdl,bedges,n_bedges,prop_name)
  implicit none
  character(len=:),allocatable      :: prop_name
  integer                           :: fhdl, n_bedges
  double precision                  :: bedges(n_bedges)     

  write(fhdl, '(A)')      "#This file has the timestep weighted freq dist of "//prop_name
  write(fhdl, '(A,I6,A)') "#The "//prop_name//" data is grouped into ", n_bedges-1, " bins."
  write(fhdl, '(A)')      "#First line of data is bin edges hence has an additional column"
  write(fhdl, '(A)')      "#Each row corresponds to a time period in pdf_time.txt "
  write(fhdl,'(500E12.4)') bedges

 end subroutine
 
 !general function that computes frequency for given bin edges and data
 function bin_data(bedges, n_bedges, data, n_data)

    integer            :: n_bedges, n_data
    double precision   :: bedges(n_bedges), data(n_data)

    !local variables
    integer            :: j, k, freq(n_bedges-1), bin_data(n_bedges-1)
    
    freq               = 0.0d0 
    do  j = 1, n_data
      do k = 1, n_bedges-1
        if ( data(j) >= bedges(k) .and. data(j) < bedges(k+1) )then
          freq(k) = freq(k) + 1
          exit 
        end if

      end do
    end do

    bin_data = freq
 end 
 
 !compute particle properties based on their type (dry aerosol, haze droplet, cloud droplet, intrument measured)
 subroutine classify_prtcl_n_prop(prtcl_type,n_prtcl_type, r_prtcl_type, ql_prtcl_type)

    !currently only count, radius, ql is returned more properties can be added 
    integer,intent(in)            :: prtcl_type
    integer,intent(inout)         :: n_prtcl_type
    double precision,intent(inout):: r_prtcl_type(n_prtcl_curr), ql_prtcl_type(n_prtcl_curr)

    !local variables
    integer                       :: j 
    double precision              :: eps

    r_prtcl_type                  = 0.0d0
    ql_prtcl_type                 = 0.0d0
    n_prtcl_type                  = 0
    eps                           = 1e-10!1e-6 * 1e-4 

    if (prtcl_type == 1) then
      !aerosol particle
      do j = 1, n_prtcl_curr
        if (r_prtcl_curr(j) <= (r_aero_curr(j) + eps)) then
          !1e-10 = 0.0001 micron as margin for numerical error
          n_prtcl_type                = n_prtcl_type+1
          r_prtcl_type(n_prtcl_type)  = r_prtcl_curr(j) 
          ql_prtcl_type(n_prtcl_type) = ql_prtcl_curr(j) 
        endif 
      enddo  

    elseif ( prtcl_type == 2) then
        !haze droplet
      do j = 1, n_prtcl_curr
        if( r_prtcl_curr(j) > (r_aero_curr(j) + eps) .and. r_prtcl_curr(j) < r_crit_curr(j) ) then
        !1e-10 = 0.0001 micron as margin for numerical error
          n_prtcl_type                = n_prtcl_type+1
          r_prtcl_type(n_prtcl_type)  = r_prtcl_curr(j) 
          ql_prtcl_type(n_prtcl_type) = ql_prtcl_curr(j) 
        endif
      enddo

    elseif ( prtcl_type == 3) then 
      !cloud droplet
      do j = 1, n_prtcl_curr
        if (r_prtcl_curr(j) >= r_crit_curr(j)) then
          n_prtcl_type                = n_prtcl_type+1
          r_prtcl_type(n_prtcl_type)  = r_prtcl_curr(j) 
          ql_prtcl_type(n_prtcl_type) = ql_prtcl_curr(j) 
        end if   
      enddo

    elseif (prtcl_type == 4) then  
      do j = 1, n_prtcl_curr
        if (r_prtcl_curr(j) >= r_min_instru)then
          n_prtcl_type                = n_prtcl_type+1
          r_prtcl_type(n_prtcl_type)  = r_prtcl_curr(j) 
          ql_prtcl_type(n_prtcl_type) = ql_prtcl_curr(j) 
        endif
      enddo

    endif  

 end subroutine
 
 !find the array index of particular particle type (ptype) 
 subroutine get_ptype_idx_in_dom(ptype, n_prtcl, r_prtcl, r_aero, r_crit, ptype_count, ptype_idx)

  !I didn't use a variable ptype_curr to store particle type because
  !a single particle can be classified into more than one type 
  !types: dry aero, wet aersol (haze droplet), both wet and dry, cloud droplet, intrument measured)
  !hence maintaining single particle type value in a variable is not effective

  integer,intent(in)            :: ptype, n_prtcl
  double precision,intent(in)   :: r_prtcl(n_prtcl),r_aero(n_prtcl), r_crit(n_prtcl)
  integer,intent(inout)         :: ptype_idx(n_prtcl), ptype_count
  

  !local variables
  integer                       :: j 
  double precision              :: eps

  eps                           = 1e-10!1e-6 * 1e-4 
  ptype_idx                     = 0
  ptype_count                   = 0

  if (ptype == iaero) then
    !aerosol particle
    do j = 1, n_prtcl
      if (r_prtcl(j) <= (r_aero(j) + eps)) then
        !1e-10 = 0.0001 micron as margin for floating point inaccuracies
        ptype_count             = ptype_count+1
        ptype_idx(ptype_count)  = j
      endif 
    enddo  

  elseif ( ptype == ihdplt) then
    !haze droplet
    do j = 1, n_prtcl
      if( r_prtcl(j) > (r_aero(j) + eps) .and. r_prtcl(j) < r_crit(j) ) then
        !1e-10 = 0.0001 micron as margin for numerical error
        ptype_count             = ptype_count+1
        ptype_idx(ptype_count)  = j
      endif
    enddo

  elseif ( ptype == icdplt) then 
    !cloud droplet
    do j = 1, n_prtcl
      if (r_prtcl(j)          >= r_crit(j)) then
        ptype_count             = ptype_count+1
        ptype_idx(ptype_count)  = j
      end if   
    enddo

  elseif (ptype == i_instru) then  
    !particle measured by an instrument
    do j = 1, n_prtcl
      if (r_prtcl(j) >= r_min_instru)then
        ptype_count             = ptype_count+1
        ptype_idx(ptype_count)  = j
      endif
    enddo

  elseif (ptype == iprtcl) then  
    !all particles
    do j = 1, n_prtcl
      ptype_count               = ptype_count+1
      ptype_idx(ptype_count)    = j
    enddo

  else
      write(fhdl_log,'(A)')  "get_prtcl_type_idx - Error: Unknown particle type"
      flush(fhdl_log)
      stop   
  endif  

 end subroutine

 !Compute pdf of scalars values accumulated over time and space
 subroutine compute_pdf(t_ndim, deltat_ndim, T_arr,qv_arr,ss_arr,N,flush_arr)
  implicit none
  integer, intent(in)           :: N
  double precision              :: t_ndim, deltat_ndim
  double precision              :: T_arr(1:N), qv_arr(1:N), ss_arr(1:N)
  integer                       :: flush_arr

  !local variables
  double precision, allocatable :: pdf(:)
  double precision              :: bulk_start, bulk_end, pdf_period, time_dim,  deltat_dim
  integer                       :: bulk_sidx, bulk_eidx, bulk_len, j, debug_lcl 
  integer                       :: sumpdfT, sumpdfqv, sumpdfss, sumpdfr
  integer                       :: aero_count, hdplt_count, cdplt_count, instru_count, ptype_count, ptype_idx(n_prtcl_curr)
  double precision              :: r_aero(n_prtcl_curr), r_hdplt(n_prtcl_curr), r_cdplt(n_prtcl_curr), r_instru(n_prtcl_curr)

  debug_lcl                     = 0 
  time_dim                      = t_ndim * tscale_nu
  deltat_dim                    = deltat_ndim * tscale_nu

  if (do_rec_pdf_n_prof  /= 1) return

  !check the start and end time for recording pdf and prof
  if (time_dim < rec_stime_pdf_n_prof)  then
    if(debug_gbl == 1 .or. debug_lcl == 1)then
      write(fhdl_log, '(A8,A)') "stats:", "time < rec_start_time for pdfs and profiles"
    endif  
    return
  endif

  if (time_dim > rec_etime_pdf_n_prof) then 
      do_rec_pdf_n_prof  = 0
      flush_arr          = 1 
  endif

  bulk_start                    = 0.2
  bulk_end                      = 0.8
  bulk_sidx                     = int(bulk_start * N) + 1
  bulk_eidx                     = int(bulk_end * N)
  bulk_len                      = bulk_eidx - bulk_sidx + 1

  bidx_pdf      = bidx_pdf + 1
  
  if (bidx_pdf == 1)then
    pdf_stime  = t_ndim - deltat_ndim
  endif  

  pdf_deltat_sum=  pdf_deltat_sum + deltat_ndim
  pdf           = bin_data( bedges_T,n_bedges_T,T_arr(1:N),N)
  pdf_dom_T     = pdf_dom_T + pdf * deltat_ndim/N
  pdf           = bin_data( bedges_T,n_bedges_T,T_arr(bulk_sidx:bulk_eidx),bulk_len)
  pdf_bulk_T    = pdf_bulk_T + pdf * deltat_ndim/bulk_len

  pdf           = bin_data( bedges_qv,n_bedges_qv,qv_arr(1:N),N)
  pdf_dom_qv    = pdf_dom_qv + pdf * deltat_ndim/N
  pdf           = bin_data( bedges_qv,n_bedges_qv,qv_arr(bulk_sidx:bulk_eidx),bulk_len)
  pdf_bulk_qv   = pdf_bulk_qv + pdf * deltat_ndim/bulk_len
  
  do j = 1,N
    !to avoid negative zeros affecting pdf 
    if (ss_arr(j) > -1.0e-8 .and. ss_arr(j) < 1.0e-8  ) ss_arr(j) = 0.0d0
  end do
  pdf           = bin_data( bedges_ss,n_bedges_ss,ss_arr(1:N),N)
  pdf_dom_ss    = pdf_dom_ss + pdf * deltat_ndim/N
  pdf           = bin_data( bedges_ss,n_bedges_ss,ss_arr(bulk_sidx:bulk_eidx),bulk_len)
  pdf_bulk_ss   = pdf_bulk_ss + pdf * deltat_ndim/bulk_len


  
  if(n_prtcl_curr > 0) then
    !pdf of diff particle types
    ! call classify_prtcl_n_prop(1,aero_count,r_aero,ql_aero)
    ! call classify_prtcl_n_prop(2,hdplt_count,r_hdplt,ql_hdplt)
    ! call classify_prtcl_n_prop(3,cdplt_count,r_cdplt,ql_cdplt)
    ! call classify_prtcl_n_prop(4,instru_count,r_instru,ql_instru)

    call get_ptype_idx_in_dom(iaero,n_prtcl_curr,r_prtcl_curr,r_aero_curr,r_crit_curr,ptype_count,ptype_idx)
    aero_count   = ptype_count
    r_aero       = r_prtcl_curr(ptype_idx(1:ptype_count))

    call get_ptype_idx_in_dom(ihdplt,n_prtcl_curr,r_prtcl_curr,r_aero_curr,r_crit_curr,ptype_count,ptype_idx)
    hdplt_count   = ptype_count
    r_hdplt       = r_prtcl_curr(ptype_idx(1:ptype_count))

    call get_ptype_idx_in_dom(icdplt,n_prtcl_curr,r_prtcl_curr,r_aero_curr,r_crit_curr,ptype_count,ptype_idx)
    cdplt_count   = ptype_count
    r_cdplt       = r_prtcl_curr(ptype_idx(1:ptype_count))

    call get_ptype_idx_in_dom(iaero,n_prtcl_curr,r_prtcl_curr,r_aero_curr,r_crit_curr,ptype_count,ptype_idx)
    instru_count   = ptype_count
    r_instru       = r_prtcl_curr(ptype_idx(1:ptype_count))
    
    !If particle or dplt  count is zero for a time step; the pdf is zero which is same as skipping the calculation
    !Dividing the pdf by total period not just period with particles
    !will bring down the average. Therfore, the area under the curve for pdfs from time averaging maynot be one
    !if particles are not present during short intervals.
    !However, it will one for typical experiments, when particles are present during entire pdf period
    if(n_prtcl_curr > 0) then
      pdf           = bin_data(bedges_r_prtcl,n_bedges_r_prtcl,r_prtcl_curr,n_prtcl_curr)
      pdf_r_prtcl   = pdf_r_prtcl + pdf * deltat_ndim/n_prtcl_curr
    endif 

    if(aero_count > 0) then
      pdf           = bin_data(bedges_r_prtcl,n_bedges_r_prtcl,r_aero,aero_count)
      pdf_r_aero    = pdf_r_aero + pdf * deltat_ndim/aero_count
    endif  

    if(hdplt_count > 0) then
      pdf           = bin_data(bedges_r_prtcl,n_bedges_r_prtcl,r_hdplt,hdplt_count)
      pdf_r_hdplt   = pdf_r_hdplt + pdf * deltat_ndim/hdplt_count
    endif  

    if(cdplt_count > 0) then
      pdf           = bin_data(bedges_r_prtcl,n_bedges_r_prtcl,r_cdplt,cdplt_count)
      pdf_r_cdplt   = pdf_r_cdplt + pdf * deltat_ndim/cdplt_count
    endif  

    if(instru_count > 0) then
      pdf            = bin_data(bedges_r_prtcl,n_bedges_r_prtcl,r_instru,instru_count)
      pdf_r_instru   = pdf_r_instru + pdf * deltat_ndim/instru_count
    endif 

  endif

  !write pdfs to file
  if (flush_arr == 1 .or. bidx_pdf == n_buff_pdf) then
    
    write(fhdl_log,"(A8,42x,A)") "data:", "************ writing PDFs ****"
    pdf_etime           = t_ndim
    pdf_period         = pdf_etime - pdf_stime 

    pdf_dom_T         = pdf_dom_T/pdf_period
    pdf_bulk_T        = pdf_bulk_T/pdf_period
    pdf_dom_qv        = pdf_dom_qv/pdf_period
    pdf_bulk_qv       = pdf_bulk_qv/pdf_period
    pdf_dom_ss        = pdf_dom_ss/pdf_period
    pdf_bulk_ss       = pdf_bulk_ss/pdf_period

    pdf_r_prtcl       = pdf_r_prtcl/pdf_period
    pdf_r_aero        = pdf_r_aero/pdf_period
    pdf_r_hdplt       = pdf_r_hdplt/pdf_period
    pdf_r_cdplt       = pdf_r_cdplt/pdf_period
    pdf_r_instru      = pdf_r_instru/pdf_period

    if (is_write_fhdr_pdf== 1)then
      write(fhdl_pdf_time,'(A,A11,5A12)')"#","stime_ndim", "etime_ndim", "intvl_ndim", "stime_dim", "etime_dim", "intvl_dim"
      write(fhdl_pdf_time,'(A,A11,5A12)') "#","ndim", "ndim", "ndim", "secs", "secs", "secs"
      is_write_fhdr_pdf = 0
    endif    

    write(fhdl_pdf_time,'(3E12.4,3F12.4)') pdf_stime,           pdf_etime,           pdf_period,&
                                            pdf_stime*tscale_nu, pdf_etime*tscale_nu, pdf_period*tscale_nu
    write(fhdl_pdf_dT,  '(500E12.4)' ) pdf_dom_T
    write(fhdl_pdf_bT,  '(500E12.4)' ) pdf_bulk_T
    write(fhdl_pdf_dqv, '(500E12.4)' ) pdf_dom_qv
    write(fhdl_pdf_bqv, '(500E12.4)' ) pdf_bulk_qv
    write(fhdl_pdf_dss, '(500E12.4)' ) pdf_dom_ss
    write(fhdl_pdf_bss, '(500E12.4)' ) pdf_bulk_ss
    
    write(fhdl_pdf_r_prtcl,   '(500E12.4)' ) pdf_r_prtcl
    write(fhdl_pdf_r_aero,    '(500E12.4)' ) pdf_r_aero
    write(fhdl_pdf_r_hdplt,   '(500E12.4)' ) pdf_r_hdplt
    write(fhdl_pdf_r_cdplt,   '(500E12.4)' ) pdf_r_cdplt
    write(fhdl_pdf_r_instru,  '(500E12.4)' ) pdf_r_instru

    flush(fhdl_pdf_time)
    flush(fhdl_pdf_dT)
    flush(fhdl_pdf_bT)
    flush(fhdl_pdf_dqv)
    flush(fhdl_pdf_bqv)
    flush(fhdl_pdf_dss)
    flush(fhdl_pdf_bss)

    flush(fhdl_pdf_r_prtcl)
    flush(fhdl_pdf_r_aero)
    flush(fhdl_pdf_r_hdplt)
    flush(fhdl_pdf_r_cdplt)
    flush(fhdl_pdf_r_instru)

    !verification
    sumpdfT           = nint(sum(pdf_dom_T))
    sumpdfqv          = nint(sum(pdf_dom_qv))
    sumpdfss          = nint(sum(pdf_dom_ss))
    if ( sumpdfT /= 1 .or. sumpdfqv /= 1  .or. sumpdfss /= 1)then
      write(fhdl_log,'(A)') 'Error: pdf sum of scalars in the domain is not equal to 1'
      write(fhdl_log,'(A,F16.4)') "Sum of dom pdf T", sum(pdf_dom_T)
      write(fhdl_log,'(A,F16.4)') "Sum of dom pdf qv", sum(pdf_dom_qv)
      write(fhdl_log,'(A,F16.4)') "Sum of dom pdf ss", sum(pdf_dom_ss)
      flush(fhdl_log)
      stop
    endif

    sumpdfT           = nint(sum(pdf_bulk_T))
    sumpdfqv          = nint(sum(pdf_bulk_qv))
    sumpdfss          = nint(sum(pdf_bulk_ss))
    if ( sumpdfT /= 1 .or. sumpdfqv /= 1  .or. sumpdfss /= 1)then
      write(fhdl_log,'(A)') 'Error: pdf sum of scalars in the bulk is not equal to 1'
      write(fhdl_log,'(A,F16.4)') "Sum of bulk pdf T", sum(pdf_bulk_T)
      write(fhdl_log,'(A,F16.4)') "Sum of bulk pdf qv", sum(pdf_bulk_qv)
      write(fhdl_log,'(A,F16.4)') "Sum of bulk pdf ss", sum(pdf_bulk_ss)
      flush(fhdl_log)
      stop
    endif

    !cld prtcl is zero for a time step; the pdf is zero
    !it bring down the average; hence sum of pdf will not be one
    sumpdfr            = nint(sum(pdf_r_prtcl))
    if ( sumpdfr /= 1 )then
      write(fhdl_log,'(A)') 'Error: pdf sum of r prtcl is not equal to 1'
      flush(fhdl_log)
      !stop
    endif
    
    pdf_stime         = 0.0d0
    pdf_etime         = 0.0d0
    bidx_pdf          = 0 
    pdf_deltat_sum    = 0.0d0

    pdf_dom_T         = 0.0d0
    pdf_bulk_T        = 0.0d0
    pdf_dom_qv        = 0.0d0
    pdf_bulk_qv       = 0.0d0
    pdf_dom_ss        = 0.0d0
    pdf_bulk_ss       = 0.0d0

    pdf_r_prtcl       = 0.0d0
    pdf_r_aero        = 0.0d0
    pdf_r_hdplt       = 0.0d0
    pdf_r_cdplt       = 0.0d0
    pdf_r_instru      = 0.0d0
  endif


 end subroutine 

 !compute time weighted profile of scalars
 subroutine compute_prof(t_ndim, deltat_ndim, T_arr, qv_arr, ss_arr,N,flush_arr)
  implicit none 
  double precision, intent(in) :: t_ndim, deltat_ndim
  double precision, intent(in) :: T_arr(N), qv_arr(N), ss_arr(N)
  integer                      :: N, flush_arr

  !local variables
  integer                      :: debug_lcl, j
  double precision             :: prof_period, rmean_verf, time_dim, deltat_dim 
  integer                      :: ptype_count_vbin(n_bins_pprof),ptype_count, ptype_idx(n_prtcl_curr)
  double precision             :: rsum_vbin(n_bins_pprof), r2sum_vbin(n_bins_pprof),r3sum_vbin(n_bins_pprof)
  logical                      :: mask(n_prtcl_curr) 

  debug_lcl        = 0 
  time_dim         = t_ndim * tscale_nu
  deltat_dim       = deltat_ndim * tscale_nu
  
  if (do_rec_pdf_n_prof  /= 1) return

  !check the start and end time for recording prtcl hist
  if (time_dim < rec_stime_pdf_n_prof)  then
    if(debug_gbl == 1 .or. debug_lcl == 1)then
      write(fhdl_log, '(A8,A)') "stats:", "time < rec_start_time for pdfs and profiles"
    endif  
    return
  endif

  if (time_dim > rec_etime_pdf_n_prof) then 
    do_rec_pdf_n_prof  = 0
    flush_arr          = 1 
  endif

  bidx_prof        = bidx_prof  + 1    
  if (bidx_prof    == 1) then
    !prof_stime is a module level variable
    prof_stime    = t_ndim - deltat_ndim
  endif  

  prof_Tsum        = prof_Tsum   + T_arr       * deltat_ndim
  prof_T2sum       = prof_T2sum  + T_arr**2.0  * deltat_ndim
  prof_qvsum       = prof_qvsum  + qv_arr      * deltat_ndim
  prof_qv2sum      = prof_qv2sum + qv_arr**2.0 * deltat_ndim
  prof_sssum       = prof_sssum  + ss_arr      * deltat_ndim
  prof_ss2sum      = prof_ss2sum + ss_arr**2.0 * deltat_ndim
  
  !**** compute particle profiles
  if (n_prtcl_curr > 0) then  

    call get_ptype_idx_in_dom(iprtcl,n_prtcl_curr,r_prtcl_curr,r_aero_curr,r_crit_curr,ptype_count,ptype_idx)
    call compute_prtcl_prof(ptype_count,ptype_idx(1:ptype_count),ptype_count_vbin, rsum_vbin,r2sum_vbin,r3sum_vbin)

    !verification
    if  ( int(sum(ptype_count_vbin)) /= n_prtcl_curr) then
      write(fhdl_log,"(A8,A)") "stats:", " error - mismatch between prtcl count from prof and variable"
      stop
    end if
    rmean_verf       = sum(r_prtcl_curr(1:n_prtcl_curr))/n_prtcl_curr
    if  ( abs(rmean_verf - sum(rsum_vbin)/sum(ptype_count_vbin)) > 1e-4) then
      write(fhdl_log,"(A8,A)") "stats:", " error - mismatch between prtcl mean radius from prof and variable"
      stop
    endif   

    if (debug_lcl == 1 .or. debug_gbl == 1) then
      write(fhdl_log,"(A8,A48, I16)")  "stats:", "sum of prtcl count from prof   :", sum(ptype_count_vbin)
      write(fhdl_log,"(A8,A48,E16.4)") "stats:", "mean of prtcl radius from prof :", sum(rsum_vbin)/sum(ptype_count_vbin)
      write(fhdl_log,"(A8,A48, I16)")  "stats:", "prtcl count from variable      :", n_prtcl_curr
      write(fhdl_log,"(A8,A48, E16.4)")"stats:", "mean prtcl radius from variable:", rmean_verf
    endif  

    
    !time weighted sum of prtcl properties
    prof_prtcl_nsum  = prof_prtcl_nsum  + ptype_count_vbin * deltat_ndim
    prof_prtcl_rsum  = prof_prtcl_rsum  + rsum_vbin   * deltat_ndim
    prof_prtcl_r2sum = prof_prtcl_r2sum + r2sum_vbin  * deltat_ndim
    prof_prtcl_r3sum = prof_prtcl_r3sum + r3sum_vbin  * deltat_ndim


    !**** compute cdplt profiles
    call get_ptype_idx_in_dom(icdplt,n_prtcl_curr,r_prtcl_curr,r_aero_curr,r_crit_curr,ptype_count,ptype_idx)
    if ( ptype_count > 0)then
      call compute_prtcl_prof(ptype_count,ptype_idx(1:ptype_count),ptype_count_vbin, rsum_vbin, r2sum_vbin,r3sum_vbin)
      prof_cdplt_nsum  = prof_cdplt_nsum  + ptype_count_vbin * deltat_ndim
      prof_cdplt_rsum  = prof_cdplt_rsum  + rsum_vbin   * deltat_ndim
      prof_cdplt_r2sum = prof_cdplt_r2sum + r2sum_vbin  * deltat_ndim
      prof_cdplt_r3sum = prof_cdplt_r3sum + r3sum_vbin  * deltat_ndim
    endif  

    !**** compute hdplt profiles
    call get_ptype_idx_in_dom(ihdplt,n_prtcl_curr,r_prtcl_curr,r_aero_curr,r_crit_curr,ptype_count,ptype_idx)
    if ( ptype_count > 0)then
      call compute_prtcl_prof(ptype_count,ptype_idx(1:ptype_count),ptype_count_vbin, rsum_vbin, r2sum_vbin,r3sum_vbin)
      prof_hdplt_nsum  = prof_hdplt_nsum  + ptype_count_vbin * deltat_ndim
      prof_hdplt_rsum  = prof_hdplt_rsum  + rsum_vbin   * deltat_ndim
      prof_hdplt_r2sum = prof_hdplt_r2sum + r2sum_vbin  * deltat_ndim
      prof_hdplt_r3sum = prof_hdplt_r3sum + r3sum_vbin  * deltat_ndim
    endif

  endif

  if (n_inj_curr > 0)then
    ptype_count         = 0 
    ptype_idx           = 0
    !assumes that particle injected in this timestep stayed in the domain
    !and did not fall out
    ptype_count         = n_inj_curr
    ptype_idx(1:ptype_count) = (/ (j, j = n_prtcl_curr - ptype_count+1, n_prtcl_curr)  /) 
    call compute_prtcl_prof(ptype_count,ptype_idx(1:ptype_count),ptype_count_vbin, rsum_vbin, r2sum_vbin,r3sum_vbin)
    prof_inj_nsum       = prof_inj_nsum  + ptype_count_vbin
    !Don't need the rsum, r2sum, r3sum profile for now

    !verification 
    if(sum(ptype_count_vbin) /= n_inj_curr) then
      write(fhdl_log,"(A8,A)") "stats:", " error-mismatch between n_inj_dom and sum(n_inj_vertical_bins)"
      stop
    endif   

  endif  

  if (n_actvd_curr > 0 )then
    ptype_count         = 0 
    ptype_idx           = 0
    ptype_count         = n_actvd_curr
    mask                = (n_actvd_curr == 1) 
    ptype_idx(1:ptype_count) = find_idx_true(n_prtcl_curr,mask)
    call compute_prtcl_prof(ptype_count,ptype_idx(1:ptype_count),ptype_count_vbin, rsum_vbin, r2sum_vbin,r3sum_vbin)
    prof_actvd_nsum     = prof_actvd_nsum  + ptype_count_vbin
    !Don't need the rsum, r2sum, r3sum profile for now

    !verification 
    if(sum(ptype_count_vbin) /= n_actvd_curr) then
      write(fhdl_log,"(A8,A)") "stats:", " error-mismatch between n_actvd_dom and sum(n_actvd_vertical_bins)"
      stop
    endif

  end if

  if (n_deactvd_curr > 0 )then
    ptype_count         = 0 
    ptype_idx           = 0
    ptype_count         = n_deactvd_curr
    mask                = (n_actvd_curr == -1) 
    ptype_idx(1:ptype_count) = find_idx_true(n_prtcl_curr,mask)
    call compute_prtcl_prof(ptype_count,ptype_idx(1:ptype_count),ptype_count_vbin, rsum_vbin, r2sum_vbin,r3sum_vbin)
    prof_deactvd_nsum   = prof_deactvd_nsum  + ptype_count_vbin
    !Don't need the rsum, r2sum, r3sum profile for now

    !verification 
    if(sum(ptype_count_vbin) /= n_deactvd_curr) then
      write(fhdl_log,"(A8,A )") "stats:", " error-mismatch between n_deactvd_dom and sum(n_deactvd_vertical_bins)"
      stop
    endif   

  end if
  

  if (bidx_prof == n_buff_prof .or. flush_arr == 1) then

    write(fhdl_log,"(A8,42x,A)") "data:", "************ writing the profiles ****"

    prof_etime      = t_ndim
    prof_period     = prof_etime - prof_stime

    if (is_write_fhdr_prof == 1)then
      write(fhdl_prof_time,'(A,A11,5A12)')"#","stime_ndim", "etime_ndim", "intvl_ndim", "stime_dim", "etime_dim", "intvl_dim"
      write(fhdl_prof_time,'(A,A11,5A12)') "#","ndim", "ndim", "ndim", "secs", "secs", "secs"
      write(fhdl_prof_Tsum, '(A,A,I8,A)') "#"," Each column corresponds to a grid cell from bot to top and there are ", N," gcells"
      write(fhdl_prof_Tsum, '(A,A)') "#"," Each row corresponds to time period in prof_time.txt " 
      write(fhdl_prof_T2sum,'(A,A,I8,A)') "#"," Each column corresponds to a grid cell from bot to top and there are ", N," gcells"
      write(fhdl_prof_T2sum,'(A,A)') "#"," Each row corresponds to time period in prof_time.txt " 
      write(fhdl_prof_qvsum,'(A,A,I8,A)') "#"," Each column corresponds to a grid cell from bot to top and there are ", N," gcells"
      write(fhdl_prof_qvsum,'(A,A)') "#"," Each row corresponds to time period in prof_time.txt " 
      write(fhdl_prof_qv2sum,'(A,A,I8,A)')"#"," Each column corresponds to a grid cell from bot to top and there are ", N," gcells"
      write(fhdl_prof_qv2sum,'(A,A)') "#"," Each row corresponds to time period in prof_time.txt " 
      write(fhdl_prof_sssum,'(A,A,I8,A)') "#"," Each column corresponds to a grid cell from bot to top and there are ", N," gcells"
      write(fhdl_prof_sssum,'(A,A)') "#"," Each row corresponds to time period in prof_time.txt " 
      write(fhdl_prof_ss2sum,'(A,A,I8,A)')"#"," Each column corresponds to a grid cell from bot to top and there are ", N," gcells"
      write(fhdl_prof_ss2sum,'(A,A)') "#"," Each row corresponds to time period in prof_time.txt " 
      is_write_fhdr_prof = 0
    endif
    
    write(fhdl_prof_time,'(3E12.4,3F12.4)') prof_stime,           prof_etime,           prof_period,&
                                            prof_stime*tscale_nu, prof_etime*tscale_nu, prof_period*tscale_nu

    write(fhdl_prof_Tsum,'(10000E12.4)')   prof_Tsum
    write(fhdl_prof_T2sum,'(10000E12.4)')  prof_T2sum
    write(fhdl_prof_qvsum,'(10000E12.4)')  prof_qvsum
    write(fhdl_prof_qv2sum,'(10000E12.4)') prof_qv2sum
    write(fhdl_prof_sssum,'(10000E12.4)')  prof_sssum
    write(fhdl_prof_ss2sum,'(10000E12.4)') prof_ss2sum

    write(fhdl_profprtcl_nsum,'(10000E12.4)') prof_prtcl_nsum
    write(fhdl_profprtcl_rsum,'(10000E12.4)') prof_prtcl_rsum
    write(fhdl_profprtcl_r2sum,'(10000E12.4)') prof_prtcl_r2sum
    write(fhdl_profprtcl_r3sum,'(10000E12.4)') prof_prtcl_r3sum
    write(fhdl_profcdplt_nsum,'(10000E12.4)') prof_cdplt_nsum
    write(fhdl_profcdplt_rsum,'(10000E12.4)') prof_cdplt_rsum
    write(fhdl_profcdplt_r2sum,'(10000E12.4)') prof_cdplt_r2sum
    write(fhdl_profcdplt_r3sum,'(10000E12.4)') prof_cdplt_r3sum
    write(fhdl_profhdplt_nsum,'(10000E12.4)') prof_hdplt_nsum
    write(fhdl_profhdplt_rsum,'(10000E12.4)') prof_hdplt_rsum
    write(fhdl_profhdplt_r2sum,'(10000E12.4)') prof_hdplt_r2sum
    write(fhdl_profhdplt_r3sum,'(10000E12.4)') prof_hdplt_r3sum
    write(fhdl_profinj_nsum,'(10000E12.4)') prof_inj_nsum
    write(fhdl_profactvd_nsum,'(10000E12.4)') prof_actvd_nsum
    write(fhdl_profdeactvd_nsum,'(10000E12.4)') prof_deactvd_nsum

    flush(fhdl_prof_time)
    flush(fhdl_prof_Tsum)
    flush(fhdl_prof_T2sum)
    flush(fhdl_prof_qvsum)
    flush(fhdl_prof_qv2sum)
    flush(fhdl_prof_sssum)
    flush(fhdl_prof_ss2sum)
    flush(fhdl_profprtcl_nsum)
    flush(fhdl_profprtcl_rsum)
    flush(fhdl_profprtcl_r2sum)
    flush(fhdl_profprtcl_r3sum)
    flush(fhdl_profcdplt_nsum)
    flush(fhdl_profcdplt_rsum)
    flush(fhdl_profcdplt_r2sum)
    flush(fhdl_profcdplt_r3sum)
    flush(fhdl_profhdplt_nsum)
    flush(fhdl_profhdplt_rsum)
    flush(fhdl_profhdplt_r2sum)
    flush(fhdl_profhdplt_r3sum)
    flush(fhdl_profinj_nsum)
    flush(fhdl_profactvd_nsum)
    flush(fhdl_profdeactvd_nsum)

    !clear array
    bidx_prof        = 0
    prof_stime       = 0.0d0
    prof_etime       = 0.0d0

    prof_Tsum        = 0.0d0
    prof_T2sum       = 0.0d0
    prof_qvsum       = 0.0d0
    prof_qv2sum      = 0.0d0
    prof_sssum       = 0.0d0
    prof_ss2sum      = 0.0d0
    prof_prtcl_nsum  = 0.0d0
    prof_prtcl_rsum  = 0.0d0
    prof_prtcl_r2sum = 0.0d0
    prof_prtcl_r3sum = 0.0d0
    prof_cdplt_nsum  = 0.0d0
    prof_cdplt_rsum  = 0.0d0
    prof_cdplt_r2sum = 0.0d0
    prof_cdplt_r3sum = 0.0d0
    prof_hdplt_nsum  = 0.0d0
    prof_hdplt_rsum  = 0.0d0
    prof_hdplt_r2sum = 0.0d0
    prof_hdplt_r3sum = 0.0d0
    prof_inj_nsum    = 0.0d0
    prof_actvd_nsum  = 0.0d0
    prof_deactvd_nsum= 0.0d0

  endif

 end subroutine  

 !compute particle profile using veritcal bins
 subroutine compute_prtcl_prof(ptype_count,ptype_idx,ptype_count_vbin,rsum_vbin,r2sum_vbin,r3sum_vbin)
  implicit none
  !Save the profile n_prtcl, rmean, r2mean for all particles

  double precision,intent(inout) :: rsum_vbin(n_bins_pprof), r2sum_vbin(n_bins_pprof), r3sum_vbin(n_bins_pprof)
  integer,intent(inout)          :: ptype_count_vbin(n_bins_pprof)
  integer,intent(in)             ::  ptype_count, ptype_idx(1:ptype_count)
  !local variables
  integer                        :: ptype_idx_vbin(n_prtcl_curr, n_bedges_pprof)
  integer                        :: j, pcount_vbin

  ptype_count_vbin     = 0
  rsum_vbin            = 0.0d0
  r2sum_vbin           = 0.0d0


  call get_ptype_idx_in_vbin(ptype_count,ptype_idx,ptype_count_vbin, ptype_idx_vbin)
  
  do j = 1, n_bins_pprof
      pcount_vbin      = ptype_count_vbin(j)
    if (pcount_vbin > 0) then
      rsum_vbin(j)  = sum(r_prtcl_curr(ptype_idx_vbin(1:pcount_vbin,j)))
      r2sum_vbin(j) = sum(r_prtcl_curr(ptype_idx_vbin(1:pcount_vbin,j))**2.0d0)
      r3sum_vbin(j) = sum(r_prtcl_curr(ptype_idx_vbin(1:pcount_vbin,j))**3.0d0)

    endif  

  end do  
  
 end subroutine
 
 !find array idx of a particle type within a vertical bin
 subroutine get_ptype_idx_in_vbin(ptype_count,ptype_idx,ptype_count_vbin, ptype_idx_vbin)
  implicit none

  integer, intent(inout)   :: ptype_count_vbin(n_bins_pprof), ptype_idx_vbin(n_prtcl_curr, n_bins_pprof)
  integer, intent(in)      :: ptype_count
  integer,intent(in)       :: ptype_idx(1:ptype_count)
  
  !local variables
  integer                  :: j, k, idx
  double precision         :: x_ptype(n_prtcl_curr)
  ptype_idx_vbin           = 0 
  ptype_count_vbin         = 0
  x_ptype                  = 0
  !get x pos values only for this particle types through the index
  x_ptype                  = x_prtcl_curr(ptype_idx(1:ptype_count))


  do k = 1, n_bins_pprof
    idx                    = 0

    do j = 1, ptype_count 
      if (x_ptype(j) >= bedges_pprof(k) .and. x_ptype(j) < bedges_pprof(k+1)) then
        idx                    = idx + 1
        ptype_idx_vbin(idx,k)  = ptype_idx(j)
      end if
    end do
    
    ptype_count_vbin(k)        = idx 
  end do  

 end subroutine 

 !List of index of array elements with value <> 0 
 function find_idx_true(n,arr)
   integer          :: n
   logical          :: arr(n)

   

   !local variables
   integer              :: j, count
   integer,allocatable  :: find_idx_true(:)
   integer              :: tidx(n) 

   count            = 0 
   do j             = 1, n
     if(arr(j)) then
      count         = count + 1
      tidx(count)   = j
     endif 
   end do
   
   allocate(find_idx_true(count))
   find_idx_true(1:count) = tidx(1:count)
   
   return
 end
 
 !evolution of particle properties - statistic computed from all particles in the domain at each time step
 subroutine evol_prtcl_prop(t_ndim, deltat_ndim, flush_arr)
    use const, only             : pi43, rho_w
    implicit none

    integer, intent(in)         :: flush_arr
    double precision, intent(in)::  t_ndim, deltat_ndim

    !local variables
    integer                     :: debug_lcl, j, k
    integer                     :: ptype_count, fhdl_evol_ptype(n_ptype)
    integer,dimension(n_prtcl_curr) ::  ptype_idx
    integer, allocatable        :: pidx(:)
    double precision            :: conv2cm3


    debug_lcl                   = 0
    conv2cm3                    = inv_vol_dom * 1e-6 
    fhdl_evol_ptype(1)          = fhdl_evol_aero
    fhdl_evol_ptype(2)          = fhdl_evol_hdplt
    fhdl_evol_ptype(3)          = fhdl_evol_cdplt
    fhdl_evol_ptype(4)          = fhdl_evol_instru
    fhdl_evol_ptype(5)          = fhdl_evol_prtcl

    bidx_prtcl_evol              = bidx_prtcl_evol + 1
    if (n_prtcl_curr > 0) then 

      !time and n_prtcl_alltime 
      time_prtcl_evol(bidx_prtcl_evol)        = t_ndim
      deltat_prtcl_evol(bidx_prtcl_evol)      = deltat_ndim
      n_inj_evol(bidx_prtcl_evol)             = n_inj_curr
      !this is same as n_inj_alltime
      n_prtcl_alltime_evol(bidx_prtcl_evol)   = n_prtcl_alltime
      n_fallout_alltime_evol(bidx_prtcl_evol) = n_fallout_alltime
      n_actvd_evol(bidx_prtcl_evol)           = n_actvd_curr
      n_deactvd_evol(bidx_prtcl_evol)         = n_deactvd_curr


      !Stats of particles currently in the domain by particle type
      do j = 1, n_ptype
        call get_ptype_idx_in_dom(j,n_prtcl_curr,r_prtcl_curr,r_aero_curr,r_crit_curr,ptype_count,ptype_idx)

        if (ptype_count > 0) then
          pidx  = ptype_idx(1:ptype_count)
          !write(fhdl_log,'(A48,I16)') "index array size:", size(pidx)
          n_ptype_evol(j,bidx_prtcl_evol)     = ptype_count
          rmin_ptype_evol(j,bidx_prtcl_evol)  = minval(r_prtcl_curr(pidx))
          rmax_ptype_evol(j,bidx_prtcl_evol)  = maxval(r_prtcl_curr(pidx))
          rmean_ptype_evol(j,bidx_prtcl_evol) = sum(r_prtcl_curr(pidx))/ptype_count
          r2mean_ptype_evol(j,bidx_prtcl_evol) = sum(r_prtcl_curr(pidx)**2.0d0)/ptype_count

          !r3mean can be obtained from ql
          if (j == 1) then
            !for dry aerosol
            ql_ptype_evol(j,bidx_prtcl_evol)  = 0.0d0
          else
            ql_ptype_evol(j,bidx_prtcl_evol)  = sum(ql_prtcl_curr(pidx))
          endif  

        endif 

      end do  

      !Stats of particles that fellout by particle type
      if (n_fallout_curr > 0) then
        do j = 1, n_ptype
          call get_ptype_idx_in_dom(j,n_fallout_curr,r_prtcl_fallout,r_aero_fallout,r_crit_fallout,ptype_count,ptype_idx)
  
          if (ptype_count > 0) then
            pidx  = ptype_idx(1:ptype_count)
            !write(fhdl_log,'(A48,I16)') "index array size:", size(pidx)
            n_fall_ptype_evol(j,bidx_prtcl_evol)     = ptype_count
            rmean_fall_ptype_evol(j,bidx_prtcl_evol) = sum(r_prtcl_fallout(pidx))/ptype_count
            r2mean_fall_ptype_evol(j,bidx_prtcl_evol)= sum(r_prtcl_fallout(pidx)**2.0d0)/ptype_count
            ql_fall_ptype_evol(j,bidx_prtcl_evol)    =  sum(ql_prtcl_fallout(pidx))/ptype_count
          endif 
        end do
      endif    


    else
      time_prtcl_evol(bidx_prtcl_evol)    = t_ndim
      deltat_prtcl_evol(bidx_prtcl_evol)  = deltat_ndim
      !all buffer array are set to 0        
      !zeros will be written for those time period without any particles in the domain
    endif

    if (flush_arr == 1 .or. bidx_prtcl_evol == n_buff_evol)then
      write(fhdl_log,"(A8,42x,A)") "data:", "************ writing evolution data ****"

      !write headers with column name and units one-time only
      if (is_fhdr_evol_prtcl  == 1 )then
        write(fhdl_evol_mphy,'(A,A15,22A16)') "#","time_nondim", "time_dim", "deltat_dim",& 
        "n_inj_alltime","n_inj",        "n_fall_alltime", "n_fall",     "n_fall_cdplt", "n_fall_hdplt", &
        "n_prtcl",      "n_cdplt",      "n_hdplt",       "n_actvd",     "n_deactvd",    "rmean_prtcl", &
        "rmean_cdplt",  "rmean_hdplt",  "r2mean_prtcl","r2mean_cdplt",  "r2mean_hdplt", "ql_prtcl", &
        "ql_cdplt",     "ql_hdplt" 
  
        write(fhdl_evol_mphy,'(A,A15,22A16)') "#","ndim", "secs", "secs",& 
               "#/cm3",   "#/cm3", "#/cm3", "#/cm3",  "#/cm3",  "#/cm3",  "#/cm3", "#/cm3",&
               "#/cm3",   "#/cm3", "#/cm3", "m",      "m",      "m",      "m2",    "m2",&
               "m2",      "kg/m3", "kg/m3", "kg/m3"      
        
        write(fhdl_evol_fallout,'(A,A15,18A16)') "#","time_nondim", "time_dim", "deltat_dim",& 
               "n_fall",      "n_fall_cdplt",     "n_fall_hdplt",     "n_fall_aero", &
               "rmean_fall",  "rmean_fall_cdp",   "rmean_fall_hdp", "rmean_fall_aero",&   
               "r2mean_fall", "r2mean_fall_cdp",  "r2mean_fall_hdp","r2mean_fall_aero",&
               "ql_fall",     "ql_fall_cdplt",    "ql_fall_hdplt",    "ql_fall_aero"

        write(fhdl_evol_fallout,'(A,A15,18A16)') "#","ndim", "secs", "secs",& 
               "#/cm3", "#/cm3", "#/cm3", "#/cm3", &
               "m",     "m",      "m",    "m", &
               "m2",    "m2",     "m2",   "m2",&
               "kg/m3", "kg/m3",  "kg/m3","kg/m3" 
  
        do k = 1, n_ptype
          write(fhdl_evol_ptype(k),'(A,A15,8A16)') "#","time_nondim", "time_dim", "deltat_dim",& 
          "num conc", "rmin", "rmax", "rmean", "r2mean", "ql"       
  
          write(fhdl_evol_ptype(k),'(A,A15,8A16)') "#","ndim", "secs", "secs",& 
          "#/cm3", "m", "m", "m", "m2", "kg/m3"       
        end do
  
          is_fhdr_evol_prtcl  = 0
      endif  

      !conv to number/cm3
      n_ptype_evol           =  n_ptype_evol * conv2cm3
      n_fall_ptype_evol      =  n_fall_ptype_evol * conv2cm3

      !write data from buffers to file
      do j = 1, bidx_prtcl_evol

        write(fhdl_evol_mphy,'(E16.8,2F16.8,11F16.4,9E16.8)') time_prtcl_evol(j), &
              time_prtcl_evol(j)*tscale_nu,       deltat_prtcl_evol(j)*tscale_nu, & 
              n_prtcl_alltime_evol(j)* conv2cm3,  n_inj_evol(j)* conv2cm3,          n_fallout_alltime_evol(j)* conv2cm3,&
              n_fall_ptype_evol(iprtcl,j),        n_fall_ptype_evol(icdplt,j),      n_fall_ptype_evol(ihdplt,j), &
              n_ptype_evol(iprtcl,j),             n_ptype_evol(icdplt,j),           n_ptype_evol(ihdplt,j), &
              n_actvd_evol(j)* conv2cm3,          n_deactvd_evol(j)*conv2cm3, &
              rmean_ptype_evol(iprtcl,j),         rmean_ptype_evol(icdplt,j),       rmean_ptype_evol(ihdplt,j),  &
              r2mean_ptype_evol(iprtcl,j),        r2mean_ptype_evol(icdplt,j),      r2mean_ptype_evol(ihdplt,j),  &
              ql_ptype_evol(iprtcl,j),            ql_ptype_evol(icdplt,j),          ql_ptype_evol(ihdplt,j)
        
        write(fhdl_evol_fallout,'(E16.8,2F16.8,4F16.4,12E16.8)') time_prtcl_evol(j),& 
        time_prtcl_evol(j)*tscale_nu,     deltat_prtcl_evol(j)*tscale_nu, &     
        n_fall_ptype_evol(iprtcl,j),      n_fall_ptype_evol(icdplt,j),&
        n_fall_ptype_evol(ihdplt,j),      n_fall_ptype_evol(iaero,j),&
        rmean_fall_ptype_evol(iprtcl,j),  rmean_fall_ptype_evol(icdplt,j),&
        rmean_fall_ptype_evol(ihdplt,j),  rmean_fall_ptype_evol(iaero,j),&
        r2mean_fall_ptype_evol(iprtcl,j), r2mean_fall_ptype_evol(icdplt,j),&
        r2mean_fall_ptype_evol(ihdplt,j), r2mean_fall_ptype_evol(iaero,j),&
        ql_fall_ptype_evol(iprtcl,j),     ql_fall_ptype_evol(icdplt,j),&
        ql_fall_ptype_evol(ihdplt,j),     ql_fall_ptype_evol(iaero,j)


        do k = 1, n_ptype
          write(fhdl_evol_ptype(k),'(E16.8,2F16.8,F16.4,5E16.8)') time_prtcl_evol(j),& 
          time_prtcl_evol(j)*tscale_nu, deltat_prtcl_evol(j)*tscale_nu, & 
          n_ptype_evol(k,j),            rmin_ptype_evol(k,j),   rmax_ptype_evol(k,j),&
          rmean_ptype_evol(k,j),        r2mean_ptype_evol(k,j), ql_ptype_evol(k,j)
        end do

      end do

      !clear the arrays
      bidx_prtcl_evol                    = 0
      time_prtcl_evol                    = 0.0d0
      deltat_prtcl_evol                  = 0.0d0
      n_inj_evol                         = 0
      n_prtcl_alltime_evol               = 0
      n_fallout_alltime_evol             = 0
      n_actvd_evol                       = 0
      n_deactvd_evol                     = 0

      n_ptype_evol                       = 0
      rmin_ptype_evol                    = 0.0d0
      rmax_ptype_evol                    = 0.0d0
      rmean_ptype_evol                   = 0.0d0
      r2mean_ptype_evol                  = 0.0d0
      ql_ptype_evol                      = 0.0d0

      n_fall_ptype_evol                 = 0      
      rmean_fall_ptype_evol             = 0.0d0
      r2mean_fall_ptype_evol            = 0.0d0
      ql_fall_ptype_evol                = 0.0d0

      !flush the file buffer
      flush(fhdl_evol_mphy)

    endif

 end subroutine

 !save particle injection information
 subroutine save_prtcl_injection(flush_arr)
  integer                   :: j, flush_arr

  !write inject rate data: inject_time,   inject_prtcl_count
  
  if (flush_arr == 1 .or. bidx_inj == n_buff_evol)then
    write(fhdl_log,"(A8,42x,A)") "data:", "************ Saving particle injection ****"
    !write to the file
    do j = 1, bidx_inj
      write(fhdl_inject,'(E12.4,I12)') inject_time(j), inject_prtcl_count(j)
    end do

    !reset the arrays
    bidx_inj          = 0
    inject_time       = 0.0d0
    inject_prtcl_count = 0
    
    !flush the file buffer
    flush(fhdl_inject)

  endif

  return

 end subroutine

 !master subroutine that calls other stats subprograms
 subroutine compute_microphy_stats(T_arr, qv_arr, ss_arr, N, t_ndim, deltat_ndim, flush_arr)
  implicit none
  integer, intent(in)               :: N, flush_arr
  double precision, intent(in)      :: T_arr(N), qv_arr(N), ss_arr(N), t_ndim, deltat_ndim

  !local variables
  integer                           :: debug_lcl

  debug_lcl   = 0
  if (debug_gbl == 1 .or. debug_lcl ==1)then
    write(fhdl_log, "(A8,42x,A)")    "stats:", "************ Calc microphysics stats ****"
    flush(fhdl_log)
  endif
  call evol_prtcl_prop(t_ndim,deltat_ndim,flush_arr)
  call compute_pdf(t_ndim,deltat_ndim,T_arr, qv_arr, ss_arr, N, flush_arr);
  call compute_prof(t_ndim,deltat_ndim,T_arr, qv_arr, ss_arr, N, flush_arr);
  call rec_prtcl_hist(T_arr, qv_arr,ss_arr,N, t_ndim, deltat_ndim, flush_arr);
  call save_prtcl_injection(flush_arr)

 end subroutine

 !Tally particle count - just a verification to check if the code is working properly
 subroutine verify_prtcl_count()

  !There are multiple variables that hold the particle count
  !This subprogram will verify if they tally
  implicit none
  double precision    temp

  write(fhdl_log,'(A8,42x,A)') "Verf:", "************ Verifying alltime particle count"
  write(fhdl_log,"(A8,A48,I16)") "Verf:", "n_prtcl_alltime              :", n_prtcl_alltime
  if (n_prtcl_alltime /= prtcl_counter) then
    write(fhdl_log,'(A8,A)') "Verf:", "Error: n_prtcl_alltime and prtcl_counter didn't match"      
  end if

  temp            =  n_inj_alltime
  write(fhdl_log,"(A8,A48,F16.1)") "Verf:", " inject_prtcl_alltime     :", temp
  if (n_prtcl_alltime /= temp) then
    write(fhdl_log,'(A8,A)') "Verf:", "Error: n_prtcl_alltime and n_inj_alltime didn't match"
  end if
  
  temp            = sum(dis_n_alltime)
  write(fhdl_log,"(A8,A48,F16.1)") "Verf:", "sum(dis_n_alltime)        :", temp
  if (n_prtcl_alltime /= temp) then
    write(fhdl_log,'(A8,A)') "Verf:", "Error: n_prtcl_alltime and sum(dis_n_alltime) didn't match"
  end if
  
  temp            = n_prtcl_curr + n_fallout_alltime
  write(fhdl_log,"(A8,A48,I16)") "Verf:", "n_prtcl_curr:", n_prtcl_curr
  write(fhdl_log,"(A8,A48,I16)") "Verf:", "n_fallout_alltime:", n_fallout_alltime
  write(fhdl_log,"(A8,A48,F16.1)") "Verf:", "n_prtcl_curr + n_fallout_alltime :", temp
  if (n_prtcl_alltime /= temp) then
    write(fhdl_log,'(A8,A)') "Verf:", "Error: n_prtcl_alltime and (n_prtcl_curr + n_fallout) didn't match"
  end if
  
  flush(fhdl_log)
 end subroutine

 !verify, close files, and cleanup
 subroutine final_step_microphy()

   call verify_prtcl_count()
   call close_files()

 end subroutine 
 
end module
  