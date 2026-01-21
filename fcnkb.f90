! -------------------------------------------------------------------------------------------------
! This subroutine calculates the derivitives of ndmax different droplet sizes and of the thermo-
! dynamic variables
! -------------------------------------------------------------------------------------------------
!
! Date:        January 2016
!
!
! Updates:
! --------
!
! January 25th, 2016
! -- the old FORTRAN77 code was rewritten in FORTRAN90
! -- full control over code was gained by declaring each variable
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------

!Mani:
!This function provides derivates for y variables in drop_radius at time t
!input:     t               -> time
!input:     drop_radius     -> y variables
!output:    drdt            -> time derivative of all y variables
SUBROUTINE fcnkb(t,drop_radius,drdt)

    USE const
    USE array

    IMPLICIT NONE

    REAL*8  :: drop_radius(8), drdt(8)
    REAL*8  :: ck, cr, denom
    REAL*8  :: e, es, ew
    REAL*8  :: falpha, fbeta, rho, rhol, sw, t

    ! -- molecular weight of the solute (NaCl, (NH4)2SO4, ...)
    REAL*8  :: Ms

    ! --
    REAL*8  :: radius, qv, temp, s, press, ver_vel, height, ql, nions
    !        REAL*8  :: cpm

    REAL*8  :: lalpha, lbeta
    REAL*8, PARAMETER :: alpha = 1
    REAL*8, PARAMETER :: beta = 0.04

    REAL*8  :: c7
    REAL*8, PARAMETER :: c7am   = 0.4363021
    REAL*8, PARAMETER :: c7nacl = 0.5381062
    REAL*8, PARAMETER :: sigma  = 7.392730e-2

    !Mani: added to array module
    !if (err_dgm == -1) return
    nions   = 0.0d0

    radius  = drop_radius(1)
    qv      = drop_radius(2)
    temp    = drop_radius(3)
    s       = drop_radius(4)
    press   = drop_radius(5)
    ver_vel = drop_radius(6)
    height  = drop_radius(7)
    ql      = drop_radius(8)

    !Mani: constanst used from module const
    cpm     = cp*((1.0+cpv/cp*qv)/(1.0+qv))

    ! -- initialize diffusivities D, Ktemp and Lv

    ! -- Equation 2 from Bolton, MWR (1980)
    Lv      = (2.501-0.00237*(temp-273.15))*1.e6

    ! -- Table 7.1 (p103): linear interpolation between -10 and 30 deg C
    ! -- Rogers and Yau (1989)
    Ktemp   = 7.7e-5*(temp-273.15)+0.02399
    D       = 1.57e-7*(temp-273.15)+2.211e-5
    D       = D*1.e5/press

    es      = ew(temp)*100.0

    ! -- calculate lalpha and lbeta
    !Mani: I found that sqrt of negative temp produces floating point error, so I am catching it before occurs
    if (temp < 0) then
    !if (temp < 273 .or. temp > 320) then
      write(*,*) ""
      write(*,'(3(A,E16.8))') "fcnb error: T: ", temp, " qv: ", qv,  " radius: " , radius
      write(*,*) "fcnb error: Ktemp: ", Ktemp
      write(*,*) "fcnb error: 2.0*pi*Md*R_uni : ", 2.0*pi*Md*R_uni
      flush(6)
      err_dgm   = -1
      !return
      ! write(*,*) "fcnb error: SQRT(2.0*pi*Md*R_uni*temp): ", SQRT(2.0*pi*Md*R_uni*temp)
      ! flush(6)
    end if

    !Mani: If the temperature oscillates to a negative value, sqrt() will return NaN and throws a floating point 
    !exception. If compiled with Floating point flags on, the exception will terminate the program. 
    !This is difficult to debug. therefore compile without floating point flags, and handle it in the program.
    !So you know, where the issue occurs.

    lalpha  = Ktemp * SQRT(2.0*pi*Md*R_uni*temp)/(alpha*press*(cv+R_uni/2.0))
    lbeta   = SQRT(2.0*pi*Mw/(Rv*temp))*D/beta

    !Mani:dmaxa is used from the module array
    !so set the value before calling odeint
    IF (dmaxa .EQ. 1) THEN
      !sodium cholride NaCl
       Ms = 58.4428e-3
       c7 = c7nacl
       nions = 2.0d0
    ELSE IF(dmaxa .EQ. 2) THEN
      !Ammonium sulphate
       Ms = 132.1395e-3
       !density of Ammonium bisulphate and sulphate are similar so using the same density factor
       c7 = c7am
       nions = 3.0d0
    ELSE IF(dmaxa .EQ. 3) THEN
      !Ammonium bisulphate
        Ms = 115.11e-3
        !density of Ammonium bisulphate and sulphate are similar so using the same density factor
        c7 = c7am
        nions = 2.0d0
    ELSE  
      write(*,*) "fcnkb: Error: aerosol type unknown"
      stop
    END IF


    falpha  = radius/(radius+lalpha)
    fbeta   = radius/(radius+lbeta)
    rhol    = (radius**3*pi43*rho_w+solute_mass*c7)/(radius**3*pi43)

    ! -- calculate drdt
    ck      = (2*sigma/(Rv*temp*rhol*radius))
    cr      = nions*(Mw/Ms)*solute_mass/(pi43*radius**3*rhol-solute_mass)
    denom   = rhol*(Rv*temp/(fbeta*D*es)+Lv**2/(falpha*Ktemp*Rv*temp**2))
    drdt(1) = 1.0/radius*(s-ck+cr)/denom

    ! -- calculate change of water vapor
    drdt(2) = -pi4*grid_scale*radius**2*drdt(1)*rho_w
    ! -- calculate change of temperature
    drdt(3) = -Lv/cpm*drdt(2)   !-ver_vel*g/cpm
    ! -- calculate change of supersaturation
    e       = qv*press/(eps+qv)
    rho     = (press-e)/(Rd*temp)*(1.0+ql)+e/(Rv*temp)
    drdt(4) = (1.0+s)*(Rd/(qv*(Rv*qv+Rd))*drdt(2)-g*rho*ver_vel/press-Lv/(Rv*temp**2)*drdt(3))
    ! -- calculate change of pressure
    drdt(5) = -1.0*rho*g*ver_vel
    ! -- calculate change of vertical velocity
    drdt(6) = 0.0
    ! -- calculate change of height
    drdt(7) = ver_vel
    ! -- calculate change of liquid water
    drdt(8) = -1.0*drdt(2)

RETURN

END
