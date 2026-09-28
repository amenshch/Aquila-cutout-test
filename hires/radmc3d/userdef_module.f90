module code_id
!__________________________________________________________________________________________________________________________________
!
! Identifier of the code version (and compilation date).
!__________________________________________________________________________________________________________________________________
!
  implicit none
!___________________________________________________
!
  character(8), public :: version = "2.2606012" !<<< !*** DON'T FORGET TO CHANGE THE VERSION ***!
!___________________________________________________
!
  character( 5), parameter, public :: cfpptime = __TIME__
  character(21), parameter, public :: compda = __DATE__//' '//cfpptime
  character(41), parameter, public :: codeID = "RADMC3D compiled on "//compda
  character(78), parameter, public :: author = "Alexander Men’shchikov, DAp IRFU CEA Saclay, France"

end module code_id

module userdef_module
use code_id
use amr_module
use amrray_module
use rtglobal_module
use constants_module
use dust_module
use lines_module
use stars_module
use quantum_module
use montecarlo_module
use namelist_module
!__________________________________________________________________________________________________________________________________
!
! Here you can define your own variables, arrays etc. 
!__________________________________________________________________________________________________________________________________
!
 character*32 :: tstamp, cexec_time_sec, cexec_time_min, cexec_time_hrs
 character*100, parameter :: outname = '+mc3d.out', outputname = '+mc3d.output'
 integer :: userdef_examplevariable, userdef_nzones, userdef_nzembc, userdef_doback, userdef_tmassav, userdef_external
 double precision, parameter :: xRsun = 6.96d10, xMsun = 1.98892d33, xLsun = 3.828d33, xJy = 1.0d-23, AtomicMU = 1.66054d-24 &
                              , muH2 = 2.8d0, GasConst = 8.314510d+07, GravC = 6.67259d-8, SBconst = 5.6704d-5 &
                              , xpc = 3.0856775814913673d18, xAU = 1.495978707d13
 double precision, parameter :: z0 = 0.0d0, z1 = 1.0d0, z2 = 2.0d0, z3 = 3.0d0, z4 = 4.0d0, z5 = 5.0d0, z6 = 6.0d0, z7 = 7.0d0 &
                              , z8 = 8.0d0, z9 = 9.0d0, z10 = 10.0d0, z05 = 0.5d0, z01 = 0.1d0, z09 = 0.9d0, z001 = 0.01d0 &
                              , zfloor30 = 1.0d-30, zfloor99 = 1.0d-99, zeps3 =  1.0d-03, zeps5 = 1.0d-05, zeps10 = 1.0d-10 &
                              , zeps14 = 1.0d-14, zeps20 = 1.0d-20, z025 = 0.25d0, almostzero = 1.0d-31
 double precision, parameter :: z24 = 24.0d0, z60 = 60.0d0, z90 = 90.0d0, z100 = 100.0d0, z180 = 180.0d0, z360 = 360.0d0 &
                              , z3600 = 3600.0d0, z1000 = 1000.0d0, z1440 = 1440.0d0, z86400 = 86400.0d0
 double precision :: userdef_Rinner, userdef_Rplawx, userdef_Rbesph, userdef_Rgauss, userdef_gfwhm, userdef_Router, userdef_fgeom &
                   , userdef_dplawx, userdef_dplout, userdef_rhofac, userdef_rhoplx, userdef_dust2gas, userdef_distpc &
                   , userdef_masstot, userdef_bestem, kappa300, userdef_besden, userdef_hbbfac, userdef_temstar &
                   , userdef_lumstar, userdef_fracaper, userdef_fraczero, userdef_mcrscale, userdef_radscale, userdef_denscale &
                   , userdef_gplaw, userdef_embden, userdef_ximax, userdef_ftrans, gplaw & !!, userdef_opascale
                   , Rmodel, Router, Rinner, Rbesph, Rgauss, Rplawx, gfwhm, dr1, sf, sfmz1, masstot, exec_time_sec, exec_time_min &
                   , exec_time_hrs, aperture(1000)
 contains

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_defaults ()
!__________________________________________________________________________________________________________________________________
!
! This subroutine allows you to specify defaults values for your variables.
!
! WARNING: In this subroutine you are not allowed to use write(stdo,*), because the stdo is not yet set. Reason: The defaults are 
!          set before the command-line options are interpreted (so that the defaults can be overwritten by command-line options),
!          but the stdo depends on whether the user calls RADMC-3D as a child or not, which is given by the command-line options.
!__________________________________________________________________________________________________________________________________
!
  implicit none
!__________________________________________________________________________________________________________________________________
!
  userdef_nzones = 100
  userdef_nzembc = 15
  userdef_fgeom = 1.0d0
  userdef_ftrans = 1.3d0
  userdef_Rinner = 0.1d0
  userdef_Rbesph = 1.0d4
  userdef_Rgauss = 1.0d4
  userdef_gfwhm = 1.0d2
  userdef_gplaw = 2.0d0
  userdef_Rplawx = 1.0d4
  userdef_Router = 3.0d4
  userdef_fracaper = 0.0d0
  userdef_fraczero = 1.0d0
  userdef_doback = 0
  userdef_tmassav = 0
  userdef_bestem = 7.0d0
  userdef_besden = 5.19d-18
  userdef_dplawx = 2.0d0
  userdef_dplout = 0.0d0
  userdef_rhofac = 1.0d0
  userdef_embden = 3.1084e-21
  userdef_ximax = 6.451d0
  userdef_dust2gas = 0.01d0
!!  userdef_opascale = 1.0d0
  userdef_distpc = 140.0d0
  userdef_external = 1
  userdef_hbbfac = 6.19E-16
  userdef_temstar = 5778.0d0
  userdef_masstot = 0.1d0
  userdef_lumstar = xLsun
  userdef_rhoplx = 1.457d-23
  userdef_mcrscale = 1.0d0
  userdef_radscale = 1.0d0
  userdef_denscale = 1.0d0
  stdo = 1111
  open ( unit=stdo, file=trim(outputname), access='append' )
  tstamp = timestmp ()
  exec_time_sec = timer ( 0.0d0 )
!__________________________________________________________________________________________________________________________________
!
end subroutine userdef_defaults

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_commandline ( buffer, numarg, iarg, fromstdi, gotit )
!__________________________________________________________________________________________________________________________________
!
! Here you can interpret your own command-line options.
!__________________________________________________________________________________________________________________________________
!
  implicit none
  character*100 :: buffer
  integer :: iarg,numarg
  logical :: gotit,fromstdi
!__________________________________________________________________________________________________________________________________
!
! Below is just an example. You can delete or replace all of it.
!__________________________________________________________________________________________________________________________________
!
  if (buffer(1:6) .eq. 'doback') then
    userdef_doback = 1
    gotit = .true.
  else if (buffer(1:7) .eq. 'tmassav') then
    userdef_tmassav = 1
    gotit = .true.
  else if (buffer(1:21) .eq. 'examplecomlinargument') then
     if (iarg .gt. numarg) then
        write (stdo,'(a)') ' USERDEF_COMMANDLINE:'
        write (stdo,'(a)') '   ERROR: Command line options: cannot read ''sizeau''.'
        write (stdo,'(a)') '          Expecting 1 integer after ''examplecomlinargument.'''
        stop
     endif
     call ggetarg ( iarg, buffer, fromstdi )
     iarg = iarg+1
     read (buffer,*) userdef_examplevariable
     gotit = .true.
  else
!
! NOTE: It is useful to keep this, as it will stop RADMC-3D if the user  accidently mistypes a command-line keyword instead
!       of simply ignoring it (and thus possibly leaving the user convinced he/she did it right). 
!
     write (*,'(a)') ' USERDEF_COMMANDLINE:'
     write (*,'(a)') '   ERROR: Unrecognized command line option: '//trim(buffer)
     stop
  endif
!!  write (*,'(a)') buffer
!
end subroutine userdef_commandline

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_commandline_postprocessing ()
!__________________________________________________________________________________________________________________________________
!
! Here you can do some postprocessing after all command line options have been read. No example given here.
!__________________________________________________________________________________________________________________________________
!
  implicit none
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
end subroutine userdef_commandline_postprocessing

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_parse_main_namelist ()
!__________________________________________________________________________________________________________________________________
!
! Here you can parse your own keyword entries from the 'radmc3d.inp' file.
!__________________________________________________________________________________________________________________________________
!
  implicit none
!__________________________________________________________________________________________________________________________________
!
! And example how to include an integer keyword. Use parse_input_double and parse_input_word for reading doubleprecision
! variables or string variables. You can delete the below example. Note that the keyword name string should always have length 30,
! hence many whitespaces.
!__________________________________________________________________________________________________________________________________
!
  call parse_input_integer ( 'nzones@                       ', userdef_nzones   )
  call parse_input_integer ( 'nzembc@                       ', userdef_nzembc   )
  call parse_input_double  ( 'bestem@                       ', userdef_bestem   )
  call parse_input_double  ( 'besden@                       ', userdef_besden   )
  call parse_input_double  ( 'Rinner@                       ', userdef_Rinner   )
  call parse_input_double  ( 'Rbesph@                       ', userdef_Rbesph   )
  call parse_input_double  ( 'Rgauss@                       ', userdef_Rgauss   )
  call parse_input_double  ( 'gfwhm@                        ', userdef_gfwhm    )
  call parse_input_double  ( 'gplaw@                        ', userdef_gplaw    )
  call parse_input_double  ( 'Rplawx@                       ', userdef_Rplawx   )
  call parse_input_double  ( 'Router@                       ', userdef_Router   )
  call parse_input_double  ( 'dplawx@                       ', userdef_dplawx   )
  call parse_input_double  ( 'dplout@                       ', userdef_dplout   )
  call parse_input_double  ( 'fgeom@                        ', userdef_fgeom    )
  call parse_input_double  ( 'ftrans@                       ', userdef_ftrans   )
  call parse_input_double  ( 'rhoplx@                       ', userdef_rhoplx   )
  call parse_input_double  ( 'rhofac@                       ', userdef_rhofac   )
  call parse_input_double  ( 'embden@                       ', userdef_embden   )
  call parse_input_double  ( 'ximax@                        ', userdef_ximax    )
  call parse_input_double  ( 'masstot@                      ', userdef_masstot  )
  call parse_input_double  ( 'dust2gas@                     ', userdef_dust2gas )
!!  call parse_input_double  ( 'opascale@                     ', userdef_opascale )
  call parse_input_double  ( 'distpc@                       ', userdef_distpc   )
  call parse_input_integer ( 'external@                     ', userdef_external )
  call parse_input_double  ( 'hbbfac@                       ', userdef_hbbfac   )
  call parse_input_double  ( 'temstar@                      ', userdef_temstar  )
  call parse_input_double  ( 'lumstar@                      ', userdef_lumstar  )
  call parse_input_double  ( 'fracaper@                     ', userdef_fracaper )
  call parse_input_double  ( 'fraczero@                     ', userdef_fraczero )
  call parse_input_double  ( 'mcrscale@                     ', userdef_mcrscale )
  call parse_input_double  ( 'radscale@                     ', userdef_radscale )
  call parse_input_double  ( 'denscale@                     ', userdef_denscale )
!__________________________________________________________________________________________________________________________________
!
end subroutine userdef_parse_main_namelist

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_main_namelist_postprocessing ()
!__________________________________________________________________________________________________________________________________
!
! Here you can do some post-processing after the radmc3d.inp namelist reading. 
!__________________________________________________________________________________________________________________________________
!
  implicit none
!__________________________________________________________________________________________________________________________________
!
  write (*,'( )')
  write (*,'(a)') ' Original screen output is redirected to '''//trim(outputname)//''''
  write (*,'(a)') ' ______________________________________________________________________________'
  write (*,'(a)') '                                                                               '
  write (*,'(a)') '                   RADMC3D: A 3D CONTINUUM AND LINE RT SOLVER                  '
  write (*,'(a)') '                  Version 2.0 (c) 2008-2023 Cornelis Dullemond                 '
  write (*,'(a)') '                     Modifications by Alexander Men’shchikov                   '
  write (*,'(a)') '                      '//codeID
  write (*,'(a)') '                                 Version '//version
  write (*,'(a)') ' ______________________________________________________________________________'
  write (*,'( )')
  write (stdo,'( )')
  write (stdo,'(a)') ' Original screen output is redirected to '''//trim(outputname)//''''
  write (stdo,'(a)') ' ______________________________________________________________________________'
  write (stdo,'(a)') '                                                                               '
  write (stdo,'(a)') '                   RADMC3D: A 3D CONTINUUM AND LINE RT SOLVER                  '
  write (stdo,'(a)') '                  Version 2.0 (c) 2008-2023 Cornelis Dullemond                 '
  write (stdo,'(a)') '                     Modifications by Alexander Men’shchikov                   '
  write (stdo,'(a)') '                      '//codeID
  write (stdo,'(a)') '                                 Version '//version
  write (stdo,'(a)') ' ______________________________________________________________________________'
  write (stdo,'( )')
  !__________________________________________________________________________________________________________________________________
!
end subroutine userdef_main_namelist_postprocessing

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_prep_model ()
!__________________________________________________________________________________________________________________________________
!
! This is the place where you can define your own (base) grid setup, read your own frequency grid or set up your own stellar
! sources. No example given here, because it would interfere with basic operations.
!__________________________________________________________________________________________________________________________________
!
  implicit none
  logical :: lexist
  integer, parameter :: itmax = 999
  integer :: iter, irad, nzonesm1, inu, i, icore, ixmax, nbesf, ir, n_cloud, n_core
  double precision, parameter :: tolerr = 1.0d-6
  double precision :: error, const, sflog, radstar, sndspeed, beam500um, rad0, T1, T2,xi1, factcloud &
                    , Rtrans, span, lamb(1000), delta, L, th, u, xi, arg, ddelta, Fderived, g, gp, s0, sh1, sh2
  double precision, allocatable :: ksiBES(:), psiBES(:), radBES(:), rhoBES(:), rgrid(:), rgridn(:)
!__________________________________________________________________________________________________________________________________
!
  nzonesm1 = userdef_nzones - 1
  
  if (userdef_Rbesph .gt. almostzero) then

    inquire ( file='besphere.inp', exist=lexist )

    if (lexist) then
      open ( unit=1, file='besphere.inp', status='old' )
      do i=1,10
        read (1,*)
      enddo
      nbesf = 16403
      allocate ( rhoBES(nbesf), radBES(nbesf) )
      allocate ( ksiBES(nbesf), psiBES(nbesf) )
      do i=1,nbesf
        read (1,*) ksiBES(i), psiBES(i)
      enddo
      close (1)

      sndspeed = sqrt ( GasConst * userdef_bestem / 2.33d0 )      !<-- Sound speed in cm/s
      rad0 = sndspeed / sqrt ( fourpi * GravC * userdef_besden )  !<~ conversion factor to physical units
      radBES = rad0 * ksiBES / xAU                                !<-- Radial grid converted from cm to AU
      radBES = userdef_radscale * userdef_mcrscale * radBES
      rhoBES = userdef_besden * exp ( -z1 * psiBES )              !<-- Density distribution of Bonnor-Ebert sphere
      rhoBES = userdef_denscale * rhoBES
      ixmax = 1
      icore = 1
      if (userdef_ximax .gt. almostzero) then
         call hunt ( ksiBES, nbesf, userdef_ximax, ixmax )
         if (ixmax .lt. 1 .or. ixmax .gt. nbesf) then
            write (*,'(a   )') ' RADMC3D: USERDEF_SETUP_MODEL (1):'
            write (*,'(a,i0)') '   ERROR: BES table index out of range: ', icore
            stop 91
         endif
         Rmodel = radBES(ixmax)
         Rbesph = Rmodel
      endif
      deallocate ( ksiBES, psiBES )
    else
      write (*,'(a)') ' RADMC3D: USERDEF_SETUP_MODEL:'
      write (*,'(a)') '   ERROR: Data file ''besphere.inp'' not found.'
      stop 99
    endif
  endif
  
! Router = R_cloud = ratio · √( (beam500 · distance)² + R_BE² )

  beam500um = 0.0d0 !! 36.3d0
  factcloud = 30.0d0
  
  if (userdef_Router .gt. almostzero) then
    Router = userdef_radscale * userdef_Router
  else
    Router = sqrt ( (beam500um * userdef_distpc)**2 + (userdef_radscale * factcloud * Rmodel)**2 )
    userdef_Router = Router
  endif
  if (userdef_Rinner .gt. almostzero) then
    Rinner = userdef_radscale * userdef_Rinner
  else
!!    Rinner = userdef_radscale * Rmodel / 100.0d0
!!    write (*,*) Rinner
    Rinner = userdef_radscale * Rmodel / (70.0d0 * userdef_ximax)
!!    write (*,*) Rinner
    userdef_Rinner = Rinner
  endif
  
  if (userdef_Rbesph .gt. almostzero) then
    Rbesph = userdef_radscale * userdef_mcrscale * Rbesph
    Rmodel = Rbesph
  endif
  if (userdef_Rgauss .gt. almostzero) then
    Rgauss = userdef_radscale * userdef_mcrscale * userdef_Rgauss
    Rmodel = Rgauss
  endif
  if (userdef_Rplawx .gt. almostzero) then
    Rplawx = userdef_radscale * userdef_mcrscale * userdef_Rplawx
    Rmodel = Rplawx
  endif

  allocate ( rgrid(userdef_nzones+1), rgridn(userdef_nzones+1) )

  if (userdef_ftrans < almostzero) then

! Geometric logarithmic grid.

    sf = z10**((log10 ( Router * xAU ) - log10 ( Rinner * xAU )) / dble ( nzonesm1 ))
    sfmz1 = sf - z1
    const = (sf**(nzonesm1) - z1) / sfmz1 / userdef_fgeom
  
    do iter=1,itmax
       sfmz1 = (const * sfmz1 + z1)**(z1 / (nzonesm1)) - z1
       error = const * sfmz1 / ((z1 + sfmz1)**(nzonesm1) - z1) - z1
       if (abs ( error ) < tolerr) exit
       if (iter == itmax) then
          write (*,'(a)') ' USERDEF_PREP_MODEL:'
          write (*,'(a,i3,a)') '   ERROR (A): No convergence in geometric grid after ', iter, ' iterations.'
          write (*,'(a,1pe9.2,a,e9.2)') '          Error:', error, ' >', tolerr
          stop 99
       endif
    enddo
  
    sf = sfmz1 + z1
    sflog = z1 / log ( sf )
    dr1 = (Router * xAU - Rinner * xAU) * sfmz1 / (sf**(nzonesm1) - z1)
    do irad=1,userdef_nzones
      rgrid(irad) = Rinner * xAU + dr1 * (sf**(irad - 1) - z1) / sfmz1
    enddo
  else

! Transition: a bit beyond the BE edge, where T has flattened.
! Split the cell budget between core and embedding cloud.
  
    Rtrans = userdef_ftrans * Rmodel * xAU        ! ftrans ~ 1.5-3, a userdef parameter
    Rtrans = min ( Rtrans, z05 * Router * xAU )   ! guard: keep a real cloud segment
    n_cloud = userdef_nzembc
    n_core = nzonesm1 - n_cloud

! Core segment: keep existing solver, but targeted at [Rinner, Rtrans].

    span = log10 ( Rtrans ) - log10 ( Rinner * xAU )
    sf = z10**(span / dble ( n_core ))
    sfmz1 = sf - z1
    const = (sf**(n_core) - z1) / sfmz1 / userdef_fgeom
    
    do iter=1,itmax
      sfmz1 = (const * sfmz1 + z1)**(z1 / (n_core)) - z1
      error = const * sfmz1 / ((z1 + sfmz1)**(n_core) - z1) - z1
      if (abs ( error ) < tolerr) exit
      if (iter == itmax) then
        write (*,'(a)') ' USERDEF_PREP_MODEL:'
        write (*,'(a,i3,a)') '   ERROR (B): No convergence in geometric grid after ', iter, ' iterations.'
        write (*,'(a,1pe9.2,a,e9.2)') '          Error:', error, ' >', tolerr
        stop 99
      endif
    enddo

    sf = sfmz1 + z1
    sflog = z1 / log ( sf )
    dr1 = (Rtrans - Rinner * xAU) * sfmz1 / (sf**(n_core) - z1)
    do irad=1,n_core+1
      rgrid(irad) = Rinner * xAU + dr1 * (sf**(irad - 1) - z1) / sfmz1
    enddo

! n_cloud fixed; solve delta so first cloud interval, core's last.

    L = log ( Router * xAU ) - log ( Rtrans )
    s0 = log ( rgrid(n_core+1) ) - log ( rgrid(n_core) )    ! core's last cell width (ln r)
    xi1 = z1 / dble ( n_cloud )
    arg = xi1 - z05                                   ! < 0

! tanh clustering only makes the ends FINER than uniform (L/n_cloud);
! a match needs s0 < L/n_cloud, else n_cloud is too large for this model's span.

    if (s0 >= L / dble ( n_cloud )) then
      delta = 1.0d-3                    ! fall back to ~uniform log
      write (*,'(a)') ' USERDEF_PREP_MODEL:'
      write (*,'(a)') '   WARNING: n_cloud too large to match core join;'
      write (*,'(a)') '            cloud is ~uniform, join not exact. Reduce n_cloud or move Rtrans in.'
    else
      delta = 2.0d0                     ! initial guess
      do iter=1,itmax
        T1 = tanh ( delta * arg )
        T2 = tanh ( z05 * delta )
        sh1 = z1 / cosh ( delta * arg )**2       ! sech^2(delta*arg)
        sh2 = z1 / cosh ( z05 * delta )**2       ! sech^2(delta/2)
        g = z05 * L * (z1 + T1 / T2) - s0
        gp = z05 * L * (arg * sh1 * T2 - z05 * T1 * sh2) / T2**2
        ddelta = g / gp
        delta = delta - ddelta
        if (delta <= z0) delta = 1.0d-3       ! keep positive
        if (delta > 60.0d0) delta = 60.0d0    ! clamp extreme clustering
        if (abs ( ddelta ) < tolerr) exit
        if (iter == itmax) then
          write(*,'(a)') ' USERDEF_PREP_MODEL:'
          write(*,'(a,1pe9.2)') '   ERROR: Cloud delta not converged, last step=', ddelta
          stop 99
        endif
      enddo
    endif

! Build cloud nodes with the solved delta.

    th = tanh ( z05 * delta )
    do i = 0, n_cloud
      xi = dble ( i ) / dble ( n_cloud )
      u  = z05 * (z1 + tanh ( delta * (xi - z05) ) / th)
      rgrid(n_core+1+i) = exp ( log ( Rtrans ) + L * u )
    enddo

    Fderived = cosh ( z05 * delta )**2     ! middle/end coarsening ratio (log it)
  
    if (Fderived > 30.0d0) then
      write(*,'(a)') ' USERDEF_PREP_MODEL:'
      write(*,'(a,1pe11.4)') '   ERROR: Cloud grid too coarse: ', Fderived
      stop 99
    endif
  endif

! After building the grid, correct boundary cell before assigning any density. 
! Then every cell will be unambiguously pure BE (rout <= Rmodel) or pure cloud (rin >= Rmodel); 
! the boundary cell ends exactly at Rmodel with no clamp and no mixing. 

!!  imin = minloc ( abs ( rgrid(1:userdef_nzones) - Rmodel * xAU ), 1 )
!!  rgrid(imin) = Rmodel * xAU

!!  ir = 0
!!  do irad=1,userdef_nzones+1
!!    ir = ir + 1
!!    if (rgrid(irad) .lt. Rmodel * xAU .and. rgrid(irad+1) .gt. Rmodel * xAU) then
!!      write (*,*) irad, ir, rgrid(irad), rgridn(ir)
!!      ir = ir + 1
!!      rgridn(ir) = Rmodel * xAU
!!    else
!!      rgridn(ir) = rgrid(irad)
!!    endif
!!    
!!    write (*,*) irad, ir, rgrid(irad), rgridn(ir)
!!  enddo

  ir = 0
  do irad=1,userdef_nzones
    ir = ir + 1
    rgridn(ir) = rgrid(irad)                       ! always copy the original wall
    if (irad .le. userdef_nzones) then             ! guard the irad+1 read
      if (rgrid(irad) .lt. Rmodel*xAU .and. rgrid(irad+1) .gt. Rmodel*xAU) then
        ir = ir + 1
        rgridn(ir) = Rmodel * xAU                  ! insert R_BE right after it
      endif
    endif
!!    write (*,*) irad, ir, rgrid(irad), rgridn(ir), Rmodel * xAU
  enddo

  userdef_nzones = ir
  nzonesm1 = userdef_nzones - 1

!!  do irad=1,userdef_nzones
!!    write (*,*) irad, rgrid(irad), rgridn(irad), Rmodel * xAU
!!  enddo

  open ( unit=1, file='amr_grid.inp', status='unknown' )
  write (1,'(a)'   ) '1               ; iformat'
  write (1,'(a)'   ) '0               ; AMR grid style  (0 = regular grid, no AMR)'
  write (1,'(a)'   ) '100             ; Coordinate system'
  write (1,'(a)'   ) '0               ; gridinfo'
  write (1,'(a)'   ) '1,0,0           ; Include x,y,z coordinate'
  write (1,'(i0,a)') nzonesm1,',1,1          ; Size of grid'

  do irad=1,userdef_nzones
     write (1,*) rgridn(irad)
  enddo

  write (1,'(a)') '0.0d0'
  write (1,'(a)') '1.0d0'
  write (1,'(a)') '0.0d0'
  write (1,'(a)') '1.0d0'
  close (1)

  deallocate ( rgrid )

  open ( unit=1, file='wavelength_micron.inp', status='unknown' )
  read (1,*) freq_nr
  do inu=1,freq_nr
    read (1,*) lamb(inu)
  enddo
  close (1)

  if (userdef_lumstar .gt. almostzero) then
    open ( unit=3, file='stars.inp', status='unknown' )
    write (3,'(a)'        ) '2               ; iformat'
    write (3,'(a,1x,i0,a)') '1', freq_nr,'          ; Number of wavelengths'
    radstar = max ( sqrt ( userdef_lumstar * xLsun / (fourpi * SBconst * userdef_temstar**4) ) &
                  , sqrt ( zeps10 * xLsun / (fourpi * SBconst * userdef_temstar**4) ) )
    write (3,'(5(1pe15.8))') radstar, xMsun, z0, z0, z0
    do inu=1,freq_nr
      write (3,*) lamb(inu)
    enddo
    write (3,'(1pe16.8)') -z1 * userdef_temstar
    close (3)
  else
    inquire ( file='stars.inp', exist=lexist )
    if (lexist) then
      open ( unit=2, file='stars.inp', status='old' )
      close ( 2, status='delete' )
    endif
  endif
  !__________________________________________________________________________________________________________________________________
!
end subroutine userdef_prep_model

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_setup_model ()
!__________________________________________________________________________________________________________________________________
!
! This is the place where you can set up your own model. By the time this subroutine is called the grid (at least a basic grid) is
! already set up. You can still modify the basic grid by adding more refinement, but to tell the AMR module to reserve space for 
! more grid points you need to take matters into your own hand and create and init the base grid yourself in the 
! userdef_prep_model() routine above. No example given here, because it would interfere with basic operations.
!__________________________________________________________________________________________________________________________________
!
  implicit none
  logical :: lexist
  character*(50) :: cbestem, cbesden, crhobou, cmasscr, cradiusau, cradiuspc, cradiusas, crhoratio, cembden, crhoemb
  character(*), parameter :: esc = achar(27)
  integer :: irad, nzonesm1, i, ic, is, idx, icore, nbesf, k1, k2, k3, inu, istar, ixmax, ia, ib
  double precision :: massgas, sndspeed, lam1, lam2, lam3, eps1, eps2, eps3, int1, int2 &
                    , ifuv, G0, dum, dust_mass, tdust_massav, dzone, opacity, odepthv, odepth100, dust_mass_embedding &
                    , masstot_emb, fracaper, bedenratio, radcell, hotbbtemp, isrflux, opacity_ext &
                    , opacity_abs, opacity_sca, drad, radcellp1, radcellm1, densRmodel, argex, gauss, plaw, dgausdr, dplawdr &
                    , fitfact, hwhm, sigm2g, sigm2g2, dlograd, radfit, expfit, plwfit, correction, Cin, Cout, r0eff, rin &
                    , rout, xin, xout, xedge, xfull, rad0, massBE, embvoldens
  double precision, allocatable :: ksiBES(:), psiBES(:), cBES(:), radBES(:), rhoBES(:), density(:), radius(:), dslope(:), lamb(:) &
                                 , odepth(:)
!__________________________________________________________________________________________________________________________________
!
  if (userdef_Rbesph .gt. almostzero) then

    inquire ( file='besphere.inp', exist=lexist )

    if (lexist) then
      open ( unit=1, file='besphere.inp', status='old' )
      do i=1,10
        read (1,*)
      enddo
      nbesf = 16403
      allocate ( rhoBES(nbesf), radBES(nbesf) )
      allocate ( ksiBES(nbesf), psiBES(nbesf), cBES(nbesf) )
      do i=1,nbesf
        read (1,*) ksiBES(i), psiBES(i)
      enddo
      close (1)

! Find the number of radial point for the BE sphere.

      call hunt ( ksiBES, nbesf, userdef_ximax, ixmax )

! ρ(ξ) = ρ_c e^{−ψ},  ρ(r) = ρ_c e^{−ψ(r/r₀)},  r₀ = c_s/√(4πGρ_c),  M = 4π ρ_c r₀³ ∫₀^{6.451} e^{−ψ} ξ² dξ
! r₀ is the Bonnor-Ebert scale length - the dimensional factor that turns the dimensionless ξ 
! into a physical radius: r₀ = c_s / √(4πGρ_c), with c_s = √(k_B T/(μm_H)) = √(GasConst·T/2.33)
! c_s² = kT/(μm_H),  c_s = √(kT/μm_H),  M = 4.43 c_s³/(G^{3/2} √ρ_c)
! Correct inversion for ρ_c at fixed T is ρ_c = (4.43 c_s³ / (G^{3/2} M))²

      sndspeed = sqrt ( GasConst * userdef_bestem / 2.33d0 )      !<~ Sound speed in cm/s
      rad0 = sndspeed / sqrt ( fourpi * GravC * userdef_besden )  !<~ conversion factor to physical units

! Integrate mass of the BE sphere:

      massBE = z0
      do i=2,ixmax
        massBE = massBE + fourpi * userdef_besden * rad0**3 &
                          * z05 * ( exp ( -psiBES(i-1) ) * (ksiBES(i-1))**2   &
                                  + exp ( -psiBES(i  ) ) * (ksiBES(i  ))**2 ) &
                                * ( ksiBES(i) - ksiBES(i-1) )
      enddo
      massBE = massBE / xMsun

! Final BE sphere radii and densities in physical units: 

      radBES = rad0 * ksiBES / xAU                              !<~ Radial grid in AU
      radBES = userdef_radscale * userdef_mcrscale * radBES
      rhoBES = userdef_besden * exp ( -psiBES )                !<~ Density distribution of Bonnor-Ebert sphere
      rhoBES = userdef_denscale * rhoBES

      embvoldens = userdef_embden * muH2 * AtomicMU / (z2 * sqrt ( (userdef_Router * xAU)**2 - (userdef_Rbesph * xAU)**2 ))
  
!!      write (*,'(1pe11.5)') userdef_Router, xAU, userdef_Router * xAU, userdef_embden, embvoldens

      open ( unit=1, file=trim(outname)//'.be' )
      write (cbestem,  '(1pe11.5)') userdef_bestem
      write (cbesden,  '(1pe11.5)') userdef_besden
      write (cembden,  '(1pe11.5)') userdef_embden
      write (crhobou,  '(1pe11.5)') rhoBES(ixmax)
      write (crhoemb,  '(1pe11.5)') embvoldens
      write (crhoratio,'(1pe11.5)') rhoBES(1) / rhoBES(ixmax)
      write (cmasscr,  '(1pe11.5)') massBE
      write (cradiuspc,'(1pe11.5)') radBES(ixmax) * xAU / parsec
      write (cradiusau,'(1pe11.5)') radBES(ixmax)
      write (cradiusas,'(1pe11.5)') radBES(ixmax) / userdef_distpc
      cbestem = trim ( adjustl ( cbestem ) )
      cbesden = trim ( adjustl ( cbesden ) )
      cembden = trim ( adjustl ( cembden ) )
      cmasscr = trim ( adjustl ( cmasscr ) )
      crhobou = trim ( adjustl ( crhobou ) )
      crhoemb = trim ( adjustl ( crhoemb ) )
      crhoratio = trim ( adjustl ( crhoratio ) )
      cradiusau = trim ( adjustl ( cradiusau ) )
      cradiuspc = trim ( adjustl ( cradiuspc ) )
      cradiusas = trim ( adjustl ( cradiusas ) )
      write (1,'(a)') '! True properties of the BE sphere:'
      write (1,'(a)') '               Mass (Msun): '//cmasscr(1:11)
      write (1,'(a)') '           Temperature (K): '//cbestem(1:11)
      write (1,'(a)') '  Central density (g/cm^3): '//cbesden(1:11)
      write (1,'(a)') ' Boundary density (g/cm^3): '//crhobou(1:11)
      write (1,'(a)') '    Cloud density (g/cm^3): '//crhoemb(1:11)
      write (1,'(a)') '    Cloud sdensity (cm^-2): '//cembden(1:11)
      write (1,'(a)') '             Density ratio: '//crhoratio(1:11)
      write (1,'(a)') '      Boundary radius (AU): '//cradiusau(1:11)
      write (1,'(a)') '      Boundary radius (pc): '//cradiuspc(1:11)
      write (1,'(a)') '      Boundary radius (as): '//cradiusas(1:11)
      close (1)
      
! Mass conservation for cell densities. Instead of assigning each radial cell the density sampled at its center, 
! assign it the analytic mass-conserving average over the cell:
! ρ_cell = (1/V_cell) ∫cell ρ_BE(r) dV = ∫ from {r_i} to {r_{i+1}} of ρ_c e^{−ψ(r)} 4πr² dr / [(4π/3)(r_{i+1}³ − r_i³)]
!
! But we don't need to do a sub-cell numerical integral at all, because the Lane-Emden equation gives the enclosed mass in closed
! form. Since ∫₀^ξ e^{−ψ(ξ′)} ξ′² dξ′ = ξ² ψ′(ξ), the mass inside any radius is M(<ξ) = 4π ρ_c r₀³ · ξ²ψ′(ξ), and the shell mass is
! just the difference at the two walls. Working in ξ = r/r₀ (so r₀ and 4π all cancel):
! ρ_cell = 3 ρ_c [ξ_{i+1}² ψ′(ξ_{i+1}) − ξ_i² ψ′(ξ_i)] / (ξ_{i+1}³ − ξ_i³)
! where ξ_i = r_i/r₀ are your cell walls in dimensionless units. Since besphere.inp table carries ψ but not ψ′, the most robust
! route avoids differentiating: precompute the cumulative dimensionless mass once:
! C(ξ) = ∫₀^ξ e^{−psiBES} ksiBES² d(ksiBES) (trapezoid over the table), 
! interpolate C onto your cell walls, and use ρ_cell = 3 ρ_c [C(ξ_{i+1}) − C(ξ_i)] / (ξ_{i+1}³ − ξ_i³).
! This uses only ψ (which you have), needs no per-cell sub-integration, and gives Σ ρ_cell V_cell = analytic mass to the accuracy
! of that one cumulative integral — exact regardless of how few zones you run. As a sanity check, C(6.451) should come out ≈ 15.7.
!
! C(ξ) is the integral of e^{−ψ(ξ′)} ξ′² dξ′, taken from ξ′ = 0 to ξ′ = ξ. In your code it's the running cumulative trapezoidal sum
! of exp(−psiBES) * ksiBES**2 along the table, and C evaluated at ξ = 6.451 should land near 15.7.

      cBES(1) = z0
      do i=2,nbesf
        cBES(i) = cBES(i-1) + z05 * ( exp ( -psiBES(i-1) ) * ksiBES(i-1)**2 + exp ( -psiBES(i) ) * ksiBES(i)**2 ) &
                                  * ( ksiBES(i) - ksiBES(i-1) )
      enddo

      icore = 1
      if (userdef_ximax .gt. almostzero) then
         if (ixmax .lt. 1 .or. ixmax .gt. nbesf) then
            write (*,'(a   )') ' RADMC3D: USERDEF_SETUP_MODEL (1):'
            write (*,'(a,i0)') '   ERROR: BES table index out of range: ', icore
            stop 91
         endif
         Rmodel = radBES(ixmax)
         Rbesph = Rmodel
         bedenratio = rhoBES(1) / rhoBES(ixmax)
      else
         call hunt ( radBES, nbesf, Rbesph, icore )
         if (icore .lt. 1 .or. icore .gt. nbesf) then
            write (*,'(a   )') ' RADMC3D: USERDEF_SETUP_MODEL (2):'
            write (*,'(a,i0)') '   ERROR: BES table index out of range: ', icore
            stop 92
         endif
      endif
    else
      write (*,'(a)') ' RADMC3D: USERDEF_SETUP_MODEL (3):'
      write (*,'(a)') '   ERROR: Data file ''besphere.inp'' not found.'
      stop 99
    endif
  endif
!__________________________________________________________________________________________________________________________________
!
  nzonesm1 = userdef_nzones - 1

  allocate ( density(userdef_nzones), dslope(userdef_nzones), radius(userdef_nzones) )

  open ( unit=1, file='amr_grid.inp', status='unknown' )
  read (1,*) ! '1               ; iformat'
  read (1,*) ! '0               ; AMR grid style  (0 = regular grid, no AMR)'
  read (1,*) ! '100             ; Coordinate system'
  read (1,*) ! '0               ; gridinfo'
  read (1,*) ! '1,0,0           ; Include x,y,z coordinate'
  read (1,*) ! nzonesm1,',1,1          ; Size of grid'
  do irad=1,userdef_nzones
     read (1,*) radius(irad)
     radius(irad) = radius(irad) / xAU
  enddo
  close (1)

  open ( unit=1, file='dust_density.inp', status='unknown' )
  write (1,'(a)'   ) '1         ; Format number'
  write (1,'(i0,a)') nzonesm1, '        ; Nr of cells'
  write (1,'(a)'   ) '1         ; Nr of dust species'

  if (userdef_Rgauss .gt. almostzero) then
    hwhm = userdef_gfwhm / 2.0d0
    sigm2g = hwhm**2 / log ( 4.0d0 )
    sigm2g2 = 2.0d0 * sigm2g

! Fit parameters for the Gaussian and power-law functions at radius where they derivatives match.
     
    radfit = sqrt ( userdef_gplaw * sigm2g )
    expfit = exp ( -min ( radfit**2 / sigm2g2, 50.0d0 ))
    plwfit = (1.0d0 / radfit)**userdef_gplaw
    fitfact = expfit / plwfit
  endif

! For mass conservation with any radial zone number: ρ_cell = 3 ρ_c [C(ξ_{i+1}) − C(ξ_i)] / (ξ_{i+1}³ − ξ_i³).
  
  r0eff = radBES(2) / ksiBES(2)   !<~ AU per unit xi (constant; includes radscale*mcrscale)
  xedge = Rmodel / r0eff
  density = z0
  
  do irad=1,nzonesm1
    drad = radius(irad+1) - radius(irad)
    radcell = (radius(irad) + radius(irad+1)) / z2
    radcellm1 = (radius(max ( irad-1, 1 )) + radius(irad)) / z2
    radcellp1 = (radius(irad+1) + radius(min ( irad+2, userdef_nzones ))) / z2
    dlograd = log10 ( radcell ) - log10 ( radcellm1 )
    dslope(irad) = 1.0d-30

    rin  = radius(irad)
    rout = radius(irad+1)

    if (userdef_Rbesph .gt. almostzero) then
      if (rin .lt. Rmodel) then                   ! cell at least partly inside the BE sphere
        xin   = rin  / r0eff                      ! inner wall (dimensionless)
        xfull = rout / r0eff                      ! true outer wall, BEFORE clamp
        if (rout .gt. Rmodel) rout = Rmodel       ! clamp the BE/cloud split point
        xout  = rout / r0eff                      ! = min(xfull, edge) — the BE edge for this cell
        call hunt ( ksiBES, nbesf, xin,  ia )
        call hunt ( ksiBES, nbesf, xout, ib )

        if (xin > 0.3d0) then
          if (ia < 1 .or. ia >= nbesf) then
             write(*,'(a,1pe12.4,a,i0)') ' USERDEF: xin=',xin,' outside cBES table, ia=',ia
             stop 99
          endif
          Cin = cBES(ia) + (xin - ksiBES(ia)) / (ksiBES(ia+1) - ksiBES(ia)) * (cBES(ia+1) - cBES(ia))
        else
          Cin = xin**3 / 3.0d0 - xin**5 / 30.0d0 + xin**7 / 840.0d0
        endif
        if (xout > 0.3d0) then
          if (ib < 1 .or. ib >= nbesf) then
             write(*,'(a,1pe12.4,a,i0)') ' USERDEF: xout=',xout,' outside cBES table, ib=',ib
             stop 99
          endif
          Cout = cBES(ib) + (xout - ksiBES(ib)) / (ksiBES(ib+1) - ksiBES(ib)) * (cBES(ib+1) - cBES(ib))
        else
          Cout = xout**3 / 3.0d0 - xout**5 / 30.0d0 + xout**7 / 840.0d0
        endif
        density(irad) = (3.0d0 * userdef_besden * (Cout - Cin) + embvoldens * (xfull**3 - xout**3)) / (xfull**3 - xin**3)
      endif
    endif

    if (userdef_Rgauss .gt. almostzero) then
      if (radcellm1 .le. Rmodel) then
        if (userdef_gplaw .gt. almostzero) then
          plaw = (1.0d0 / radcell)**userdef_gplaw
          dplawdr = userdef_gplaw
        endif
        argex = min ( radcell**2 / sigm2g2, 50.0d0 )
        gauss = exp ( -argex )
        dgausdr = radcell**2 / sigm2g
        density(irad) = userdef_denscale * gauss
        if (userdef_gplaw .gt. almostzero) then
          if (radcell .ge. radfit) then
            density(irad) = userdef_denscale * plaw * fitfact
          endif
        endif
        correction = 2.0d0**(2.0d0 / userdef_gplaw) - 1.0d0
        density(irad) = (1.0d0 + correction * (2.0d0 * radcell / userdef_gfwhm)**2)**(-userdef_gplaw / 2.0d0)
      endif
    endif
    
    if (userdef_Rplawx .gt. almostzero) then
      if (radcellm1 .le. Rmodel) then
        density(irad) = userdef_denscale * userdef_rhoplx * (radius(1) / radcell)**userdef_dplawx
      endif
    endif

    if (irad .gt. 1) then
      if (radius(irad) .lt. Rmodel) then
        densRmodel = density(irad)
        dslope(irad) = (log10 ( density(irad) ) - log10 ( density(irad-1) )) / dlograd
      endif
    endif
    
    write (1,*) density(irad), z0
  enddo
  
  close (1)

  if (userdef_Rbesph .gt. almostzero) then
    densRmodel = rhoBES(ixmax)
  endif

  if (allocated(ksiBES)) deallocate ( ksiBES )
  if (allocated(psiBES)) deallocate ( psiBES )
  if (allocated(radBES)) deallocate ( radBES )
  if (allocated(rhoBES)) deallocate ( rhoBES )
!__________________________________________________________________________________________________________________________________
!
  call read_dustdata ( 1 )
  call read_dust_density ( 1 )

  masstot = z0
  do irad=1,nzonesm1
    masstot = masstot + density(irad) * (fourpi / 3.0d0) * ((radius(irad+1) * xAU)**3 - (radius(irad) * xAU)**3)
  enddo  
    
! Create model density input file.
  
  open ( unit=1, file='dust_density.inp', status='unknown' )
  write (1,'(a)') '1         ; Format number'
  write (1,'(i0,a)') nzonesm1, '     ; Nr of cells'
  write (1,'(a)') '1         ; Nr of dust species'

  do irad=1,nzonesm1
    radcell = (radius(irad) + radius(irad+1)) / z2
    density(irad) = density(irad) * userdef_dust2gas
    if (radius(irad) .ge. Rmodel) then
      if (userdef_rhofac .gt. almostzero) then
        density(irad) = userdef_rhofac * densRmodel * (Rmodel / radcell)**userdef_dplout
      else
        density(irad) = embvoldens * userdef_dust2gas
      endif
    endif
  enddo

  do irad=1,nzonesm1
    radcell = (radius(irad) + radius(irad+1)) / z2
    if (userdef_fraczero .gt. almostzero) then
      if (radcell .gt. userdef_fraczero * Router) then
        density(irad) = z0
      endif
    endif

    write (1,*) density(irad), z0
  enddo

  close (1)

  allocate ( lamb(freq_nr), odepth(freq_nr) )

  open ( unit=3, file='wavelength_micron.inp', status='unknown' )
  read (3,*) freq_nr
  do inu=1,freq_nr
    read (3,*) lamb(inu)
  enddo
  close (3)
  
  fracaper = 1.0d0
  if (userdef_fracaper .gt. almostzero) then
    fracaper = userdef_fracaper
  endif

  open ( unit=2, file='aperture_info.inp', status='unknown' )
  write (2,'(a)'   ) '1             ; iformat'
  write (2,'(i0,a)') freq_nr, '           ; Number of wavelengths'
  do inu=1,freq_nr
    if (userdef_fracaper .gt. almostzero) then
      aperture(inu) = fracaper * Rmodel / userdef_distpc
    else
      aperture(inu) = Rmodel / userdef_distpc
    endif
    write (2,*) lamb(inu), aperture(inu)
  enddo
  close (2)
!__________________________________________________________________________________________________________________________________
!
  call read_dust_density ( 2 )

  inquire ( file='dust_temperature.dat', exist=lexist )

  if (.not.do_montecarlo_therm .and. lexist) then
    call read_dust_temperature (1)
  else
    if (.not.(do_montecarlo_therm .or. lexist)) then
      write (*,'(a)') ' RADMC3D: USERDEF_SETUP_MODEL (6):'
      write (*,'(a)') '   ERROR: Data file ''dust_temperature.dat'' not found.'
      stop 99
    endif
  endif

  dust_mass = z0
  dust_mass_embedding = z0
  tdust_massav = z0
  do is=1,dust_nr_species
    dust_massdust(is) = z0
  enddo

  do ic=1,nrcells
    idx = cellindex(ic)
    radcell = (radius(ic) + radius(ic+1)) / z2
    if (radcell .lt. Rmodel) then
      do is=1,dust_nr_species
        dust_massdust(is) = dust_massdust(is) + cellvolume(idx) * dustdens(is,idx)
        dust_mass = dust_mass + cellvolume(idx) * dustdens(is,idx)
        if (.not.do_montecarlo_therm .and. lexist) then
          tdust_massav = tdust_massav + cellvolume(idx) * dustdens(is,idx) * dusttemp(is,idx)
        endif
      enddo
    endif
    if (radcell .gt. Rmodel) then
      do is=1,dust_nr_species
        dust_mass_embedding = dust_mass_embedding + cellvolume(idx) * dustdens(is,idx)
      enddo
    endif
  enddo
  if (.not.do_montecarlo_therm .and. lexist) then
    tdust_massav = tdust_massav / dust_mass
  endif

  masstot_emb = dust_mass_embedding / userdef_dust2gas

! Replace temperature distribution with a mass-averaged temperature.

  if (userdef_tmassav .eq. 1 .and. .not.do_montecarlo_therm .and. lexist) then
    do ic=1,nrcells
      idx = cellindex(ic)
      radcell = (radius(ic) + radius(ic+1)) / z2
      if (radcell .lt. Rmodel) then
        do is=1,dust_nr_species
          dusttemp(is,idx) = tdust_massav
        enddo
      endif
    enddo
    open ( unit=1, file='dust_temperature.dat' )
    write (1,*) 1 ! Format number
    write (1,*) nrcells
    write (1,*) dust_nr_species
    do is=1,dust_nr_species
      do ic=1,nrcells
        idx = cellindex(ic)
        write (1,*) dusttemp(is,idx)
      enddo
    enddo
    close ( 1 )
  endif

! Re-create model density input file for limb-brightened background calculation.

  if (do_raytrace_image .and. userdef_doback .eq. 1) then
    open ( unit=1, file='dust_density.inp', status='unknown' )
    write (1,'(a)') '1         ; Format number'
    write (1,'(i0,a)') nzonesm1, '     ; Nr of cells'
    write (1,'(a)') '1         ; Nr of dust species'
  
    do irad=1,nzonesm1
      radcell = (radius(irad) + radius(irad+1)) / z2
      if (do_raytrace_image .and. userdef_doback .eq. 1 .and. radcell .lt. Rmodel) then
        density(irad) = z0
      endif
      write (1,*) density(irad), z0
    enddo
    close (1)
    massgas = z0
    masstot = z0
    call read_dust_density ( 2 )
  endif
!__________________________________________________________________________________________________________________________________
!
  write (stdo,'()')
  write (*   ,'(a,1pe10.3)'         ) '                 Number of photons:', dble(rt_mcparams%nphot_therm) 
  write (stdo,'(a,1pe10.3)'         ) '                 Number of photons:', dble(rt_mcparams%nphot_therm) 
  write (*   ,'()')
  write (stdo,'()')
  write (*   ,'(a,1pe10.3,a)'       ) '                 Distance to model:', userdef_distpc,' pc'
  write (stdo,'(a,1pe10.3,a)'       ) '                 Distance to model:', userdef_distpc,' pc'
  write (*   ,'(a,1pe9.2 ,a)'       ) '            Dust-to-gas mass ratio:', userdef_dust2gas
  write (stdo,'(a,1pe9.2 ,a)'       ) '            Dust-to-gas mass ratio:', userdef_dust2gas
!!  write (*   ,'(a,1pe11.4,a)'       ) '            Opacity scaling factor:', userdef_opascale
!!  write (stdo,'(a,1pe11.4,a)'       ) '            Opacity scaling factor:', userdef_opascale
  write (*   ,'(a,1x,i0,a)'         ) '                     External ISRF:', userdef_external,' (1|0)'
  write (stdo,'(a,1x,i0,a)'         ) '                     External ISRF:', userdef_external,' (1|0)'
  write (*   ,'(a,1pe10.3,a)'       ) '        ISRF hot BB scaling factor:', userdef_hbbfac
  write (stdo,'(a,1pe10.3,a)'       ) '        ISRF hot BB scaling factor:', userdef_hbbfac
  write (*   ,'()')
  write (stdo,'()')
  write (*   ,'(a,2(1pe15.8,a))'    ) '       Inner radius of model space:', userdef_Rinner,' AU =' &
                                                                           , userdef_Rinner / userdef_distpc,' as'
  write (stdo,'(a,2(1pe15.8,a))'    ) '       Inner radius of model space:', userdef_Rinner,' AU =' &
                                                                           , userdef_Rinner / userdef_distpc,' as'
  write (*   ,'(a,2(1pe15.8,a))'    ) '       Outer radius of model space:', userdef_Router,' AU =' &
                                                                           , userdef_Router / userdef_distpc,' as'
  write (stdo,'(a,2(1pe15.8,a))'    ) '       Outer radius of model space:', userdef_Router,' AU =' &
                                                                           , userdef_Router / userdef_distpc,' as'
  write (*   ,'(a,1x,i0)'           ) '      Total number of radial zones:', userdef_nzones
  write (stdo,'(a,1x,i0)'           ) '      Total number of radial zones:', userdef_nzones
  write (*   ,'(a,1x,i0)'           ) '          Zones in embedding cloud:', userdef_nzembc
  write (stdo,'(a,1x,i0)'           ) '          Zones in embedding cloud:', userdef_nzembc
  write (*   ,'(a,1pe9.2 ,a)'       ) '             Geometric grid factor:', userdef_fgeom
  write (stdo,'(a,1pe9.2 ,a)'       ) '             Geometric grid factor:', userdef_fgeom
  write (*   ,'(a,1pe9.2 ,a)'       ) '      Fractional transition radius:', userdef_ftrans
  write (stdo,'(a,1pe9.2 ,a)'       ) '      Fractional transition radius:', userdef_ftrans
  write (*   ,'()')
  write (stdo,'()')
  if (userdef_Rgauss .gt. almostzero) then
  write (*   ,'(a,2(1pe15.8,a))'    ) '            Radius of Gauss sphere:', userdef_Rgauss,' AU =' &
                                                                           , userdef_Rgauss / userdef_distpc,' as'
  write (stdo,'(a,2(1pe15.8,a))'    ) '            Radius of Gauss sphere:', userdef_Rgauss,' AU =' &
                                                                           , userdef_Rgauss / userdef_distpc,' as'
  write (*   ,'(a,2(1pe15.8,a))'    ) '              FWHM of Gauss sphere:', userdef_gfwhm,' AU =' &
                                                                           , userdef_gfwhm / userdef_distpc,' as'
  write (stdo,'(a,2(1pe15.8,a))'    ) '              FWHM of Gauss sphere:', userdef_gfwhm,' AU =' &
                                                                           , userdef_gfwhm / userdef_distpc,' as'
  write (*   ,'(a,2(1pe15.8,a))'    ) '             Exponent of power law:', userdef_gplaw
  write (stdo,'(a,2(1pe15.8,a))'    ) '             Exponent of power law:', userdef_gplaw
  endif
  if (userdef_Rbesph .gt. almostzero) then
  write (*   ,'(a,2(1pe15.8,a))'    ) '               Radius of BE sphere:', userdef_Rbesph,' AU =' &
                                                                           , userdef_Rbesph / userdef_distpc,' as'
  write (stdo,'(a,2(1pe15.8,a))'    ) '               Radius of BE sphere:', userdef_Rbesph,' AU =' &
                                                                           , userdef_Rbesph / userdef_distpc,' as'
  write (*   ,'(a,1pe11.4,a)'       ) '      Central density of BE sphere:', userdef_besden,' g/cm^3'
  write (stdo,'(a,1pe11.4,a)'       ) '      Central density of BE sphere:', userdef_besden,' g/cm^3'
  write (*   ,'(a,1pe11.4,a)'       ) '          Temperature of BE sphere:', userdef_bestem,' K'
  write (stdo,'(a,1pe11.4,a)'       ) '          Temperature of BE sphere:', userdef_bestem,' K'
  write (*   ,'(a,1pe11.4,a)'       ) '               xi_max of BE sphere:', userdef_ximax
  write (stdo,'(a,1pe11.4,a)'       ) '               xi_max of BE sphere:', userdef_ximax
  endif
  if (userdef_Rplawx .gt. almostzero) then
  write (*   ,'(a,2(1pe15.8,a))'    ) '      Radius of power-law envelope:', userdef_Rplawx,' AU =' &
                                                                           , userdef_Rplawx / userdef_distpc,' as'
  write (stdo,'(a,2(1pe15.8,a))'    ) '      Radius of power-law envelope:', userdef_Rplawx,' AU =' &
                                                                           , userdef_Rplawx / userdef_distpc,' as'
  write (*   ,'(a,1pe15.8,a)'       ) '   Power-law exponent of densities:', userdef_dplawx
  write (stdo,'(a,1pe15.8,a)'       ) '   Power-law exponent of densities:', userdef_dplawx
  endif
  write (*   ,'()')
  write (stdo,'()')
  write (*   ,'(a,1pe15.8)'         ) ' Density factor of embedding cloud:', userdef_rhofac
  write (stdo,'(a,1pe15.8)'         ) ' Density factor of embedding cloud:', userdef_rhofac
  write (*   ,'(a,1pe15.8)'         ) '   Power-law exponent of densities:', userdef_dplout
  write (stdo,'(a,1pe15.8)'         ) '   Power-law exponent of densities:', userdef_dplout
  write (*   ,'(a,1pe11.4,a)'       ) '   Embedding cloud surface density:', userdef_embden,' cm^-2'
  write (stdo,'(a,1pe11.4,a)'       ) '   Embedding cloud surface density:', userdef_embden,' cm^-2'
  write (*   ,'()')
  write (stdo,'()')
  write (*   ,'(a,1pe10.3,a)'       ) '                        Model mass:', userdef_masstot,' Msun'
  write (stdo,'(a,1pe10.3,a)'       ) '                        Model mass:', userdef_masstot,' Msun'
  write (*   ,'(a,1pe10.3,a)'       ) '                Stellar luminosity:', userdef_lumstar,' Lsun'
  write (stdo,'(a,1pe10.3,a)'       ) '                Stellar luminosity:', userdef_lumstar,' Lsun'
  write (*   ,'(a,1pe10.3,a)'       ) '             Effective temperature:', userdef_temstar,' K'
  write (stdo,'(a,1pe10.3,a)'       ) '             Effective temperature:', userdef_temstar,' K'
  write (*   ,'()')
  write (stdo,'()')
  write (*   ,'(a,1pe15.8)'         ) '      Model scaling factor (radii):', userdef_mcrscale
  write (stdo,'(a,1pe15.8)'         ) '      Model scaling factor (radii):', userdef_mcrscale
  write (*   ,'(a,1pe15.8)'         ) '     Global scaling factor (radii):', userdef_radscale
  write (stdo,'(a,1pe15.8)'         ) '     Global scaling factor (radii):', userdef_radscale
  write (*   ,'(a,1pe15.8)'         ) ' Global scaling factor (densities):', userdef_denscale
  write (stdo,'(a,1pe15.8)'         ) ' Global scaling factor (densities):', userdef_denscale
  write (*   ,'(a,2(1pe15.8,a))'    ) '  Radius fraction for obs aperture:', fracaper
  write (stdo,'(a,2(1pe15.8,a))'    ) '  Radius fraction for obs aperture:', fracaper
  write (*   ,'(a,2(1pe15.8,a))'    ) ' Fraction of Router of empty space:', userdef_fraczero
  write (stdo,'(a,2(1pe15.8,a))'    ) ' Fraction of Router of empty space:', userdef_fraczero
  write (*   ,'()')
  write (stdo,'()')
  write (*   ,'(a,2(1pe15.8,a))'    ) '               Radius of the model:', Rmodel,' AU =', Rmodel / userdef_distpc,' as'
  write (stdo,'(a,2(1pe15.8,a))'    ) '               Radius of the model:', Rmodel,' AU =', Rmodel / userdef_distpc,' as'
  write (*   ,'(a,2(1pe15.8,a))'    ) '  Radius of observational aperture:', aperture(1) * userdef_distpc,' AU =', aperture(1),' as'
  write (stdo,'(a,2(1pe15.8,a))'    ) '  Radius of observational aperture:', aperture(1) * userdef_distpc,' AU =', aperture(1),' as'
  write (*   ,'(a,2(1pe15.8,a))'    ) '       Radius of outer empty space:', userdef_fraczero * Router,' AU =' &
                                                                          , userdef_fraczero * Router/ userdef_distpc,' as'
  write (stdo,'(a,2(1pe15.8,a))'    ) '       Radius of outer empty space:', userdef_fraczero * Router,' AU =' &
                                                                           , userdef_fraczero * Router/ userdef_distpc,' as'
  write (*   ,'()')
  write (stdo,'()')
  write (*   ,'(a,1pe15.8,a)'       ) '     Density at the inner boundary:', density(1) / userdef_dust2gas,' g/cm^3'
  write (stdo,'(a,1pe15.8,a)'       ) '     Density at the inner boundary:', density(1) / userdef_dust2gas,' g/cm^3'
  if (userdef_Rbesph .gt. almostzero) then
  if (userdef_ximax .gt. almostzero) then         ! BES critical density ratio must be 14.04
  write (*   ,'(a,1pe15.8,a,e10.3)') ' Density at the BE sphere boundary:', densRmodel,' g/cm^3  ratio:', bedenratio
  write (stdo,'(a,1pe15.8,a,e10.3)') ' Density at the BE sphere boundary:', densRmodel,' g/cm^3  ratio:', bedenratio
  else
  write (*   ,'(a,1pe15.8,a)'       ) ' Density at the BE sphere boundary:', densRmodel,' g/cm^3'
  write (stdo,'(a,1pe15.8,a)'       ) ' Density at the BE sphere boundary:', densRmodel,' g/cm^3'
  endif
  endif
  if (userdef_Rplawx .gt. almostzero) then
  write (*   ,'(a,1pe15.8,a)'       ) '  Density at the envelope boundary:', densRmodel,' g/cm^3'
  write (stdo,'(a,1pe15.8,a)'       ) '  Density at the envelope boundary:', densRmodel,' g/cm^3'
  endif
  write (*   ,'(a,1pe15.8,a)'       ) ' Density at computational boundary:', density(nzonesm1) / userdef_dust2gas,' g/cm^3'
  write (stdo,'(a,1pe15.8,a)'       ) ' Density at computational boundary:', density(nzonesm1) / userdef_dust2gas,' g/cm^3'

  write (*   ,'(a)') ' ______________________________________________________________________________'
  write (stdo,'(a)') ' ______________________________________________________________________________'
!__________________________________________________________________________________________________________________________________
!
  if (scattering_mode .eq. 0) then
    write (*   ,'()')
    write (stdo,'()')
    write (*   ,'(a)') ' No scattering requested for this run...'
    write (stdo,'(a)') ' No scattering requested for this run...'
  endif
  if (scattering_mode .eq. 1) then
    write (*   ,'()')
    write (stdo,'()')
    write (*   ,'(a)') ' Isotropic scattering requested for this run...'
    write (stdo,'(a)') ' Isotropic scattering requested for this run...'
  endif
  if (scattering_mode .eq. 2) then
    write (*   ,'()')
    write (stdo,'()')
    write (*   ,'(a)') ' Anisotropic scattering (Henyey-Greenstein) requested for this run...'
    write (stdo,'(a)') ' Anisotropic scattering (Henyey-Greenstein) requested for this run...'
  endif
  if (scattering_mode .eq. 3) then
    write (*   ,'()')
    write (stdo,'()')
    write (*   ,'(a)') ' Anisotropic scattering (tabulated phase function) requested for this run...'
    write (stdo,'(a)') ' Anisotropic scattering (tabulated phase function) requested for this run...'
  endif

  write (*   ,'()')
  write (stdo,'()')
  write (*   ,'(a,1pe15.8,a)') ' Total mass (gas+dust):', masstot / xMsun, ' Msun'
  write (stdo,'(a,1pe15.8,a)') ' Total mass (gas+dust):', masstot / xMsun, ' Msun'
  write (*   ,'(a,1pe15.8,a)') ' of embedding envelope:', masstot_emb / xMsun, ' Msun'
  write (stdo,'(a,1pe15.8,a)') ' of embedding envelope:', masstot_emb / xMsun, ' Msun'
!__________________________________________________________________________________________________________________________________
!
  starlumtot = z0
  do istar=1,nstars
    star_lum(istar) = z0
    do inu=1,freq_nr
      dum = pi * star_spec(inu,istar) * fourpi * star_r(istar)**2
      star_lum(istar) = star_lum(istar) + dum * freq_dnu(inu)
    enddo
    starlumtot = starlumtot + star_lum(istar)
  enddo

  if (userdef_external .eq. 1) then
    incl_extlum = 1
  else
    incl_extlum = 0
  endif

! Antonio Parravano, David J. Hollenbach, and Christopher F. McKee (2002)
! Time Dependence of the Ultraviolet Radiation Field in the Local Interstellar Medium:
!
! Far Ultraviolet (FUV, 6 eV< hν <13.6 eV) radiation has been recognized as the main source of heating of the neutral interstellar 
! gas, and, as a consequence, it determines whether the thermal balance of the neutral gas results in cold (T ∼ 50 − 100K) clouds 
! (CNM), warm (T ∼ 10^4K) clouds (WNM), or a combination of the two.
!
! The band [912 − 2070  ̊A], hereafter the “FUV band”. The mean energy density in this band is related to the parameter G0 that 
! is often used to parameterize the photoelectric dust heating rate (Tielens & Hollenbach 1985). Let Uband ≡ U(band)/∆band be 
! the average radiation energy density per  ̊A in the wavelength band ∆band. The energy density of FUV radiation is generally 
! expressed in units of the Habing (1968) typical energy density at the solar circle averaged between 6 eV≤ hν ≤ 13.6 eV 
! (i.e. [912 − 2070  ̊A]), G0 ≡ UFUV / UHFUV = UFUV / 4.6×10−17 erg cm^−3 ̊ A^−1. If the radiation field is one-dimensional, it 
! is convenient to express the Habing field as a flux, c × UHFUV × ∆FUV = 1.6 × 10−3 erg cm^−2 s^−1.
!
! Factor G0 computed as defined above (FUV: 912−2070 A):

  G0 = z0
  extlumtot = z0

  call read_externalsource ( 1 )

  do inu=1,freq_nr
    lamb(inu) = cc / freq_nu(inu) 
  enddo
  
  hotbbtemp = 22000.0d0

  if (incl_extlum .eq. 1) then
    isrflux = z0
    do inu=1,freq_nr
      extlum_intens(inu) = extlum_intens(inu) + userdef_hbbfac * bplanck ( hotbbtemp, freq_nu(inu) )
      dum = pi * extlum_intens(inu) * fourpi * grid_contsph_r**2
      extlumtot = extlumtot + dum * freq_dnu(inu)
      if (lamb(inu) .gt. 1.0d-5) then
        isrflux = isrflux + fourpi * extlum_intens(inu) * freq_dnu(inu)
      endif
    enddo
  endif

  open ( unit=1, file='isrf_intensity.tab' )
  do inu=1,freq_nr
    write (1,*) lamb(inu) * 1.0d4, extlum_intens(inu), userdef_hbbfac * bplanck ( hotbbtemp, freq_nu(inu) )
  enddo
  close ( 1 )

  lam1 = 0.0912d0 * 1.0d-4
  lam2 = 0.2070d0 * 1.0d-4

  call hunt ( lamb, freq_nr, lam1, k1 )
  call hunt ( lamb, freq_nr, lam2, k2 )

  if(k1 .lt. 1 .or. k1 .gt. freq_nr .or. k2 .lt. 1 .or. k2 .gt. freq_nr) then
     write (*,'(a         )') ' RADMC3D: USERDEF_SETUP_MODEL (7):'
     write (*,'(a,i0,1x,i0)') '   ERROR: Wave index out of range: ', k1, k2
     stop 99
  endif

  ifuv = z0
  if (incl_extlum .eq. 1) then
    eps1 = (lam1 - lamb(k1)) / (lamb(k1+1) - lamb(k1))
    eps2 = (lam2 - lamb(k2)) / (lamb(k2+1) - lamb(k2))
    int1 = (z1 - eps1) * extlum_intens(k1) + eps1 * extlum_intens(k1+1)
    int2 = (z1 - eps2) * extlum_intens(k2) + eps2 * extlum_intens(k2+1)
    ifuv = z05 * (int1 + extlum_intens(k1+1)) * (cc / lam1 - cc / lamb(k1+1))
    if (k2 - k1 .gt. 1) then
      do inu=k1+2,k2-1
        ifuv = ifuv + z05 * (extlum_intens(inu-1) + extlum_intens(inu)) * (cc / lamb(inu-1) - cc / lamb(inu))
      enddo
    endif
    ifuv = ifuv + z05 * (extlum_intens(k2) + int2) * (cc / lamb(k2) - cc / lam2)
  endif
  G0 = (fourpi * ifuv) / 1.6d-3

  write (*   ,'(a,1pe14.8,a)'      ) '    Central luminosity: ', starlumtot / xLsun, ' Lsun'
  write (stdo,'(a,1pe14.8,a)'      ) '    Central luminosity: ', starlumtot / xLsun, ' Lsun'
  write (*   ,'(a,1pe14.8,a,e10.3)') '   External luminosity: ', extlumtot  / xLsun, ' Lsun'
  write (stdo,'(a,1pe14.8,a,e10.3)') '   External luminosity: ', extlumtot  / xLsun, ' Lsun'
  write (*   ,'(a,1pe14.8,a,e10.3)') '  ISRF bolometric flux: ', isrflux, ' erg/s/cm^2  G0:', G0
  write (stdo,'(a,1pe14.8,a,e10.3)') '  ISRF bolometric flux: ', isrflux, ' erg/s/cm^2  G0:', G0

  write (*,'()')
  write (*,'(a)') ' Writing ''mc.odepths.tab''...'
  write (stdo,'()')
  write (stdo,'(a)') ' Writing ''mc.odepths.tab''...'
  tstamp = timestmp ()

  open ( unit=1, file='mc.odepths.tab' )
  write (1,'(a)') '! '//tstamp//"  Created by "//codeID
  write (1,'(a)') '! '//author
  write (1,'(a)') '!'
  write (1,'(a)') '! MODEL RADIAL OPTICAL DEPTHS AND OPACITIES'
  write (1,'(a)') '!'
  write (1,'(2(a,1pe15.8))') '! Distance (pc):', userdef_distpc, '  Mass-averaged Tdust (K):', tdust_massav
  write (1,'(a)') '!'
  write (1,'(a)') '!   n     WAVE(um)       ODEPTH       OPACITY(abs)   OPACITY(sca)   OPACITY(ext)'
  write (1,'(a)') '!'

!!  do inu=1,freq_nr
!!    do is=1,dust_nr_species
!!      dust_kappa_abs(inu,is) = userdef_opascale * dust_kappa_abs(inu,is)
!!      dust_kappa_scat(inu,is) = userdef_opascale * dust_kappa_scat(inu,is)
!!    enddo
!!  enddo

  odepth(:) = z0
  do inu=1,freq_nr
    do ic=1,nrcells
      idx = cellindex(ic)
      if (radius(ic) .lt. Rmodel) then
        dzone = (radius(ic+1) - radius(ic)) * xAU
        do is=1,dust_nr_species
          opacity = dust_kappa_abs(inu,is) + dust_kappa_scat(inu,is)
          odepth(inu) = odepth(inu) + z2 * opacity * dustdens(is,idx) * dzone
        enddo
      endif
    enddo
    opacity_abs = 0.0d0
    opacity_sca = 0.0d0
    opacity_ext = 0.0d0
    do is=1,dust_nr_species
      opacity_abs = opacity_abs + dust_kappa_abs(inu,is)
      opacity_sca = opacity_sca + dust_kappa_scat(inu,is)
      opacity_ext = opacity_ext + dust_kappa_abs(inu,is) + dust_kappa_scat(inu,is)
      odepth(inu) = odepth(inu) + z2 * opacity_ext * dustdens(is,idx) * dzone
    enddo
    write (1,'(i5,5(1pe15.7))') inu, lamb(inu) * 1.0d4, odepth(inu), opacity_abs, opacity_sca, opacity_ext
  enddo
  close ( 1 )

  lam1 = 0.55d0 * 1.0d-4
  lam2 = 100.0d0 * 1.0d-4

  call hunt ( lamb, freq_nr, lam1, k1 )
  call hunt ( lamb, freq_nr, lam2, k2 )

  eps1 = (lam1 - lamb(k1)) / (lamb(k1+1) - lamb(k1))
  eps2 = (lam2 - lamb(k2)) / (lamb(k2+1) - lamb(k2))
  odepthv = (z1 - eps1) * odepth(k1) + eps1 * odepth(k1+1)
  odepth100 = (z1 - eps2) * odepth(k2) + eps2 * odepth(k2+1)

  write (*   ,'(/a,1pe11.4)') ' Model radial optical depth (V band):', odepthv
  write (stdo,'(/a,1pe11.4)') ' Model radial optical depth (V band):', odepthv
  write (*   ,'( a,1pe11.4)') ' Model radial optical depth (100 µm):', odepth100
  write (stdo,'( a,1pe11.4)') ' Model radial optical depth (100 µm):', odepth100

  lam3 = 300.0d0 * 1.0d-4
  call hunt ( lamb, freq_nr, lam3, k3 )
  eps3 = (lam3 - lamb(k3)) / (lamb(k3+1) - lamb(k3))
  
  do is=1,dust_nr_species
    kappa300 = (z1 - eps3) * dust_kappa_abs(k3,is) + eps3 * dust_kappa_abs(k3+1,is)
    write (*   ,'(a,1pe11.4)') '     Reference dust opacity (300 µm):', kappa300
    write (stdo,'(a,1pe11.4)') '     Reference dust opacity (300 µm):', kappa300
  enddo

  if (.not.do_montecarlo_therm) then
  write (*   ,'(/a,1pe15.8,a)'     ) ' Mass-averaged dust temperature:', tdust_massav, ' K'
  write (stdo,'(/a,1pe15.8,a)'     ) ' Mass-averaged dust temperature:', tdust_massav, ' K'
  endif

  deallocate ( lamb, odepth )
  deallocate ( density, dslope, radius )

  write (*,'(/1x,a)') trim(tstamp)
  write (stdo,'(/1x,a/)') trim(tstamp)
  
  if (userdef_doback .eq. 1) then
    write (*   ,'()')
!!    write (stdo,'()')
    write (*   ,'(a)') esc//'[1;31m'//' Computing limb-brightened backgrounds'//esc//'[0m'//'...'
    write (stdo,'(a)') ' Computing limb-brightened backgrounds'//'...'
  endif
  if (userdef_tmassav .eq. 1) then
    if (userdef_doback .eq. 0) then
      write (*   ,'()')
      write (stdo,'()')
    endif
    write (*   ,'(a)') esc//'[1;31m'//' Replacing temperatures with mass-averaged value'//esc//'[0m'//'...'
    write (stdo,'(a)') ' Replacing temperatures with mass-averaged value'//'...'
  endif

  if (rt_mcparams%nphot_therm < 10) then
    write (*,'(/a )') ' RADMC3D: USERDEF_SETUP_MODEL (8):'
    write (*,'(/a/)') '   ERROR: Number of photons NPHOT < 10'
    write (stdo,'(/a )') ' RADMC3D: USERDEF_SETUP_MODEL (8):'
    write (stdo,'(/a/)') '   ERROR: Number of photons NPHOT < 10'
    call sleep (1) 
    stop 1
  endif
!__________________________________________________________________________________________________________________________________
!
end subroutine userdef_setup_model

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_dostuff ()
!__________________________________________________________________________________________________________________________________
!
! If you want to do some calculation for the model after the main calculation, you can do it here. 
! Here you can also write stuff to file. 
!__________________________________________________________________________________________________________________________________
!
  implicit none
  logical :: lexist, lbackexist, lunix = .true., lfitsimage = .true., limage, lsurfdens, ltausurf, ltau
  character*1 :: cunit
  character*2 :: cfr, cnu
  character*3 :: cnan
  character*4 :: cny
  character*5 :: cwave, fwave
  character*6 :: cfitsversion
  character*7 :: clibname
  character*8 :: ctime
  character*8 :: object = "RT Model", creator = "RADMC3D", ctype1 = "RA---TAN", ctype2 = "DEC--TAN", bunit = "MJy/sr", history = ""
  character*10 :: cdate
  character*16 :: cdust2gas, ckappa300, ctdust_massav, cmassval, cradiusau, cradiuspc, cradiusas, crspaceau, crspacepc, crspaceas
  character*100 :: ctsurf
  character*128 :: header, footer(2)
  character*500 :: imgname, imgnamebs, headline
  integer :: is, ic, idx, inu, nfreq, ix, iy, nx, ny, nxo, nyo, istat, reclen, ndate, inmx, ix1, iy1, istar, nfq, bgpix, k, idot &
           , ix0, iy0, inuap, nfreq_cam, npts, i1, j1, j, irad, iw, nm, nwin, blank = 0
  real :: fitsvers
  double precision :: fact, freq, cam_flux, cam_flux_Jy, dx, dy, farcs, dxas, dyas, isrfbg &
                    , denstot, dennumb, relsize, cam_back, cam_back_Jy, cam_stars, cam_stars_Jy, cam_flux_bgs, cam_flux_bgs_Jy &
                    , areax, fstars, flux, fluxbgs, backflux, backintens, eps, lambda(1000), lambap(1000), lamb(1000) &
                    , lamb_cam(1000), flux_cam(1000), intensbs(1000), fluxesbs(1000), dust_mass, tdust_massav, isrfflux &
                    , cam_isrf_Jy, rfact, Has, mxco, myco, amom, bmom, atheta, equivrad, equivsize, elongation
  double precision :: beam = 0.0d0, crpix1 = 0.0d0, crpix2 = 0.0d0, crval1 = 0.0d0, crval2 = 0.0d0, crota1 = 0.0d0, funtmp &
                    , crota2 = 0.0d0, epoch = 2.0d3, equinox = 2.0d3, bzero = 0.0d0, bscale = 1.0d0, datamin, datamax, dataminb &
                    , datamaxb, cd11 = 0.0d0, cd12 = 0.0d0, cd21 = 0.0d0, cd22 = 0.0d0, rann1, rann2, rdist, rx2, ry2, tsurf &
                    , tmpmax, massrad, dzone, dustflux, tmavflux, solidangle, kabs, temp, transparent_borderline, radcellcm &
                    , radcellas, radcellAU, radcellpc, massval, convmass, dlograd, densityfwhm, halfmax, selectk &
                    , masstot, radius
                    
  double precision, allocatable :: image(:,:), imagebs(:,:), arg(:), prof(:,:,:), odepth(:) &
                                 , odepthi(:,:), dtau(:,:), opacity(:,:), radiuscm(:), slope(:,:,:), slom1(:) &
                                 , slom2(:), slo1(:), slo2(:)
!__________________________________________________________________________________________________________________________________
!
  allocate ( radiuscm(userdef_nzones) )

  open ( unit=1, file='amr_grid.inp', status='unknown' )
  read (1,*) ! '1               ; iformat'
  read (1,*) ! '0               ; AMR grid style  (0 = regular grid, no AMR)'
  read (1,*) ! '100             ; Coordinate system'
  read (1,*) ! '0               ; gridinfo'
  read (1,*) ! '1,0,0           ; Include x,y,z coordinate'
  read (1,*) !nzonesm1,',1,1    ; Size of grid'
  do irad=1,userdef_nzones
     read (1,*) radiuscm(irad)
  enddo
  close (1)

  dust_mass = z0
  tdust_massav = z0
  do ic=1,nrcells
    idx = cellindex(ic)
    if (radiuscm(ic) .lt. Rmodel * xAU) then
      do is=1,dust_nr_species
        dust_mass = dust_mass + cellvolume(idx) * dustdens(is,idx)
        tdust_massav = tdust_massav + cellvolume(idx) * dustdens(is,idx) * dusttemp(is,idx)
      enddo
    endif
  enddo

  tdust_massav = tdust_massav / dust_mass
  massval = dust_mass / userdef_dust2gas
  if(amrray_mirror_equator) then
    massval = massval * z2
  endif
  
  if (do_montecarlo_therm) then
    open ( unit=1, file=trim(outname)//'.mc' )
    write (cdust2gas,'(1pe11.5)') userdef_dust2gas
    write (ckappa300,'(1pe11.5)') kappa300
    write (ctdust_massav,'(1pe11.5)') tdust_massav
    write (cmassval, '(1pe11.5)') massval / xMsun
    write (cradiusau,'(1pe11.5)') Rmodel
    write (cradiuspc,'(1pe11.5)') Rmodel * xAU / parsec
    write (cradiusas,'(1pe11.5)') Rmodel / userdef_distpc
    write (crspaceau,'(1pe11.5)') Router
    write (crspacepc,'(1pe11.5)') Router * xAU / parsec
    write (crspaceas,'(1pe11.5)') Router / userdef_distpc
    cdust2gas = trim ( adjustl ( cdust2gas ) )
    ckappa300 = trim ( adjustl ( ckappa300 ) )
    ctdust_massav = trim ( adjustl ( ctdust_massav ) )
    cmassval = trim ( adjustl ( cmassval ) )
    cradiusau = trim ( adjustl ( cradiusau ) )
    cradiuspc = trim ( adjustl ( cradiuspc ) )
    cradiusas = trim ( adjustl ( cradiusas ) )
    crspaceau = trim ( adjustl ( crspaceau ) )
    crspacepc = trim ( adjustl ( crspacepc ) )
    crspaceas = trim ( adjustl ( crspaceas ) )
    write (1,'(a)') '! Model properties:'
    write (1,'(a)') '               Dust-to-gas mass ratio: '//cdust2gas(1:11)
    write (1,'(a)') ' Reference opacity at 300 µm (cm^2/g): '//ckappa300(1:11)
    write (1,'(a)') '  Masss-averaged dust temperature (K): '//ctdust_massav(1:11)
    write (1,'(a)') '                    Model mass (Msun): '//cmassval(1:11)
    write (1,'(a)') '                    Model radius (AU): '//cradiusau(1:11)
    write (1,'(a)') '                    Model radius (pc): '//cradiuspc(1:11)
    write (1,'(a)') '                    Model radius (as): '//cradiusas(1:11)
    write (1,'(a)') '      Computational space radius (AU): '//crspaceau(1:11)
    write (1,'(a)') '      Computational space radius (pc): '//crspacepc(1:11)
    write (1,'(a)') '      Computational space radius (as): '//crspaceas(1:11)
    close (1)
    write (*,'(/a,1pe15.8,a)') ' Mass-averaged dust temperature:', tdust_massav, ' K'
    write (*,'()')
    write (*,'(a)') ' Writing ''mc.denstemps.tab''...'
    write (stdo,'()')
    write (stdo,'(a)') ' Writing ''mc.denstemps.tab''...'
  endif

  if (userdef_doback .eq. 0) then
    open ( unit=1, file='mc.denstemps.tab' )
  else
    open ( unit=1, file='mc.denstemps.bg.tab' )
  endif
  write (1,'(a)') '! '//tstamp//"  Created by "//codeID
  write (1,'(a)') '! '//author
  write (1,'(a)') '!'
  write (1,'(a)') '! PROFILES OF DENSITIES & TEMPERATURES'
  write (1,'(a)') '!'
  write (1,'(a)') '! Note: radii in this table are cell-centered'
  write (1,'(a)') '!'
  write (1,'(2(a,1pe15.8))') '! Distance (pc):', userdef_distpc, '  Mass-averaged Tdust (K):', tdust_massav
  write (1,'(a)') '!'
  write (1,'(a)') '!   N    RADIUS(cm)    RELZONESIZE     RADIUS(au)     RADIUS(pc)     RADIUS(as)    DENS(g/cm^3)' &
                //'  NUMDEN(1/cm^3)   DENSDUST       NUMDENDUST     MASSFRAC      TEMPERATURE'
  write (1,'(a)') '!'

  dust_mass = z0
  massval = z0
  densityfwhm = z0
  do ic=1,nrcells
    idx = cellindex(ic)
    radcellcm = (radiuscm(ic) + radiuscm(ic+1)) / z2
    radcellAU = radcellcm / xAU
    radcellpc = radcellcm / parsec
    radcellas = radcellAU / userdef_distpc
    denstot = sum ( dustdens(:,idx) ) / userdef_dust2gas
    dennumb = denstot / (muH2 * AtomicMU)
    relsize = (radiuscm(ic+1) - radiuscm(ic)) / radcellcm
    do is=1,dust_nr_species
      dust_mass = dust_mass + cellvolume(idx) * dustdens(is,idx)
    enddo
    massrad = dust_mass / userdef_dust2gas
    if(amrray_mirror_equator) then
      massrad = massrad * z2
    endif
    if (radcellcm .le. Rmodel * xAU) then
      massval = massrad
    endif
    write (1,'(i5,12(1pe15.7))') idx, radcellcm, relsize, radcellAU, radcellpc, radcellas, &
                                 max ( denstot, almostzero ), max ( dennumb, almostzero ), &
                                (max ( dustdens(is,idx), almostzero ), &
                                 max ( dustdens(is,idx) / (muH2 * AtomicMU), almostzero ), &
                                 massrad / masstot, max ( dusttemp(is,idx), almostzero ), is=1,dust_nr_species)
                                 
! Find actual FWHM of volume densities.

    halfmax = dustdens(is,cellindex(1)) / 2.0d0

    if (dustdens(is,idx) .lt. halfmax .and. dustdens(is,idx-1) .ge. halfmax) then
      eps = (halfmax - dustdens(is,idx-1)) / (dustdens(is,idx) - dustdens(is,idx-1))
      densityfwhm = 2.0d0 * ((z1 - eps) * radiuscm(ic-1) + eps * radiuscm(ic)) / xAU / userdef_distpc
    endif
  enddo

  close ( 1 )
!__________________________________________________________________________________________________________________________________
!
  if (do_raytrace_image .or. do_raytrace_tausurf .or. do_raytrace_tau) then

    inquire ( file='image.out', exist=lexist )

    if (lexist) then
      open ( unit=2, file='image.out', status='old' )
      read (2,*)
      read (2,*) nxo, nyo
      read (2,*) nfreq, tsurf
      read (2,*) dx, dy
      read (2,*) (lambda(inu), inu=1,nfreq)
      read (2,*)

      if (do_raytrace_tausurf) then
        ltau = .false.
        limage = .false.
        lsurfdens = .false.
        ltausurf = .true.
      elseif (do_raytrace_tau) then
        ltau = .true.
        limage = .false.
        lsurfdens = .false.
        ltausurf = .false.
      elseif (nfreq .eq. 1 .and. nint ( lambda(1) ) .eq. 999) then
        ltau = .false.
        limage = .false.
        lsurfdens = .true.
        ltausurf = .false.
      else
        ltau = .false.
        limage = .true.
        lsurfdens = .false.
        ltausurf = .false.
      endif

      nx = nxo
      ny = nyo
      ix1 = (nxo - nx) / 2
      iy1 = (nyo - ny) / 2
      write (cny,'(i4)') ny
      write (cfr,'(i2)') nfreq
      farcs = z1 / userdef_distpc
      dxas = dx / xAU * farcs
      dyas = dy / xAU * farcs
      cd11 =-dxas
      cd22 = dyas
      crpix1 = (nx - 1) / 2 + 1
      crpix2 = (ny - 1) / 2 + 1
      inmx = max ( nx, ny )
      npts = min ( nx, ny ) / 2
      i1 = nint ( crpix1 )
      j1 = nint ( crpix2 )

      inquire ( file='aperture_info.inp', exist=lexist )
      if (lexist) then
        open ( unit=3, file='aperture_info.inp', status='old' )
        read (3,*)
        read (3,*) nfq
        do inu=1,nfq
          read (3,*) lambap(inu), aperture(inu)
        enddo
        close ( 3 )
      endif

      ix0 = (nxo - 1) / 2 + 1
      iy0 = (nyo - 1) / 2 + 1
      intensbs(:) = z0
      fluxesbs(:) = z0

      allocate ( prof(nfreq,npts+1+10,2), arg(npts+1+10), slope(nfreq,npts+1,2), slo1(npts+1), slo2(npts+1), slom1(npts+1), &
                 slom2(npts+1), stat=istat )

      if (istat > 0) then
         write (*,'(a)') ' RADMC3D: USERDEF_DOSTUFF (1):'
         write (*,'(a)') '   ERROR: Unable to allocate memory.'
         stop
      endif    
      prof(:,:,:) = 0.0d0

      if (lfitsimage) then
        call osystem ( lunix )
        call ftvers ( fitsvers )
        write (cfitsversion,'(f6.3)') fitsvers
        clibname = 'CFITSIO'
      else
        write (*,'(a)') ' Writing ''mc.images.tab''...'
        reclen = nxo * 11
        open ( unit=1, file='mc.images.tab', recl=reclen )
        write (1,'(a)') '! '//tstamp//"  Created by "//codeID
        write (1,'(a)') '! '//author
        write (1,'(a)') '!'
        write (1,'(a)') '! IMAGES AT SELECTED WAVELENGTHS (MJy/sr)'
        write (1,'(a)') '!'
        write (1,'(2(a,1pe15.8))') '! Distance (pc):', userdef_distpc, '  Mass-averaged Tdust (K):', tdust_massav
        write (1,'(a)') '!'
        write (1,'(a, i2,a)') '! number of wavelengths: ', nfreq
        write (1,'(a,2i4,a)') '! numbers of pixels: ', nxo, nyo
        write (1,'(a,2(1pe15.8),a)') '! pixel sizes (cm): ', dx, dy
        write (1,'(a,'//cfr//'(1pe15.8),a)') '! wavelengths (µm): ', (lambda(inu), inu=1,nfreq)
        write (1,'(a)') '!'
      endif

      do inu=1,nfreq

        allocate ( image(nxo,nyo), imagebs(nxo,nyo), stat=istat )
  
        if (istat > 0) then
           write (*,'(a)') ' RADMC3D: USERDEF_DOSTUFF (1):'
           write (*,'(a)') '   ERROR: Unable to allocate memory.'
           stop
        endif    

        do iy=1,nyo
          do ix=1,nxo
            read (2,*) image(ix,iy)
          enddo
        enddo
        read (2,*)

        isrfbg = image(1,1)

        call hunt ( lambap, nfq, lambda(inu), inuap )

        if(inuap .lt. 1 .or. inuap .gt. nfq) then
           write (*,'(a   )') ' RADMC3D: USERDEF_DOSTUFF (2):'
           write (*,'(a,i0)') '   ERROR: Wave index out of range: ', inuap
           stop 99
        endif

! Circular annulus for background evaluation is defined by the aperture radius with thickness of one pixel.

        eps = (log(lambda(inu)) - log(lambap(inuap))) / (log(lambap(inuap+1)) - log(lambap(inuap)))
        rann1 = exp ( (z1 - eps) * log ( aperture(inuap) ) + eps * log ( aperture(inuap+1) ) ) + dxas
        rann2 = rann1 + dxas

        backintens = z0
        bgpix = z0
        do iy=1,nyo
          ry2 = (dyas * dble (iy - iy0))**2
          do ix=1,nxo
            rx2 = (dxas * dble (ix - ix0))**2
            rdist = sqrt ( rx2 + ry2 )

            if (rdist .gt. rann1 .and. rdist .le. rann2) then
              if (.not.(ltausurf .and. image(ix,iy) / xAU .le. -1.0d10)) then
                bgpix = bgpix + 1
                backintens = backintens + image(ix,iy)
              endif
            endif
          enddo
        enddo
        if (bgpix .gt. z0) backintens = backintens / bgpix

        fact = dx * dy / (userdef_distpc * parsec)**2
        convmass = dxas * dyas * (userdef_distpc * xAU)**2 * AtomicMU * muH2 / xMsun

        do iy=1,nyo
          ry2 = (dyas * dble (iy - iy0))**2
          do ix=1,nxo
            rx2 = (dxas * dble (ix - ix0))**2
            rdist = sqrt ( rx2 + ry2 )

            imagebs(ix,iy) = image(ix,iy) - backintens  ! <- backintens contains ISRF background
            if (rdist .gt. aperture(inuap)) then
              imagebs(ix,iy) = z0
            endif

            image(ix,iy) = max ( image(ix,iy) - isrfbg, almostzero )  ! <- ISRF background is to be subtracted

            if (limage .and. userdef_doback .eq. 0) then
              fluxesbs(inu) = fluxesbs(inu) + imagebs(ix,iy) * fact / xJy
            endif
            if (lsurfdens .and. userdef_doback .eq. 0) then
              fluxesbs(inu) = fluxesbs(inu) + imagebs(ix,iy) / userdef_dust2gas / (muH2 * AtomicMU) * convmass
            endif
          enddo
        enddo
        if (limage .and. userdef_doback .eq. 0) then
          intensbs(inu) = imagebs(ix0,iy0) / xJy / 1.0d6
          fluxesbs(inu) = max ( fluxesbs(inu), 1.0d-31 )
        endif
        if (lsurfdens .and. userdef_doback .eq. 0) then
          intensbs(inu) = imagebs(ix0,iy0) / userdef_dust2gas / (muH2 * AtomicMU) / 1.0d20
          fluxesbs(inu) = max ( fluxesbs(inu), 1.0d-31 )
        endif
      
        if (inu .eq. 1) write (*,'()')

!!      if (limage .and. userdef_doback .eq. 0) then
!!
!!        open ( unit=4, file='mc.fluxes.sed.cat', recl=1000 )
!!        write (*,'(a)') ' Writing ''mc.fluxes.sed.cat''...'
!!        write (cnu,'(i2)') nfreq
!!        write (4,'(a)') '! '//tstamp//"  Created by "//codeID
!!        write (4,'(a)') '! '//author
!!        write (4,'(a)') '!'
!!        write (4,'(a)') '! CATALOG OF TOTAL FLUXES FOR SED FITTING'
!!        write (4,'(a)') '!'
!!        write (4,'(a)') '! NO = 1: flux measured by ray tracing for BS images (surfdens mask not used)'
!!        write (4,'(a)') '! NO = 2: flux measured by ray tracing for SEDs'
!!        write (4,'(a)') '!'
!!        write (4,'(2(a,1pe15.8))') '! Distance (pc):', userdef_distpc, '  Mass-averaged Tdust (K):', tdust_massav
!!        write (4,'(a)') '!'
!!        write (4,'(a,'//cnu//'(1x,a1,i2.2,a1,1pe10.3))') '! Wavelengths:', ('(',inu,')', lambda(inu), inu=1,nfreq)
!!        write (4,'(a)') '!'
!!        headline = '!   NO     XCO_P   YCO_P    SIGN    GOOD'
!!        do inu=1,nfreq
!!           write (cnu,'(i2.2)') inu
!!           headline = trim(headline)//' SIG_MONO'//cnu//'  FXP_BST'//cnu//'  FXP_ERR'//cnu//'  FXT_BST'//cnu//'  FXT_ERR'//cnu
!!        enddo
!!        write (4,'(a)') trim(headline)
!!        write (4,'(a)') '!'
!!
!!        write (4,'(i6,1x,2f8.1,1pe11.3,0pf6.2,'//cnu//'(1pe11.3,4(1pe11.3)))') &
!!               1,z0,z0,z100,z10, (z100, intensbs(inu), intensbs(inu) / z100, fluxesbs(inu), fluxesbs(inu) / z100, inu=1,nfreq)
!!        close (4)
!!      endif

        if (lfitsimage) then
        
          do k=1,npts+1
            arg(k) = min ( radiuscm(1) / xAU * farcs, z025 * dxas ) + dxas * dble ( k - 1 )
          enddo
        
!!          if (limage .and. userdef_doback .eq. 0) then
!!            open ( unit=1, file=trim(outname)//'.im' )
!!            write (1,'(a)') '! Model sizes (arcsec), peaks (MJy/sr), and fluxes (Jy):'
!!            write (1,'(a)') '!              halfmax     moments    footprint     peakint     totfluxes'
!!          endif
        
!!          do inu=1,nfreq
        
            call hunt ( lambap, nfq, lambda(inu), inuap )
        
            if(inuap .lt. 1 .or. inuap .gt. nfq) then
               write (*,'(a   )') ' RADMC3D: USERDEF_DOSTUFF (3):'
               write (*,'(a,i0)') '   ERROR: Wave index out of range: ', inuap
               stop 99
            endif
        
            if (lambda(inu) .lt. 1.0d3) then
              write (fwave,'(f5.1)') lambda(inu)
              fwave(4:4) = 'p'
              cunit = 'u'
              if (fwave(1:1) .eq. ' ') fwave(1:1) = '0'
              if (fwave(2:2) .eq. ' ') fwave(2:2) = '0'
            else
              if (lambda(inu) .lt. 1.0d4) then
                write (cwave,'(f5.2)') lambda(inu) / 1.0d3
                cunit = 'm'
              else
                write (cwave,'(f5.2)') lambda(inu) / 1.0d4
                cunit = 'c'
              endif
              fwave(1:1) = '_'
              fwave(2:5) = cwave(2:5)
              fwave(3:3) = 'p'
            endif
            iw = index ( fwave, 'p' ) - 1
        
            datamin =  1.0d90
            datamax = -1.0d90
            dataminb =  1.0d90
            datamaxb = -1.0d90
            tmpmax = -1.0d90
            cnan = 'NAN'
        
            do iy=1,ny
              do ix=1,nx
                if (limage) then
                  image(ix,iy) = image(ix,iy) / xJy / 1.0d6
                  imagebs(ix,iy) = imagebs(ix,iy) / xJy / 1.0d6
                endif
                if (lsurfdens) then
                  image(ix,iy) = image(ix,iy) / userdef_dust2gas / (muH2 * AtomicMU)
                  imagebs(ix,iy) = imagebs(ix,iy) / userdef_dust2gas / (muH2 * AtomicMU)
                endif
                if (ltausurf) then
                  funtmp = image(ix,iy) / xAU
                  image(ix,iy) = max ( funtmp, -1.0d10 )
                  if (funtmp .lt. -1.0d10) read(cnan,*) image(ix,iy)
                  funtmp = imagebs(ix,iy) / xAU
                  imagebs(ix,iy) = max ( funtmp, -1.0d10 )
                  if (funtmp .lt. -1.0d10) read(cnan,*) imagebs(ix,iy)
                  if (funtmp .gt. -1.0d10 .and. funtmp > tmpmax) tmpmax = funtmp
                endif
                if (image(ix,iy) < datamin) datamin = image(ix,iy)
                if (image(ix,iy) > datamax) datamax = image(ix,iy)
                if (imagebs(ix,iy) < dataminb) dataminb = imagebs(ix,iy)
                if (imagebs(ix,iy) > datamaxb) datamaxb = imagebs(ix,iy)
              enddo
            enddo
        
            if (userdef_doback .eq. 0) then
              call halfmaxsizes ( nx, ny, imagebs, dxas, dyas, Has, crpix1, crpix2 )
              call momentsizes ( 1, nx, 1, ny, nx, ny, imagebs, dxas, dyas, mxco, myco, amom, bmom, atheta, equivrad, elongation )
              equivsize = z2 * equivrad
            endif
        
            if (limage) then
              if (userdef_doback .eq. 1) then
                imgname = 'mc.'//fwave(1:iw)//cunit//'m.bgl.fits'
              else
                imgname = 'mc.'//fwave(1:iw)//cunit//'m.fits'
!!                write (cHas,'(f16.9)') Has
!!                write (csizeamom,'(f16.9)') sqrt ( amom * bmom)
!!                write (cequivsize,'(f16.9)') equivsize
!!                write (cintens,'(f16.9)') intensbs(inu)
!!                write (cflux,'(f16.9)') fluxesbs(inu)
!!                cHas = trim ( adjustl ( cHas ) )
!!                csizeamom = trim ( adjustl ( csizeamom ) )
!!                cequivsize = trim ( adjustl ( cequivsize ) )
!!                cintens = trim ( adjustl ( cintens ) )
!!                cflux = trim ( adjustl ( cflux ) )
!!                write (1,'(a,3(f12.8),a)') ' '//fwave(1:iw)//' µm bs   '//cHas(1:11)//' '//csizeamom(1:11) &
!!                                         //' '//cequivsize(1:11)//'  '//cintens(1:11)//'  '//cflux(1:11)
              endif
              imgnamebs = 'mc.'//fwave(1:iw)//cunit//'m.bs.fits'
              header = 'Model intensity at '//trim ( fwave(1:iw) )//' µm'
              footer(1) = 'Model intensity at '//trim ( fwave(1:iw) )//' µm (MJy|S| |N|sr|S|-1|N|)'
              bunit = "MJy/sr"
            endif
            if (lsurfdens) then
              if (userdef_doback .eq. 1) then
                imgname = 'mc.surfdens.bgl.fits'
              else
                imgname = 'mc.surfdens.fits'
!!                write (cHas,'(f16.9)') Has
!!                write (csizeamom,'(f16.9)') sqrt ( amom * bmom)
!!                write (cequivsize,'(f16.9)') equivsize
!!                write (cintens,'(f16.9)') intensbs(inu)
!!                write (cmass,'(f16.9)') fluxesbs(inu)
!!                cHas = trim ( adjustl ( cHas ) )
!!                csizeamom = trim ( adjustl ( csizeamom ) )
!!                cequivsize = trim ( adjustl ( cequivsize ) )
!!                cintens = trim ( adjustl ( cintens ) )
!!                cmass = trim ( adjustl ( cmass ) )
!!                open ( unit=1, file=trim(outname)//'.sd' )
!!                write (1,'(a)') '! Model sizes (arcsec), peaks (H2/cm^2/1e20), and masses (Msun):'
!!                write (1,'(a)') '!              halfmax     moments     footprint   peaksdens     sdmass'
!!                write (1,'(a)') '  sdens bs   '//cHas(1:11)//' '//csizeamom(1:11)//' '//cequivsize(1:11) &
!!                              //'  '//cintens(1:11)//'  '//cmass(1:11)
!!                close (1)
              endif
              imgnamebs = 'mc.surfdens.bs.fits'
              header = 'Model surface densities N_H2'
              footer(1) = 'Model surface densities N_H2 (g|S| |N|cm|S|-3|N|)'
              bunit = "g/cm^3"
            endif
            if (ltausurf) then
              if (tsurf .ge. z1) write (ctsurf,'(f0.1)') tsurf
              if (tsurf .lt. z1) write (ctsurf,'(f6.4)') tsurf
              idot = index ( ctsurf, '.' )
              if (idot .gt. 0) ctsurf(idot:idot) = 'p'
              imgname = 'mc.tau.'//trim ( ctsurf )//'.surf.'//fwave(1:iw)//cunit//'m.fits'
              imgnamebs = 'mc.tau.'//trim ( ctsurf )//'.surf.'//fwave(1:iw)//cunit//'.bs.fits'
              header = 'Model surface of tau='//trim ( ctsurf )//' at '//trim ( fwave(1:iw) )//' µm'
              footer(1) = 'Model surface of tau='//trim ( ctsurf )//' at '//trim ( fwave(1:iw) )//' µm'
              bunit = "AU"
            endif
            if (ltau) then
              imgname = 'mc.odepth.'//fwave(1:iw)//cunit//'m.fits'
              imgnamebs = 'mc.odepth.'//fwave(1:iw)//cunit//'m.bs.fits'
              header = 'Model optical depths at '//trim ( fwave(1:iw) )//' µm'
              footer(1) = 'Model optical depths at '//trim ( fwave(1:iw) )//' µm'
              bunit = ""
            endif
        
            call when ( lunix, ctime, cdate, ndate, 4 )
        
!!            write (*,'(a)') ' Writing ''image.fits''...'
!!            call wfits ( cfitsversion, nxo, nxo, nyo, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, image, dxas, dyas,  &
!!                         object, crval1, crval2, 'image.fits', header, footer, cdate, ctime, creator, beam, blank, crota1, & 
!!                         crota2, cd11, cd12, cd21, cd22, epoch, equinox, bzero, bscale, lambda(inu), datamin, datamax, history )
        
            if (.not.ltau .and. .not.(ltausurf .and. datamax .le. almostzero)) then
              write (*,'(a)') ' Writing '''//trim(imgname)//'''...'
              call wfits ( cfitsversion, inmx, nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, image, dxas, dyas,  &
                           object, crval1, crval2, trim(imgname), header, footer, cdate, ctime, creator, beam, blank, crota1, & 
                           crota2, cd11, cd12, cd21, cd22, epoch, equinox, bzero, bscale, lambda(inu), datamin, datamax, history )
            else
              if (ltausurf .and. datamax .le. almostzero) then
                write (*,'(a,1pe8.1)') ' Writing '''//trim(imgname)//''' ~> skipping: tau <', tsurf
              endif
            endif
        
!!            if (.not.ltausurf) then
!!              write (*,'(a)') ' Writing '''//trim(imgnamebs)//'''...'
!!              call wfits ( cfitsversion, inmx, nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, imagebs, dxas, dyas,  &
!!                           object, crval1, crval2, trim(imgnamebs), header, footer, cdate, ctime, creator, beam, blank, crota1, &
!!                           crota2, cd11, cd12, cd21, cd22, epoch, equinox, bzero, bscale, lambda(inu), dataminb, datamaxb, history )
!!            else
!!              if (.not.ltau .and. .not.(ltausurf .and. datamax .le. almostzero)) then
!!                write (*,'(a)')
!!              endif
!!            endif
            
            iy = j1
            do k=1,npts
              ix = i1 + k - 1
              if (isnan ( image(ix,iy) )) then
                prof(inu,k,1) = z0
              else
                prof(inu,k,1) = image(ix,iy)
              endif
              if (isnan ( imagebs(ix,iy) )) then
                prof(inu,k,2) = z0
              else
                if (userdef_doback .eq. 1) then
                  prof(inu,k,2) = (image(ix,iy) - imagebs(ix,iy))
                else
                  prof(inu,k,2) = max ( imagebs(ix,iy), 1.0d-99 )
                endif
              endif
            enddo
            
! Comp  ute logarithmic slopes of the profiles.
        
            slope(inu,1,1) = 1.0d-99
            slope(inu,1,2) = 1.0d-99
            if (npts .gt. 1) then
              do k=2,npts-1
                dlograd = log10 ( arg(k) ) - log10 ( arg(k-1) )
                if (prof(inu,k-1,1) .gt. almostzero .and. prof(inu,k,1) .gt. almostzero) then
                  slope(inu,k,1) = abs ( (log10 ( prof(inu,k,1) ) - log10 ( prof(inu,k-1,1) )) / dlograd )
                else
                  slope(inu,k,1) = 1.0d-99
                endif
                if (abs ( slope(inu,k,1) ) .lt. 1.0d-99) slope(inu,k,1) = 1.0d-99
        
                if (prof(inu,k-1,2) .gt. almostzero .and. prof(inu,k,2) .gt. almostzero) then
                  slope(inu,k,2) = abs ( (log10 ( prof(inu,k,2) ) - log10 ( prof(inu,k-1,2) )) / dlograd )
                else
                  slope(inu,k,2) = 1.0d-99
                endif
                if (abs ( slope(inu,k,2) ) .lt. 1.0d-99) slope(inu,k,2) = 1.0d-99
              enddo
              slope(inu,1,1) = slope(inu,2,1)
              slope(inu,1,2) = slope(inu,2,2)
              slope(inu,npts,1) = slope(inu,npts-1,1)
              slope(inu,npts,2) = slope(inu,npts-1,2)
        
! Smoo  th the slopes by median filtering.
        
              nwin = 3
              do k=1,npts
                nm = 0
                do j=max(k-nwin,1),k
                  nm = nm + 1
                  slo1(nm) = slope(inu,j,1)
                  slo2(nm) = slope(inu,j,2)
                enddo
                if (nm .ge. 1) then
                  slom1(k) = selectk ( (nm + 1) / 2, nm, slo1 )
                  slom2(k) = selectk ( (nm + 1) / 2, nm, slo2 )
                else
                  slom1(k) = slope(inu,k,1)
                  slom2(k) = slope(inu,k,2)
                endif
              enddo
              do k=1,npts
                slope(inu,k,1) = slom1(k)
                slope(inu,k,2) = slom2(k)
              enddo
            endif
        
!!          enddo
!!          if (limage .and. userdef_doback .eq. 0) then
!!            close ( 1 )
!!          endif
        else
          write (1,'(1pe15.8,a)') lambda(inu), ' ! wavelength (µm)'
          do ix=1,nx
            write (1,'('//cny//'(1pe11.4))') (image(ix,iy), iy=1,ny)
          enddo
          close (1)
        endif

        deallocate ( image, imagebs )
      enddo

      close (2)

      tstamp = timestmp ()
      if (limage) then
        if (userdef_doback .eq. 1) then
          open ( unit=1, file='mc.intens.pro.bg.tab' )
          write (*,'(a)') ' Writing ''mc.intens.pro.bg.tab''...'
        else
          open ( unit=1, file='mc.intens.pro.tab' )
          write (*,'(a)') ' Writing ''mc.intens.pro.tab''...'
        endif
      endif
      
      if (lsurfdens) then
        if (userdef_doback .eq. 1) then
          open ( unit=1, file='mc.surfdens.pro.bg.tab' )
          write (*,'(a)') ' Writing ''mc.surfdens.pro.bg.tab''...'
        else
          open ( unit=1, file='mc.surfdens.pro.tab' )
          write (*,'(a)') ' Writing ''mc.surfdens.pro.tab''...'
        endif
      endif
      if (ltausurf .and. datamax .gt. almostzero) then
        open ( unit=1, file='mc.tau'//trim(ctsurf)//'.surf.pro.tab' )
        write (*,'(a)') ' Writing ''mc.tau'//trim(ctsurf)//'.surf.pro.tab''...'
      else 
        if (ltausurf) then
          write (*,'(a,1pe8.1)') ' Writing ''mc.tau'//trim(ctsurf)//'.surf.pro.tab'' ~> skipping: tau <', tsurf
        endif
      endif
      if (ltau) then
        open ( unit=1, file='mc.odepth.images.pro.tab' )
        write (*,'(a)') ' Writing ''mc.odepth.images.pro.tab''...'
      endif

      if (.not.ltausurf .or. (ltausurf .and. datamax .gt. almostzero)) then
        write (1,'(a)') '! '//tstamp//"  Created by "//codeID
        write (1,'(a)') '! '//author
        write (1,'(a)') '!'
        if (limage) then
          if (userdef_doback .eq. 1) then
            write (1,'(a)') '! RADIAL PROFILES OF BACKGROUND INTENSITIES (MJy/sr)'
          else
            write (1,'(a)') '! RADIAL PROFILES OF INTENSITIES (MJy/sr)'
          endif
        endif
        if (lsurfdens) then
          if (userdef_doback .eq. 1) then
            write (1,'(a)') '! RADIAL PROFILES OF BACKGROUND SURFACE DENSITIES (cm^-2)'
          else
            write (1,'(a)') '! RADIAL PROFILES OF SURFACE DENSITIES (cm^-2)'
          endif
        endif
        if (ltausurf) then
          if (tsurf .ge. z1) write (ctsurf,'(f0.1)') tsurf
          if (tsurf .lt. z1) write (ctsurf,'(f6.4)') tsurf
          write (1,'(a)') '! RADIAL PROFILES OF A TAU = '//trim(ctsurf)//' SURFACE'
        endif
        if (ltau) then
          write (1,'(a)') '! RADIAL PROFILES OF OPTICAL DEPTHS'
        endif
        write (1,'(a)') '!'
        write (1,'(2(a,1pe15.8))') '! Distance (pc):', userdef_distpc, '  Mass-averaged Tdust (K):', tdust_massav
        write (1,'(a)') '!'
        if (.not.lsurfdens) then
        write (cnu,'(i2)') nfreq
        write (1,'(a,'//cnu//'(1x,a1,i2.2,a1,1pe10.3))') '! Wavelengths:', ('(',inu,')', lambda(inu), inu=1,nfreq)
        write (1,'(a)') '!'
        endif
        headline = '!   N   RADIUS(au)     RADIUS(pc)     RADIUS(as)'
        do inu=1,nfreq
          write (cnu,'(i2.2)') inu
          if (.not.lsurfdens) then
            if (userdef_doback .eq. 1) then
              headline = trim(headline)//'    PROFILEBG('//cnu//')   SLOPEBG('//cnu//')   PROFILEB('//cnu//')    SLOPEB('//cnu//')'
            else
              if (inu .eq. 1) then
             headline = trim(headline)//'     PROFILE('//cnu//')     SLOPE('//cnu//')    PROFILEBS('//cnu//')   SLOPEBS('//cnu//')'
              else
             headline = trim(headline)//'    PROFILE('//cnu//')     SLOPE('//cnu//')    PROFILEBS('//cnu//')   SLOPEBS('//cnu//')'
              endif
            endif
          else
            if (userdef_doback .eq. 1) then
              headline = trim(headline)//'      PROFILEBG      SLOPEBG        PROFILEB       SLOPEB'
            else
              headline = trim(headline)//'      PROFILE        SLOPE          PROFILEBS      SLOPEBS'
              !!'        SDENS' &
              !!                                        //'          SDENSFG        SDENSX         SDENSXFGX      XSLOPE'
            endif
          endif
        enddo
        write (1,'(a)') trim(headline)
        write (1,'(a)') '!'

        do k=1,npts
          if (k .gt. npts) arg(k) = arg(k-1) * 1.08
          radius = arg(k) * userdef_distpc
          write (1,'(i5,200(1pe15.7))') k, radius, radius * xAU / parsec, arg(k) &
                                      , (prof(inu,k,1), slope(inu,k,1), prof(inu,k,2), slope(inu,k,2), inu=1,nfreq)
        enddo

        close ( 1 )
      endif
      
      deallocate ( prof, arg, slope, slo1, slo2, slom1, slom2 )
    endif
  endif  !<~ do_raytrace_image .or. do_raytrace_tausurf .or. do_raytrace_tau
!__________________________________________________________________________________________________________________________________
!
  if (do_raytrace_spectrum) then

    inquire ( file='spectrum.out', exist=lexist )

    if (lexist) then
      open ( unit=2, file='spectrum.out', status='old' )
      inquire ( file='specback.out', exist=lbackexist )
      if (lbackexist) then
        open ( unit=3, file='specback.out', status='old' )
      endif
      read (2,*)
      read (2,*) nfreq
      read (2,*)
      if (lbackexist) then
        read (3,*)
      endif

      write (*,'()')
      write (*,'(a)') ' Writing ''mc.all.fluxes.cat''...'

      open ( unit=1, file='mc.all.fluxes.cat' )
      write (1,'(a)') '! '//tstamp//"  Created by "//codeID
      write (1,'(a)') '! '//author
      write (1,'(a)') '!'
      write (1,'(a)') '! TOTAL FLUXES & SPECTRAL ENERGY DISTRIBUTION'
      write (1,'(a)') '!'
      write (1,'(2(a,1pe15.8))') '! Distance (pc):', userdef_distpc, '  Mass-averaged Tdust (K):', tdust_massav
      write (1,'(a)') '!'
      write (1,'(a)') '!   N     WAVE(um)       FREQ(Hz)       FNU(Jy)       FNU(bs)        FNU(bg)        FNU(star)' &
                    //'      FNU(isrf)     NUFNU(HzJy)     NUFNU(bs)     NUFNU(bg)      NUFNU(star)    NUFNU(isrf)' &
                    //'      FNU(*)       NUFNU(Hz*)     EMIS(Td)mav    EMIS(Tdmav)   NUEMIS(Td)mav  NUEMIS(Tdmav)'
      write (1,'(a)') '!'

      fact = (z1 / userdef_distpc)**2
      solidangle = z1 / (userdef_distpc * parsec)**2
      rfact = Router / Rmodel
      backflux = z0
      isrfflux = z0
      do inu=1,nfreq
        read (2,*) lamb(inu), flux
        if (lbackexist) then
          read (3,*) backflux, isrfflux
        endif
        freq = 1.0d4 * cc / lamb(inu)
        cam_back = max ( fact * backflux, almostzero )
        cam_back_Jy = max ( fact * backflux / xJy, almostzero )
        cam_isrf_Jy = max ( fact * isrfflux * rfact**2 / xJy, almostzero )
        fluxbgs = fact * (flux - backflux)
        cam_flux_bgs = max ( fluxbgs, almostzero )
        cam_flux_bgs_Jy = max ( fluxbgs / xJy, almostzero )
        cam_flux = max ( fact * flux, almostzero )
        cam_flux_Jy = max ( fact * flux / xJy, almostzero )

        fstars = z0
        do istar=1,nstars
          areax = pi * (star_r(istar) / parsec)**2
          fstars = fstars + find_starlight_interpol(freq,istar) * fact * areax
        enddo
        cam_stars = max ( fstars, almostzero )
        cam_stars_Jy = max ( fstars / xJy, almostzero )

        temp = 100
        dustflux = z0
        do ic=1,nrcells
          idx = cellindex(ic)
          radcellcm = (radiuscm(ic) + radiuscm(ic+1)) / z2
          transparent_borderline = 175.0d0 * sqrt ( masstot / xMsun ) * (starlumtot / xLsun)**(z1/5.0d0)
          if (radcellcm .gt. transparent_borderline .and. radiuscm(ic) .lt. Rmodel * xAU) then
            do is=1,dust_nr_species
              kabs = find_dust_kappa_interpol ( freq, is, temp, 1, 0, 0 )
              dustflux = dustflux + kabs * bplanck ( dusttemp(is,idx), freq ) * dustdens(is,idx) * cellvolume(idx)
              tmavflux = tmavflux + kabs * bplanck ( tdust_massav, freq )
            enddo
          endif
        enddo
        dustflux = max ( dustflux * solidangle / xJy, almostzero )
        tmavflux = z0
        do is=1,dust_nr_species
          kabs = find_dust_kappa_interpol ( freq, is, temp, 1, 0, 0 )
          tmavflux = tmavflux + kabs * bplanck ( tdust_massav, freq )
        enddo
        tmavflux = max ( tmavflux * dust_mass * solidangle / xJy, almostzero )

        write (1,'(i5,30(1pe15.7))') inu, lamb(inu), freq, cam_flux_Jy, cam_flux_bgs_Jy, cam_back_Jy, cam_stars_Jy, cam_isrf_Jy, &
                                     freq * cam_flux_Jy, freq * cam_flux_bgs_Jy, freq * cam_back_Jy, freq * cam_stars_Jy, &
                                     freq * cam_isrf_Jy, cam_flux, freq * cam_flux, dustflux, tmavflux, freq * dustflux, & 
                                     freq * tmavflux
        fluxesbs(inu) = cam_flux_bgs_Jy
      enddo

      close ( 1 )
      close ( 2 )
      if (lbackexist) then
        close ( 3 )
      endif

      inquire ( file='camera_wavelength_micron.inp', exist=lexist )

      if (lexist) then
        open ( unit=1, file='camera_wavelength_micron.inp', status='old' )
        read (1,*) nfreq_cam

        lamb(1) = lamb(1) / (z1 + zeps14)
        lamb(nfreq) = lamb(nfreq) * (z1 + zeps14)
        do inu=1,nfreq_cam
          read (1,*) lamb_cam(inu)

          call hunt ( lamb, nfreq, lamb_cam(inu), k )

          if (k .lt. 1 .or. k .gt. nfreq) then
             write (*,'(a   )') ' RADMC3D: USERDEF_DOSTUFF (4):'
             write (*,'(a,i0)') '   ERROR: Wave index out of range: ', k
             stop 99
          endif
          eps = (lamb_cam(inu) - lamb(k)) / (lamb(k+1) - lamb(k))
          flux_cam(inu) = (z1 - eps) * fluxesbs(k) + eps * fluxesbs(k+1)
        enddo
        close ( 1 )
      endif

!!      inquire ( file='mc.fluxes.sed.cat', exist=lexist )
!!
!!      if (lexist) then
!!        write (*,'(a)') ' Appending ''mc.fluxes.sed.cat''...'
!!        open ( unit=4, file='mc.fluxes.sed.cat', recl=1000, access='append' )
!!        backspace ( 4 )
!!        read (4,'(i6)') numb
!!        if (numb .eq. 2) then
!!          backspace (4)
!!          backspace (4)
!!          read (4,'(i6)') numb
!!        endif
!!        if (numb .eq. 1) then
!!          write (cnu,'(i2)') nfreq_cam
!!          write (4,'(i6,1x,2f8.1,1pe11.3,0pf6.2,'//cnu//'(1pe11.3,4(1pe11.3)))') &
!!                 2,z0,z0,z100,z10, (z100, flux_cam(inu), flux_cam(inu) / z100, flux_cam(inu), flux_cam(inu)/z100,inu=1,nfreq_cam)
!!        endif
!!        close ( 4 )
!!      endif
    endif
  endif !<~ do_raytrace_spectrum
!__________________________________________________________________________________________________________________________________
!
  inquire ( file='camera_wavelength_micron.inp', exist=lexist )

  if (lexist) then

    open ( unit=1, file='wavelength_micron.inp', status='old' )
    open ( unit=2, file='camera_wavelength_micron.inp', status='old' )
    read (1,*) nfreq
    read (2,*) nfreq_cam
    do inu=1,nfreq
      read (1,*) lamb(inu)
    enddo
    allocate ( odepth(nfreq_cam), odepthi(nrcells,nfreq_cam), dtau(nrcells,nfreq_cam), opacity(nfreq_cam,dust_nr_species) )

    lamb(1) = lamb(1) / (z1 + zeps14)
    lamb(nfreq) = lamb(nfreq) * (z1 + zeps14)
    do inu=1,nfreq_cam
      read (2,*) lamb_cam(inu)

      call hunt ( lamb, nfreq, lamb_cam(inu), k )

      if (k .lt. 1 .or. k .gt. nfreq) then
        write (*,'(a   )') ' RADMC3D: USERDEF_DOSTUFF (5):'
        write (*,'(a,i0)') '   ERROR: Wave index out of range: ', k
        stop 99
      endif
      eps = (lamb_cam(inu) - lamb(k)) / (lamb(k+1) - lamb(k))
      do is=1,dust_nr_species
        opacity(inu,is) = (z1 - eps) * dust_kappa_abs (k,is) + eps * dust_kappa_abs (k+1,is) + &
                          (z1 - eps) * dust_kappa_scat(k,is) + eps * dust_kappa_scat(k+1,is)
      enddo
    enddo
    close ( 2 )
    close ( 1 )

    tstamp = timestmp ()
    open ( unit=4, file='mc.odepth.pro.tab', recl=1000 )
    write (*,'(a)') ' Writing ''mc.odepth.pro.tab''...'
    write (cnu,'(i2)') nfreq_cam
    write (4,'(a)') '! '//tstamp//"  Created by "//codeID
    write (4,'(a)') '! '//author
    write (4,'(a)') '!'
    write (4,'(a)') '! RADIAL PROFILES OF OPTICAL DEPTHS'
    write (4,'(a)') '!'
    write (4,'(2(a,1pe15.8))') '! Distance (pc):', userdef_distpc, '  Mass-averaged Tdust (K):', tdust_massav
    write (4,'(a)') '!'
    write (4,'(a,'//cnu//'(1x,a1,i2.2,a1,1pe10.3))') '! Wavelengths:', ('(',inu,')', lamb_cam(inu), inu=1,nfreq_cam)
    write (4,'(a)') '!'
    headline = '!   N    RADIUS(au)     RADIUS(pc)     RADIUS(as)'
    do inu=1,nfreq_cam
       write (cnu,'(i2.2)') inu
       headline = trim(headline)//'      DTAU('//cnu//')'//'     IN->OUT('//cnu//')'//'    OUT->IN('//cnu//')'
    enddo
    write (4,'(a)') trim(headline)
    write (4,'(a)') '!'
    dtau(:,:) = z0
    odepthi(:,:) = z0
    do ic=nrcells,1,-1
      idx = cellindex(ic)
      dzone = radiuscm(ic+1) - radiuscm(ic)
      if (radiuscm(ic) .lt. Rmodel * xAU) then 
        do inu=1,nfreq_cam
          do is=1,dust_nr_species
            dtau(ic,inu) = dtau(ic,inu) + z2 * opacity(inu,is) * dustdens(is,idx) * dzone
          enddo
          odepthi(max(ic-1,1),inu) = odepthi(ic,inu) + dtau(ic,inu)
        enddo
      endif
    enddo
    odepth(:) = z0
    do ic=1,nrcells
      idx = cellindex(ic)
      radcellcm = (radiuscm(ic) + radiuscm(ic+1)) / z2
      radcellAU = radcellcm / xAU
      radcellpc = radcellcm / parsec
      radcellas = radcellAU / userdef_distpc
      dzone = radiuscm(ic+1) - radiuscm(ic)
      if (radiuscm(ic) .lt. Rmodel * xAU) then 
        do inu=1,nfreq_cam
          do is=1,dust_nr_species
            odepth(inu) = odepth(inu) + z2 * opacity(inu,is) * dustdens(is,idx ) * dzone
          enddo
        enddo
      endif
      write (4,'(i5,31(1pe15.7))') idx, radcellAU, radcellpc, radcellas, (dtau(ic,inu), odepth(inu), odepthi(ic,inu) &
                                 , inu=1,nfreq_cam)
    enddo
    close ( 4 )
    deallocate ( odepth, odepthi, dtau, opacity )
  endif

  inquire ( file='image.out', exist=lexist )
  if (lexist) then
    open ( unit=2, file='image.out', status='old' )
    close ( 2, status='delete' )
  endif

  deallocate ( radiuscm )

  exec_time_sec = timer ( exec_time_sec )
  exec_time_min = exec_time_sec / 60.0d0
  exec_time_hrs = exec_time_min / 60.0d0
  write (cexec_time_sec,'(f9.2)') exec_time_sec
  write (cexec_time_min,'(f9.2)') exec_time_min
  write (cexec_time_hrs,'(f9.2)') exec_time_hrs
  write (*,'(/a)') ' Execution time: '//trim ( adjustl ( cexec_time_sec ) )//' sec = ' &
                                      //trim ( adjustl ( cexec_time_min ) )//' min = ' &
                                      //trim ( adjustl ( cexec_time_hrs ) )//' hrs'
  write (stdo,'(/a)') ' Execution time: '//trim ( adjustl ( cexec_time_sec ) )//' sec = ' &
                                         //trim ( adjustl ( cexec_time_min ) )//' min = ' &
                                         //trim ( adjustl ( cexec_time_hrs ) )//' hrs'
  tstamp = timestmp ()
  write (*,'(/1x,a)') trim(tstamp)
  write (stdo,'(/1x,a/)') trim(tstamp)
  write (*,'(/a)') ' Done.'

!__________________________________________________________________________________________________________________________________
!
end subroutine userdef_dostuff

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_action ()
!------------------------------------------------------------------------
! If you want to do design your own action (like 'mctherm' or 'image' but
! now designed by you entirely, and activated with 'radmc3d myaction')
! then here is your chance!
!------------------------------------------------------------------------
  implicit none
end subroutine userdef_action

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_compute_lorentz_delta ( ray_index )
!------------------------------------------------------------------------
! If you want to use a Voigt line profile instead of a Gaussian line
! profile, you must initialize the variable lines_ray_lorentz_delta
! here. (Added by Thomas Peters 2011)
!------------------------------------------------------------------------
  implicit none
  integer :: ray_index
end subroutine userdef_compute_lorentz_delta

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_compute_levelpop ( ispec, nlevels, index, x, y, z, numberdens, levelpop )
!------------------------------------------------------------------------
! If you have a good idea how to calculate the level populations 
! of a molecule on-the-fly, you can do it here. But to activate it,
! you must put the line mode to -2.
! IMPORTANT NOTE: If you use the method of selecting a subset of the
!                 levels of a molecule, then you must be very careful
!                 in this subroutine to do it right. You must then use
!                 the "active_***" arrays and variables in the line
!                 module to figure out which levels are "active" and 
!                 which are not. 
!------------------------------------------------------------------------
  implicit none
  integer :: ispec,nlevels,index
  double precision :: x,y,z,numberdens,levelpop(nlevels)
end subroutine userdef_compute_levelpop

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_general_compute_levelpop ( ray_index, levelpop )
!------------------------------------------------------------------------
! This routine lets you calculate the level populations entirely from
! scratch. Use lines_mode = -10 to use this routine.
!------------------------------------------------------------------------
  implicit none
  integer :: ray_index
  double precision :: levelpop(1:lines_nrlevels_subset_max,1:lines_nr_species)
end subroutine userdef_general_compute_levelpop

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_srcalp ( index, nrfreq, inu0, inu1, freq, src, alp )
!------------------------------------------------------------------------
! This subroutine allows you to specify exactly according to your own
! recipes/ideas the emissivity coefficient j_nu [erg/s/cm^3/Hz/ster]
! and extinction coefficient alpha_nu [1/cm] at the wavelengths given
! and at the location given. 
!
! ARGUMENTS:
!  index        The array index of the cell. This allows you to find
!               e.g. the gas density gasdens(index) or the gas
!               temperature gastemp(index) or any other quantity,
!               provided it is read in into the code. 
!  nrfreq       Nr of frequencies of the freq(:) array
!  freq(:)      Array of frequencies [Hz]
!  inu0,inu1    Starting/ending index: Calculate only src(inu0:inu1) and
!               alp(inu0:inu1) belonging to freq(inu0:inu1).
!
! RESULTS:
!  src(:)       Emissivity [erg/s/cm^3/Hz/ster]
!  alp(:)       Extinction [1/cm]
!
! Note: To activate this, you must set incl_userdef_srcalp = 1 in the
!       radmc3d.inp input file (in the code this is the logical 
!       rt_incl_userdef_srcalp from rtglobal_module.f90).
!
! Note: By the time RADMC-3D call this code, it has already computed
!       its own src(:) and alp(:). So just ADD your own values by e.g.
!       src(:) = src(:) + yourstuff and alp(:) = alp(:) + yourstuff.
!       In this way you won't delete the standard stuff. But if you
!       want to replace RADMC-3D's own stuff, you can also just write
!       src(:) = yourstuff and alp(:) = yourstuff. This is up to you.
!
! Note: Below you find some example, but commented-out, stuff.
!------------------------------------------------------------------------
  implicit none
  integer :: index,nrfreq,inu0,inu1
  double precision :: freq(1:nrfreq),src(1:nrfreq,1:4),alp(1:nrfreq)
  ! !
  ! ! If the index.lt.1 then we are not in a cell
  ! !
  ! if(index.lt.1) return
  ! ! 
  ! ! Replace src with dummy emission 
  ! !  
  ! src(:,1) = gasdens(index)
  ! alp(:) = 1d-40
  ! !
end subroutine userdef_srcalp

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_writemodel ()
!------------------------------------------------------------------------
! Here you can write model setup arrays (the stuff you have set up here
! in the userdef_module.f90) to standard RADMC-3D-readable files.
! This will only be done if radmc3d receives the 'writemodel' command
! on the command line. 
!------------------------------------------------------------------------
  implicit none
  call write_grid_file ()
  call write_dust_density ()
end subroutine userdef_writemodel

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine userdef_reset_flags ()
!------------------------------------------------------------------------
! Reset some action flags for next command?
!------------------------------------------------------------------------
  implicit none
end subroutine userdef_reset_flags

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  subroutine halfmaxsizes ( nx, ny, image, dx, dy, fwhm, xcentr, ycentr )
!__________________________________________________________________________________________________________________________________
!
! Direct measurements of the minimum and maximum sizes (at peak half-maximum) for an arbitrary intensity distribution.
!__________________________________________________________________________________________________________________________________
!
  implicit      none
  logical       lcond1, lcond2, lcond3, lcond4, lcond5, lcond6, lcond7, lcond8
  integer       nn, i, j, ii, jj, imax, jmax, nx, ny, nedge, im1, jm1, ip1, jp1
  real*8        dx, dy, dblii2, dbljj2, fwhm, imhalf, funmax, almostzero, radius2, drad2, radhalf, hwhm, sigm2 &
              , dlogim, sqpixel, xhalf, yhalf, xcentr, ycentr, imageijabs, image(nx,ny) !!, mask(nx,ny)
!__________________________________________________________________________________________________________________________________
!                  
  funmax = -1.0d99
  imax = 0
  jmax = 0
  do j=1,ny
    do i=1,nx
!!      if (mask(i,j) .gt. almostzero) then
        if (abs ( image(i,j) ) .ge. funmax) then
          imax = i
          jmax = j
          funmax = abs ( image(i,j) )
        endif
!!      endif
    enddo
  enddo
  
  hwhm = 0.0d0
  xcentr = 0.0d0
  ycentr = 0.0d0
  imhalf = abs ( funmax ) / 2.0d0
  nedge = 0
  
  if (imax .gt. 0 .and. jmax .gt. 0) then
    
    do j=1,ny
      jj = j - jmax
      do i=1,nx
        ii = i - imax
        imageijabs = abs ( image(i,j) )
        
        if (imageijabs .ge. imhalf) then     !! .and. mask(i,j) .gt. 0.9d0
          im1 = max ( i - 1, 1 )
          jm1 = max ( j - 1, 1 )
          ip1 = min ( i + 1, nx )
          jp1 = min ( j + 1, ny )
          
          lcond1 = abs ( image(im1,j  ) ) .lt. imhalf .and. abs ( image(im1,j  ) ) .ge. 0.0d0    !! .gt. almostzero
          lcond2 = abs ( image(ip1,j  ) ) .lt. imhalf .and. abs ( image(ip1,j  ) ) .ge. 0.0d0    !! .gt. almostzero
          lcond3 = abs ( image(i  ,jm1) ) .lt. imhalf .and. abs ( image(i  ,jm1) ) .ge. 0.0d0    !! .gt. almostzero
          lcond4 = abs ( image(i  ,jp1) ) .lt. imhalf .and. abs ( image(i  ,jp1) ) .ge. 0.0d0    !! .gt. almostzero
          lcond5 = abs ( image(im1,jm1) ) .lt. imhalf .and. abs ( image(im1,jm1) ) .ge. 0.0d0    !! .gt. almostzero
          lcond6 = abs ( image(im1,jp1) ) .lt. imhalf .and. abs ( image(im1,jp1) ) .ge. 0.0d0    !! .gt. almostzero
          lcond7 = abs ( image(ip1,jm1) ) .lt. imhalf .and. abs ( image(ip1,jm1) ) .ge. 0.0d0    !! .gt. almostzero
          lcond8 = abs ( image(ip1,jp1) ) .lt. imhalf .and. abs ( image(ip1,jp1) ) .ge. 0.0d0    !! .gt. almostzero
          
          if (lcond1 .or. lcond2 .or. lcond3 .or. lcond4 .or. lcond5 .or. lcond6 .or. lcond7 .or. lcond8) then
            nedge = nedge + 1
            sqpixel = dx * dy
            dblii2 = dble ( ii )**2
            dbljj2 = dble ( jj )**2
            radius2 = (dblii2 + dbljj2) * sqpixel
            dlogim = log ( imageijabs ) - log ( abs ( imhalf ) )
            radhalf = 0.0d0
            xhalf = 0.0d0
            yhalf = 0.0d0
            nn = 0
                     
! Gaussian interpolation to obtain accurate values of full width at half-maximum for Gaussian-like shapes.
! From two equations: ln(G1) - ln(G2) = -x1^2 + x2^2 / (2 sigma^2) and ln(G1) - ln(Gh) = -x1^2 + xh^2 / (2 sigma^2):
! (2 sigma^2) = (x2^2 - x1^2) / (ln(G1) - ln(G2)) and xh = sqrt ( x1^2 + (2 sigma^2) * (ln(G1) - ln(Gh)) )
                 
            if (lcond1) then
              nn = nn + 1
              drad2 = (dble ( im1 - imax )**2 + dbljj2) * sqpixel - radius2
              sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(im1,j) ) ))
              radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
              xhalf = xhalf + dble ( im1 )
              yhalf = yhalf + dble ( j )
            endif
            if (lcond2) then
              nn = nn + 1
              drad2 = (dble ( ip1 - imax )**2 + dbljj2) * sqpixel - radius2
              sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(ip1,j) ) ))
              radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
              xhalf = xhalf + dble ( ip1 )
              yhalf = yhalf + dble ( j )
            endif
            if (lcond3) then
              nn = nn + 1
              drad2 = (dblii2 + dble ( jm1 - jmax )**2) * sqpixel - radius2
              sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(i,jm1) ) ))
              radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
              xhalf = xhalf + dble ( i )
              yhalf = yhalf + dble ( jm1 )
            endif
            if (lcond4) then
              nn = nn + 1
              drad2 = (dblii2 + dble ( jp1 - jmax )**2) * sqpixel - radius2
              sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(i,jp1) ) ))
              radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
              xhalf = xhalf + dble ( i )
              yhalf = yhalf + dble ( jp1 )
            endif
            if (lcond5) then
              nn = nn + 1
              drad2 = (dble ( im1 - imax )**2 + dble ( jm1 - jmax )**2) * sqpixel - radius2
              sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(im1,jm1) ) ))
              radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
              xhalf = xhalf + dble ( im1 )
              yhalf = yhalf + dble ( jm1 )
            endif
            if (lcond6) then
              nn = nn + 1
              drad2 = (dble ( im1 - imax )**2 + dble ( jp1 - jmax )**2) * sqpixel - radius2
              sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(im1,jp1) ) ))
              radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
              xhalf = xhalf + dble ( im1 )
              yhalf = yhalf + dble ( jp1 )
            endif
            if (lcond7) then
              nn = nn + 1
              drad2 = (dble ( ip1 - imax )**2 + dble ( jm1 - jmax )**2) * sqpixel - radius2
              sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(ip1,jm1) ) ))
              radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
              xhalf = xhalf + dble ( ip1 )
              yhalf = yhalf + dble ( jm1 )
            endif
            if (lcond8) then
              nn = nn + 1
              drad2 = (dble ( ip1 - imax )**2 + dble ( jp1 - jmax )**2) * sqpixel - radius2
              sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(ip1,jp1) ) ))
              radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
              xhalf = xhalf + dble ( ip1 )
              yhalf = yhalf + dble ( jp1 )
            endif
            
            if (nn .gt. 0) then
              radhalf = radhalf / dble ( nn )
              xhalf = xhalf / dble ( nn )
              yhalf = yhalf / dble ( nn )
            endif
            hwhm = hwhm + radhalf
            xcentr = xcentr + xhalf
            ycentr = ycentr + yhalf
          endif
        endif
      enddo
    enddo
         
    if (nedge .gt. 0) then
      hwhm = hwhm / dble ( nedge )
      xcentr = xcentr / dble ( nedge )
      ycentr = ycentr / dble ( nedge )
      fwhm = max ( 2.0d0 * hwhm, sqrt ( sqpixel ) )
    endif
  endif
       
  if (xcentr .lt. almostzero) xcentr = dble ( nx ) / 2.0d0       
  if (ycentr .lt. almostzero) ycentr = dble ( ny ) / 2.0d0       

  return
  end subroutine halfmaxsizes

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  subroutine momentsizes ( nxmin, nxmax, nymin, nymax, nx, ny, image, dx, dy, momxco, momyco, afwhm, bfwhm, atheta, equivrad &
                         , elongation )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
  implicit      none
  integer       nx, ny, l, m, nxmin, nxmax, nymin, nymax
  real*8        dx, dy, pi, sqpixel, sourcearea, dblm, dbll, convtodegrees, sqreightlogtwo, intsxx, intsyy, intsxy &
              , xcorsc, ycorsc, momp2, momm2, sqrmo, asig, bsig, intsd, ints0m, ints0, almostzero, intx, inty, int0, int0m &
              , intd, nonzero
  real*8        image(nx,ny), momxco, momyco, afwhm, bfwhm, atheta, equivrad, elongation, momxx, momyy, momxy
  parameter   ( pi = 3.14159265358979d0, convtodegrees = 180.0d0 / pi, sqreightlogtwo = sqrt ( 8.0d0 * log ( 2.0d0 ) ) )
!__________________________________________________________________________________________________________________________________
!             
  sqpixel = dx * dy
  intx = 0.0d0                        
  inty = 0.0d0                        
  int0 = 0.0d0
  nonzero = 0.0d0
  
  do m=nymin,nymax
    dblm = dble ( m )                 
    do l=nxmin,nxmax            
      dbll = dble ( l )               
      intd = image(l,m) * sqpixel
      int0 = int0 + intd
      intx = intx + dbll * intd
      inty = inty + dblm * intd
      if (image(l,m) .gt. almostzero) then
        nonzero = nonzero + 1.0d0
      endif
    enddo
  enddo
  int0m = int0 + almostzero
  momxco = intx / int0m
  momyco = inty / int0m
  
  ints0 = 0.0d0
  intsxx = 0.0d0
  intsyy = 0.0d0
  intsxy = 0.0d0           
  
  do m=nymin,nymax
    dblm = dble ( m )
    do l=nxmin,nxmax
      dbll = dble ( l )       
      intsd = image(l,m) * sqpixel
      ints0 = ints0 + intsd
      xcorsc = dbll - momxco
      ycorsc = dblm - momyco
      intsxx = intsxx + xcorsc * xcorsc * intsd
      intsyy = intsyy + ycorsc * ycorsc * intsd
      intsxy = intsxy + xcorsc * ycorsc * intsd
    enddo
  enddo
  ints0m = ints0 + almostzero
  momxx = intsxx / ints0m
  momyy = intsyy / ints0m
  momxy = intsxy / ints0m

  momp2 = (momxx + momyy) / 2.0d0
  momm2 = (momxx - momyy) / 2.0d0
  
  if (abs ( momxy ) .lt. almostzero) then
    if (momxx .gt. momyy) then
      momp2 = momxx
      momm2 = momyy
      atheta = 90.0d0
    else
      momp2 = momyy
      momm2 = momxx
      atheta = 0.0d0
    endif
    sqrmo = sqrt ( momm2**2 + momxy**2 )
    asig  = sqrt ( (momp2 + sqrmo) * sqpixel ) 
    if (momp2 .gt. sqrmo) then 
      bsig = sqrt ( (momp2 - sqrmo) * sqpixel )
    else
      asig = 0.0d0
      bsig = 0.0d0
    endif
  else       
    sqrmo = sqrt ( momm2**2 + momxy**2 )
    if (momp2 + sqrmo .ge. 0.0d0) then 
      asig = sqrt ( (momp2 + sqrmo) * sqpixel )
    endif
    if (momp2 - sqrmo .ge. 0.0d0) then 
      bsig = sqrt ( (momp2 - sqrmo) * sqpixel )
    endif
    if (momp2 + sqrmo .lt. 0.0d0 .or. momp2 - sqrmo .lt. 0.0d0) then 
      asig = 0.0d0
      bsig = 0.0d0
    endif
    atheta = 180.0d0 - atan2 ( sqrmo + momm2, momxy ) * convtodegrees
  endif                                        

  afwhm = max ( sqreightlogtwo * asig, sqrt ( dx * dy ) )
  bfwhm = max ( sqreightlogtwo * bsig, sqrt ( dx * dy ) )

  elongation = afwhm / bfwhm

  sourcearea = nonzero * sqpixel
  equivrad = sqrt ( sourcearea / pi )

  return
  end subroutine momentsizes

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  subroutine wfits ( cfitsversion, inmx, nx, ny, bunit, ctype1, ctype2, rp1, rp2, xr, yr, fun, dx, dy, object, ra, dec, fname &
                   , header, footer, cdate, ctime, creator, beam, blank, rot1, rot2, cd11, cd12, cd21, cd22, epoch, equinox &
                   , bzero, bscale, wave, datamn, datamx, history )
!__________________________________________________________________________________________________________________________________
!
! Write an image in FITS format.
!__________________________________________________________________________________________________________________________________
!
  implicit      none
  logical       simple, extend
  character*(*) fname, ctype1, ctype2, object, creator, header, footer(2), ctime, cdate, bunit, cfitsversion, history
  character*19  cdatetime
  integer       inmx, i, j, nx, ny, status, blocksze, group, naxis, pcount, gcount, bitpix, unit, naxes(2), iat, iflen, blank
  real*8        xr, yr, dx, dy, rp1, rp2, ra, dec, beam, bzero, bscale, as2deg, dfun, rot1, rot2, epoch, equinox, datamn     &
              , cosrota, sinrota, rotangle, datamx, dxdeg, dydeg, wave, cd11, cd12, cd21, cd22, almostzero, pi
  real*8        fun(inmx,inmx)
  parameter   ( pi = 3.14159265358979d0, as2deg = 1.0d0 / 3600.0d0, almostzero = 1.0d-30 )
!__________________________________________________________________________________________________________________________________
!
! Create a FITS primary array containing a 2-D image. Initialize parameters about the FITS image.

  status   = 0
  blocksze = 1
  simple   = .true.
  extend   = .false.
  group    = 1
  naxis    = 2
  pcount   = 0
  gcount   = 1
  naxes(1) = nx
  naxes(2) = ny

! Compute offset and scaling factor so as to represent real numbers most efficiently and accurately with given 'bitpix'.

  datamn = 1.0d30
  datamx =-1.0d30
  do j=1,ny
    do i=1,nx
      dfun = fun(i,j)
      if (dfun .lt. datamn ) datamn = dfun
      if (dfun .gt. datamx ) datamx = dfun
    enddo
  enddo

! When packing numbers using the formulas below, I got values of the function different in the 5-th digit from DATAMIN and
!  DATAMAX placed in the header when writing the data. From now on I will use bzero = 0 bscale = 1 and bitpix = -32, as these
! seem to cure the problem.

!!!   bitpix = 16
!!!   bzero  = 0.5d0 * (funmin + funmax)
!!!   bscale = dabs ( funmax - bzero ) / 32767.0d0

  bitpix = -32     !<-- single precision floating-point values
!!    bzero  = 0.0d0
!!    bscale = 1.0d0

  cdatetime = cdate(7:10)//'-'//cdate(4:5)//'-'//cdate(1:2)//'T'//ctime

! Delete the file if it already exists.

  if (cfitsversion .eq. '5.030') then
    call fdelete ( fname, status )
  endif

! Get an unused Logical Unit Number to use to open the FITS file.

  call ftgiou ( unit, status )

! Create the new empty FITS file.

  if (cfitsversion .eq. '5.030') then
    call ftinit ( unit, fname, blocksze, status )
  else
    call ftinit ( unit, '!'//fname, blocksze, status )
  endif

! Write the required header keywords.

  call ftphpr ( unit, simple, bitpix, naxis, naxes, pcount, gcount, extend, status )

! Write an optional comment.

  call ftpkys ( unit, 'CREATOR', creator,'Alexander Men’shchikov, DAp IRFU CEA Saclay', status )
  call ftpcom ( unit,' ',status)
  call ftpkys ( unit, 'HEADER' , header, 'title', status )

  iat = 0
  iat = index ( footer(1), ' |F' )
  if (iat .eq. 0) iat = index ( footer(1), ' mm ' )
  iflen = len ( trim(footer(1)) )

  if (iat .gt. 0) then
    call ftpkys ( unit, 'FOOTER1', footer(1)(1:iat-1), 'part 1', status )
    call ftpkys ( unit, 'FOOTER2', footer(1)(iat:iflen), 'part 2', status )
  else
    call ftpkys ( unit, 'FOOTER1', footer(1), 'part 1', status )
    call ftpkys ( unit, 'FOOTER2', ' ', 'part 2', status )
  endif
  call ftpcom ( unit,' ',status)
  call ftpkys ( unit, 'DATE'  , cdatetime, 'creation date and time', status )

! Write all parameters.

  dxdeg = dx * as2deg
  dydeg = dy * as2deg
       
! Correct if the input images had cd11 or cd22 set to zero.

  if (abs ( cd11 ) .lt. almostzero .or. abs ( cd11 ) .gt. dxdeg) then
    cd11 = -dxdeg
    cd12 = 0.0d0
  endif
  if (abs ( cd22 ) .lt. almostzero .or. abs ( cd22 ) .gt. dxdeg) then
    cd21 = 0.0d0
    cd22 = dydeg
  endif
  
  cosrota = - cd11 / dxdeg
  sinrota = - cd12 / dxdeg
  rotangle = acos ( cosrota ) * 180.0d0 / pi
  
  if (cosrota .ge. 0.0d0 .and. sinrota .ge. 0.0d0) then
    rot1 = rotangle
  elseif (cosrota .lt. 0.0d0 .and. sinrota .ge. 0.0d0) then
    rot1 = 180.0d0 - rotangle
  elseif (cosrota .lt. 0.0d0 .and. sinrota .lt. 0.0d0) then
    rot1 = 180.0d0 + rotangle
  elseif (cosrota .ge. 0.0d0 .and. sinrota .lt. 0.0d0) then
    rot1 = - rotangle
  endif
  rot2 = rot1

  call ftpkyd ( unit, 'BZERO'  , bzero , 13 , 'zero point in scaling equation', status )
  call ftpkyd ( unit, 'BSCALE' , bscale, 13 , 'linear factor in scaling equation' , status )
  call ftpkyd ( unit, 'DATAMAX', datamx, 13 , 'maximum data value', status )
  call ftpkyd ( unit, 'DATAMIN', datamn, 13 , 'minimum data value', status )
  if  (blank .gt. 0) &
  call ftpkyj ( unit, 'BLANK'  , blank      , 'value used for undefined array elements', status )
  call ftpkys ( unit, 'BUNIT'  , bunit      , 'physical units of the array values', status )
  call ftpkys ( unit, 'CTYPE1' , ctype1     , 'name of the coordinate axis', status )
  call ftpkys ( unit, 'CTYPE2' , ctype2     , 'name of the coordinate axis', status )
  call ftpkyd ( unit, 'CRPIX1' , rp1   , 13 , 'coordinate system reference pixel', status )
  call ftpkyd ( unit, 'CRPIX2' , rp2   , 13 , 'coordinate system reference pixel', status )
  call ftpkyd ( unit, 'CROTA1' , rot1  , 13 , 'coordinate system rotation angle', status )
  call ftpkyd ( unit, 'CROTA2' , rot2  , 13 , 'coordinate system rotation angle', status )
  call ftpkyd ( unit, 'CD1_1'  , cd11  , 13 , 'linear projection matrix', status )
  call ftpkyd ( unit, 'CD1_2'  , cd12  , 13 , 'linear projection matrix', status )
  call ftpkyd ( unit, 'CD2_1'  , cd21  , 13 , 'linear projection matrix', status )
  call ftpkyd ( unit, 'CD2_2'  , cd22  , 13 , 'linear projection matrix', status )
  call ftpkyd ( unit, 'CRVAL1' , xr    , 13 , 'coordinate value at reference pixel', status )
  call ftpkyd ( unit, 'CRVAL2' , yr    , 13 , 'coordinate value at reference pixel', status )
  call ftpkyd ( unit, 'CDELT1' , -dxdeg, 13 , 'coordinate increment along axis', status )
  call ftpkyd ( unit, 'CDELT2' , dydeg , 13 , 'coordinate increment along axis', status )
  call ftpkys ( unit, 'OBJECT' , object     , 'object identifier', status )
  call ftpkyd ( unit, 'RA'     , ra     , 13, 'right ascension', status )
  call ftpkyd ( unit, 'DEC'    , dec    , 13, 'declination', status )
  call ftpkyd ( unit, 'WAVE'   , wave   , 13, 'wavelength (microns)', status )
  call ftpkyd ( unit, 'EPOCH'  , epoch  , 13, ' ', status )
  call ftpkyd ( unit, 'EQUINOX', equinox, 13, ' ', status )
  call ftphis ( unit, history, status )

  if (beam .ne. 0.0d0) call ftpkyd ( unit, 'BEAM', beam, 13, 'image was convolved with the beam (arcsec)', status )

! Write the array to the FITS file.

  call ftp2dd ( unit, group, inmx, nx, ny, fun, status )

! Close the file.

  call ftclos ( unit, status )

! Free the unit number.

  call ftfiou ( unit, status )

! Check for any error, and if so print out error messages

  if (status .gt. 0) call printerr ( status )

  return
  end subroutine wfits

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 function timer ( oldtime ) result ( time_sec )
!__________________________________________________________________________________________________________________________________
!
!                                     MCM-3D, Alexander Men'shchikov, 20/08/2007
!
!  PURPOSE: Timer for measuring execution time and CPU time.
!
!  INPUT VARIABLES:
!    oldtime     previously taken measurement that will be subtracted
!                from the new system timer value, giving time interval
!                (otime can be zero).
!
!  OUTPUT VARIABLES:
!    timer       time interval or initial timer value (for oldtime=0.0).
!__________________________________________________________________________________________________________________________________
!
 double precision, intent(in) :: oldtime
 double precision :: timers, time_sec, hrs, mins, secs, msec
 integer, dimension(8) :: time
!__________________________________________________________________________________________________________________________________
!
 call date_and_time ( values = time )

 hrs  = dble ( time(5) )
 mins = dble ( time(6) )
 secs = dble ( time(7) )
 msec = dble ( time(8) )

 timers = hrs * z3600 + mins * z60 + secs + msec / z1000
 if (timers < oldtime) timers = timers + 24.0d0 * z3600
 time_sec = timers - oldtime

 return
 end function timer

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 function timestmp () result ( tstamp )
!__________________________________________________________________________________________________________________________________
!
!                                     MCM-3D, Alexander Men'shchikov, 20/08/2007
!
!  PURPOSE: Time stamping
!
!  OUTPUT VARIABLES:
!    tstamp        time stamp string, e.g. "Date: 03/07/2007 Time: 13:00:47"
!__________________________________________________________________________________________________________________________________
!
 character( 2) :: ct2, ct3, ct5, ct6, ct7
 character(32) :: tstamp
 integer, dimension(8) :: time
!__________________________________________________________________________________________________________________________________
!
 call date_and_time ( values = time )

 write (ct2,'(i2)') time(2)
 write (ct3,'(i2)') time(3)
 write (ct5,'(i2)') time(5)
 write (ct6,'(i2)') time(6)
 write (ct7,'(i2)') time(7)
 if (ct2(1:1) == ' ') ct2(1:1) = '0'
 if (ct3(1:1) == ' ') ct3(1:1) = '0'
 if (ct5(1:1) == ' ') ct5(1:1) = '0'
 if (ct6(1:1) == ' ') ct6(1:1) = '0'
 if (ct7(1:1) == ' ') ct7(1:1) = '0'

 write (tstamp,'(a,i4,a)') "Date: "//ct3//"/"//ct2//"/", time(1), "  Time: "//ct5//":"//ct6//":"//ct7

 return
 end function timestmp

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 subroutine timedhms ( timesec, cdays, chours, cmins, csecs, chsec )
!__________________________________________________________________________________________________________________________________
!
!                                     MCM-3D, Alexander Men'shchikov, 20/08/2007
!
!  PURPOSE: Conversion from seconds into days, hours, minutes, seconds, and fraction of a second
!
!  INPUT VARIABLES:
!    timesec        time in seconds
!
!  OUTPUT VARIABLES:
!    cdays, chours, cmins, csecs, chsec - integer days, hours, minutes, seconds, and fraction of a second
!__________________________________________________________________________________________________________________________________
!
 integer, intent(out) :: cdays, chours, cmins, csecs, chsec
 double precision, intent(in ) :: timesec
!__________________________________________________________________________________________________________________________________
!
 cdays  = int ( timesec / z86400 )
 chours = int ( timesec /  z3600 ) -    z24 * cdays
 cmins  = int ( timesec /    z60 ) -  z1440 * cdays -   z60 * chours
 csecs  = int ( timesec          ) - z86400 * cdays - z3600 * chours - z60 * cmins
 chsec  = int ( (timesec - int ( timesec )) * z100 )

 return
 end subroutine timedhms

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module userdef_module

