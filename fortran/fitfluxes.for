
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       program fitfluxes
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       logical       lfname1, lfname2, lunix, lgood, lconv, ldebug, lfnameo, limages, lplotting, ltiming, lgoodfound, lconverged
     &             , lconvergedo, ldone, lnochange, lconvgood, lnotconvo, lsing

       integer       nparams, npts, nfun, ntries

       parameter   ( nparams = 5, npts = 200, nfun = 5, ntries = 99 )
                                            
       character*1   csnband(100)
       character*2   cwav, cno
       character*3   cjx, ckx, clx, cmx, cl, cpx, cnx, cpjx, cpkx, cplx, cpmx, cpnx, cnperc, cnan, dot, cokbad
       character*4   clam1, clam2, cx, cy
       character*6   cfitsversion, cverbose
       character*7   clibname, ccx, ccy, cwave, punits(nparams)
       character*8   ctime, cfpptime, cmodel, cnumber, yesno(2)
       character*9   cfts, ccpu, cwal, cccpu, ccwal
       character*10  cdate, cnsrc, cnsend, cnst
       character*21  compda
       character*80  object, creator, ctype1, ctype2, bunit, history
       character*500 arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, argA, argB, argC, argD, argE, outname, prefix, excat
     &             , callwaves, callsnerr, callparams, postscript, paramsfile, filename, surfdensit, temperatur, chisquarnu
     &             , rtemperror, rsurferror, paraspace, commandline
       character*50000 cline
       character*50000, allocatable :: headline(:)

       integer       firstb, lastc, i, j, k, l, m, nsrc, ilp, ion, ndate, iotty, iolog, lhl, lb, nbands, icp, ia1, ia2, ia3, ia4
     &             , ia5, ia6, ia7, ia8, ia9, iaA, iaB, iaC, iaD, iaE, nhead, nextr, ncolxcoo, ncolycoo, nhcols, iexc, ne, psunit
     &             , nn, n1, n2, n3, irc, icm, nptstofit, ncolgood, ncolsigg, nmodparams, npar(nparams), nplt, ils, ns, nfreedom
     &             , njx, nkx, nlx, nmx, nl, ics, ip, nl1, nl2, nbande, nfitparams, lam1, lam2, number, npjx, npkx, nplx, npmx, ny1
     &             , numfits, nfts, ncolnumb, nx, ny, nb, blank, ilc, ilt, ilq, ncx, ncy, na, icl, npercent, icw, nx1, nband(100)
     &             , ivarornot(nparams), npart, nparts, nsbeg, nsend, nx0o, ny0o, ile, wchours, wcmins, ipr, wcsecs, cpuhours, nm
     &             , cpumins, cpusecs, ndc, ndw, nbn, nbx, itref, npx, maxitref, iter, itrefmin, itrefmax, itermin, itermax, nparo
     &             , nparmx(nparams), iterdone, npnx, n, nnx, nst, nt, ilm, ilr, idum, newran, nskp, nskip, iverbose
     
       real          fitsvers, xn(npts,nfun), fun(npts,nfun), yerrp(npts,nfun), yerrm(npts,nfun)

       real*8        dx, dy, crpix1, crpix2, crval1, crval2, ra, dec, funmin, funmax, beamx, crota1, crota2, cd11, cd12, chisqfree 
     &             , cd21, cd22, equinox, bzero, bscale, wave, datamin, datamax, pi, pio2, pi2, derivedmass, pspace1o, beam1
     &             , speedolight, distancepc, dust2gas, opacity0, chisqx, opaslope, wavepeak, frequency0, wavelength0, goodmin
     &             , sq_arcsecs_per_sterad, sigmin, amu, muH2, almostzero, Msun, pc, chisqmin, opacit0relerr, dustgasrelerr, dlogpo
     &             , distancrelerr, fluxpeak, modparam1(nparams), modparam2(nparams), dlogp(nparams), waveband(100), adderro(100)
     &             , snrmin(100), wctot, modparams1(nparams), modparams2(nparams), quantities(100), parmin(nparams), parmax(nparams)
     &             , modp1, modp2, delta, cpu, wal, cpu_sum, wal_sum, cputot, timer, offset, modelflux, cdelt1, cdelt2, chisqfreeo
     &             , parmn(nparams), parmx(nparams), numpart, numtotal, isumpct, itermean, itrefmean, nitmean, nrefmean, parmnlog
     &             , parmxlog, chisqfmax, pspace2o, dlogpmax(nparams), solidapixel, okbad, sigma, q1, ran1, gasdev
     &             , simtruemass(ntries)

       logical, allocatable :: ltakeit(:)
                                         
       integer, allocatable :: ncolsig(:), ncolfxp(:), ncolfxt(:), ncolsnr(:), ncolfpe(:), ncolfte(:), ipfitornot(:), nxmin(:,:)
     &                       , uplimit(:)
       real*8, allocatable :: sigmono(:), snrpeak(:), snrtotal(:), frequency(:), fxpbest(:), fxperro(:), fxtbest(:), fxterro(:)
     &                      , sfrequency(:), sfxtbest(:), sfxterro(:), opacity(:), sedslope(:), wavelength(:), fitparameters(:)
     &                      , dyda(:), covar(:,:), alpha(:,:), chisq(:,:,:,:,:), paramspace(:,:), fun2d(:,:), image(:,:,:), beam(:)
     &                      , paramspacemin(:), paramspacemax(:), x0o(:), y0o(:), sigglobal(:), goodness(:), space(:,:)
     &                      , surfden(:,:), tempers(:,:), chisqnu(:,:), temperr(:,:), surferr(:,:), fxterror(:)  !!, fxperror(:)
     
       parameter   ( sq_arcsecs_per_sterad = 3282.80635d0 * 3600.0d0**2, pi = 3.14159265358979d0, pio2 = pi / 2.0d0, muH2 = 2.8d0
     &             , pi2 = pi * 2.0d0, speedolight = 2.99792458d10, almostzero = 1.0d-20, amu = 1.6605402D-24, Msun = 1.9891d33
     &             , pc = 3.0856775814913673d18, ldebug = .false., ltiming = .true., maxitref = 5, chisqfmax = 0.001d0, dot = '•'
     &             , psunit = 60 )

       common / copacity / frequency0, wavelength0, opacity0, dust2gas, opacit0relerr, dustgasrelerr, iverbose

       external      thinbody, firstb, lastc, lastbs, timer, modbody, when, fitfun, osystem, ftvers, skipcomm, linescolumns
     &             , getfitshead, rfits, wfits, showprogress, fitsed, plotsed, plotnd, ps2pdf, ran1, gasdev
!__________________________________________________________________________________________________________________________________
!                    
       iotty = 6
       yesno(1) = 'constant'
       yesno(2) = 'variable'
       cnan = 'NAN'

! Determine operating system

       call osystem ( lunix )
       call ftvers ( fitsvers )
       lb = 1
       write (cfitsversion,'(f6.3)') fitsvers
       if (cfitsversion .eq. ' 5.030') then
         lb = 2
       endif
       clibname = 'CFITSIO'
       cfpptime = __TIME__
       cfpptime = cfpptime(1:5)
       compda = __DATE__//' '//cfpptime

! Get the current date and time.

       call when ( lunix, ctime, cdate, ndate, 4 )

! Get command line parameters (file names).

       call getarg (  1, arg1 )
       call getarg (  2, arg2 )
       call getarg (  3, arg3 )
       call getarg (  4, arg4 )
       call getarg (  5, arg5 )
       call getarg (  6, arg6 )
       call getarg (  7, arg7 )
       call getarg (  8, arg8 )
       call getarg (  9, arg9 )
       call getarg ( 10, argA )
       call getarg ( 11, argB )
       call getarg ( 12, argC )
       call getarg ( 13, argD )
       call getarg ( 14, argE )

       ia1 = lastc ( arg1 )
       ia2 = lastc ( arg2 )
       ia3 = lastc ( arg3 )
       ia4 = lastc ( arg4 )
       ia5 = lastc ( arg5 )
       ia6 = lastc ( arg6 )
       ia7 = lastc ( arg7 )
       ia8 = lastc ( arg8 )
       ia9 = lastc ( arg9 )
       iaA = lastc ( argA )
       iaB = lastc ( argB )
       iaC = lastc ( argC )
       iaD = lastc ( argD )
       iaE = lastc ( argE )
       
       cverbose = '-verb2'
       if (iaE .gt. 0 .and. (argE(1:iaE-1) .eq. '-verb' .or. argE(1:iaE-1) .eq. '+verb')) cverbose = argE(1:iaE)
       if (iaD .gt. 0 .and. (argD(1:iaD-1) .eq. '-verb' .or. argD(1:iaD-1) .eq. '+verb')) cverbose = argD(1:iaD)
       if (iaC .gt. 0 .and. (argC(1:iaC-1) .eq. '-verb' .or. argC(1:iaC-1) .eq. '+verb')) cverbose = argC(1:iaC)
       if (iaB .gt. 0 .and. (argB(1:iaB-1) .eq. '-verb' .or. argB(1:iaB-1) .eq. '+verb')) cverbose = argB(1:iaB)
       if (iaA .gt. 0 .and. (argA(1:iaA-1) .eq. '-verb' .or. argA(1:iaA-1) .eq. '+verb')) cverbose = argA(1:iaA)
       
       if (cverbose(1:1) .eq. '+') then
         lplotting = .true.
         cverbose(1:1) = '-'
       else
         lplotting = .false.
       endif
       read (cverbose,'(5x,i1)') iverbose
       nhcols = 0
       chisqfree = 0.0d0
       excat = argB(1:iaB)
       iexc = iaB
       
       if (iotty .gt. 0) then         !! .and. cverbose .eq. '-verb2'
!!         write (iotty,'( )')
         if (arg1(1:ia1) .eq. ':') write (iotty,'(a$)') ' '
         if (iverbose .gt. 0 .or. excat(iexc-3:iexc) .eq. '.cfg') then
           write (iotty,'(/a$)') ' FITFLUXES '//dot//' Fit Source Fluxes or Image Intensities '//dot//' '//compda
         endif
       endif
       if (arg1(1:ia1) .eq. ':') stop
       if (iotty .gt. 0) then
         if (iverbose .gt. 0 .or. excat(iexc-3:iexc) .eq. '.cfg') write (iotty,'()')
         if (iverbose .gt. 0 .or. excat(iexc-3:iexc) .eq. '.cfg') then
           write (iotty,'( a)') ' Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
           write (iotty,'( a)') ' Using MRQMIN from Numerical Recipes for F77 by William H Press et al.'
         endif
         if (iaA .eq. 0) then
           write (iotty,'(a)')
           write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USAGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
           write (iotty,'(a)')
           write (iotty,'(a)') ' fitfluxes {<number>|<xcoo>} {<sigmin>|<npart>|<ycoo>} {<goodmin>|<nparts>}'
           write (iotty,'(a)') '           [+|-]<n1>:<wave1>,[+|-]<n2>:<wave2>,...,[+|-]<nN>:<waveN>'
           write (iotty,'(a)') '           <snr1>:<adderr1>,<snr2>:<adderr2>,...,<snrN>:<adderrN>'
           write (iotty,'(a)') '           {thinbody|modbody}'
           write (iotty,'(a)') '           {<pval1>|0}:{<dlp1>|1|0},{<pval2>|0}:{<dlp2>|1|0},'
           write (iotty,'(a)') '           {<pval3>|0}:{<dlp3>|1|0},{<pval4>|0}:{<dlp4>|1|0}'
           write (iotty,'(a)') '           {<distpc>|1}[:<derr>] <kappa0>[:<kerr>] <dus2gas>[:<dgerr>]'
           write (iotty,'(a)') '           <catalog> [-o <outname>] [{-|+}verb{0|1|2}]'
           write (iotty,'(a)')
           write (iotty,'(a)') ' This utility fits fluxes of sources from a catalog or of pixels from two or'
           write (iotty,'(a)') ' more images. The source catalogs must have column headers consistent with'
           write (iotty,'(a)') ' those found in catalogs produced by GETSF. The images must first be'
           write (iotty,'(a)') ' convolved to various band resolutions, preferably by using a Bash script'
           write (iotty,'(a)') ' HIGHRES to prepare and create higher-resolution surface density images.'
           write (iotty,'(a)')
           write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
           write (iotty,'(a)')
           write (iotty,'(a)') ' {<number>|<xcoo>} .......... numb of a source/pixel to fit (0 - all sources/pixels)'
           write (iotty,'(a)') ' {<sigmin>|<npart>|<ycoo>} .. min source significance | part of all pixels to fit'
           write (iotty,'(a)') ' {<goodmin>|<nparts>} ....... min global goodness | total number of equal parts'
           write (iotty,'(a)') ' [+|-]<n1>:<wave1> .......... numb & value of waveband1 (um); [+|-] to fit or skip'
           write (iotty,'(a)') ' [+|-]<n2>:<wave2> .......... numb & value of waveband2 (um); [+|-] to fit or skip'
           write (iotty,'(a)') ' [+|-]<nN>:<waveN> .......... numb & value of wavebandN (um); [+|-] to fit or skip'
           write (iotty,'(a)') ' <snr1>:<adderr1> ........... min S/N ratio; additional flux uncertainty (relative)'
           write (iotty,'(a)') ' <snr2>:<adderr2> ........... min S/N ratio; additional flux uncertainty (relative)'
           write (iotty,'(a)') ' <snrN>:<adderrN> ........... min S/N ratio; additional flux uncertainty (relative)'
           write (iotty,'(a)') ' {thinbody|modbody} ......... name of the fitting model known to the code'
           write (iotty,'(a)') ' {<pval1>|0}:{<dlp1>|1|0} ... param1 value or 0; dlogp1max or 1 or 0 (vari/const)'
           write (iotty,'(a)') ' {<pval2>|0}:{<dlp2>|1|0} ... param2 value or 0; dlogp2max or 1 or 0 (vari/const)'
           write (iotty,'(a)') ' {<pval3>|0}:{<dlp3>|1|0} ... param3 value or 0; dlogp3max or 1 or 0 (vari/const)'
           write (iotty,'(a)') ' {<pval4>|0}:{<dlp4>|1|0} ... param4 value or 0; dlogp4max or 1 or 0 (vari/const)'
           write (iotty,'(a)') ' {<distpc>|1}[:<derr>] ...... distance to sources (pc); use <distpc> = 1 for images'
           write (iotty,'(a)') ' <kappa0>[:<kerr>] .......... ref opacity value (per gram of *dust*) at 1000 GHz'
           write (iotty,'(a)') ' <dus2gas>[:<dgerr>] ........ dust-to-gas ratio and its relative error (optional)'
           write (iotty,'(a)') ' <catalog> .................. catalog of sources or images with proper column headers'
           write (iotty,'(a)') ' [-o <outname>] ............. optional base name of output files produced by the code'
           write (iotty,'(a)') ' [{-|+}verb{0|1|2}] ......... verbosity level; "+" to produce postscript (pdf) plots'
           write (iotty,'(a)')
           write (iotty,'(a)') ' Columns with headers NO, XCO_P, YCO_P, SIG_GLOB, GOOD, FXP_BEST##, FXP_ERRO##'
           write (iotty,'(a)') ' FXT_BEST##, FXT_ERRO## must be present in the catalog (## means the number of'
           write (iotty,'(a)') ' the waveband, with a leading zero). Columns SIG_MONO## (if present) will be'
           write (iotty,'(a)') ' used for selecting sources to be fitted.'
           write (iotty,'(a)')
           write (iotty,'(a)') ' The code can use upper limits (instead of fluxes) to constrain acceptable'
           write (iotty,'(a)') ' fits by requiring them to always pass *below* upper limits. Users can add a'
           write (iotty,'(a)') ' "+" sign in front of flux values in the source catalog to indicate that they'
           write (iotty,'(a)') ' are upper limits (i.e., 35.647 is flux, but +35.647 is an upper limit).'
           write (iotty,'(a)')
           write (iotty,'(a)') ' Parameters of the thinbody model (for fitting sources): dust temperature (K),'
           write (iotty,'(a)') ' source mass (Msun), opacity slope beta (> 0).'
           write (iotty,'(a)')
           write (iotty,'(a)') ' Parameters of the modbody model (for fitting sources): dust temperature (K),'
           write (iotty,'(a)') ' source mass (Msun), opacity slope beta (> 0), solid angle diameter (arcsec).'
           write (iotty,'(a)')
           write (iotty,'(a)') ' Parameters of both thinbody and modbody models (for fitting images): dust'
           write (iotty,'(a)') ' temperature (K), surface density (cm^-2), opacity slope beta (> 0).'
           write (iotty,'(a)')
           write (iotty,'(a)') ' If a parameter value (pval1-pval4) is 0, then its optimal initial value for'
           write (iotty,'(a)') ' fitting is found by automatically exploring the multi-dimensional parameter'
           write (iotty,'(a)') ' space. It is possible to specify a fixed value for parameters or a range of'
           write (iotty,'(a)') ' parameter values, i.e., temperature <pval1> may be entered as 10 or 10~30.'
           write (iotty,'(a)') ' The opacity slope is generally to be kept fixed and its actual value must'
           write (iotty,'(a)') ' be supplied in pval3. It is suggested that one always uses the default'
           write (iotty,'(a)') ' values by setting the parameters 0:1,0:1,2:0 (thinbody) or 0:1,0:1,2:0,0:1'
           write (iotty,'(a)') ' (modbody).'
           write (iotty,'(a)')
           write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXAMPLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
           write (iotty,'(a)')
           write (iotty,'(a)') ' DERIVE MASSES AND TEMPERATURES OF SOURCES FROM AN EXTRACTION CATALOG:'
           write (iotty,'(a)')
           write (iotty,'(a)') ' fitfluxes 0 1 1 -1:070,-2:160,3:250,4:350,5:500 1:0.2,1:0.2,1:0.2,1:0.2,1:0.2 thinbody'
     &                       //' 0:1,0:1,2:0 140 10 0.01 <fieldname>.sw.sources.ok.cat +verb2'
           write (iotty,'(a)') ' fitfluxes 0 1 1 -1:070,-2:160,3:250,4:350,5:500 1:0.2,1:0.2,1:0.2,1:0.2,1:0.2 modbody'
     &                       //' 0:1,0:1,2:0,0:1 140 10 0.01 <fieldname>.sw.sources.ok.cat +verb2'
           write (iotty,'(a)')
           write (iotty,'(a)') ' DERIVE SURFACE DENSITY IMAGE FROM PIXEL-BY-PIXEL FITTING OF IMAGES:'
           write (iotty,'(a)')
           write (iotty,'(a)') ' fitfluxes 0 0 0 -1:070,-2:160,3:250,4:350,5:500 1:0.2,1:0.2,1:0.2,1:0.2,1:0.2 thinbody'
     &                       //' 0:1,0:1,2:0 1 10 0.01 surfden.cfg -verb0'
           write (iotty,'(a)')
           write (iotty,'(a)') ' DERIVE SURFACE DENSITY VALUE BY FITTING A SINGLE PIXEL (51,73) OF IMAGES:'
           write (iotty,'(a)')
           write (iotty,'(a)') ' fitfluxes 51 73 0 -1:070,-2:160,3:250,4:350,5:500 1:0.2,1:0.2,1:0.2,1:0.2,1:0.2 thinbody'
     &                       //' 0:1,0:1,2:0 1 10 0.01 surfden.cfg +verb2'
           write (iotty,'(a)')
           write (iotty,'(a)') ' DERIVE SURFACE DENSITY IMAGE FROM PIXEL-BY-PIXEL FITTING OF LARGE IMAGES,'
           write (iotty,'(a)') ' ACCELERATING THE ENTIRE PROCESS BY SPLITTING IT IN 3 INDEPENDENT PARTS'
           write (iotty,'(a)') ' (PARALLEL JOBS):'
           write (iotty,'(a)')
           write (iotty,'(a)') ' fitfluxes 0 1 3 -1:070,-2:160,3:250,4:350,5:500 1:0.2,1:0.2,1:0.2,1:0.2,1:0.2 thinbody'
     &                       //' 0:1,0:1,2:0 1 10 0.01 surfden.cfg -verb0'
           write (iotty,'(a)') ' fitfluxes 0 2 3 -1:070,-2:160,3:250,4:350,5:500 1:0.2,1:0.2,1:0.2,1:0.2,1:0.2 thinbody'
     &                       //' 0:1,0:1,2:0 1 10 0.01 surfden.cfg -verb0'
           write (iotty,'(a)') ' fitfluxes 0 3 3 -1:070,-2:160,3:250,4:350,5:500 1:0.2,1:0.2,1:0.2,1:0.2,1:0.2 thinbody'
     &                       //' 0:1,0:1,2:0 1 10 0.01 surfden.cfg -verb0'
           write (iotty,'(a)')
           write (iotty,'(a)') ' PARTIAL IMAGES RESULTING FROM THE ABOVE SPLITTING CAN EASILY BE COMBINED'
           write (iotty,'(a)') ' INTO COMPLETE IMAGES USING THE OPERATE UTILITY OR THE HIRES SCRIPT.'
           write (iotty,'(a)')
           write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
           stop 99
         endif
       endif

       dust2gas = 0.01d0
       opacity0 = 10.0d0  !<-- Default value for opacity per gram of dust (cm^2/g).
       opaslope = 2.0d0   !<-- Default value for long-wavelengh opacity slope. 

       read (arg1(1:ia1),*,err=10) number
       read (arg2(1:ia2),*,err=12) sigmin
       read (arg3(1:ia3),*,err=14) goodmin

       n1 = index ( arg8(1:ia8), ':' )
       if (n1 .eq. 0) then
         read (arg8(1:ia8),*,err=16) distancepc
         distancrelerr = 0.0d0  !<- 0.0 assuming that we are interested in sources of a certain region at the same distance.
       else
         read (arg8(1:n1-1),*,err=16) distancepc
         read (arg8(n1+1:ia8),*,err=16) distancrelerr
       endif

       n1 = index ( arg9(1:ia9), ':' )
       if (n1 .eq. 0) then
         read (arg9(1:ia9),*,err=18) opacity0
         opacit0relerr = 0.2d0
       else
         read (arg9(1:n1-1),*,err=18) opacity0
         read (arg9(n1+1:ia9),*,err=18) opacit0relerr
       endif
       
       n1 = index ( argA(1:iaA), ':' )
       if (n1 .eq. 0) then
         read (argA(1:iaA),*,err=20) dust2gas
         dustgasrelerr = 0.2d0
       else
         read (argA(1:n1-1),*,err=20) dust2gas
         read (argA(n1+1:iaA),*,err=20) dustgasrelerr
       endif

       callwaves = arg4(1:ia4)
       icw = ia4
       callsnerr = arg5(1:ia5)
       ics = ia5
       cmodel = arg6(1:ia6)
       icm = ia6
       callparams = arg7(1:ia7)
       icp = ia7
       excat = argB(1:iaB)
       iexc = iaB

! Check if the input file exists and open it.

       inquire ( file=excat(1:iexc), exist=lfname1 )
       
       if (excat(iexc-3:iexc) .eq. '.cfg') then
         limages = .true.
         if (iverbose .eq. 0) then
           iolog = 3
         else
           iolog = 0
         endif
         npart = nint ( sigmin)
         nparts = nint ( goodmin )
         na = 1
         if (number .eq. 0 .and. npart .gt. nparts .and. iotty .gt. 0) then
           write (iotty,'(/a)') ' FITFLUXES: ERROR: number = 0 & npart > nparts'
           stop 99
         endif
         if (number .ne. 0 .and. npart .eq. 0 .and. nparts .eq. 0) then
           iolog = 3
           write (cnumber,'(i8.8)') number
         else
           if (number .ne. 0 .and. npart .ne. 0 .and. nparts .eq. 0) then
             iolog = 3
             write (cnumber(1:4),'(i4.4)') number
             write (cnumber(5:8),'(i4.4)') npart
           endif
         endif
       else
         limages = .false.
         iolog = 3
         npart = 0
         nparts = 0
         na = 3
         write (cnumber,'(i8.8)') number
       endif

       if (limages .and. distancepc .gt. 1.001d0) then
         if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: For fitting images, the distance value must be 1.0'
         if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: For fitting images, the distance value must be 1.0'
         if (iolog .gt. 0) close ( iolog )
         stop 99
       endif

       goto 100
  10   continue
         if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: Trouble getting NUMBER'
         stop 99
  12   continue
         if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: Trouble getting SIGMIN'
         stop 99
  14   continue
         if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: Trouble getting GOODMIN'
         stop 99
  16   continue
         if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: Trouble getting DISTPC or DERR'
         stop 99
  18   continue
         if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: Trouble getting KAPPA0 or KERR'
         stop 99
  20   continue
         if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: Trouble getting DUST2GAS or DGERR'
         stop 99
 100   continue
       
       if (argC(1:2) .eq. '-o') then
         if (argD(iaD-3:iaD) .ne. '.cat') then
           outname = argD(1:iaC)//'.00'
         else
           outname = argD(1:iaD-4)//'.00'
         endif
       else
         if (excat(iexc-3:iexc) .ne. '.cat') then
           outname = excat(1:iexc)//'.00'
         else
           outname = excat(1:iexc-4)//'.00'
         endif
       endif
       ion = lastc ( outname )
       
!!       if (limages) then
!!         prefix = outname(1:ion)
!!       else
       prefix = cmodel(1:4)//'.'//outname(1:ion)
!!       endif
       ipr = lastc ( prefix )

       do i=0,99
         write (cno,'(i2.2)') i
         prefix = prefix(1:ipr-3)//'.'//cno
         inquire ( file=prefix(1:ipr)//'.log', exist=lfnameo )
         if (.not.lfnameo) then
           if (limages .and. npart .eq. 0) exit
           inquire ( file=prefix(1:ipr)//'.'//cnumber(na:)//'.log', exist=lfnameo )
           if (.not.lfnameo) exit
           do j=i+1,99
             write (cno,'(i2.2)') j
             prefix = prefix(1:ipr-3)//'.'//cno
             inquire ( file=prefix(1:ipr)//'.'//cnumber(na:)//'.log', exist=lfnameo )
             if (.not.lfnameo) goto 1 
           enddo
         endif
       enddo
  1    continue

       ipr = lastc ( prefix )
       if (number .gt. 0) then
         prefix = prefix(1:ipr)//'.'//cnumber(na:)
       endif
       ipr = lastc ( prefix )
       postscript = prefix(1:ipr)//'.ps'
       paramsfile = prefix(1:ipr)//'.spa'
       chisquarnu = prefix(1:ipr)//'.chisqnu.fits'
       surfdensit = prefix(1:ipr)//'.surfden.fits'
       rsurferror = prefix(1:ipr)//'.surferr.fits'
       temperatur = prefix(1:ipr)//'.tempers.fits'
       rtemperror = prefix(1:ipr)//'.temperr.fits'
       paraspace  = prefix(1:ipr)//'.space.fits'
       ilp = lastc ( postscript )
       ilc = lastc ( surfdensit )
       ilq = lastc ( chisquarnu )
       ilt = lastc ( temperatur )
       ilr = lastc ( rtemperror )
       ilm = lastc ( rsurferror )
       ils = lastc ( paramsfile )
       ile = lastc ( paraspace  )
       ncx = 1
       ncy = 1
       
       if (iolog .gt. 0) then
         open ( iolog, file=postscript(1:ilp-3)//'.log', status='unknown' )
       endif

       commandline = 'fitfluxes '//arg1(1:ia1)//' '//arg2(1:ia2)//' '//arg3(1:ia3)//' '//arg4(1:ia4)//' '//arg5(1:ia5)//' '
     &             //arg6(1:ia6)//' '//arg7(1:ia7)//' '//arg8(1:ia8)//' '//arg9(1:ia9)//' '//argA(1:iaA)//' '//argB(1:iaB)//' '
     &             //argC(1:iaC)//' '//argD(1:iaD)//' '//argE(1:iaE)
       icl = lastc ( commandline )

       if (limages) then
         write (iotty,'(/a)') ' '//commandline(1:icl)
       else
         if (iverbose .eq. 0) then
           write (iotty,'(a)') ' FITFLUXES: '//commandline(11:icl)
         endif
       endif

       if (iolog .gt. 0) then
         write (iolog,'(a)') ' FITFLUXES '//dot//' Fit Source Fluxes or Image Intensities '//dot//' '//compda
         write (iolog,'(a)') ' Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
         write (iolog,'(a)') ' Using MRQMIN from Numerical Recipes for F77 by William H Press et al.'
         write (iolog,'( )')
         write (iolog,'(a)') ' '//commandline(1:icl)
       endif

       if (cmodel .ne. 'modbody' .and. cmodel .ne. 'thinbody') then
         if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: Unknown name of a fitting model: '//cmodel
         if (iolog .gt. 0) write (iolog,'(/a)') ' FITFLUXES: ERROR: Unknown name of a fitting model: '//cmodel
         if (iolog .gt. 0) close ( iolog )
         stop 99
       endif
       if (.not.lfname1) then
         if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: File '''//excat(1:iexc)//''' not found.'
         if (iolog .gt. 0) write (iolog,'(/a)') ' FITFLUXES: ERROR: File '''//excat(1:iexc)//''' not found.'
         if (iolog .gt. 0) close ( iolog )
         stop 1
       endif

       allocate ( headline (1000), stat=irc )

       if (irc .ne. 0) then
         if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory.'
         if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory.'
         if (iolog .gt. 0) close ( iolog )
         stop
       endif
       
       open ( 15, file=excat(1:iexc), status='old' )

       call skipcomm ( 15, excat(1:iexc), nhead )
       call linescolumns ( 15, excat(1:iexc), nhead, nextr )

       rewind ( 15 )
       do k=1,nhead
         read (15,'(a)') headline(k)
         lhl = lastc ( headline(k) )
         if (lhl .gt. 5) nhcols = k
       enddo
       lhl = lastc ( headline(nhcols) )
       solidapixel = 0.0d0

       if (limages) then
         read (15,*,err=333) wave, beam1, offset, filename
         goto 334
 333     continue
         if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Trouble reading '''//excat(1:iexc)//''''
         if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Trouble reading '''//excat(1:iexc)//''''
         if (iolog .gt. 0) close ( iolog )
         stop 333
 334     continue
         backspace ( 15 )

         lhl = lastc ( filename )

! Check if the input file exists and open it.

         inquire ( file=filename(1:lhl), exist=lfname2 )
         
         if (.not.lfname2) then
           if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: File '''//filename(1:lhl)//''' not found.'
           if (iolog .gt. 0) write (iolog,'(/a)') ' FITFLUXES: ERROR: File '''//filename(1:lhl)//''' not found.'
           if (iolog .gt. 0) close ( iolog )
           stop 2
         endif

! Determine numbers of pixels in the first FITS image.
  
         call getfitshead ( filename(1:lhl), nx1, ny1, dx, dy, bunit )

         solidapixel = dx * dy / sq_arcsecs_per_sterad
         nx0o = 1
         ny0o = 1
         if (number .ne. 0 .and. npart .eq. 0 .and. nparts .eq. 0) then
           ny0o = int ( dble ( (number - 1) / nx1 ) ) + 1
           nx0o = int ( dble ( number ) ) - int ( dble ( nx1 ) * (ny0o - 1) )
         endif
         if (number .ne. 0 .and. npart .ne. 0 .and. nparts .eq. 0) then
           nx0o = number
           ny0o = npart
           number = (ny0o - 1) * nx1 + nx0o
           npart = 0
           arg2 = '0'
           ia2 = 1
         endif
         write (cx,'(i4)') nx0o
         write (cy,'(i4)') ny0o
         ncx = int ( log10 ( dble ( nx0o ) ) ) + 1
         ncy = int ( log10 ( dble ( ny0o ) ) ) + 1
       endif
       write (cnsrc,'(i10)') number
       if (number .ne. 0) then
         nn = int ( log10 ( dble ( number ) ) ) + 1
       else
         nn = 1
       endif

       if (ltiming) then
         cpu = 0.0d0
         wal = 0.0d0
         cpu_sum = 0.0d0
         wal_sum = 0.0d0
         if (iotty .gt. 0 .and. iverbose .eq. 2)
     &                     write (iotty,'(/a)') ' STARTED: '//cdate//' '//ctime
         if (iolog .gt. 0) write (iolog,'(/a)') ' STARTED: '//cdate//' '//ctime
       endif
       if (iotty .gt. 0 .and. iverbose .gt. 0)
     &     write (iotty,'(a)') ' ____________________________________________________________________________'
       if (iolog .gt. 0)
     &     write (iolog,'(a)') ' ____________________________________________________________________________'

! Report accepted values of the command-line parameters.

       if (number .ne. 0) then
         if (.not.limages) then
           if (iotty .gt. 0 .and. iverbose .gt. 0) 
     &                       write (iotty,'(/a)') ' Fitting flux distribution for single source # '//cnsrc(10-nn+1:10)
           if (iolog .gt. 0) write (iolog,'(/a)') ' Fitting flux distribution for single source # '//cnsrc(10-nn+1:10)
         else
           if (iotty .gt. 0 .and. iverbose .gt. 0) 
     &                       write (iotty,'(/a)') ' Fitting flux distribution for single pixel # '//cnsrc(10-nn+1:10)//' ('
     &                                          //cx(4-ncx+1:4)//','//cy(4-ncy+1:4)//')'
           if (iolog .gt. 0) write (iolog,'(/a)') ' Fitting flux distribution for single pixel # '//cnsrc(10-nn+1:10)//' ('
     &                                          //cx(4-ncx+1:4)//','//cy(4-ncy+1:4)//')'
         endif
       else
         if (.not.limages) then
           if (iotty .gt. 0 .and. iverbose .gt. 0)
     &                       write (iotty,'(/a)') ' Fitting flux distributions for all selected sources'
           if (iolog .gt. 0) write (iolog,'(/a)') ' Fitting flux distributions for all selected sources'
         else
           if (npart .eq. 0) then
             if (iotty .gt. 0) !! .and. iverbose .gt. 0
     &                         write (iotty,'(/a)') ' Fitting intensities for all pixels of images'
             if (iolog .gt. 0) write (iolog,'(/a)') ' Fitting intensities for all pixels of images'
           else
             if (iotty .gt. 0) !! .and. iverbose .gt. 0
     &                         write (iotty,'(/a)') ' Fitting intensities for selected pixels of images'
             if (iolog .gt. 0) write (iolog,'(/a)') ' Fitting intensities for selected pixels of images'
           endif
         endif
       endif
       if (.not.limages) then
         if (iotty .gt. 0 .and. iverbose .gt. 0)
     &                     write (iotty,'(/a)') ' Selection of sources: ... global significance > '//arg2(1:ia2)
         if (iolog .gt. 0) write (iolog,'(/a)') ' Selection of sources: ... global significance > '//arg2(1:ia2)
         if (iotty .gt. 0 .and. iverbose .gt. 0)
     &                     write (iotty,'(a )') ' Selection of sources: ....... global goodness > '//arg3(1:ia3)
         if (iolog .gt. 0) write (iolog,'(a )') ' Selection of sources: ....... global goodness > '//arg3(1:ia3)
       else
         if (iotty .gt. 0) !! .and. iverbose .gt. 0
     &                     write (iotty,'(/a)') ' Selection of pixels: .... selected part number: '//arg2(1:ia2)
         if (iolog .gt. 0) write (iolog,'(/a)') ' Selection of pixels: .... selected part number: '//arg2(1:ia2)
         if (iotty .gt. 0) !! .and. iverbose .gt. 0
     &                     write (iotty,'(a )') ' Selection of pixels: ... total number of parts: '//arg3(1:ia3)
         if (iolog .gt. 0) write (iolog,'(a )') ' Selection of pixels: ... total number of parts: '//arg3(1:ia3)
       endif

! Specify all wavelengths for fitting.

       n1 = 1
       n2 = 1
       n3 = 0
       nbands = 0
       do i=1,100
         n2 = index ( callwaves(n1:icw), ':' )
         n3 = index ( callwaves(n1:icw), ',' )
         if (n2 .gt. 0 .and. n3 .gt. 0) then
           n2 = n1 - 1 + n2
           n3 = n1 - 1 + n3
           read (callwaves(n1:n2-1  ),*,err=30,end=30) nband(i)
           read (callwaves(n2+1:n3-1),*,err=30,end=30) waveband(i)
           csnband(i) = ''
           if (callwaves(n1:n1) .eq. '+') csnband(i) = '+'
           if (callwaves(n1:n1) .eq. '-') csnband(i) = '-'
           nband(i) = abs ( nband(i) )
           n1 = n3 + 1
         else
           n2 = n1 - 1 + n2
           read (callwaves(n1:n2-1 ),*,err=30,end=30) nband(i)
           read (callwaves(n2+1:icw),*,err=30,end=30) waveband(i)
           csnband(i) = ''
           if (callwaves(n1:n1) .eq. '+') csnband(i) = '+'
           if (callwaves(n1:n1) .eq. '-') csnband(i) = '-'
           nband(i) = abs ( nband(i) )
           nbands = i
         endif
         if (nbands .eq. i) exit
       enddo

       goto 33
 30    continue
         if (iotty .gt. 0) write (iotty,'(a)')  ' FITFLUXES: ERROR: Trouble getting wavelengths.'
         if (iolog .gt. 0) write (iolog,'(a)')  ' FITFLUXES: ERROR: Trouble getting wavelengths.'
         if (iolog .gt. 0) close ( iolog )
         stop 99
 33    continue  
       if (nbands .le. 1) then
         if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Too few wavebands for fitting.'
         if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Too few wavebands for fitting.'
         if (iolog .gt. 0) close ( iolog )
         stop 99
       endif
       
! Specify selection S/N and additional errors, such as calibration uncertainties, for each band.

       n1 = 1
       n2 = 1
       n3 = 0
       nbande = 0
       if (iotty .gt. 0 .and. (iverbose .gt. 0 .or. limages)) write (iotty,'()')
       if (iolog .gt. 0) write (iolog,'()')
       do i=1,100
         n2 = index ( callsnerr(n1:ics), ':' )
         n3 = index ( callsnerr(n1:ics), ',' )
         if (n2 .gt. 0 .and. n3 .gt. 0) then
           n2 = n1 - 1 + n2
           n3 = n1 - 1 + n3
           read (callsnerr(n1:  n2-1),*,err=40,end=40) snrmin(i)
           read (callsnerr(n2+1:n3-1),*,err=40,end=40) adderro(i)
           n1 = n3 + 1
         else
           n2 = n1 - 1 + n2
           read (callsnerr(n1: n2-1),*,err=40,end=40) snrmin(i)
           read (callsnerr(n2+1:ics),*,err=40,end=40) adderro(i)
           nbande = i
         endif
         if (iverbose .gt. 0 .or. limages) then
           if (iotty .gt. 0) write (iotty,'(a,i2.2,i4,a,f5.2,a,f4.1)')  ' Waveband'//csnband(i), nband(i), nint( waveband(i) )
     &                                   , ' µm, additional relative error:', adderro(i), ', peak S/N >', snrmin(i)
         endif
         if (iolog .gt. 0) write (iolog,'(a,i2.2,i4,a,f5.2,a,f4.1)')  ' Waveband'//csnband(i), nband(i), nint( waveband(i) )
     &                                 , ' µm, additional relative error:', adderro(i), ', peak S/N >', snrmin(i)
         if (nbande .eq. i) exit                                                                 
       enddo
       
       goto 44
 40    continue
       if (i .lt. nbands) then
         if (iotty .gt. 0) write (iotty,'(a)')  ' FITFLUXES: ERROR: Trouble getting additional errors.'
         if (iolog .gt. 0) write (iolog,'(a)')  ' FITFLUXES: ERROR: Trouble getting additional errors.'
         if (iolog .gt. 0) close ( iolog )
         stop 99
       endif
 44    continue
       if (nbande .ne. nbands) then
         if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Invalid number of additional errors.'
         if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Invalid number of additional errors.'
         if (iolog .gt. 0) close ( iolog )
         stop 99
       else
         if (iotty .gt. 0 .and. (iverbose .gt. 0 .or. limages)) write (iotty,'()')
         if (iolog .gt. 0) write (iolog,'()')
       endif
       
       if (iverbose .gt. 0 .or. limages) then
         if (iotty .gt. 0) then
           write (iotty,'(a,f6.3,f6.2)') ' Adopted dust-to-gas mass ratio and relative error: ...', dust2gas, dustgasrelerr
           write (iotty,'(a,f6.2,f6.2)') ' Adopted reference dust opacity and relative error: ...', opacity0, opacit0relerr
           write (iotty,'(a,f6.1,f6.2)') ' Adopted distance (pc) and relative error: ............', distancepc, distancrelerr
         endif
       endif
       if (iolog .gt. 0) then
         write (iolog,'(a,f6.3,f6.2)') ' Adopted dust-to-gas mass ratio and relative error: ...', dust2gas, dustgasrelerr
         write (iolog,'(a,f6.2,f6.2)') ' Adopted reference dust opacity and relative error: ...', opacity0, opacit0relerr
         write (iolog,'(a,f6.1,f6.2)') ' Adopted distance (pc) and relative error: ............', distancepc, distancrelerr
       endif

       allocate ( ncolsig(nbands), ncolfxp(nbands), ncolfxt(nbands), ncolsnr(nbands), ncolfpe(nbands), ncolfte(nbands)
     &          , snrpeak(nbands), snrtotal(nbands), frequency(nbands), fxpbest(nbands), fxperro(nbands), fxtbest(nbands)
     &          , fxterro(nbands), opacity(nbands), sedslope(nbands), wavelength(nbands), sfrequency(nbands), sfxtbest(nbands)
     &          , sfxterro(nbands), sigmono(nbands), beam(nbands), fxterror(nbands), uplimit(nbands), stat=irc ) !, fxperror(nbands)

       if (irc .ne. 0) then
         if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory (5).'
         if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory (5).'
         if (iolog .gt. 0) close ( iolog )
         stop 5
       endif
       
       do i=1,nbands
         wavelength(i) = waveband(i)
       enddo
       do i=1,nparams
         ivarornot(i) = 0
         modparams1(i) = 0.0d0
         modparams2(i) = 0.0d0
         dlogp(i) = 0.0d0
         dlogpmax(i) = 0.0d0
       enddo
       do i=1,100
         quantities(i) = 0.0d0
       enddo

! Number of parameters NMODPARAMS in the fitting model.
! Specify if the corresponding fitting parameter should be varied during the fitting (=1) or held fixed (=0).

       n1 = 1
       n2 = 1
       n3 = 0
       nm = 0
       nmodparams = 0
       do i=1,nparams
         n2 = index ( callparams(n1:icp), ':' )
         n3 = index ( callparams(n1:icp), ',' )
         nm = index ( callparams(n1:icp), '~' )
         if (n2 .gt. 0 .and. n3 .gt. 0) then
           n2 = n1 - 1 + n2
           n3 = n1 - 1 + n3
           if (nm .gt. 0) then
             nm = n1 - 1 + nm
           else
             nm = n2
           endif
           read (callparams(n1:nm-1  ),*) modparams1(i)
           read (callparams(n2+1:n3-1),*) dlogpmax(i)
           if (dlogpmax(i) .gt. almostzero) then
             ivarornot(i) = 1
           else
             ivarornot(i) = 0
           endif
           n1 = n3 + 1
         else
           n2 = n1 - 1 + n2
           if (nm .gt. 0) then
             nm = n1 - 1 + nm
           else
             nm = n2
           endif
           read (callparams(n1:nm-1 ),*) modparams1(i)
           read (callparams(n2+1:icp),*) dlogpmax(i)
           if (dlogpmax(i) .gt. almostzero) then
             ivarornot(i) = 1
           else
             ivarornot(i) = 0
           endif
           nmodparams = i
         endif
         if (nm .lt. n2) then
           read (callparams(nm+1:n2-1),*) modparams2(i)
         else
           modparams2(i) = modparams1(i)
         endif
         if (nmodparams .eq. i) exit
       enddo

       nmodparams = max ( nmodparams, nparams )

       if (limages .and. cmodel .eq. 'modbody') then
         ivarornot(4) = 0
       endif       
       
       if (nmodparams .gt. nparams) then
         if (iotty .gt. 0) write (iotty,'(/a,i1,a)')  ' FITFLUXES: ERROR: Too many parameters (must be <= ', nparams, ').'
         if (iolog .gt. 0) write (iolog,'(/a,i1,a)')  ' FITFLUXES: ERROR: Too many parameters (must be <= ', nparams, ').'
         if (iolog .gt. 0) close ( iolog )
         stop 99
       endif
       
! Prepare the header of the output parameter space file for all selected sources.

       if (ldebug) then
         open ( 55, file=paramsfile(1:ils), status='unknown' )

         write (55,'(a)') '# FITFLUXES '//dot//' Fit Source Fluxes or Image Intensities '//dot//' '//compda
         write (55,'(a)') '# Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
         write (55,'(a)') '# Using MRQMIN from Numerical Recipes for F77 by William H Press et al.'
         write (55,'(a)') '#'
         write (55,'(a)') '# '//commandline(1:icl)
       endif

! Prepare the header of the output data file with results of final fitting for all selected sources.

       if (.not.limages .and. ldebug .and. lplotting) then
         open ( 66, file=postscript(1:ilp-3)//'.dat', status='unknown' )
 
         write (66,'(a)') '# FITFLUXES '//dot//' Fit Source Fluxes or Image Intensities '//dot//' '//compda
         write (66,'(a)') '# Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
         write (66,'(a)') '# Using MRQMIN from Numerical Recipes for F77 by William H Press et al.'
         write (66,'(a)') '#'
         write (66,'(a)') '# OBSERVED FLUXES AND FITTED FUNCTIONS FOR SELECTED SOURCES OR IMAGE PIXELS'
         write (66,'(a)') '#'
         write (66,'(a)') '# '//commandline(1:icl)
         write (66,'(a)') '#'
         write (66,'(a)') '# Input catalog: '//excat(1:iexc)
         write (66,'(a)') '#  Output files: '//prefix(1:ipr)//'.*'
 
         write (66,'(a)') '#'
         if (.not.limages) then
         write (66,'(a)') '# Selected sources: ... global significance > '//arg2(1:ia2)
         write (66,'(a)') '# Selected sources: ....... global goodness > '//arg3(1:ia3)
         else
         write (66,'(a)') '# Selected pixels: ......... part of all pixels to fit: '//arg2(1:ia2)
         write (66,'(a)') '# Selected pixels: ... number of equal parts of pixels: '//arg3(1:ia3)
         endif
         write (66,'(a)') '#'
         write (66,'(a,50(1pe9.2))') '#  Wavebands for SED fitting:', (waveband(i),i=1,nbands)
         write (66,'(a,50(1pe9.2))') '# Additional relative errors:', (adderro(i),i=1,nbands)
         write (66,'(a)') '#'
         write (66,'(a,2(1pe9.2))') '# Adopted dust-to-gas mass ratio and relative error: ...', dust2gas, dustgasrelerr
         write (66,'(a,2(1pe9.2))') '# Adopted reference dust opacity and relative error: ...', opacity0, opacit0relerr
         write (66,'(a,2(1pe9.2))') '# Adopted distance (pc) and relative error: ............', distancepc, distancrelerr
         write (66,'(a)') '#'
         if (cmodel .eq. 'modbody') then
         write (66,'(a)') '# Fitting model: '//cmodel//' ~> Fnu = Bnu(T) (1-exp(-tau)) Omega'
         write (66,'(a)') '#'
         write (66,'(a)') '# 1st parameter ....... temperature: '//yesno(ivarornot(1)+1)
         if (.not.limages) then
         write (66,'(a)') '# 2nd parameter ........ total mass: '//yesno(ivarornot(2)+1)
         else
         write (66,'(a)') '# 2nd parameter ... surface density: '//yesno(ivarornot(2)+1)
         endif
         write (66,'(a)') '# 3rd parameter ..... opacity slope: '//yesno(ivarornot(3)+1)
         write (66,'(a)') '# 4th parameter ....... solid angle: '//yesno(ivarornot(4)+1)
         write (66,'(a)') '# 5th parameter .......... distance: '//yesno(ivarornot(5)+1)
         endif
         if (cmodel .eq. 'thinbody') then
         write (66,'(a)') '# Fitting model: '//cmodel//' ~> Fnu = Bnu(T) kappa M D^-2'
         write (66,'(a)') '#'
         write (66,'(a)') '# 1st parameter ....... temperature: '//yesno(ivarornot(1)+1)
         if (.not.limages) then
         write (66,'(a)') '# 2nd parameter ........ total mass: '//yesno(ivarornot(2)+1)
         else
         write (66,'(a)') '# 2nd parameter ... surface density: '//yesno(ivarornot(2)+1)
         endif
         write (66,'(a)') '# 3rd parameter ..... opacity slope: '//yesno(ivarornot(3)+1)
         write (66,'(a)') '# 4th parameter .......... distance: '//yesno(ivarornot(4)+1)
         endif
         write (66,'(a)') '#'
         write (66,'(a)') '# TABULATED QUANTITIES:'
         write (66,'(a)') '#'
         write (66,'(a)') '#   NO ................ Point number'
         write (66,'(a)') '#   WAVEO ............. Observed wavelength (um)'
         write (66,'(a)') '#   OBSPOINT .......... Observed fluxes or intensities (Jy or Jy/sr)'
         write (66,'(a)') '#   WAVEOF ............ Observed wavelength to fit (um)'
         write (66,'(a)') '#   OBSPOINTF ......... Observed fluxes or intensities to fit (Jy or Jy/sr)'
         if (cmodel .eq. 'modbody') then
         write (66,'(a)') '#   WAVEG ............. Wavelength for modbody (um)'
         write (66,'(a)') '#   MODBODY ........... Fitted modbody fluxes or intensities (Jy or Jy/sr)'
         endif
         if (cmodel .eq. 'thinbody') then
         write (66,'(a)') '#   WAVEM ............. Wavelength for thinbody (um)'
         write (66,'(a)') '#   THINBODY .......... Fitted thinbody fluxes or intensities (Jy or Jy/sr)'
         endif
         write (66,'(a)') '#   WAVEB ............. Wavelength for blackbody (um)'
         write (66,'(a)') '#   BLACKBODY ......... Blackbody at derived temperature (scaled)'
         write (66,'(a)') '#   WAVERJ ............ Wavelength for Rayleigh-Jeans (um)'
         write (66,'(a)') '#   RAY-JEANS ......... Rayleigh-Jeans approximation (scaled)'
         write (66,'(a)') '#   OBSERRP ........... Measurement error+ in OBSPOINT (Jy or Jy/sr)'
         write (66,'(a)') '#   OBSERRM ........... Measurement error- in OBSPOINT (Jy or Jy/sr)'
         write (66,'(a)') '#   OBSERRPF .......... Measurement error+ in OBSPOINTF (Jy or Jy/sr)'
         write (66,'(a)') '#   OBSERRMF .......... Measurement error- in OBSPOINTF (Jy or Jy/sr)'
       endif

! Prepare the header of the output catalog with derived parameters for all selected sources.

       if (.not.limages) then
         open ( 77, file=postscript(1:ilp-3)//'.cat', status='unknown' )

         write (77,'(a)') '#___________________________________________________________________________________________________'
     &                  //'_________________________________________________________________'
         write (77,'(a)') '#'
         write (77,'(a)') '# FITFLUXES '//dot//' Fit Source Fluxes or Image Intensities '//dot//' '//compda
         write (77,'(a)') '# Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
         write (77,'(a)') '# Using MRQMIN from Numerical Recipes for F77 by William H Press et al.'
         write (77,'(a)') '#___________________________________________________________________________________________________'
     &                  //'_________________________________________________________________'
         write (77,'(a)') '#'
         write (77,'(a)') '# CATALOG OF DERIVED QUANTITIES FOR SELECTED SOURCES OR IMAGE PIXELS'
         write (77,'(a)') '#'
         write (77,'(a)') '# '//commandline(1:icl)
         write (77,'(a)') '#'
         write (77,'(a)') '# Input catalog: '//excat(1:iexc)
         write (77,'(a)') '#  Output files: '//prefix(1:ipr)//'.*'
         write (77,'(a)') '#'
         if (.not.limages) then
         write (77,'(a)') '# Selected sources: ... global significance > '//arg2(1:ia2)
         write (77,'(a)') '# Selected sources: ....... global goodness > '//arg3(1:ia3)
         else
         write (77,'(a)') '# Selected pixels: ......... part of all pixels to fit: '//arg2(1:ia2)
         write (77,'(a)') '# Selected pixels: ... number of equal parts of pixels: '//arg3(1:ia3)
         endif
         write (77,'(a)') '#'
         write (77,'(a,50(1pe9.2))') '#  Wavebands for SED fitting:', (waveband(i),i=1,nbands)
         write (77,'(a,50(1pe9.2))') '# Additional relative errors:', (adderro(i),i=1,nbands)
         write (77,'(a)') '#'
         write (77,'(a,2(1pe9.2))') '# Adopted dust-to-gas mass ratio and relative error: ...', dust2gas, dustgasrelerr
         write (77,'(a,2(1pe9.2))') '# Adopted reference dust opacity and relative error: ...', opacity0, opacit0relerr
         write (77,'(a,2(1pe9.2))') '# Adopted distance (pc) and relative error: ............', distancepc, distancrelerr
         write (77,'(a)') '#'
         if (cmodel .eq. 'modbody') then
         write (77,'(a)') '# Fitting model: '//cmodel//' ~> Fnu = Bnu(T) (1-exp(-tau)) Omega'
         write (77,'(a)') '#'
         write (77,'(a)') '# 1st parameter ....... temperature: '//yesno(ivarornot(1)+1)
         if (.not.limages) then
         write (77,'(a)') '# 2nd parameter ........ total mass: '//yesno(ivarornot(2)+1)
         else
         write (77,'(a)') '# 2nd parameter ... surface density: '//yesno(ivarornot(2)+1)
         endif
         write (77,'(a)') '# 3rd parameter ..... opacity slope: '//yesno(ivarornot(3)+1)
         write (77,'(a)') '# 4th parameter ....... solid angle: '//yesno(ivarornot(4)+1)
         write (77,'(a)') '# 5th parameter .......... distance: '//yesno(ivarornot(5)+1)
         endif
         if (cmodel .eq. 'thinbody') then
         write (77,'(a)') '# Fitting model: '//cmodel//' ~> Fnu = Bnu(T) kappa M D^-2'
         write (77,'(a)') '#'
         write (77,'(a)') '# 1st parameter ....... temperature: '//yesno(ivarornot(1)+1)
         if (.not.limages) then
         write (77,'(a)') '# 2nd parameter ........ total mass: '//yesno(ivarornot(2)+1)
         else
         write (77,'(a)') '# 2nd parameter ... surface density: '//yesno(ivarornot(2)+1)
         endif
         write (77,'(a)') '# 3rd parameter ..... opacity slope: '//yesno(ivarornot(3)+1)
         write (77,'(a)') '# 4th parameter .......... distance: '//yesno(ivarornot(4)+1)
         endif
         write (77,'(a)') '#___________________________________________________________________________________________________'
     &                  //'_________________________________________________________________'
         write (77,'(a)') '#'
         write (77,'(a)') '# TABULATED QUANTITIES:'
         write (77,'(a)') '#'
         write (77,'(a)') '#  1 NO .......... Input catalog source or image pixel number ...................'
         write (77,'(a)') '#  2 XCO_P ....... X-coordinate of current source or pixel ...................... pixels'
         write (77,'(a)') '#  3 YCO_P ....... Y-coordinate of current source or pixel ...................... pixels'
         write (77,'(a)') '#  4 SIG_GLOB .... Global significance or part of all pixels to fit .............'
         write (77,'(a)') '#  5 GOOD ........ Global goodness or number of equal parts of pixels ...........'
         write (77,'(a)') '#  6 DUST_TEMP ... Derived dust temperature of source or pixel .................. K'
         write (77,'(a)') '#  7 TEMP_ERRO ... Uncertainty of derived temperature ........................... K'
         if (cmodel .eq. 'modbody') then
         write (77,'(a)') '#  8 TOTL_MASS ... Derived total mass (gas+dust) of source ...................... Msun'
         write (77,'(a)') '#  9 MASS_ERRO ... Uncertainty of derived mass from SED fitting ................. Msun'
         write (77,'(a)') '# 10 MASS_ERRT ... Total mass uncertainty (incl. distance opacity dust/gas) ..... Msun'
         endif
         if (cmodel .eq. 'thinbody') then
         if (.not.limages) then
         write (77,'(a)') '#  8 TOTL_MASS ... Derived total mass (gas+dust) of source ...................... Msun'
         write (77,'(a)') '#  9 MASS_ERRO ... Uncertainty of derived mass from SED fitting ................. Msun'
         write (77,'(a)') '# 10 MASS_ERRT ... Total mass uncertainty (incl. distance opacity dust/gas) ..... Msun'
         else
         write (77,'(a)') '#  8 SURFACE_D ... Derived surface density of source or pixel ................... cm^-2'
         write (77,'(a)') '#  9 SURFD_ERR ... Uncertainty of derived surface density from SED fitting ...... cm^-2'
         write (77,'(a)') '# 10 SURFD_ERT ... Total surface density uncertainty (incl. opacity dust/gas) ... cm^-2'
         endif
         endif
         write (77,'(a)') '# 11 BOLO_LUMI ... Derived bolometric luminosity ................................ Lsun'
         write (77,'(a)') '# 12 LUMI_ERRO ... Uncertainty of derived luminosity ............................ Lsun'
         write (77,'(a)') '# 13 SLD_ANGLE ... Derived solid angle (diameter) ............................... arcsec'
         write (77,'(a)') '# 14 ANGL_ERRO ... Uncertainty of derived solid angle ........................... arcsec'
         write (77,'(a)') '# 15 SURFACE_D ... Derived surface density of source or pixel ................... cm^-2'
         write (77,'(a)') '# 16 SURFD_ERR ... Uncertainty of derived surface density ....................... cm^-2'
         write (77,'(a)') '# 17 OPA_BETA .... Derived (or assumed) opacity slope ...........................'
         write (77,'(a)') '# 18 BETA_ERRO ... Uncertainty of derived slope .................................'
         write (77,'(a)') '# 19 DISTANCE .... Derived (or assumed) distance ................................ pc'
         write (77,'(a)') '# 20 DIST_ERRO ... Uncertainty of derived distance .............................. pc'
         write (77,'(a)') '# 21 WAVE_TAU1 ... Wavelength where optical depth tau=1 ......................... µm'
         write (77,'(a)') '# 22 CHISQUARE ... Actual chi^2 value of the fit ................................'
         write (77,'(a)') '# 23 DEG_FREED ... Degrees of freedom of the fit ................................'
         write (77,'(a)') '# 24 QUALITY ..... Quality of the fit (ok or bad) ...............................'
         write (77,'(a)') '#    SIGMA ....... Standard deviation of the random numbers (true std = 1.0) ....'
         write (77,'(a)') '#    RANDMASSxx .. Simulated random "true" masses within the MASS_ERRT range .... Msun'
         write (77,'(a)') '#___________________________________________________________________________________________________'
     &                  //'_________________________________________________________________'
         write (77,'(a)') '#'
         write (77,'(a)') '#    1      2       3        4        5          6         7          8         9        10         11'
     &                  //'        12         13        14         15        16         17        18         19        20'
     &                  //'         21        22        23       24'

         cline = ''
         icl = 0
         if (.not.limages) then
           write (ccx,'(i7)') ntries
           ncx = int ( log10 ( dble ( ntries ) ) ) + 1
           cline = '  SIGMA :'
           icl = lastc ( cline ) - 1
           do i=1,ntries
             write (ccx,'(i7.7)') i
             cline = cline(1:icl)//' RANDMASS'//ccx(7-ncx+1:7)
             icl = lastc ( cline )
           enddo
         endif

         if (cmodel .eq. 'modbody') then
           write (77,'(a)') '#   NO    XCO_P   YCO_P  SIG_GLOB   GOOD    DUST_TEMP TEMP_ERRO  TOTL_MASS MASS_ERRO MASS_ERRT'
     &                    //'  BOLO_LUMI LUMI_ERRO  SLD_ANGLE ANGL_ERRO  SURFACE_D SURFD_ERR   OPA_BETA BETA_ERRO   DISTANCE'
     &                    //' DIST_ERRO  WAVE_TAU1 CHISQUARE DEG_FREED QUALITY'//cline(1:icl)
         endif
         
         if (cmodel .eq. 'thinbody') then
           if (.not.limages) then
             write (77,'(a)') '#   NO    XCO_P   YCO_P  SIG_GLOB   GOOD    DUST_TEMP TEMP_ERRO  TOTL_MASS MASS_ERRO MASS_ERRT'
     &                      //'  BOLO_LUMI LUMI_ERRO  SLD_ANGLE ANGL_ERRO  SURFACE_D SURFD_ERR   OPA_BETA BETA_ERRO   DISTANCE'
     &                      //' DIST_ERRO  WAVE_TAU1 CHISQUARE DEG_FREED QUALITY'//cline(1:icl)
           else
             write (77,'(a)') '#   NO    XCO_P   YCO_P  SIG_GLOB   GOOD    DUST_TEMP TEMP_ERRO  SURFACE_D SURFD_ERR SURFD_ERT'
     &                      //'  BOLO_LUMI LUMI_ERRO  SLD_ANGLE ANGL_ERRO  SURFACE_D SURFD_ERR   OPA_BETA BETA_ERRO   DISTANCE'
     &                      //' DIST_ERRO  WAVE_TAU1 CHISQUARE DEG_FREED QUALITY'
           endif
         endif
         write (77,'(a)') '#'
       endif

       allocate ( fitparameters(nmodparams), ipfitornot(nmodparams), dyda(nmodparams), covar(nmodparams,nmodparams)
     &          , alpha(nmodparams,nmodparams), paramspacemin(nmodparams), paramspacemax(nmodparams), stat=irc )

       if (irc .ne. 0) then
         if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory (6).'
         if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory (6).'
         if (iolog .gt. 0) close ( iolog )
         stop 6
       endif

       frequency0 = 1.0d12
       wavelength0 = 1.0d4 * speedolight / frequency0
       
! Actual number of fitting parameters.

       nfitparams = 0
       do i=1,nmodparams
         ipfitornot(i) = ivarornot(i)
         nfitparams = nfitparams + ipfitornot(i)
       enddo

! Create grid over the entire multidimensional parameter space.

       punits(1) = ' K'
       punits(2) = ' Msun'
       punits(3) = ''
       if (cmodel .eq. 'modbody') then
         punits(4) = ' arcsec'
         punits(5) = ' cm'
       endif
       if (cmodel .eq. 'thinbody') then
         punits(4) = ' cm'
         punits(5) = ''
       endif
       
! Default and maximum steps over parameter space.
       
       if (abs ( dlogpmax(1) - 1.0d0 ) .lt. 0.00001d0) then
         if (     limages) dlogpmax(1) = 0.5d0
         if (.not.limages) dlogpmax(1) = 0.1d0
       else
         dlogpmax(1) = min ( dlogpmax(1), 1.0d0 )
       endif
       if (abs ( dlogpmax(2) - 1.0d0 ) .lt. 0.00001d0) then
         if (     limages) dlogpmax(2) = 0.5d0
         if (.not.limages) dlogpmax(2) = 0.1d0
       else
         dlogpmax(2) = min ( dlogpmax(2), 1.0d0 )
       endif
       if (abs ( dlogpmax(3) - 1.0d0 ) .lt. 0.00001d0) then
         if (     limages) dlogpmax(3) = 0.5d0
         if (.not.limages) dlogpmax(3) = 0.1d0
       else
         dlogpmax(3) = min ( dlogpmax(3), 1.0d0 )
       endif
       if (abs ( dlogpmax(4) - 1.0d0 ) .lt. 0.00001d0) then
         if (     limages) dlogpmax(4) = 0.5d0
         if (.not.limages) dlogpmax(4) = 0.1d0
       else
         dlogpmax(4) = min ( dlogpmax(4), 1.0d0 )
       endif
       if (abs ( dlogpmax(5) - 1.0d0 ) .lt. 0.00001d0) then
         if (     limages) dlogpmax(5) = 0.5d0
         if (.not.limages) dlogpmax(5) = 0.1d0
       else
         dlogpmax(5) = min ( dlogpmax(5), 1.0d0 )
       endif

! Dust temperature.

       npar(1) = 1
       nparmx(1) = 1
       if (modparams1(1) .lt. almostzero .and. modparams2(1) .lt. almostzero) then
         modparam1(1) = 4.0d0
         modparam2(1) = 1000.0d0
       else
         modparam1(1) = modparams1(1)
         modparam2(1) = modparams2(1)
       endif

! Source mass or image surface density.
       
       npar(2) = 1
       nparmx(2) = 1
       if (modparams1(2) .lt. almostzero .and. modparams2(2) .lt. almostzero) then
         if (distancepc .gt. 1.0001d0) then
           punits(2) = ' g'
           modparam1(2) = 3.0d-3 * Msun
           modparam2(2) = 3.0d+6 * Msun
         else
           punits(2) = ' cm^-2'
           modparam1(2) = 1.0d17
           modparam2(2) = 1.0d27
         endif
       else
         if (distancepc .gt. 1.0001d0) then
           punits(2) = ' g'
           modparam1(2) = modparams1(2)
           modparam2(2) = modparams2(2)
         else
           punits(2) = ' cm^-2'
           modparam1(2) = modparams1(2)
           modparam2(2) = modparams2(2)
         endif
       endif

! Opacity slope beta.

       npar(3) = 1
       nparmx(3) = 1
       if (modparams1(3) .lt. almostzero .and. modparams2(3) .lt. almostzero) then
         modparam1(3) = 1.0d0
         modparam2(3) = 3.0d0
       else
         modparam1(3) = modparams1(3)
         modparam2(3) = modparams2(3)
       endif

! Solid angle or distance.

       npar(4) = 1
       nparmx(4) = 1
       if (modparams1(4) .lt. almostzero .and. modparams2(4) .lt. almostzero) then
         if (cmodel .eq. 'modbody') then
           punits(4) = ' ster'
           if (limages) then
             modparam1(4) = 2.0d0 * sqrt ( solidapixel * sq_arcsecs_per_sterad / pi )  !<- equivalent pixel diameter in arcseconds
             modparam2(4) = modparam1(4)
             dlogpmax(4) = 0.0d0
           else
             modparam1(4) = 1.0d-1  !<- diameter in arcseconds
             modparam2(4) = 1.0d+3  !<- diameter in arcseconds
           endif
           modparam1(4) = pi * (modparam1(4) / 2.0d0)**2 / sq_arcsecs_per_sterad  !<- solid angle in steradians
           modparam2(4) = pi * (modparam2(4) / 2.0d0)**2 / sq_arcsecs_per_sterad  !<- solid angle in steradians
         endif
         if (cmodel .eq. 'thinbody') then
           modparam1(4) = distancepc * pc
           modparam2(4) = modparam1(4)
         endif
       else
         if (cmodel .eq. 'modbody') then
           punits(4) = ' ster'
           modparam1(4) = pi * (modparams1(4) / 2.0d0)**2 / sq_arcsecs_per_sterad
           modparam2(4) = pi * (modparams2(4) / 2.0d0)**2 / sq_arcsecs_per_sterad
         endif
         if (cmodel .eq. 'thinbody') then
           modparam1(4) = modparams1(4) * pc
           modparam2(4) = modparams2(4) * pc
           if (abs ( (modparams1(4) - distancepc) * 2.0d0 / (modparams1(4) + distancepc) ) .gt. 0.1d0) then
             if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Incompatible command-line distance values.'
             if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Incompatible command-line distance values.'
             if (iolog .gt. 0) close ( iolog )
             stop 99
           endif
         endif
       endif

! Distance for modbody.

       npar(5) = 1
       nparmx(5) = 1
       if (modparams1(5) .lt. almostzero .and. modparams2(5) .lt. almostzero) then
         if (cmodel .eq. 'modbody') then
           modparam1(5) = distancepc * pc
           modparam2(5) = modparam1(5)
         else
           modparam1(5) = 1.0d0
           modparam2(5) = 1.0d0
         endif
       else
         if (cmodel .eq. 'modbody') then
           modparam1(5) = modparams1(5) * pc
           modparam2(5) = modparams2(5) * pc
           if (abs ( (modparams1(5) - distancepc) * 2.0d0 / (modparams1(5) + distancepc) ) .gt. 0.1d0) then
             if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Incompatible command-line distance values.'
             if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Incompatible command-line distance values.'
             if (iolog .gt. 0) close ( iolog )
             stop 99
           endif
         else
           modparam1(5) = 1.0d0
           modparam2(5) = 1.0d0
         endif
       endif
       
       do i=1,nmodparams
         if (ipfitornot(i) .eq. 1) then
           nparmx(i) = 2**(maxitref)
         endif
       enddo

       do i=1,nmodparams
         if (modparam1(i) .lt. almostzero .and. ipfitornot(i) .eq. 0) then
           if (cmodel .eq. 'modbody' .or. (cmodel .eq. 'thinbody' .and. i .lt. nmodparams)) then
             if (iotty .gt. 0) write (iotty,'(/a,1pe9.2,2(a,i1))')  ' FITFLUXES: ERROR: Incompatible initial values'
     &                                                           , modparam1(i), ' and ', ivarornot(i),' for parameter ', i
             if (iolog .gt. 0) write (iolog,'(/a,1pe9.2,2(a,i1))')  ' FITFLUXES: ERROR: Incompatible initial values'
     &                                                           , modparam1(i), ' and ', ivarornot(i),' for parameter ', i
             if (iolog .gt. 0) close ( iolog )
             stop 99
           endif
         endif
       enddo
       if (iotty .gt. 0 .and. (iverbose .gt. 0 .or. limages)) write (iotty,'()')
       if (iolog .gt. 0) write (iolog,'()')
       do i=1,nmodparams
         if (cmodel .eq. 'modbody' .or. i .lt. nparams) then
           if (iotty .gt. 0 .and. (iverbose .gt. 0 .or. limages)) then
             write (iotty,'(a,i1,a,4(1pe9.2,a))') ' Parameter ', i, ' ('//yesno(ivarornot(i)+1)//'):', modparams1(i), ' ~>'
     &                                          , modparam1(i), ' to', modparam2(i), punits(i), dlogpmax(i)
           endif
           if (iolog .gt. 0) then
             write (iolog,'(a,i1,a,4(1pe9.2,a))') ' Parameter ', i, ' ('//yesno(ivarornot(i)+1)//'):', modparams1(i), ' ~>'
     &                                          , modparam1(i), ' to', modparam2(i), punits(i), dlogpmax(i)
           endif
         endif
       enddo
       if (iotty .gt. 0 .and. (iverbose .gt. 0 .or. limages))
     &                   write (iotty,'(/a,i2)') ' Maximum number of adaptive refinement levels: ......', maxitref
       if (iolog .gt. 0) write (iolog,'(/a,i2)') ' Maximum number of adaptive refinement levels: ......', maxitref
       
       write (cpjx,'(i3)') nparmx(1)
       npjx = int ( log10 ( dble ( nparmx(1) ) ) ) + 1
       if (iotty .gt. 0 .and. (iverbose .gt. 0 .or. limages))
     &                   write (iotty,'(a)') ' Maximum number of zones in adaptive refinement: .... '//cpjx(3-npjx+1:3)
       if (iolog .gt. 0) write (iolog,'(a)') ' Maximum number of zones in adaptive refinement: .... '//cpjx(3-npjx+1:3)
       
       if (iotty .gt. 0 .and. (iverbose .gt. 0 .or. limages))
     &                   write (iotty,'(a,f6.3)') ' Maximum desired value of chi^2/nu in refinement: ...', chisqfmax
       if (iolog .gt. 0) write (iolog,'(a,f6.3)') ' Maximum desired value of chi^2/nu in refinement: ...', chisqfmax
       
       if (iotty .gt. 0 .and. (iverbose .gt. 0 .or. limages))
     &                   write (iotty,'(/a)') ' Names of output files: '//prefix(1:ipr)//'.*'
       if (iolog .gt. 0) write (iolog,'(/a)') ' Names of output files: '//prefix(1:ipr)//'.*'
       
       if (.not.limages) then
         if (iotty .gt. 0 .and. iverbose .gt. 0) 
     &                     write (iotty,'(/a)') ' Reading catalog of sources: '''//excat(1:iexc)//''''
         if (iolog .gt. 0) write (iolog,'(/a)') ' Reading catalog of sources: '''//excat(1:iexc)//''''
       else
         if (iotty .gt. 0 .and. (iverbose .gt. 0 .or. limages))
     &                     write (iotty,'(/a)') ' Reading catalog of images: '''//excat(1:iexc)//''''
         if (iolog .gt. 0) write (iolog,'(/a)') ' Reading catalog of images: '''//excat(1:iexc)//''''
       endif
!!       if (cmodel .eq. 'modbody') then
!!         if (iotty .gt. 0 .and. iverbose .gt. 0)
!!     &                     write (iotty,'(/a)') ' Fitting model ~> '//cmodel//': Fnu = Bnu(T) (1-exp(-tau)) Omega'
!!         if (iolog .gt. 0) write (iolog,'(/a)') ' Fitting model ~> '//cmodel//': Fnu = Bnu(T) (1-exp(-tau)) Omega'
!!       endif
       
! Start processing.

       if (ltiming) then
         cpu = timer ( 'cpu', 0.0d0 )
         wal = timer ( 'wal', 0.0d0 )
       endif
       nsbeg = 1

       ncolnumb = 0
       ncolxcoo = 0
       ncolycoo = 0
       ncolsigg = 0
       ncolgood = 0
       do i=1,nbands
         ncolsig(i) = 0
         ncolfxp(i) = 0
         ncolfxt(i) = 0
         ncolsnr(i) = 0
         ncolfpe(i) = 0
         ncolfte(i) = 0
         fxperro(i) = 0.0d0
         fxterro(i) = 0.0d0
       enddo

       if (.not.limages) then
         ncolnumb = index ( headline(nhcols)(1:lhl), '   NO ' ) - 1
         ncolxcoo = index ( headline(nhcols)(1:lhl), ' XCO_P' )
         ncolycoo = index ( headline(nhcols)(1:lhl), ' YCO_P' )
         ncolsigg = index ( headline(nhcols)(1:lhl), ' SIGN'  )
         ncolgood = index ( headline(nhcols)(1:lhl), ' GOOD'  )
         if (ncolxcoo .gt. 0 .and. ncolxcoo .lt. 5) ncolxcoo = 0        
         if (ncolycoo .gt. 0 .and. ncolycoo .lt. 5) ncolycoo = 0        
               
         do i=1,nbands
           write (cwav,'(i2.2)') nband(i)

           if (ncolsig(i) .eq. 0) ncolsig(i) = index ( headline(nhcols)(1:lhl), ' SIGNM'//cwav )
           if (ncolfxp(i) .eq. 0) ncolfxp(i) = index ( headline(nhcols)(1:lhl), 'FXP_BST'//cwav )
           if (ncolfpe(i) .eq. 0) ncolfpe(i) = index ( headline(nhcols)(1:lhl), 'FXP_ERR'//cwav )
           if (ncolfxt(i) .eq. 0) ncolfxt(i) = index ( headline(nhcols)(1:lhl), 'FXT_BST'//cwav )
           if (ncolfte(i) .eq. 0) ncolfte(i) = index ( headline(nhcols)(1:lhl), 'FXT_ERR'//cwav )

           if (ncolnumb .eq. 0 .or. ncolxcoo .eq. 0 .or. ncolycoo .eq. 0 .or. ncolsigg .eq. 0 .or. 
     &         ncolfxp(i) .eq. 0 .or. ncolfpe(i) .eq. 0 .or. ncolfxt(i) .eq. 0 .or. ncolfte(i) .eq. 0) then
             if (iotty .gt. 0) then
               write (iotty,'(/a)')  ' FITFLUXES: ERROR: Some columns not found for wavelength # '//cwav
               write (iotty,'(9i5)') i, ncolnumb, ncolxcoo, ncolycoo, ncolsigg, ncolfxp(i), ncolfpe(i), ncolfxt(i), ncolfte(i)
               write (iotty,'(a)') headline(nhcols)(1:lhl)
             endif
             if (iolog .gt. 0) then
               write (iolog,'(/a)')  ' FITFLUXES: ERROR: Some columns not found for wavelength # '//cwav
               write (iolog,'(9i5)') i, ncolnumb, ncolxcoo, ncolycoo, ncolsigg, ncolfxp(i), ncolfpe(i), ncolfxt(i), ncolfte(i)
               write (iolog,'(a)') headline(nhcols)(1:lhl)
             endif
             if (iolog .gt. 0) close ( iolog )
             stop 99
           endif
           if (number .eq. 0) then
             nsbeg = 1
             nsend = nextr
           else
             nsbeg = number
             nsend = number
           endif
         enddo
       else
         read (15,*,err=111) wave, beam(1), offset, filename
         goto 112
 111     continue
         if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Trouble reading '''//excat(1:iexc)//''''
         if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Trouble reading '''//excat(1:iexc)//''''
         if (iolog .gt. 0) close ( iolog )
         stop 111
 112     continue
 
         deallocate ( headline )

         lhl = lastc ( filename )

! Check if the input file exists and open it.

         inquire ( file=filename(1:lhl), exist=lfname2 )
         
         if (.not.lfname2) then
           if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: File '''//filename(1:lhl)//''' not found.'
           if (iolog .gt. 0) write (iolog,'(/a)') ' FITFLUXES: ERROR: File '''//filename(1:lhl)//''' not found.'
           if (iolog .gt. 0) close ( iolog )
           stop 2
         endif

! Determine numbers of pixels in the FITS image.
  
         call getfitshead ( filename(1:lhl), nx, ny, dx, dy, bunit )
  
! Allocate the image arrays.
  
         allocate ( fun2d(nx,ny), image(nx,ny,nbands), surfden(nx,ny), tempers(nx,ny), chisqnu(nx,ny), temperr(nx,ny)
     &            , surferr(nx,ny), stat=irc )
  
         if (irc .ne. 0) then
           if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory (8).'
           if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory (8).'
           if (iolog .gt. 0) close ( iolog )
           stop 8
         endif
         backspace ( 15 )
         nbx = 0
         if (iotty .gt. 0) write (iotty,'()')
         if (iolog .gt. 0) write (iolog,'()') 
         
         do nb=1,nbands
           write (cwave,'(i7)') nint ( waveband(nb) )
           ncx = int ( log10 ( dble ( nint ( waveband(nb) ) ) ) ) + 1
           nbx = nbx + 1
           do nbn=nbx,99
             read (15,*,end=114,err=113) wave, beam(nb), offset, filename
             goto 115
 113         continue
             if (iotty .gt. 0) write (iotty,'(a)')  ' FITFLUXES: ERROR: Trouble reading '''//excat(1:iexc)//''''
             if (iolog .gt. 0) write (iolog,'(a)')  ' FITFLUXES: ERROR: Trouble reading '''//excat(1:iexc)//''''
             if (iolog .gt. 0) close ( iolog )
             stop 113
 114         continue
             if (iotty .gt. 0) write (iotty,'(a)')  ' FITFLUXES: ERROR: Waveband "'//cwave(7-ncx+1:7)
     &                                           //'" not found in '''//excat(1:iexc)//''''
             if (iolog .gt. 0) write (iolog,'(a)')  ' FITFLUXES: ERROR: Waveband "'//cwave(7-ncx+1:7)
     &                                           //'" not found in '''//excat(1:iexc)//''''
             if (iolog .gt. 0) close ( iolog )
             stop 114
 115         continue
             if (nint ( waveband(nb) ) .eq. nint ( wave )) exit
           enddo
           lhl = lastc ( filename )
           
! Check if the input file exists and open it.

           inquire ( file=filename(1:lhl), exist=lfname2 )

           if (.not.lfname2) then
             if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: File '''//filename(1:lhl)//''' not found.'
             if (iolog .gt. 0) write (iolog,'(/a)') ' FITFLUXES: ERROR: File '''//filename(1:lhl)//''' not found.'
             if (iolog .gt. 0) close ( iolog )
             stop 3
           endif

! Determine numbers of pixels in the FITS image.
  
           call getfitshead ( filename(1:lhl), nx, ny, dx, dy, bunit )

           write (ccx,'(i7)') nx
           write (ccy,'(i7)') ny
           ncx = int ( log10 ( dble ( nx ) ) ) + 1
           ncy = int ( log10 ( dble ( ny ) ) ) + 1

           if (iotty .gt. 0) write (iotty,'(a)') ' Reading ('//ccx(7-ncx+1:7)//' x '//ccy(7-ncy+1:7)//') '''
     &                                         //filename(1:lhl)//''''
           if (iolog .gt. 0) write (iolog,'(a)') ' Reading ('//ccx(7-ncx+1:7)//' x '//ccy(7-ncy+1:7)//') '''
     &                                         //filename(1:lhl)//''''
           if (nx .ne. nx1 .or. ny .ne. ny1) then
             if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Image size differs between wavebands.'
             if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Image size differs between wavebands.'
             if (iolog .gt. 0) close ( iolog )
             stop 99
           endif
  
           call rfits ( nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, fun2d, dx, dy, object, ra, dec
     &                , filename(1:lhl), iotty, 0, creator, beamx, funmin, funmax, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &                , equinox, bzero, bscale, wave, datamin, datamax, history, iverbose )
  
           if (iverbose .eq. 2) then
             if (iotty .gt. 0) write (iotty,'(a,1x,2(1pe14.7))') ' Minmax values in the image:', funmin, funmax
             if (iolog .gt. 0) write (iolog,'(a,1x,2(1pe14.7))') ' Minmax values in the image:', funmin, funmax
           endif
  
           do j=1,ny
             do i=1,nx
               image(i,j,nb) = fun2d(i,j)
             enddo
           enddo
         enddo
         
         nextr = nx * ny
         if (number .eq. 0 .and. npart .gt. 0 .and. nparts .gt. 0) then
           nsbeg = (npart - 1) * nextr / nparts + 1
           nsend = npart * nextr / nparts
         else
           if (number .eq. 0) then
             nsbeg = 1
             nsend = nextr
           else
             nsbeg = number
             nsend = number
           endif
         endif
       endif

       if (cmodel .eq. 'thinbody') then
         if (.not.limages) then
           if (iotty .gt. 0 .and. iverbose .gt. 0)
     &                       write (iotty,'(/a)')  ' Fitting model ~> '//cmodel//': Fnu = Bnu(T) kappa M D^-2'
           if (iolog .gt. 0) write (iolog,'(/a)')  ' Fitting model ~> '//cmodel//': Fnu = Bnu(T) kappa M D^-2'
         else
           if (iotty .gt. 0) !! .and. iverbose .gt. 0
     &                       write (iotty,'(/a)')  ' Fitting model ~> '//cmodel//': Fnu = Bnu(T) kappa Sigma Omega'
           if (iolog .gt. 0) write (iolog,'(/a)')  ' Fitting model ~> '//cmodel//': Fnu = Bnu(T) kappa Sigma Omega'
         endif
       endif

       allocate ( chisq(nparmx(1),nparmx(2),nparmx(3),nparmx(4),nparmx(5))
     &          , paramspace(nmodparams,max(nparmx(1),nparmx(2),nparmx(3),nparmx(4),nparmx(5))), stat=irc )

       if (irc .ne. 0) then
         if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory (7).'
         if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory (7).'
         if (iolog .gt. 0) close ( iolog )
         stop 7
       endif
       
       allocate ( ltakeit(nextr), x0o(nextr), y0o(nextr), sigglobal(nextr), goodness(nextr), nxmin(nextr,nmodparams), stat=irc )

       if (irc .ne. 0) then
         if (iotty .gt. 0) write (iotty,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory (9).'
         if (iolog .gt. 0) write (iolog,'(/a)')  ' FITFLUXES: ERROR: Trouble allocating memory (9).'
         if (iolog .gt. 0) close ( iolog )
         stop 9
       endif

       do ip=1,nmodparams
         paramspacemin(ip) = 1.0d99
         paramspacemax(ip) = 0.0d0
       enddo
       nplt = 0
       nsrc = 0
       itermean = 0.0d0
       itermin = 99
       itermax = 0
       itrefmean = 0.0d0
       itrefmin = 9
       itrefmax = 0
       nitmean = 0.0d0
       nrefmean = 0.0d0
       isumpct = 0.0d0
!!       numtotal = dble ( nsend - nsbeg + 1 )
       numtotal = dble ( nextr )
       numpart = 0.0d0
       
       if (iverbose .eq. 0 .and. limages) then
         if (iotty .gt. 0) then
           write (iotty,'(11x,a)') '__________________________________________________ 100%'
           write (iotty,'(a)', advance='no') ' Progress: '
           if (iotty .gt. 0) endfile   ( iotty, err=999 )
 999       continue 
           if (iotty .gt. 0) backspace ( iotty )
         endif       
         if (iolog .gt. 0) then
           write (iolog,'(11x,a)') '__________________________________________________ 100%'
           write (iolog,'(a)', advance='no') ' Progress: '
         endif
       endif
       
!!       if (iverbose .eq. 0 .and. .not.limages) then
!!         if (iotty .gt. 0) write (iotty,'()')
!!         if (iolog .gt. 0) write (iolog,'()')
!!       endif

! Prepare arrays and explore fitting model parameter space.

       nst = 0 
       nskip = 0
       
       do ns=1,nextr
       
         nst = nst + 1
       
         if (.not.limages) then
           read (15,'(a)') cline
           do i=1,10
             if (cline(ncolnumb:ncolnumb) .ne. ' ') ncolnumb = ncolnumb - 1
           enddo
           read (cline(ncolnumb:),*) nsrc
         else
           cline = ' '
           nsrc = ns
         endif
         write (cnsrc,'(i10)') nsrc
         nn = int ( log10 ( dble ( nsrc ) ) ) + 1
         write (cnst,'(i10)') nst
         nt = int ( log10 ( dble ( nst ) ) ) + 1
         
!!         write (*,*) nsrc, cnsrc, nn, nst

         if (.not.limages) then
           do i=1,10
             if (cline(ncolxcoo:ncolxcoo) .ne. ' ') ncolxcoo = ncolxcoo - 1
             if (cline(ncolycoo:ncolycoo) .ne. ' ') ncolycoo = ncolycoo - 1
             if (cline(ncolsigg:ncolsigg) .ne. ' ') ncolsigg = ncolsigg - 1
             if (ncolgood .gt. 0 .and. cline(ncolgood:ncolgood) .ne. ' ') ncolgood = ncolgood - 1
           enddo
           read (cline(ncolxcoo:),*) x0o(nst)
           read (cline(ncolycoo:),*) y0o(nst)
           read (cline(ncolsigg:),*) sigglobal(nst)
           if (ncolgood .gt. 0) read (cline(ncolgood:),*) goodness(nst)
         else
           y0o(nst) = dble ( (nst - 1) / nx ) + 1.0d0
           x0o(nst) = dble ( nst ) - dble ( nx ) * (y0o(nst) - 1.0d0)
           sigglobal(nst) = 1.0d30
           goodness(nst) = 1.0d30
           if (iverbose .eq. 0) then
             numpart = dble ( nst )
             call showprogress ( iotty, iolog, lunix, numpart, numtotal, isumpct )
           endif
         endif
         nx0o = nint ( x0o(nst) )
         ny0o = nint ( y0o(nst) )

         ltakeit(ns) = .false.
         
         if (nst .ge. nsbeg .and. nst .le. nsend .and. cline(1:1) .ne. '!' .and. cline(1:1) .ne. '#') then
           if (sigglobal(nst) .lt. sigmin .or. (ncolgood .gt. 0 .and. abs ( goodness(nst) ) .lt. goodmin)) then
             if (.not.limages) then
               if (iotty .gt. 0 .and. iverbose .gt. 0)
     &                           write (iotty,'(/a)') ' Skipping source # '//cnsrc(10-nn+1:10)//': too low SIGGLOB or GOODNESS'
               if (iolog .gt. 0) write (iolog,'(/a)') ' Skipping source # '//cnsrc(10-nn+1:10)//': too low SIGGLOB or GOODNESS'
             endif
             goto 777
           endif

           do i=1,nbands
             sigmono(i) = sigglobal(nst)
             uplimit(i) = 0
             if (.not.limages) then
               if (ncolsig(i) .gt. 0) then
                 read (cline(ncolsig(i):),*) sigmono(i)
               endif
               do m=1,5
                 if (cline(ncolfxp(i):ncolfxp(i)) .ne. ' ' .and. cline(ncolfxp(i):ncolfxp(i)) .ne. '-') ncolfxp(i) = ncolfxp(i) - 1
                 if (cline(ncolfxt(i):ncolfxt(i)) .ne. ' ' .and. cline(ncolfxt(i):ncolfxt(i)) .ne. '-') ncolfxt(i) = ncolfxt(i) - 1
               enddo
               read (cline(ncolfxp(i):),*) fxpbest(i)
               read (cline(ncolfpe(i):),*) fxperro(i)
               read (cline(ncolfxt(i):),*) fxtbest(i)
               read (cline(ncolfte(i):),*) fxterro(i)
               if (index ( cline(ncolfxt(i):ncolfxt(i)+5), ' +' ) .gt. 0) uplimit(i) = 1
             else
               fxtbest(i) = max ( image(nx0o,ny0o,i), 0.0d0 ) * 1.0d6  !<-- Convert intensity from MJy/sr to Jy/sr
               if (cmodel .eq. 'modbody' .or. distancepc .gt. 1.0001d0) then
                 fxtbest(i) = fxtbest(i) * solidapixel                 !<-- Convert from intensity in Jy/sr to flux in Jy
               endif
               fxterro(i) = 0.0d0
             endif
             if (csnband(i) .eq. '-' .and. uplimit(i) .eq. 1) uplimit(i) = -1
             frequency(i) = speedolight * 1.0d4 / wavelength(i)
             fxterror(i) = sqrt ( fxterro(i)**2 + (fxtbest(i) * adderro(i))**2 )
             snrpeak(i) = fxpbest(i) / (fxperro(i) + almostzero)
             snrtotal(i) = fxtbest(i) / (fxterro(i) + almostzero)
           enddo
           
           if (.not.limages .and. iverbose .gt. 0) then
             if (iotty .gt. 0) write (iotty,'(a)') ' ____________________________________________________________________________'
             if (iolog .gt. 0) write (iolog,'(a)') ' ____________________________________________________________________________'
           endif

! Here we can choose to use fewer data points for fitting. Make sure that the number is at least NMODPARAMS.

           nptstofit = 0
           wavepeak = 0.0d0
           fluxpeak = 0.0d0

           if (.not.limages .and. iotty .gt. 0 .and. iverbose .gt. 1) then
             write (iotty,'(/a)') ' Analyzing flux qualities for source # '//cnsrc(10-nn+1:10)
           endif
           if (.not.limages .and. iolog .gt. 0) then
             write (iolog,'(/a)') ' Analyzing flux qualities for source # '//cnsrc(10-nn+1:10)
           endif
           nskp = 0

           do i=1,nbands
             write (cwave,'(i7)') nint ( wavelength(i) )
             ncx = int ( log10 ( dble ( nint ( wavelength(i) ) ) ) ) + 1

             if ((.not.limages .and. uplimit(i) .eq. 0 .and. csnband(i) .ne. '-' .and.
     &            (csnband(i) .eq. '+' .or. (sigmono(i) .ge. sigmin .and. 
     &            snrpeak(i) .ge. snrmin(i) .and. snrtotal(i) .gt. snrmin(i))) .and. 
     &            fxpbest(i) .gt. almostzero .and. fxperro(i) .gt. almostzero .and.
     &            fxtbest(i) .gt. almostzero .and. fxterro(i) .gt. almostzero) .or.
     &           (limages .and. csnband(i) .ne. '-' .and. 
     &            (csnband(i) .eq. '+' .or. snrtotal(i) .gt. 0.01d0) .and.               
     &            fxtbest(i) .gt. almostzero)) then

               nptstofit = nptstofit + 1
               sfrequency(nptstofit) = frequency(i)
               sfxtbest(nptstofit) = fxtbest(i)
               sfxterro(nptstofit) = fxterror(i)

               if (fxtbest(i) .gt. fluxpeak) then
                 fluxpeak = fxtbest(i)
                 wavepeak = wavelength(i)
               endif
!!               if (.not.limages) then
!!                 if (iotty .gt. 0) write (iotty,'(a)') ' Accepting '//cwave(7-ncx+1:7)//' µm fluxes for fitting'
!!                 if (iolog .gt. 0) write (iolog,'(a)') ' Accepting '//cwave(7-ncx+1:7)//' µm fluxes for fitting'
!!               endif
             else
               if (.not.limages .and. csnband(i) .ne. '-') then
                 nskp = nskp + 1
                 if (iotty .gt. 0) then
                   if (nskp .eq. 1 .and. iverbose .gt. 0) write (iotty,'()')
                   if (iverbose .gt. 0)
     &             write (iotty,'(a,i1,12(1pe9.2))') ' ~> Skipping flux at '//cwave(7-ncx+1:7)//' µm: '
     &                                             , uplimit(i), sigmono(i), sigmin, snrpeak(i), snrmin(i), snrtotal(i), snrmin(i)
                 endif
                 if (iolog .gt. 0) then
                   if (nskp .eq. 1) write (iolog,'()')
                   write (iolog,'(a,i1,12(1pe9.2))') ' ~> Skipping flux at '//cwave(7-ncx+1:7)//' µm: '
     &                                             , uplimit(i), sigmono(i), sigmin, snrpeak(i), snrmin(i), snrtotal(i), snrmin(i)
                 endif
               endif
             endif
           enddo
           if (wavepeak .lt. almostzero) then
             if (.not.limages) then
               if (iotty .gt. 0 .and. iverbose .gt. 0)
     &                           write (iotty,'(/a)') ' ~> Skipping source # '//cnsrc(10-nn+1:10)//': acceptable fluxes not found'
               if (iolog .gt. 0) write (iolog,'(/a)') ' ~> Skipping source # '//cnsrc(10-nn+1:10)//': acceptable fluxes not found'
             endif
             nskip = nskip + 1
             goto 777
           endif
           if (nptstofit .lt. nfitparams) then
             if (.not.limages) then
               if (iotty .gt. 0 .and. iverbose .gt. 0) then
                 write (iotty,'(/a,i1,a)') ' ~> Skipping source # '//cnsrc(10-nn+1:10)//': too few (', nptstofit, ') fluxes to fit'
               endif
               if (iolog .gt. 0) then
                 write (iolog,'(/a,i1,a)') ' ~> Skipping source # '//cnsrc(10-nn+1:10)//': too few (', nptstofit, ') fluxes to fit'
               endif
             endif
             nskip = nskip + 1
             goto 777
           endif
           if (limages .and. iverbose .gt. 0) then
             if (iotty .gt. 0) write (iotty,'(a)') ' ____________________________________________________________________________'
             if (iolog .gt. 0) write (iolog,'(a)') ' ____________________________________________________________________________'
           endif
           lam1 = nint ( 1.0d4 * speedolight / sfrequency(1) )
           lam2 = nint ( 1.0d4 * speedolight / sfrequency(nptstofit) )
           write (clam1,'(i4)') lam1
           write (clam2,'(i4)') lam2
           nl1 = int ( log10 ( dble ( lam1 ) ) ) + 1
           nl2 = int ( log10 ( dble ( lam2 ) ) ) + 1
           write (cx,'(i4)') nx0o
           write (cy,'(i4)') ny0o
           ncx = int ( log10 ( dble ( max ( nx0o, 1 ) ) ) ) + 1
           ncy = int ( log10 ( dble ( max ( ny0o, 1 ) ) ) ) + 1
           write (cnsend,'(i10)') nsend
           ne = int ( log10 ( dble ( nsend ) ) ) + 1
           npercent = nint ( dble ( nst ) / dble ( nsend ) * 100.0d0 )
           write (cnperc,'(i3)') npercent
           
           if (iverbose .ne. 0 .or. .not.limages) then
             if (iotty .gt. 0) then
               if (iverbose .ne. 0) then
                 write (iotty,'()')
                 if (.not.limages) then
                   write (iotty,'(a)') ' Fitting '//clam1(4-nl1+1:4)//'~'//clam2(4-nl2+1:4)//' µm fluxes: source # '
     &                               //cnst(10-nt+1:10)//' of '//cnsend(10-ne+1:10)//': '//cnperc//'%'
                 else
                   write (iotty,'(a)') ' Fitting '//clam1(4-nl1+1:4)//'~'//clam2(4-nl2+1:4)//' µm fluxes: pixel ('
     &                               //cx(4-ncx+1:4)//','//cy(4-ncy+1:4)//') # '//cnst(10-nt+1:10)
     &                               //' of '//cnsend(10-ne+1:10)//': '//cnperc//'%'
                 endif
               endif
             endif
             if (iolog .gt. 0) then
               write (iolog,'()')
               if (.not.limages) then
                 write (iolog,'(a)') ' Fitting '//clam1(4-nl1+1:4)//'~'//clam2(4-nl2+1:4)//' µm fluxes: source # '
     &                             //cnst(10-nt+1:10)//' of '//cnsend(10-ne+1:10)//': '//cnperc//'%'
               else
                 write (iolog,'(a)') ' Fitting '//clam1(4-nl1+1:4)//'~'//clam2(4-nl2+1:4)//' µm fluxes: pixel ('
     &                             //cx(4-ncx+1:4)//','//cy(4-ncy+1:4)//') # '//cnst(10-nt+1:10)
     &                             //' of '//cnsend(10-ne+1:10)//': '//cnperc//'%'
               endif
             endif
           endif

           do ip=1,nmodparams
             parmin(ip) = modparam1(ip)
             parmax(ip) = modparam2(ip)
           enddo
           ltakeit(nst) = .true.
           chisqfreeo = 1.0d0
           lconvergedo = .false.
           ldone = .false.

           itref = 0
 1111      continue
           lnochange = .true.
           
           do ip=1,nmodparams
             nparo = npar(ip)
             dlogpo = dlogp(ip)
             if (ipfitornot(ip) .eq. 1) then
               npar(ip) = 2**(itref)
             endif
             if (npar(ip) .gt. nparmx(ip)) then
               if (iotty .gt. 0) then
                 write (iotty,'(/a,3i4)') ' FITFLUXES: ERROR: Too many zones in parameter space:', itref, npar(ip), nparmx(ip)
               endif
               if (iolog .gt. 0) then
                 write (iolog,'(/a,3i4)') ' FITFLUXES: ERROR: Too many zones in parameter space:', itref, npar(ip), nparmx(ip)
               endif
               if (iolog .gt. 0) close ( iolog )
               stop 99
             endif
             modp1 = log ( parmin(ip) )
             modp2 = log ( parmax(ip) )
             pspace1o = paramspace(ip,1)
             pspace2o = paramspace(ip,nparo)
             dlogp(ip) = (modp2 - modp1) / dble ( npar(ip) )
             if (dlogp(ip) .gt. almostzero) then
               dlogp(ip) = max ( dlogp(ip), dlogpmax(ip) )
               npar(ip) = nint ( (modp2 - modp1) / dlogp(ip) )
               dlogp(ip) = (modp2 - modp1) / dble ( npar(ip) )
             endif
             modp1 = modp1 + dlogp(ip) / 2.0d0
             do j=1,npar(ip)
               paramspace(ip,j) = exp ( modp1 + dble ( j - 1 ) * dlogp(ip) )
             enddo
             if (dlogp(ip) .gt. almostzero .and. (npar(ip) .ne. nparo .or. 
     &           abs ( dlogp(ip) - dlogpo ) .gt. almostzero .or. abs ( paramspace(ip,1) - pspace1o ) .gt. almostzero .or.
     &           abs ( paramspace(ip,npar(ip)) - pspace2o ) .gt. almostzero)) then
               lnochange = .false.
             endif
           enddo

! Exit refinement loop if the above discretization of parameter space did not change since the last iteration.

           if (itref .gt. 0 .and. lnochange) then
             itref = itref - 1
             goto 3333
           endif
           
           write (cpjx,'(i3)') npar(1)
           npjx = int ( log10 ( dble ( npar(1) ) ) ) + 1
           write (cpkx,'(i3)') npar(2)
           npkx = int ( log10 ( dble ( npar(2) ) ) ) + 1
           write (cplx,'(i3)') npar(3)
           nplx = int ( log10 ( dble ( npar(3) ) ) ) + 1
           write (cpmx,'(i3)') npar(4)
           npmx = int ( log10 ( dble ( npar(4) ) ) ) + 1
           write (cpnx,'(i3)') npar(5)
           npnx = int ( log10 ( dble ( npar(5) ) ) ) + 1
           numfits = npar(1) * npar(2) * npar(3) * npar(4) * npar(5)
           write (cfts,'(i9)') numfits
           nfts = int ( log10 ( dble ( numfits ) ) ) + 1

           if (iverbose .eq. 2) then
             if (iotty .gt. 0) then
               if (cmodel .eq. 'modbody') then
                 write (iotty,'(/a)') ' Exploring '//cmodel//' parameter space ('//cpjx(3-npjx+1:3)//','//cpkx(3-npkx+1:3)//','
     &                              //cplx(3-nplx+1:3)//','//cpmx(3-npmx+1:3)//','//cpnx(3-npnx+1:3)//') in '
     &                              //cfts(9-nfts+1:9)//' fits'
               else
                 write (iotty,'(/a)') ' Exploring '//cmodel//' parameter space ('//cpjx(3-npjx+1:3)//','//cpkx(3-npkx+1:3)//','
     &                              //cplx(3-nplx+1:3)//','//cpmx(3-npmx+1:3)//') in '
     &                              //cfts(9-nfts+1:9)//' fits'
               endif
             endif
             if (iolog .gt. 0) then
               if (cmodel .eq. 'modbody') then
                 write (iolog,'(/a)') ' Exploring '//cmodel//' parameter space ('//cpjx(3-npjx+1:3)//','//cpkx(3-npkx+1:3)//','
     &                              //cplx(3-nplx+1:3)//','//cpmx(3-npmx+1:3)//','//cpnx(3-npnx+1:3)//') in '
     &                              //cfts(9-nfts+1:9)//' fits'
               else
                 write (iolog,'(/a)') ' Exploring '//cmodel//' parameter space ('//cpjx(3-npjx+1:3)//','//cpkx(3-npkx+1:3)//','
     &                              //cplx(3-nplx+1:3)//','//cpmx(3-npmx+1:3)//') in '
     &                              //cfts(9-nfts+1:9)//' fits'
               endif
             endif
             if (iotty .gt. 0) write (iotty,'()')
             if (iolog .gt. 0) write (iolog,'()')
           endif
           
           do i=1,nmodparams
             write (cpx,'(i3)') npar(i)
             npx = int ( log10 ( dble ( npar(i) ) ) ) + 1
             if (iverbose .eq. 2) then
               if (cmodel .eq. 'modbody' .or. i .lt. nparams) then
                 if (iotty .gt. 0) then
                   write (iotty,'(a,i1,a,4(1pe9.2,a))') ' Parameter ', i, ' ('//yesno(ivarornot(i)+1)//'):', modparams1(i), ' ~>'
     &                                                , paramspace(i,1), ' to', paramspace(i,npar(i)), punits(i), dlogp(i)
                 endif
                 if (iolog .gt. 0) then
                   write (iolog,'(a,i1,a,4(1pe9.2,a))') ' Parameter ', i, ' ('//yesno(ivarornot(i)+1)//'):', modparams1(i), ' ~>'
     &                                                , paramspace(i,1), ' to', paramspace(i,npar(i)), punits(i), dlogp(i)
                 endif
               endif
             endif
             if (modparam1(i) .lt. almostzero .and. ipfitornot(i) .eq. 0) then
               if (iotty .gt. 0) write (iotty,'(/a,1pe9.2,2(a,i1))')  ' FITFLUXES: ERROR: Incompatible initial values'
     &                                                             , modparam1(i), ' and ', ivarornot(i),' for parameter ', i
               if (iolog .gt. 0) write (iolog,'(/a,1pe9.2,2(a,i1))')  ' FITFLUXES: ERROR: Incompatible initial values'
     &                                                             , modparam1(i), ' and ', ivarornot(i),' for parameter ', i
               if (iolog .gt. 0) close ( iolog )
               stop 99
             endif
           enddo

           nxmin(nst,1) = 1
           nxmin(nst,2) = 1
           nxmin(nst,3) = 1
           nxmin(nst,4) = 1
           nxmin(nst,5) = 1
           iterdone = 0
           chisqmin = 1.0d55
           chisqfree = 1.0d55
           lgoodfound = .false.
           lconverged = .false.

           do j=1,npar(1)
             do k=1,npar(2)
               do l=1,npar(3)
                 do m=1,npar(4)
                   do n=1,npar(5)
                   fitparameters(1) = paramspace(1,j)
                   fitparameters(2) = paramspace(2,k)
                   fitparameters(3) = paramspace(3,l)
                   fitparameters(4) = paramspace(4,m)
                   fitparameters(5) = paramspace(5,n)
                   do i=1,nbands
                     opacity(i) = opacity0 * (frequency(i) / frequency0)**fitparameters(3)
                   enddo
                   
                   if (cmodel .eq. 'modbody') then
                     call fitsed ( iotty, iolog, nsrc, nbands, nmodparams, nfitparams, nptstofit, wavelength, frequency, fxtbest
     &                           , fxterror, sfrequency, sfxtbest, sfxterro, sedslope, distancepc, derivedmass, fitparameters
     &                           , ipfitornot, dyda, covar, alpha, chisq(j,k,l,m,n), iter, lgood, lconv, lsing, limages, cmodel
     &                           , modbody, 'parspace', quantities, distancrelerr, uplimit, lplotting, solidapixel, npts, nfun
     &                           , xn, fun, yerrp, yerrm, ldebug )
                   endif
                   if (cmodel .eq. 'thinbody') then
                     call fitsed ( iotty, iolog, nsrc, nbands, nmodparams, nfitparams, nptstofit, wavelength, frequency, fxtbest
     &                           , fxterror, sfrequency, sfxtbest, sfxterro, sedslope, distancepc, derivedmass, fitparameters
     &                           , ipfitornot, dyda, covar, alpha, chisq(j,k,l,m,n), iter, lgood, lconv, lsing, limages, cmodel
     &                           , thinbody, 'parspace', quantities, distancrelerr, uplimit, lplotting, solidapixel, npts, nfun
     &                           , xn, fun, yerrp, yerrm, ldebug )
                   endif

! Discard the source fit if it goes above an upper limit value.

                   if (.not.limages) then
                     do i=1,nbands
                       if (uplimit(i) .eq. 1) then
                         if (cmodel .eq. 'modbody') then
                           call modbody ( frequency(i), fitparameters, modelflux, dyda, nmodparams )
                         endif
                         if (cmodel .eq. 'thinbody') then
                           call thinbody ( frequency(i), fitparameters, modelflux, dyda, nmodparams )
                         endif
                         if (modelflux .gt. fxtbest(i)) lgood = .false.
                       endif
                     enddo
                   endif

! Identify the minimum values of the parameter space.

                   if (chisq(j,k,l,m,n) .lt. chisqmin .and. .not.lsing) then
                     if (lgood) lgoodfound = .true.
                     if (lconv) lconverged = .true.
                     chisqmin = chisq(j,k,l,m,n)
                     iterdone = iter
                     nxmin(nst,1) = j
                     nxmin(nst,2) = k
                     nxmin(nst,3) = l
                     nxmin(nst,4) = m
                     nxmin(nst,5) = n
                     nfreedom = nptstofit - nfitparams
                     chisqfree = chisqmin / (dble ( nfreedom ) + 1.0d0 )
                   endif
                   if (lconv) then
                     itermin = min ( itermin, iter )
                     itermax = max ( itermax, iter )
                     itermean = itermean + dble ( iter )
                     nitmean = nitmean + 1.0d0
                   endif
                   enddo
                 enddo
               enddo
             enddo
           enddo
           npjx = 1
           write (cpjx,'(i3)') iterdone
           if (iterdone .gt. 0) npjx = int ( log10 ( dble ( iterdone ) ) ) + 1
           if (iverbose .eq. 2 .and. iotty .gt. 0) then
             write (iotty,'(/a,1pe21.14,a,3l)') 
     &                    ' Minimum chi^2 value', chisqmin, ' found in '//cpjx(3-npjx+1:3)//' iterations', lgood, lconv, lsing
           endif
           if (iverbose .ne. 0) then
             write (iolog,'(/a,1pe21.14,a,3l)')
     &                    ' Minimum chi^2 value', chisqmin, ' found in '//cpjx(3-npjx+1:3)//' iterations', lgood, lconv, lsing
           endif

           do ip=1,nmodparams
             parmn(ip) = 1.0d99
             parmx(ip) = 0.0d0
           enddo
           do j=1,npar(1)
             do k=1,npar(2)
               do l=1,npar(3)
                 do m=1,npar(4)
                   do n=1,npar(5)
                   if (chisq(j,k,l,m,n) .le. 2.0d0 * chisqmin) then
                     parmn(1) = min ( parmn(1), paramspace(1,j) )
                     parmx(1) = max ( parmx(1), paramspace(1,j) )
                     parmn(2) = min ( parmn(2), paramspace(2,k) )
                     parmx(2) = max ( parmx(2), paramspace(2,k) )
                     parmn(3) = min ( parmn(3), paramspace(3,l) )
                     parmx(3) = max ( parmx(3), paramspace(3,l) )
                     parmn(4) = min ( parmn(4), paramspace(4,m) )
                     parmx(4) = max ( parmx(4), paramspace(4,m) )
                     parmn(5) = min ( parmn(5), paramspace(5,n) )
                     parmx(5) = max ( parmx(5), paramspace(5,n) )
                   endif
                   enddo
                 enddo
               enddo
             enddo
           enddo
           do ip=1,nmodparams
             parmxlog = log ( parmx(ip) )
             parmnlog = log ( parmn(ip) )
             delta = max ( (parmxlog - parmnlog) / 10.0d0, 2.0d0 * dlogp(ip) )
             if (parmn(ip) .lt. 1.0d99) parmin(ip) = max ( exp ( parmnlog - delta ), modparam1(ip) )
             if (parmx(ip) .gt. 0.0d0 ) parmax(ip) = min ( exp ( parmxlog + delta ), modparam2(ip) )
           enddo

           if (iverbose .eq. 2 .and. iotty .gt. 0) then
             write (iotty,'()')
             do ip=1,nmodparams
               if (cmodel .eq. 'modbody' .or. ip .lt. nparams) then
                 write (iotty,'(a,i1,a,2(1pe9.2,a))') ' Minimum chi^2 domain in parameter ', ip, ':'
     &                                              , parmin(ip), ' to', parmax(ip), punits(ip)
               endif
             enddo
           endif
           if (iverbose .ne. 0 .and. iolog .gt. 0) then
             write (iolog,'()')
             do ip=1,nmodparams
               if (cmodel .eq. 'modbody' .or. ip .lt. nparams) then
                 write (iolog,'(a,i1,a,2(1pe9.2,a))') ' Minimum chi^2 domain in parameter ', ip, ':'
     &                                              , parmin(ip), ' to', parmax(ip), punits(ip)
               endif
             enddo
           endif

! Check if we need to continue parameter space refinement iterations.

           lconvgood = lconverged .and. lgoodfound
           lnotconvo = .not.lconverged .and. .not.lconvergedo
           ldone = .false.
           if (itref .gt. 0) then
             if ((lconvgood  .and. chisqfree .lt. chisqfmax) .or.
     &           (lconverged .and. abs ( chisqfree - chisqfreeo ) .lt. 1.0d-3 * chisqfreeo) .or.
     &           (lnotconvo  .and. abs ( chisqfree - chisqfreeo ) .gt. 1.0d-3 * chisqfreeo .and. 
     &            chisqfree .gt. chisqfreeo / 2.0d0 .and. itref .ge. 3)) then
!!!     &           (lnotconvo  .and. chisqfree .gt. chisqfreeo / 2.0d0 .and. itref .ge. 2)) then
               ldone = .true.
             endif
           endif
           if (itref .lt. maxitref .and. .not.ldone) then 
             chisqfreeo = chisqfree
             lconvergedo = lconverged
             itref = itref + 1
             goto 1111
           endif
 3333      continue
           itrefmin = min ( itrefmin, itref )
           itrefmax = max ( itrefmax, itref )
           itrefmean = itrefmean + dble ( itref )
           nrefmean = nrefmean + 1.0d0

! Find minimum and maximum values of the optimal initial parameters for all best fits.

           if (ldone .and. chisqfree .le. chisqfmax) then
             do ip=1,nmodparams
               paramspacemin(ip) = min ( paramspacemin(ip), parmin(ip) )
               paramspacemax(ip) = max ( paramspacemax(ip), parmax(ip) )
             enddo
           endif

! Save chi^2 values in the model parameter space.

           if (ldebug) then
             if (cmodel .eq. 'modbody') then
               write (55,'(/a)') '# CHI^2 VALUES IN A modbody PARAMETER SPACE ('//cpjx(3-npjx+1:3)//','//cpkx(3-npkx+1:3)//','
     &                         //cplx(3-nplx+1:3)//','//cpmx(3-npmx+1:3)//') FOR SOURCE # '//cnsrc(10-nn+1:10)
             endif
             if (cmodel .eq. 'thinbody') then
               write (55,'(/a)') '# CHI^2 VALUES IN THE THINBODY PARAMETER SPACE ('//cpjx(2-npjx+1:2)//','//cpkx(2-npkx+1:2)//','
     &                         //cplx(2-nplx+1:2)//','//cpmx(2-npmx+1:2)//') FOR SOURCE # '//cnsrc(10-nn+1:10)
             endif
             write (cjx,'(i3)') nxmin(nst,1)
             njx = int ( log10 ( dble ( nxmin(nst,1) ) ) ) + 1
             write (ckx,'(i3)') nxmin(nst,2)
             nkx = int ( log10 ( dble ( nxmin(nst,2) ) ) ) + 1
             write (clx,'(i3)') nxmin(nst,3)
             nlx = int ( log10 ( dble ( nxmin(nst,3) ) ) ) + 1
             write (cmx,'(i3)') nxmin(nst,4)
             nmx = int ( log10 ( dble ( nxmin(nst,4) ) ) ) + 1
             write (cnx,'(i3)') nxmin(nst,5)
             nnx = int ( log10 ( dble ( nxmin(nst,5) ) ) ) + 1
             write (55,'(/a,1pe10.3)') '# BEST FIT FOUND AT POINT ('//cjx(3-njx+1:3)//','//ckx(3-nkx+1:3)//','//clx(3-nlx+1:3)
     &                               //','//cmx(3-nmx+1:3)//'): CHI^2 =', chisq(nxmin(nst,1),nxmin(nst,2),nxmin(nst,3),nxmin(nst,4)
     &                               ,nxmin(nst,5))
            
             write (55,'(/a,1pe10.3)') '# Slice through the 4th dimension for source # '//cnsrc(10-nn+1:10)//' at point # '
     &                               //cmx(3-nmx+1:3)//', value:', paramspace(4,nxmin(nst,4))
             do l=1,npar(3)
               write (cl,'(i3)') l
               nl = int ( log10 ( dble ( l ) ) ) + 1
               write (55,'(/a,1pe10.3)') '# Slice through the 3rd dimension for source # '//cnsrc(10-nn+1:10)//' at point # '
     &                                 //cl(3-nl+1:3)//', value:', paramspace(3,l)
               write (55,'(/15x,1000(2x,i3,4x))') (k, k=1,npar(2))
               write (55,'( 14x,1000(1pe9.2))') (paramspace(2,k), k=1,npar(2))
               write (55,'()')
               do j=1,npar(1)
                 write (55,'(i3,1x,1pe9.2,1x,1000(e9.2))') j, paramspace(1,j), (chisq(j,k,nxmin(nst,3),nxmin(nst,4),nxmin(nst,5))
     &                                                   , k=1,npar(2))
               enddo
             enddo
           endif
           do ip=1,nmodparams
             fitparameters(ip) = paramspace(ip,nxmin(nst,ip))
           enddo
           do i=1,nbands
             opacity(i) = opacity0 * (frequency(i) / frequency0)**fitparameters(3)
           enddo

           if (cmodel .eq. 'modbody') then
             call fitsed ( iotty, iolog, nsrc, nbands, nmodparams, nfitparams, nptstofit, wavelength, frequency, fxtbest, fxterror
     &                   , sfrequency, sfxtbest, sfxterro, sedslope, distancepc, derivedmass, fitparameters, ipfitornot, dyda
     &                   , covar, alpha, chisqx, iter, lgood, lconv, lsing, limages, cmodel, modbody, 'finalfit', quantities
     &                   , distancrelerr, uplimit, lplotting, solidapixel, npts, nfun, xn, fun, yerrp, yerrm, ldebug )
           endif
           if (cmodel .eq. 'thinbody') then
             call fitsed ( iotty, iolog, nsrc, nbands, nmodparams, nfitparams, nptstofit, wavelength, frequency, fxtbest, fxterror
     &                   , sfrequency, sfxtbest, sfxterro, sedslope, distancepc, derivedmass, fitparameters, ipfitornot, dyda
     &                   , covar, alpha, chisqx, iter, lgood, lconv, lsing, limages, cmodel, thinbody, 'finalfit', quantities
     &                   , distancrelerr, uplimit, lplotting, solidapixel, npts, nfun, xn, fun, yerrp, yerrm, ldebug )

           endif
           
!ps           call plotsed ( iotty, iolog, psunit, nsrc, nplt, npts, nfun, nx0o, ny0o, nmodparams, ctime, cdate
!ps     &                  , postscript(1:ilp), fitparameters, chisqx, limages, cmodel, quantities, lplotting
!ps     &                  , nfitparams, nptstofit, xn, fun, yerrp, yerrm )

           if (.not.limages) then
             okbad = dble ( quantities(17) ) / dble ( quantities(18) + 1.0 )
             if (okbad .le. 1.0d0) then
               cokbad = 'ok'
             else
               cokbad = 'bad'
             endif

! Define random noise at each pixel using a random number generator 'gasdev' with a normal (Gaussian) distribution.
! Re-initialize the random sequence first (idum = -1).

             idum = -1
             sigma = 0.0d0

! Inner loop is to make a new independent random sequence for every wavelength, to avoid correlation between the quasi-random
! noise at any wavelength. The initialization consists of multiplying the integer number equal to the wavelength by a uniform 
! random number in the range [0,1]. So many cycles will be skept out of the Gaussian random sequence between the rows of pixels.

             newran = int ( 10 * nsrc * ran1 ( idum ) )
             do i=1,newran
               q1 = gasdev ( idum )
             enddo
             do i=1,ntries
               q1 = gasdev ( idum )
               simtruemass(i) = quantities(3) + q1 * quantities(5) !<- (5) is the total mass error, not just from SED fitting.
               sigma = sigma + q1**2
             enddo
             sigma = sqrt ( sigma / real ( ntries ) )
!!             if (iverbose .eq. 2) write (iotty,'(''  Standard deviation: '',$)')
!!             if (iverbose .eq. 2) write (iotty,*) nsrc, sigma
!!             if (iverbose .eq. 2) write (iotty,*) (randmasserr(i),i=1,100)

             write (77,'(i7,2f8.1,1x,2(1pe9.2),1x,2(1pe10.3),1x,3(1pe10.3),6(1x,2(1pe10.3)),1pe10.3,3x,a,0pf9.3,1x,1000(1pe11.3))')
     &                  nsrc, x0o(nst), y0o(nst), sigglobal(nst), goodness(nst), (quantities(ip), ip=1,18), cokbad, sigma
     &                , (simtruemass(i),i=1,ntries)
           endif

! chisq / dble ( nu + 1 ) .le. 1.0d0
     
           if (limages) then
             chisqnu(nx0o,ny0o) = quantities(17) / (quantities(18) + 1.0d0)
             if (chisqnu(nx0o,ny0o) .le. 10.0d0) then
               tempers(nx0o,ny0o) = quantities(1)  !<-- temperatures
               surfden(nx0o,ny0o) = quantities(10)  !<-- surface densities
               if (chisqnu(nx0o,ny0o) .le. 1.0d0) then
                 temperr(nx0o,ny0o) = quantities(2)  !<-- error in tempers
                 surferr(nx0o,ny0o) = quantities(11) !<-- error in surfden
               else
                 temperr(nx0o,ny0o) = 0.0d0  !<-- error in tempers
                 surferr(nx0o,ny0o) = 0.0d0  !<-- error in surfden
               endif
             else
               tempers(nx0o,ny0o) = 0.0d0  !<-- temperatures
               surfden(nx0o,ny0o) = 0.0d0  !<-- surface densities
               temperr(nx0o,ny0o) = 0.0d0  !<-- error in tempers
               surferr(nx0o,ny0o) = 0.0d0  !<-- error in surfden
             endif
           endif
           goto 888
         endif
 777     continue
         if (limages) then
           chisqnu(nx0o,ny0o) = 0.0d0  !<-- chi^2 values
           tempers(nx0o,ny0o) = 0.0d0  !<-- temperatures
           surfden(nx0o,ny0o) = 0.0d0  !<-- surface densities
           temperr(nx0o,ny0o) = 0.0d0  !<-- error in tempers
           surferr(nx0o,ny0o) = 0.0d0  !<-- error in surfden
         endif
 888     continue
       enddo

!ps! Convert plots from PostScript into PDF using macos 'ps2pdf'.
!ps
!ps       if (lplotting .and. nextr .gt. 0 .and. nextr .gt. nskip) then
!ps          call plotnd ( psunit )
!ps          if (iotty .gt. 0 .and. iverbose .gt. 0)
!ps     &       write (iotty,'(a)') ' ____________________________________________________________________________'
!ps          if (iolog .gt. 0)
!ps     &       write (iolog,'(a)') ' ____________________________________________________________________________'
!ps          if (iotty .gt. 0 .and. iverbose .gt. 0)
!ps     &       write (iotty,'(/a)') ' Converting '''//postscript(1:ilp)//''' to PDF with ''ps2pdf'''
!ps          if (iolog .gt. 0) write (iolog,'(/a)') ' Converting '''//postscript(1:ilp)//''' to PDF with ''ps2pdf'''
!ps          call ps2pdf ( postscript(1:ilp) )
!ps       endif

       close ( 15 )
       if (ldebug) close ( 55 )
       if (.not.limages) then
          if (ldebug) close ( 66 )
         close ( 77 )
       endif
       
       if (nitmean .gt. 0) itermean = itermean / nitmean
       if (nrefmean .gt. 0) itrefmean = itrefmean / nrefmean

       if (ldebug .and. lplotting .and. npar(1) .ge. 2 .and. npar(2) .ge. 2) then
         
         allocate ( space(npar(1),npar(2)), stat=irc )
        
         if (irc .ne. 0) then
           if (iotty .gt. 0) write (iotty,'(/a)') ' FITFLUXES: ERROR: Trouble allocating memory (10).'
           if (iolog .gt. 0) write (iolog,'(/a)') ' FITFLUXES: ERROR: Trouble allocating memory (10).'
           if (iolog .gt. 0) close ( iolog )
           stop 10
         endif

         bunit = ' '
         creator = 'FITFLUXES'
         ctype1 = 'Temper'
         if (cmodel .eq. 'modbody' .or. limages) then
           ctype2 = 'SurfDens'
         else
           ctype2 = 'Mass'
         endif
         object = 'None'
         history = ' '
         blank = 0
         datamax = 0.0d0
         datamin = 0.0d0
         cdelt1 = -dlogp(1) * 3600.0d0
         cdelt2 = dlogp(2) * 3600.0d0
         crpix1 = dble ( npar(1) ) / 2.0d0
         crpix2 = dble ( npar(2) ) / 2.0d0
         crota1 = 0.0d0
         crota2 = 0.0d0
         cd11 = -dlogp(1) * 3600.0d0
         cd12 = 0.0d0
         cd21 = 0.0d0
         cd22 = dlogp(2) * 3600.0d0
         crval1 = (log ( paramspace(1,1) ) + log ( paramspace(1,npar(1)) )) / 2.0d0 - dlogp(1) / 2
         crval2 = (log ( paramspace(2,1) ) + log ( paramspace(2,npar(2)) )) / 2.0d0 - dlogp(2) / 2
         equinox = 0.0d0
         bzero = 0.0d0
         bscale = 1.0d0
         wave = 0.0d0
         beamx = 0.0d0
         
         do j=1,npar(1)
           do k=1,npar(2)
             space(j,k) = chisq(j,k,nxmin(nsbeg,3),nxmin(nsbeg,4),nxmin(nsbeg,5))
             if (space(j,k) .gt. 10.0d0) read(cnan,*) space(j,k)
           enddo
         enddo
        
         if (iotty .gt. 0) write (iotty,'(/a)') ' Writing output: '''//paraspace(1:ile)//''''
         if (iolog .gt. 0) write (iolog,'(/a)') ' Writing output: '''//paraspace(1:ile)//''''
        
         call wfits ( cfitsversion, npar(1), npar(2), bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, space, cdelt1
     &              , cdelt2, object, crval1, crval2, paraspace(1:ile), cdate, ctime, creator, beamx, blank, crota1, crota2
     &              , cd11, cd12, cd21, cd22, equinox, bzero, bscale, wave, datamin, datamax, history )
        
         deallocate ( space )
       endif

       if (iotty .gt. 0 .and. iverbose .eq. 2) 
     &    write (iotty,'(a)') ' ____________________________________________________________________________'
       if (iolog .gt. 0)
     &    write (iolog,'(a)') ' ____________________________________________________________________________'
       npjx = 1
       npkx = 1
       nplx = 1
       write (cpjx,'(i3)') itermin
       if (itermin .gt. 0) npjx = int ( log10 ( dble ( itermin ) ) ) + 1
       write (cpkx,'(i3)') nint ( itermean )
       if (nint ( itermean ) .gt. 0) npkx = int ( log10 ( dble ( nint ( itermean ) ) ) ) + 1
       write (cplx,'(i3)') itermax
       if (itermax .gt. 0) nplx = int ( log10 ( dble ( itermax ) ) ) + 1
       if (iotty .gt. 0 .and. iverbose .eq. 2) then
         write (iotty,'(/a)')    ' Min avg max numbers of iterations in fitting: '
     &                      //cpjx(3-npjx+1:3)//' '//cpkx(3-npkx+1:3)//' '//cplx(3-nplx+1:3)
         write (iotty,'(a,3i2)') ' Min avg max numbers of refinement iterations:', itrefmin, nint ( itrefmean ), itrefmax
       endif
       if (iolog .gt. 0) then
         write (iolog,'(/a)')    ' Min avg max numbers of iterations in fitting: '
     &                      //cpjx(3-npjx+1:3)//' '//cpkx(3-npkx+1:3)//' '//cplx(3-nplx+1:3)
         write (iolog,'(a,3i2)') ' Min avg max numbers of refinement iterations:', itrefmin, nint ( itrefmean ), itrefmax
       endif
       
       deallocate ( ltakeit, x0o, y0o, sigglobal, goodness, nxmin )

       deallocate ( chisq, paramspace )

       deallocate ( ncolsig, ncolfxp, ncolfxt, ncolsnr, ncolfpe, ncolfte, snrpeak, snrtotal, frequency, fxpbest, fxperro, fxtbest
     &            , fxterro, opacity, sedslope, wavelength, sfrequency, sfxtbest, sfxterro, sigmono, beam, fxterror  !!, fxperror
     &            , uplimit )
       
       deallocate ( fitparameters, ipfitornot, dyda, covar, alpha, paramspacemin, paramspacemax )
       
       if (limages) then
         if (iotty .gt. 0) write (iotty,'(/a)') ' Writing output: '''//temperatur(1:ilt)//''''
         if (iolog .gt. 0) write (iolog,'(/a)') ' Writing output: '''//temperatur(1:ilt)//''''
        
         call wfits ( cfitsversion, nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, tempers, dx, dy, object, crval1
     &              , crval2, temperatur(1:ilt), cdate, ctime, creator, beamx, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history )
        
         if (iotty .gt. 0) write (iotty,'(a)') ' Writing output: '''//surfdensit(1:ilc)//''''
         if (iolog .gt. 0) write (iolog,'(a)') ' Writing output: '''//surfdensit(1:ilc)//''''
        
         call wfits ( cfitsversion, nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, surfden, dx, dy, object, crval1
     &              , crval2, surfdensit(1:ilc), cdate, ctime, creator, beamx, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history )

         if (iotty .gt. 0) write (iotty,'(a)') ' Writing output: '''//rsurferror(1:ilm)//''''
         if (iolog .gt. 0) write (iolog,'(a)') ' Writing output: '''//rsurferror(1:ilm)//''''
        
         call wfits ( cfitsversion, nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, surferr, dx, dy, object, crval1
     &              , crval2, rsurferror(1:ilm), cdate, ctime, creator, beamx, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history )

         if (iotty .gt. 0) write (iotty,'(a)') ' Writing output: '''//rtemperror(1:ilr)//''''
         if (iolog .gt. 0) write (iolog,'(a)') ' Writing output: '''//rtemperror(1:ilr)//''''
        
         call wfits ( cfitsversion, nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, temperr, dx, dy, object, crval1
     &              , crval2, rtemperror(1:ilr), cdate, ctime, creator, beamx, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history )
        
         if (iotty .gt. 0) write (iotty,'(a)') ' Writing output: '''//chisquarnu(1:ilq)//''''
         if (iolog .gt. 0) write (iolog,'(a)') ' Writing output: '''//chisquarnu(1:ilq)//''''

         do j=1,ny
           do i=1,nx
             if (chisqnu(i,j) .gt. 10.0d0) read(cnan,*) chisqnu(i,j)
           enddo
         enddo
        
         call wfits ( cfitsversion, nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, chisqnu, dx, dy, object, crval1
     &              , crval2, chisquarnu(1:ilq), cdate, ctime, creator, beamx, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history )
    
         deallocate ( fun2d, image, surfden, tempers, chisqnu, temperr, surferr )
       endif

       if (iotty .gt. 0 .and. iverbose .gt. 0) 
     &     write (iotty,'(a)') ' ____________________________________________________________________________'
       if (iolog .gt. 0)
     &     write (iolog,'(a)') ' ____________________________________________________________________________'

       if (ltiming) then
         cpu_sum = cpu_sum + timer ( 'cpu', cpu )
         wal_sum = wal_sum + timer ( 'wal', wal )
       endif

! Total cpu usage, expressed in hours, minutes, and seconds.

       cputot = cpu_sum
       cpuhours = int ( cputot / 3600.0d0 )
       cpumins  = int ( cputot / 60.0d0 ) -  60 * cpuhours
       cpusecs  = int ( cputot + 0.5d0 ) - 3600 * cpuhours - 60 * cpumins

! Final wall clock time, expressed in hours, minutes, and seconds.

       wctot = wal_sum
       wchours = int ( wctot / 3600.0d0 )
       wcmins  = int ( wctot / 60.0d0 ) -  60 * wchours
       wcsecs  = int ( wctot + 0.5d0 ) - 3600 * wchours - 60 * wcmins

       call when ( lunix, ctime, cdate, ndate, 4 )

       write (ccpu,'(1pe9.2)') cpu_sum
       write (cwal,'(1pe9.2)') wal_sum
       write (cccpu,'(i3.2,a,i2.2,a,i2.2)') cpuhours,':',cpumins,':',cpusecs
       write (ccwal,'(i3.2,a,i2.2,a,i2.2)') wchours,':',wcmins,':',wcsecs
       ndc = 1
       ndw = 1
       if (cccpu(1:1) .eq. ' ') ndc = 2
       if (ccwal(1:1) .eq. ' ') ndw = 2

       if (iotty .gt. 0 .and. iverbose .eq. 2) write (iotty,'(/a)') ' DONE. '//cdate//' '//ctime//' CPU:'//ccpu//' s '//cccpu(ndc:)
     &                                                            //' WALL:'//cwal//' s '//ccwal(ndw:)
       if (iolog .gt. 0) write (iolog,'(/a)') ' DONE. '//cdate//' '//ctime//' CPU:'//ccpu//' s '//cccpu(ndw:)
     &                                      //' WALL:'//cwal//' s '//ccwal(ndw:)
       if (iolog .gt. 0) close ( iolog )

!!       stop !<-- commented out because it would lead to run-time messages about denormalized values, when using gfortran.
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine fitsed 
       
     &            ( iotty, iolog, nsrc, nbands, nmodparams, nfitparams, nptstofit, wavelength, frequency, fluxvalue, fluxerror
     &            , sfrequency, sfluxvalue, sfluxerror, sedslope, distancepc, derivedmass, fitparameters, ipfitornot, dyda, covar
     &            , alpha, chisq, iter, lgood, lconv, lsing, limages, cmodel, modelfun, cwhat, quantities
     &            , distancrelerr, uplimit, lplotting, solidapixel, npts, nfun, xn, fun, yerrp, yerrm, ldebug )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       logical       lgood, lconv, lsing, limages, lplotting, ldebug

       character*7   cnu
       character*10  cnsrc
       character*(*) cmodel, cwhat

       integer       npts, nfun, iotty, iolog
       integer       nsrc, nmodparams, nbands, nfitparams, i, j, nn, lastc, nm, nptstofit, npf, iter, nrj, nu, iverbose
       integer       ibe(nfun), ien(nfun), uplimit(nbands), ipfitornot(nmodparams)

       real          funmx, funmn, xnmn, xnmx, xnlmn, xnlmx

       real          fun(npts,nfun), xn(npts,nfun), xerrp(npts,nfun), xerrm(npts,nfun), yerrp(npts,nfun), yerrm(npts,nfun)
     &             , modslope(nfun)                                             
                                                                                                     
       real*8        distance2, pc, fq, derivedmass, derivedmassx, speedolight, fitflux, opacity0, angulardiam, hh, k_B, Msun
     &             , modintens1, modintens2, chisq, wave4mass, distancepc, dust2gas, freq4mass, scalingfactor, fitflux4mass
     &             , rayleighjeans, frequency0, wavelength0, bbody, almostzero, fluxpeak, fluxmin, dcdenrelerr
     &             , sq_arcsecs_per_sterad, fun2max, fun3max, freqmeanhi, freqmeanlo, Lsun, pi, pi4, sedluminosity, fitluminosity
     &             , sedbolflux, fitbolflux, amu, muH2, tauwave4mass, wavedepthone, dlog10freq, fitlumerror, bbodyZmass
     &             , dmasserror, sedbolfluxerr, fitbolfrelerr, fitfluxrelerr, opacit0relerr, dustgasrelerr
     &             , distancrelerr, dmassrelerr, opacdepthone, fitflux4massup, derivedcden, dcdenerror, solidapixel

       real*8        frequency(nbands), wavelength(nbands), fluxvalue(nbands), fluxerror(nbands), sedslope(nbands), fitparamsup(5)
     &             , sfrequency(nptstofit), sfluxvalue(nptstofit), sfluxerror(nptstofit), waves(npts), sedfit(npts)
     &             , fitparameters(nmodparams), dyda(nmodparams), covar(nmodparams,nmodparams), alpha(nmodparams,nmodparams)
     &             , quantities(*)

       parameter   ( pc = 3.085678d18, speedolight = 2.99792458d10, k_B = 1.380658d-16, hh = 6.626d-27, Msun = 1.9891d33
     &             , Lsun = 3.826D+33, sq_arcsecs_per_sterad = 3282.80635d0 * 3600.0d0**2, pi = 3.14159265358979d0
     &             , pi4 = pi * 4.0d0, almostzero = 1.0d-20, amu = 1.6605402D-24, muH2 = 2.8d0, nrj = 50 )

       common / copacity / frequency0, wavelength0, opacity0, dust2gas, opacit0relerr, dustgasrelerr, iverbose
  
       external      modelfun, fitfun, planck, lastc
!__________________________________________________________________________________________________________________________________
!
       write (cnsrc,'(i10)') nsrc
       nn = int ( log10 ( dble ( nsrc ) ) ) + 1
       
! Number of the degrees of freedom NU: for reliable fits, CHI^2 is similar to NU.

       nu = nptstofit - nfitparams
       write (cnu,'(i7)') nu
       if (nu .gt. 0) then
         nm = int ( log10 ( dble ( nu ) ) ) + 1
       else
         nm = 1
       endif

       distance2 = (distancepc * pc)**2

       call fitfun ( iotty, iolog, nmodparams, nptstofit, sfrequency, sfluxvalue, sfluxerror, fitparameters, ipfitornot, covar
     &             , alpha, chisq, iter, lgood, lconv, lsing, cmodel, modelfun, cwhat )
     
! Integrate rough bolometric flux under the input fluxes.

       sedbolflux = 0.0d0
       freqmeanhi = sfrequency(1) + (sfrequency(1) - 0.5d0 * (sfrequency(1) + sfrequency(2)))
       do i=1,nptstofit-1
         freqmeanlo = 0.5d0 * (sfrequency(i) + sfrequency(i+1))
         sedbolflux = sedbolflux + 1.0d-23 * sfluxvalue(i) * (freqmeanhi - freqmeanlo)
         freqmeanhi = freqmeanlo
       enddo
       
       sedluminosity = sedbolflux * pi4 * distance2 / Lsun

! Integrate bolometric flux under the fitted modbody fluxes.

       xnmn = real ( min ( wavelength(1), 30.0d0 ) )
       xnmx = real ( max ( wavelength(nptstofit), 1000.0d0 ) )
       xnlmn = log10 ( xnmn )
       xnlmx = log10 ( xnmx )

       do i=1,npts
         waves(i) = 10.0d0**( xnlmn + (i - 1) * (xnlmx - xnlmn) / (npts - 1) )

         call modelfun ( 1.0d4 * speedolight / waves(i), fitparameters, sedfit(i), dyda, nmodparams )

       enddo

       fitbolflux = 0.0d0
       freqmeanhi = 1.0d4 * speedolight * (1.0d0 / waves(1) + (1.0d0 / waves(1) - 2.0d0 / (waves(1) + waves(2))))
       do i=1,npts-1
         freqmeanlo = 1.0d4 * speedolight * 2.0d0 / (waves(i) + waves(i+1))
         fitbolflux = fitbolflux + 1.0d-23 * sedfit(i) * (freqmeanhi - freqmeanlo)
         freqmeanhi = freqmeanlo
       enddo

! Error propagation for sums and differences for random and independent quantities: 
! total uncertainty of sums and differences is the quadratic sum of the individual uncertainties.
! Error propagation for products and quotients for random and independent quantities: 
! total uncertainty of products and quotients is the quadratic sum of the individual *relative* uncertainties.

       sedbolfluxerr = 0.0d0
       freqmeanhi = sfrequency(1) + (sfrequency(1) - 0.5d0 * (sfrequency(1) + sfrequency(2)))
       do i=1,nptstofit-1
         freqmeanlo = 0.5d0 * (sfrequency(i) + sfrequency(i+1))
         sedbolfluxerr = sedbolfluxerr + (1.0d-23 * sfluxerror(i) * (freqmeanhi - freqmeanlo))**2
         freqmeanhi = freqmeanlo
       enddo
       fitbolfrelerr = sqrt ( sedbolfluxerr ) / sedbolflux

       fitparamsup(1) = fitparameters(1) + sqrt ( abs ( covar(1,1) ) )
       fitparamsup(2) = fitparameters(2) + sqrt ( abs ( covar(2,2) ) )
       fitparamsup(3) = fitparameters(3) + sqrt ( abs ( covar(3,3) ) )
       fitparamsup(4) = fitparameters(4) + sqrt ( abs ( covar(4,4) ) )
       if (nmodparams .eq. 5) fitparamsup(5) = fitparameters(5) + sqrt ( abs ( covar(5,5) ) )

       call planck ( fitparameters(1), 1.0d4 / wavelength(nptstofit), bbodyZmass )
       
       wave4mass = 100000.0d0 !<-- 100 millimeters
       freq4mass = 1.0d4 * speedolight / wave4mass

       call modelfun ( freq4mass, fitparameters, fitflux4mass, dyda, nmodparams )
       call modelfun ( freq4mass, fitparamsup, fitflux4massup, dyda, nmodparams )

       fitfluxrelerr = 0.0d0
       if (fitflux4mass .gt. almostzero) then
         fitfluxrelerr = (fitflux4massup - fitflux4mass) / fitflux4mass
       endif
       fitluminosity = fitbolflux * pi4 * distance2 / Lsun
       fitlumerror = fitluminosity * sqrt ( fitbolfrelerr**2 + 2.0d0 * distancrelerr**2 )

       dmasserror = 0.0d0
       derivedcden = 0.0d0
       dcdenerror = 0.0d0
       angulardiam = 0.0d0
       dmassrelerr = 0.0d0

! Derive quantities at fitted temperature.

       if (cmodel .eq. 'modbody') then

! Standard formula for mass derivation: M=Fnu*D^2/(kappanu*Bnu(T)). It can be converted to the alternative formula (see below) 
! by simply using the definition: Fnu=Bnu(T)*Omega, assuming a blackbody source with uniform brightness at temperature T.
! Simpler formula for mass derivation: M=Omega*D^2*Sigma (mass = area * column_number_density * muH2 * mH).

         if (distancepc .gt. 1.0001d0) then
           derivedmass = fitparameters(2) / Msun
           dmassrelerr = sqrt ( abs ( covar(2,2) )) / fitparameters(2)
           dmasserror = derivedmass * sqrt ( dmassrelerr**2 + 2.0d0 * distancrelerr**2 + opacit0relerr**2 + dustgasrelerr**2 )

           derivedcden = fitparameters(2) / (amu * muH2) / fitparameters(4) / distance2
           dcdenerror = derivedcden * sqrt ( dmassrelerr**2 + 2.0d0 * distancrelerr**2 + opacit0relerr**2 + dustgasrelerr**2 )
         else
           derivedcden = fitparameters(2)
           dcdenrelerr = sqrt ( abs ( covar(2,2) )) / fitparameters(2)
           dcdenerror = derivedcden * sqrt ( dcdenrelerr**2 + opacit0relerr**2 + dustgasrelerr**2 )

           derivedmass = fitparameters(2) * distance2 * solidapixel * (amu * muH2) / Msun
           dmasserror = derivedmass * sqrt ( dcdenrelerr**2 + 2.0d0 * distancrelerr**2 + opacit0relerr**2 + dustgasrelerr**2 )
         endif
     
! Another mass formula (Pezzutto et al. 2012). From tau=(frq/frq1)^beta or tau=(lam1/lam)^beta; kappa=kappa0*(lam0/lam)^beta;
! Omega=pi*R^2/D^2; tau=kappa0*(lam0/lam)^beta*M/(pi*R^2), one can get: M=Omega*D^2*/kappa1, where kappa1=kappa0*(lam0/lam1)^beta
! is the opacity of gas+dust at wavelength where the optical depth is unity.

         tauwave4mass = derivedcden * (amu * muH2) * dust2gas * opacity0 * (wavelength0 / wave4mass)**fitparameters(3)
         wavedepthone = (wave4mass**fitparameters(3) * tauwave4mass)**(1.0d0 / fitparameters(3))
         opacdepthone = opacity0 * (wavelength0 / wavedepthone)**fitparameters(3)
         derivedmassx = fitparameters(4) * distance2 / dust2gas / opacdepthone / Msun
         angulardiam = 2.0d0 * sqrt ( fitparameters(4) * sq_arcsecs_per_sterad / pi )
       endif
       if (cmodel .eq. 'thinbody') then
         if (distancepc .gt. 1.0001d0) then
           derivedmass = fitparameters(2) / Msun
           dmassrelerr = sqrt ( abs ( covar(2,2) )) / fitparameters(2)
           dmasserror = derivedmass * sqrt ( dmassrelerr**2 + 2.0d0 * distancrelerr**2 + opacit0relerr**2 + dustgasrelerr**2 )
           if (limages) then
             derivedcden = fitparameters(2) / (amu * muH2) / solidapixel / distance2
             dcdenerror = derivedcden * sqrt ( dmassrelerr**2 + 2.0d0 * distancrelerr**2 + opacit0relerr**2 + dustgasrelerr**2 )
           else
             derivedcden = 0.0d0
             dcdenerror = 0.0d0
           endif
         else
           derivedcden = fitparameters(2)
           dcdenrelerr = sqrt ( abs ( covar(2,2) )) / fitparameters(2)
           dcdenerror = derivedcden * sqrt ( dcdenrelerr**2 + opacit0relerr**2 + dustgasrelerr**2 )
           derivedmass = 0.0d0
           dmasserror = 0.0d0
         endif
         tauwave4mass = derivedcden * (amu * muH2) * dust2gas * opacity0 * (wavelength0 / wave4mass)**fitparameters(3)
         wavedepthone = (wave4mass**fitparameters(3) * tauwave4mass)**(1.0d0 / fitparameters(3))
         angulardiam = 2.0d0 * sqrt ( solidapixel * sq_arcsecs_per_sterad / pi )
       endif

! Collect derived parameters.

       quantities( 1) = fitparameters(1)
       quantities( 2) = sqrt ( abs ( covar(1,1) ) )
       quantities( 3) = derivedmass
       quantities( 4) = derivedmass * dmassrelerr
       quantities( 5) = dmasserror
       quantities( 6) = fitluminosity
       quantities( 7) = fitlumerror
       quantities(10) = derivedcden
       quantities(11) = dcdenerror
       quantities(12) = fitparameters(3)
       quantities(13) = sqrt ( abs ( covar(3,3) ) )
         
       if (cmodel .eq. 'modbody') then
         quantities( 8) = angulardiam
         quantities( 9) = 2.0d0 * sqrt ( sqrt ( abs ( covar(4,4) ) ) * sq_arcsecs_per_sterad / pi )
         quantities(14) = fitparameters(5) / pc
         quantities(15) = sqrt ( abs ( covar(5,5) ) ) / pc
       endif
       if (cmodel .eq. 'thinbody') then
         quantities( 8) = angulardiam
         quantities( 9) = 0.0d0
         quantities(14) = fitparameters(4) / pc
         quantities(15) = sqrt ( abs ( covar(4,4) ) ) / pc
       endif

       quantities(16) = wavedepthone
       quantities(17) = chisq
       quantities(18) = dble ( nu )
       
! Show modbody slopes between the PACS and SPIRE wavelengths.

       if (cwhat .eq. 'finalfit' .and. iotty .gt. 0 .and. iverbose .ne. 0) then
         if (iverbose .eq. 2) then

! Compute SED slopes between the PACS and SPIRE wavelengths and slopes of the fitted SED.

           call modelfun ( sfrequency(1), fitparameters, modintens1, dyda, nmodparams )
    
           do i=1,nptstofit-1
    
             call modelfun ( sfrequency(i+1), fitparameters, modintens2, dyda, nmodparams )
    
             dlog10freq = -(log10 ( sfrequency(i+1) ) - log10 ( sfrequency(i) ))
             sedslope(i) = real ( (log10 ( sfluxvalue(i+1) ) - log10 ( sfluxvalue(i) )) / dlog10freq )
             modslope(i) = real ( (log10 ( modintens2 ) - log10 ( modintens1 )) / dlog10freq )
             modintens1 = modintens2
           enddo

           write (iotty,'(/a,f7.3,10f8.3)' ) '   FITSED: Original SED slopes:', (sedslope(i), i=1,nptstofit-1)
           write (iotty,'( a,f7.3,10f8.3)' ) '   FITSED:  Modbody fit slopes:', (modslope(i), i=1,nptstofit-1)
           write (iotty,'()')
           if (nmodparams .eq. 5 .and. wavedepthone .gt. almostzero) then  !!! .and. distancepc .lt. 1.0001d0   
             write (iotty,'(a,2(1pe9.2,a))') '   FITSED: Optical depth is 1.00 at', wavedepthone, ' µm'
     &                                     , 1.0d4 * speedolight / wavedepthone / 1.0d9, ' GHz'
           endif
           write (iotty,'(a,f6.2,a,2(1pe9.2,a))') '   FITSED:  Dust opacity is', opacity0,' at', wavelength0, ' µm'
     &                                          , frequency0 / 1.0d9, ' GHz'
         endif
         write (iotty,'(/a,1pe10.3,a,2(e9.2,a))') '   FITSED:     Temperature:', quantities(1), ' K     ±'
     &                                          , quantities(2), ' K     ±', quantities(2) / (quantities(1) + almostzero)
         write (iotty,'(a,1pe10.3,2(a,e9.2))'   ) '   FITSED:      Total mass:', quantities(3), ' Msun  ±'
     &                                          , quantities(5), ' Msun  ±', quantities(5) / (quantities(3) + almostzero)
         write (iotty,'(a,1pe10.3,2(a,e9.2))'   ) '   FITSED: Surface density:', quantities(10), ' cm^-2 ±'
     &                                          , quantities(11), ' cm^-2 ±', quantities(11) / (quantities(10) + almostzero)
         write (iotty,'(a,1pe10.3,2(a,e9.2))'   ) '   FITSED:   Opacity slope:', quantities(12), '       ±'
     &                                          , quantities(13), '       ±', quantities(13) / (quantities(12) + almostzero)
         if (cmodel .eq. 'modbody') then
         write (iotty,'(a,1pe10.3,a,2(e9.2,a))' ) '   FITSED:     Solid angle:', fitparameters(4), ' ster  ±'
     &                                          , sqrt ( abs ( covar(4,4) ) ), ' ster  ±'
     &                                          , sqrt ( abs ( covar(4,4) ) ) / (fitparameters(4) + almostzero)
         endif
         write (iotty,'(a,1pe10.3,2(a,e9.2))') '   FITSED:        Distance:', quantities(14), ' pc    ±'
     &                                       , quantities(15), ' pc    ±', quantities(15) / (quantities(14) + almostzero)
         if (iverbose .eq. 2) then
           if (cmodel .eq. 'modbody') then
             write (iotty,'(/a,3(1pe10.3,a))'  )'   FITSED: Derived luminosity:', quantities(6), ' Lsun ±', quantities(7)
     &                                         , ' Lsun ±', quantities(7) / (quantities(6) + almostzero)
             write (iotty,'(a,2(1pe10.3,a))'   ) '   FITSED: Approx. luminosity:', sedluminosity, ' Lsun (from observed fluxes)'
             write (iotty,'(a,2(1pe10.3,a))'   ) '   FITSED:   Angular diameter:', quantities(8), ' arcsec'
           endif
           if (cmodel .eq. 'thinbody' .and. iverbose .ne. 0) then    !!! .and. distancepc .gt. 1.0001d0
             write (iotty,'(/a,3(1pe10.3,a))'  ) '   FITSED: Derived luminosity:', quantities(6), ' Lsun ±', quantities(7)
     &                                         , ' Lsun ±', quantities(7) / (quantities(6) + almostzero)
             write (iotty,'(a,2(1pe10.3,a))'   ) '   FITSED: Approx. luminosity:', sedluminosity, ' Lsun (from observed fluxes)'
             write (iotty,'(a,2(1pe10.3,a))'   ) '   FITSED:   Angular diameter:', quantities(8), ' arcsec'
           endif
         endif
       endif
       
       if (cwhat .eq. 'finalfit' .and. iolog .gt. 0 .and. iverbose .ne. 0) then
         if (iverbose .eq. 2) then
           write (iolog,'(/a,f7.3,10f8.3)' ) '   FITSED: Original SED slopes:', (sedslope(i), i=1,nptstofit-1)
           write (iolog,'( a,f7.3,10f8.3)' ) '   FITSED:  Modbody fit slopes:', (modslope(i), i=1,nptstofit-1)
           write (iolog,'()')
           if (nmodparams .eq. 5 .and. wavedepthone .gt. almostzero) then        !!! .and. distancepc .lt. 1.0001d0
             write (iolog,'(a,2(1pe9.2,a))') '   FITSED: Optical depth is 1.00 at', wavedepthone, ' µm'
     &                                      , 1.0d4 * speedolight / wavedepthone / 1.0d9, ' GHz'
           endif
           write (iolog,'(a,f6.2,a,2(1pe9.2,a))') '   FITSED:  Dust opacity is', opacity0,' at', wavelength0, ' µm'
     &                                          , frequency0 / 1.0d9, ' GHz'
         endif
         write (iolog,'(/a,1pe10.3,a,2(e9.2,a))') '   FITSED:     Temperature:', quantities(1), ' K     ±'
     &                                          , quantities(2), ' K     ±', quantities(2) / (quantities(1) + almostzero)
         write (iolog,'(a,1pe10.3,2(a,e9.2))'   ) '   FITSED:      Total mass:', quantities(3), ' Msun  ±'
     &                                          , quantities(5), ' Msun  ±', quantities(5) / (quantities(3) + almostzero)
         write (iolog,'(a,1pe10.3,2(a,e9.2))'   ) '   FITSED: Surface density:', quantities(10), ' cm^-2 ±'
     &                                          , quantities(11), ' cm^-2 ±', quantities(11) / (quantities(10) + almostzero)
         write (iolog,'(a,1pe10.3,2(a,e9.2))'   ) '   FITSED:   Opacity slope:', quantities(12), '       ±'
     &                                          , quantities(13), '       ±', quantities(13) / (quantities(12) + almostzero)
         if (cmodel .eq. 'modbody') then
         write (iolog,'(a,1pe10.3,a,2(e9.2,a))' ) '   FITSED:     Solid angle:', fitparameters(4), ' ster  ±'
     &                                          , sqrt ( abs ( covar(4,4) ) ), ' ster  ±'
     &                                          , sqrt ( abs ( covar(4,4) ) ) / (fitparameters(4) + almostzero)
         endif
         write (iolog,'(a,1pe10.3,2(a,e9.2))') '   FITSED:        Distance:', quantities(14), ' pc    ±'
     &                                       , quantities(15), ' pc    ±', quantities(15) / (quantities(14) + almostzero)
         if (cmodel .eq. 'modbody') then
           write (iolog,'(/a,3(1pe10.3,a))'  )'   FITSED: Derived luminosity:', quantities(6), ' Lsun ±', quantities(7)
     &                                       , ' Lsun ±', quantities(7) / (quantities(6) + almostzero)
           write (iolog,'(a,2(1pe10.3,a))'   ) '   FITSED: Approx. luminosity:', sedluminosity, ' Lsun (from observed fluxes)'
           write (iolog,'(a,2(1pe10.3,a))'   ) '   FITSED:   Angular diameter:', quantities(8), ' arcsec'
         endif
         if (cmodel .eq. 'thinbody' .and. iverbose .ne. 0) then        !!! .and. distancepc .gt. 1.0001d0
           write (iolog,'(/a,3(1pe10.3,a))'  ) '   FITSED: Derived luminosity:', quantities(6), ' Lsun ±', quantities(7)
     &                                       , ' Lsun ±', quantities(7) / (quantities(6) + almostzero)
           write (iolog,'(a,2(1pe10.3,a))'   ) '   FITSED: Approx. luminosity:', sedluminosity, ' Lsun (from observed fluxes)'
           write (iolog,'(a,2(1pe10.3,a))'   ) '   FITSED:   Angular diameter:', quantities(8), ' arcsec'
         endif
       endif
       
! Plot fitted SEDs.

       if (cwhat .eq. 'finalfit' .and. lplotting) then
        
         do j=1,nfun
           ibe(j) = 1
           do i=1,npts
             xerrp(i,j) = 1.0e-31
             xerrm(i,j) = 1.0e-31
             yerrp(i,j) = 1.0e-31
             yerrm(i,j) = 1.0e-31
           enddo
         enddo

! Wavelengths.

         do i=1,nbands
           xn(i,1) = real ( 1.0e4 * speedolight / frequency(i) )
           fun(i,1) = real ( max ( fluxvalue(i), 0.000001d0 ) )
           if (uplimit(i) .eq. 0) then
             yerrp(i,1) = real ( fluxerror(i) )
           else
             yerrp(i,1) = 0.0d0
           endif
           yerrm(i,1) = real ( min ( fluxerror(i), 0.999999 * fun(i,1) ) )
         enddo
         ien(1) = nbands

         do i=1,nptstofit
           xn(i,2) = real ( 1.0e4 * speedolight / sfrequency(i) )
           fun(i,2) = real ( max ( sfluxvalue(i), 0.000001d0 ) )
           yerrp(i,2) = real ( max ( sfluxerror(i), 1.0e-31 ) )
           yerrm(i,2) = real ( max ( min ( sfluxerror(i), 0.999999 * fun(i,2) ), 1.0e-31 ) )
         enddo
         ien(2) = nptstofit

! If there are upper limits, make them the same color as the fitted points.

         npf = nptstofit
         do i=1,nbands
           if (uplimit(i) .eq. 1) then
             npf = npf + 1
             xn(npf,2) = xn(i,1)
             fun(npf,2) = real ( max ( fluxvalue(i), 0.000001d0 ) )
             yerrp(npf,2) = 1.0e-31
             yerrm(npf,2) = real ( max ( min ( fluxerror(i), 0.999999 * fun(i,1) ), 1.0e-31 ) )
           endif
         enddo
         ien(2) = npf

! Modbody fitted to the input fluxes.

         do i=1,npts
           xn(i,3) = 10.0**( xnlmn + (i - 1) * (xnlmx - xnlmn) / (npts - 1) )
           call modelfun ( 1.0d4 * speedolight / dble ( xn(i,3) ), fitparameters, fitflux, dyda, nmodparams )
           fun(i,3) = real ( fitflux )
         enddo
         ien(3) = npts

! Blackbody normalized to the maximum of the modbody fit.

         fun2max = -1.0e+30
         fun3max = -1.0e+30
         
         do i=1,npts
           xn(i,4) = 10.0**( xnlmn + (i - 1) * (xnlmx - xnlmn) / (npts - 1) )
           call planck ( fitparameters(1), 1.0d4 / dble ( xn(i,4) ), bbody )
           fun(i,4) = real ( bbody )
           fun2max = max ( fun(i,3), fun2max )
           fun3max = max ( fun(i,4), fun3max )
         enddo
         scalingfactor = fun2max / fun3max
         do i=1,npts
           call planck ( fitparameters(1), 1.0d4 / dble ( xn(i,4) ), bbody )
           fun(i,4) = real ( bbody * scalingfactor )
         enddo
         ien(4) = npts
        
! Rayleigh-Jeans curve, normalized to blackbody at 100 millimeters.

         fq = 1.0d4 * speedolight / 1.0d5
         call planck ( fitparameters(1), 1.0d4 / 1.0d5, bbody )
         rayleighjeans = 2.0 * fq**2 * k_B * fitparameters(1) / speedolight**2
         scalingfactor = scalingfactor * (bbody / rayleighjeans)
        
         do i=1,nrj
           xn(i,5) = 500.0 + (i - 1) * (xnmx - 500.0) / (nrj - 1)
           fq = 1.0d4 * speedolight / xn(i,5)
           rayleighjeans = 2.0 * fq**2 * k_B * fitparameters(1) / speedolight**2
           fun(i,5) = real ( rayleighjeans * scalingfactor )
         enddo
         ien(5) = nrj

! Limit displayed values of fluxes.

         fluxpeak = 0.0d0
         fluxmin = 1.0d30
         do i=1,nbands
           if (fluxvalue(i) .gt. fluxpeak) then
             fluxpeak = fluxvalue(i)
           endif
           if (fluxvalue(i) .le. fluxmin .and. fluxvalue(i) .gt. 0.000001) then
             fluxmin = fluxvalue(i)
           endif
         enddo
         funmx = -1.0e+30
         funmn =  1.0e+30
         do j=1,nfun
           do i=1,ien(j)
             funmx = max ( fun(i,j), funmx )
             if (fun(i,j) .gt. 0.000001) funmn = min ( fun(i,j), funmn )
           enddo
         enddo
         funmn = real ( max ( 0.5 * fluxmin, funmn ) )  !!!, 0.0001
         funmx = real ( min ( 1.5 * fluxpeak, 1.1 * funmx ) )

! Limit function values to avoid problems when creating PostScript plots.

         do j=1,nfun
           do i=1,npts
             if (i .le. ien(j)) then
               fun(i,j) = max ( fun(i,j), funmn / 100.0 )
               fun(i,j) = min ( fun(i,j), funmx * 10.0 )
             else
               xn(i,j) = xn(ien(j),j) * 1000.0
               fun(i,j) = fun(ien(j),j) / 100.0
             endif
           enddo
         enddo
           
         if (.not.limages .and. ldebug) then
           write (66,'(a)') '#___________________________________________________________________________________________'
           write (66,'(a)') '#'
           write (66,'(a)') '# DATA FOR SOURCE # '//cnsrc(10-nn+1:10)
           write (66,'(a)') '#'
           write (66,'(a)') '#  1      2          3          4          5          6          7          8          9'
     &                    //'          10         11         12         13         14         15'
           if (cmodel .eq. 'modbody') then
             write (66,'(a)') '# NO    WAVEO     OBSPOINT    WAVEOF   OBSPOINTF    WAVEG     MODBODY     WAVEB    BLACKBODY'
     &                      //'    WAVER    RAY-JEANS   OBSERRP    OBSERRM    OBSERRPF   OBSERRMF'
           endif
           if (cmodel .eq. 'thinbody') then
             write (66,'(a)') '# NO    WAVEO     OBSPOINT    WAVEOF   OBSPOINTF    WAVEM     THINBODY    WAVEB    BLACKBODY'
     &                      //'    WAVERJ   RAY-JEANS   OBSERRP    OBSERRM    OBSERRPF   OBSERRMF'
           endif
           write (66,'(a)') '#'
          
           do i=1,npts
             write (66,'(i4,100(1pe11.3))') i, (xn(i,j), fun(i,j), j=1,nfun), (yerrp(i,j), yerrm(i,j), j=1,2)
           enddo
         endif
       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!ps       subroutine plotsed 
!ps       
!ps     &            ( iotty, iolog, psunit, nsrc, nplt, npts, nfun, nx0o, ny0o, nmodparams, ctime, cdate, postscript
!ps     &            , fitparameters, chisq, limages, cmodel, quantities, lplotting, nfitparams, nptstofit, xn, fun, yerrp, yerrm)
!ps!________________________________________________________________________________________________________________________________
!ps!
!ps!________________________________________________________________________________________________________________________________
!ps!
!ps       implicit      none
!ps       logical       limages, lplotting, lportrait
!ps       integer       npts, nfun, iotty, iolog
!ps       character*4   cx, cy
!ps       character*5   xhlp, yhlp
!ps       character*6   clam1, clam2, ctdf
!ps       character*7   cnu
!ps       character*8   cchi
!ps       character*9   cso, clo
!ps       character*10  cnsrc
!ps       character*80  symbol(nfun), color(nfun)
!ps       character*128 xlab(2), ylab(2), phdr, varnam(nfun), pftr(4)
!ps       character*(*) ctime, cdate, cmodel, postscript
!ps       integer       lastc, nmodparams, i, j, nu, nn, nl1, nl2, nm, nt, nx0o, ny0o, ncx, ncy, allxl, allyl, iverbose
!ps     &             , nsrc, psunit, nplt, nfitparams, nptstofit
!ps       integer       ndsh(nfun), ibe(nfun), ien(nfun), lwid(nfun)
!ps       real          funmx, funmn, xnmn, xnmx, plt1x0, plt1y0, xoffset, yoffset, xsize, ysize, font, delta, xlam1, xlam2
!ps     &             , annsizefact
!ps       real          fun(npts,nfun), xn(npts,nfun), ssze(nfun), uppr(3), rght(3), xerrp(npts,nfun), xerrm(npts,nfun)
!ps     &             , yerrp(npts,nfun), yerrm(npts,nfun)
!ps       real*8        chisq, almostzero, wavelength0, frequency0, opacity0, dust2gas, opacit0relerr, dustgasrelerr
!ps       real*8        fitparameters(nmodparams), quantities(*)
!ps       parameter   ( almostzero = 1.0d-20, lportrait = .true. )
!ps       external      lastc, plotnd, ps2pdf, newdev, psinit, chopit, gsav, psmv1d, grest, pschars
!ps       common / copacity / frequency0, wavelength0, opacity0, dust2gas, opacit0relerr, dustgasrelerr, iverbose
!ps!________________________________________________________________________________________________________________________________
!ps!
!ps! Plot fitted SEDs.
!ps
!ps       if (lplotting) then
!ps
!ps         nu = nptstofit - nfitparams
!ps         write (cnu,'(i7)') nu
!ps         if (nu .gt. 0) then
!ps           nm = int ( log10 ( dble ( nu ) ) ) + 1
!ps         else
!ps           nm = 1
!ps         endif
!ps         write (cnsrc,'(i10)') nsrc
!ps         nn = int ( log10 ( dble ( nsrc ) ) ) + 1
!ps         write (ctdf,'(f6.1)') fitparameters(1)
!ps         nt = int ( log10 ( fitparameters(1) ) ) + 1
!ps         if (iverbose .eq. 2) then
!ps           if (iotty .gt. 0) write (iotty,'(/a)') '  PLOTSED: Plotting: '''//postscript//''''
!ps           if (iolog .gt. 0) write (iolog,'(/a)') '  PLOTSED: Plotting: '''//postscript//''''
!ps         endif
!ps         nplt = nplt + 1
!ps        
!ps! Open PostScript file for plotting (if required).
!ps! NOTE: Don't modify the  below: it is a preprocessing directive in my compilation script.
!ps
!ps         if (nplt .eq. 1) then
!ps           call newdev ( postscript )
!ps           call psinit ( psunit, lportrait, ctime, cdate )
!ps         else
!ps           call chopit ( 0., 0. )
!ps         endif
!ps         phdr = ''
!ps         pftr(1) = ''
!ps         pftr(2) = ''
!ps         pftr(3) = ''
!ps         pftr(4) = ''
!ps         write (pftr(1),'(a)') postscript
!ps         write (pftr(2),'(a10,''  '',a8)') cdate, ctime
!ps        
!ps!!         xlab(1) = 'Wavelength |F5|l|BF0| |N|(um)'
!ps         xlab(1) = 'Wavelength |F5|l|BF0| |N|(|F5|u|F0|m)'
!ps         ylab(1) = '  Total flux |F10|F|BF5|n |NF0|(Jy)'
!ps         xlab(2) = 'Frequency |F5|n|F0| (GHz)'
!ps         ylab(2) = ' '
!ps         varnam(1) = ' '
!ps         varnam(2) = '|F10|F|BF5|n|BF0|,obs'
!ps         if (cmodel .eq. 'modbody') then
!ps           varnam(3) = '|F10|F|BF5|n|BF0|,GB ('//ctdf(6-nt-1:6)//'|B|K)'
!ps         else
!ps           varnam(3) = '|F10|F|BF5|n|BF0|,TB ('//ctdf(6-nt-1:6)//'|B|K)'
!ps         endif
!ps         varnam(4) = '|F10|F|BF5|n|BF0|,BB ('//ctdf(6-nt-1:6)//'|B|K)'
!ps         varnam(5) = '|F10|F|BF5|n|BF0|,BBRJ'
!ps         xhlp = 'all'
!ps         yhlp = 'all'
!ps
!ps         xnmx = 1000.0
!ps         xnmn = 1000.0
!ps         funmx = -1000.0
!ps         funmn = 0.1
!ps         xlam2 = xn(1,2)
!ps         do j=1,nfun
!ps           ibe(j) = 1
!ps           ien(j) = npts
!ps           do i=1,npts
!ps             xnmn = min ( 0.9 * xn(i,1), xnmn )
!ps             if (j .le. 2 .and. nint ( xn(i,1) ) .ne. 255) funmx = max ( fun(i,1), funmx )
!ps             if (xn(i,2) .lt. 1.0e5) xlam2 = xn(i,2)
!ps           enddo
!ps         enddo
!ps         xnmn = 1.25 * xnmn
!ps         funmx = 1.25 * funmx
!ps         funmn = funmx / 100.0
!ps
!ps         do j=1,nfun
!ps           do i=1,npts
!ps             xerrp(i,j) = 1.0e-31
!ps             xerrm(i,j) = 1.0e-31
!ps             if (fun(i,j) .lt. funmn) then
!ps               yerrp(i,j) = 0.0
!ps               yerrm(i,j) = 0.0
!ps             endif
!ps           enddo
!ps         enddo
!ps
!ps! NOTE: Don't modify the "ps" below: it is a preprocessing directive in my compilation script.
!ps
!ps         uppr(1) = 3.0e5
!ps         uppr(2) =-1.0
!ps         uppr(3) = 0.0
!ps         rght(1) = 0.0
!ps         rght(2) = 1.0
!ps         rght(3) = 0.0
!ps        
!ps         xlam1 = xn(1,2)
!ps         write (clam1,'(i6)') nint ( xlam1 )
!ps         write (clam2,'(i6)') nint ( xlam2 )
!ps         nl1 = int ( log10 ( dble ( nint ( xlam1 ) ) ) ) + 1
!ps         nl2 = int ( log10 ( dble ( nint ( xlam2 ) ) ) ) + 1
!ps        
!ps         if (clam1(6-nl1+1:6) .eq. '071' .or. clam1(6-nl1+1:6) .eq. '075' .or. clam1(6-nl1+1:6) .eq. '101' .or. 
!ps     &       clam1(6-nl1+1:6) .eq. '105' .or. clam1(6-nl1+1:6) .eq. '161' .or. clam1(6-nl1+1:6) .eq. '165' .or. 
!ps     &       clam1(6-nl1+1:6) .eq. '251' .or. clam1(6-nl1+1:6) .eq. '255') then
!ps           color (1) = 'bro'
!ps          else
!ps           color (1) = 'gre'
!ps         endif
!ps         symbol(1) = 'circle'
!ps         ssze  (1) = 1.20
!ps         ndsh  (1) = 0
!ps         lwid  (1) = 0
!ps         color (2) = 'blu'
!ps         symbol(2) = 'circle'
!ps         ssze  (2) = 1.20
!ps         ndsh  (2) = 0
!ps         lwid  (2) = 0
!ps         color (3) = 'red'
!ps         symbol(3) = 'circle'
!ps         ssze  (3) = 0.0
!ps         ndsh  (3) = 0
!ps         lwid  (3) = 7
!ps         color (4) = 'bla'
!ps         symbol(4) = 'circle'
!ps         ssze  (4) = 0.0
!ps         ndsh  (4) = 4
!ps         lwid  (4) = 3
!ps         color (5) = 'bla'
!ps         symbol(5) = 'circle'
!ps         ssze  (5) = 0.0
!ps         ndsh  (5) = 8
!ps         lwid  (5) = 3
!ps
!ps! Position plot on page.
!ps
!ps         plt1x0 =  4.5
!ps         plt1y0 = 17.6
!ps         xsize = 10.0
!ps         ysize = 10.0
!ps         allxl = 2
!ps         allyl = 2
!ps         
!ps! Produce the plot.         
!ps
!ps         call gsav ()
!ps
!ps         annsizefact = 0.8
!ps
!ps         call psmv1d ( xn, fun, ibe, ien, npts, nfun, varnam, xlab, ylab, phdr, pftr, plt1x0, plt1y0, 'loglog', 'nono', xnmn
!ps     &               , xnmx, uppr, funmn, funmx, rght, xhlp, yhlp, allxl, allyl, color, ndsh, lwid, ssze, symbol, xerrp, xerrm
!ps     &               , yerrp, yerrm, xsize, ysize, annsizefact )
!ps     
!ps         call grest ()
!ps
!ps         xoffset = 10.3
!ps         yoffset = 5.0
!ps         if (nmodparams .lt. 4) yoffset = yoffset - 2.4
!ps         font = 0.22
!ps         delta = 0.45
!ps        
!ps         if (.not.limages) then
!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset + 5.4, font
!ps     &                  , 'Source # '//cnsrc(10-nn+1:10), 0., 0 )
!ps         else
!ps           write (cx,'(i4)') nx0o
!ps           write (cy,'(i4)') ny0o
!ps           ncx = int ( log10 ( dble ( nx0o ) ) ) + 1
!ps           ncy = int ( log10 ( dble ( ny0o ) ) ) + 1
!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset + 5.4, font
!ps     &                  , 'Pixel # '//cnsrc(10-nn+1:10), 0., 0 )
!ps           call pschars ( plt1x0 + xoffset + 1.3, plt1y0 + yoffset + 5.4, 0.2
!ps     &                  , '('//cx(6-ncx+1:6)//','//cy(6-ncy+1:6)//')', 0., 0 )
!ps         endif
!ps        
!ps         call pschars ( plt1x0 + xoffset, plt1y0 + yoffset + 0.5, font
!ps     &                , '|F5|l|F0| '//clam1(6-nl1+1:6)//'|F5|-|F0|'//clam2(6-nl2+1:6)//'|B| |F5|u|F0|m', 0., 0 )
!ps!!     &                , '|F5|l|F0| '//clam1(6-nl1+1:6)//'|F5|-|F0|'//clam2(6-nl2+1:6)//'|B| |F0|um', 0., 0 )
!ps        
!ps         yoffset = yoffset - 0.15
!ps         write (cso,'(f6.2)') quantities(1)
!ps         write (clo,'(f6.3)') quantities(2) / (quantities(1) + 1.0e-30)
!ps         ncx = 1
!ps         if (cso(1:1) .eq. ' ') ncx = 1
!ps         if (cso(2:2) .eq. ' ') ncx = 2
!ps         ncy = 1
!ps         if (clo(1:1) .eq. ' ') ncy = 2
!ps         if (clo(2:2) .eq. ' ') ncy = 3
!ps         call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!ps     &                , '|F10|T|B|  |NF0|='//cso(ncx:6)//'|B| |N|K |BF0| |F0|', 0., 0 )
!ps         call plsmin ( 999., 999., 0.9*font, 0., 0 )
!ps         call pschars ( 999., 999., 0.8*font, ' '//clo(ncy:6), 0., 0 )
!ps        
!ps         yoffset = yoffset - delta
!ps         if (cmodel .eq. 'modbody') then
!ps           write (cso,'(f8.4)') quantities(3)
!ps           write (clo,'(f6.3)') quantities(5) / (quantities(3) + 1.0e-30)
!ps           ncx = 1
!ps           if (cso(1:1) .eq. ' ') ncx = 1
!ps           if (cso(2:2) .eq. ' ') ncx = 2
!ps           ncy = 1
!ps           if (clo(1:1) .eq. ' ') ncy = 2
!ps           if (clo(2:2) .eq. ' ') ncy = 3
!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!ps     &                  , '|F10|M|B| |NF0|='//cso(ncx:8+ncx-2)//'|B| |NF10|M|BF18|*|NBF0|  |F0|'//clo(1:6), 0., 0 )
!ps           call plsmin ( 999., 999., 0.9*font, 0., 0 )
!ps           call pschars ( 999., 999., 0.8*font, ' '//clo(ncy:6), 0., 0 )
!ps         endif
!ps         if (cmodel .eq. 'thinbody') then
!ps           write (cso,'(f8.4)') quantities(3)
!ps           write (clo,'(f6.3)') quantities(5) / (quantities(3) + 1.0e-30)
!ps           ncx = 1
!ps           if (cso(1:1) .eq. ' ') ncx = 1
!ps           if (cso(2:2) .eq. ' ') ncx = 2
!ps           ncy = 1
!ps           if (clo(1:1) .eq. ' ') ncy = 2
!ps           if (clo(2:2) .eq. ' ') ncy = 3
!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!ps     &                  , '|F10|M|B| |NF0|='//cso(ncx:8+ncx-2)//'|B| |NF10|M|BF18|*|NBF0|  |F0|', 0., 0 )
!ps           call plsmin ( 999., 999., 0.9*font, 0., 0 )
!ps           call pschars ( 999., 999., 0.8*font, ' '//clo(ncy:6), 0., 0 )
!ps         endif
!ps        
!ps         if (nmodparams .ge. 4) then
!ps           yoffset = yoffset - delta
!ps           write (cso,'(f5.2)') quantities(12)
!ps           write (clo,'(f6.3)') quantities(13) / (quantities(12) + 1.0e-30)
!ps           ncy = 1
!ps           if (clo(1:1) .eq. ' ') ncy = 2
!ps           if (clo(2:2) .eq. ' ') ncy = 3
!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!ps     &                  , '|F5|b|B|  |NF0|='//cso(1:5)//'  |F0|', 0., 0 )
!ps           call plsmin ( 999., 999., 0.9*font, 0., 0 )
!ps           call pschars ( 999., 999., 0.8*font, ' '//clo(ncy:6), 0., 0 )
!ps           if (cmodel .eq. 'modbody') then
!ps             yoffset = yoffset - delta
!ps             write (cso,'(1pe9.2)') fitparameters(4)
!ps             cso(6:6) = 'e'
!ps             write (clo,'(f6.3)') quantities(9) / (quantities(8) + 1.0e-30)
!ps             ncy = 1
!ps             if (clo(1:1) .eq. ' ') ncy = 2
!ps             if (clo(2:2) .eq. ' ') ncy = 3
!ps             call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!ps     &                    , '|F5|W|B| |NF0|='//cso//'|B| |N|ster  |F0|', 0., 0 )
!ps             call plsmin ( 999., 999., 0.9*font, 0., 0 )
!ps             call pschars ( 999., 999., 0.8*font, ' '//clo(ncy:6), 0., 0 )
!ps           endif
!ps           yoffset = yoffset - delta
!ps           write (cso,'(i5)') nint ( quantities(14) )
!ps           write (clo,'(f6.3)') quantities(15) / (quantities(14) + 1.0e-30)
!ps           ncx = int ( log10 ( quantities(14) ) ) + 1
!ps           ncy = 1
!ps           if (clo(1:1) .eq. ' ') ncy = 2
!ps           if (clo(2:2) .eq. ' ') ncy = 3
!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!ps     &                  , '|F10|D|B|  |NF0|='//cso(5-ncx:5)//'|B| |N|pc  |F0|', 0., 0 )
!ps           call plsmin ( 999., 999., 0.9*font, 0., 0 )
!ps           call pschars ( 999., 999., 0.8*font, ' '//clo(ncy:6), 0., 0 )
!ps         endif
!ps
!ps         yoffset = yoffset - delta - 0.15
!ps         write (cchi,'(f8.4)') chisq
!ps         ncx = 1
!ps         if (cchi(1:1) .eq. ' ') ncx = 1
!ps         if (cchi(2:2) .eq. ' ') ncx = 2
!ps         call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!ps     &                , '|F5|c|SF0|2 |NF0|='//cchi(ncx:8+ncx-2)//' ', 0.,0 )
!ps         if (chisq / dble ( nu + 1 ) .le. 1.0d0) then
!ps           call pschars ( 999., 999., 0.8*font, '  ('//cnu(7-nm+1:7)//') ok', 0.,0 )
!ps         else
!ps           call pschars ( 999., 999., 0.8*font, '  ('//cnu(7-nm+1:7)//') bad', 0.,0 )
!ps         endif
!ps        
!ps         yoffset = yoffset - delta - 0.15
!ps         write (cso,'(f8.4)') quantities(6)
!ps         write (clo,'(f6.3)') quantities(7) / (quantities(6) + almostzero)
!ps         ncx = 1
!ps         if (cso(1:1) .eq. ' ') ncx = 1
!ps         if (cso(2:2) .eq. ' ') ncx = 2
!ps         ncy = 1
!ps         if (clo(1:1) .eq. ' ') ncy = 2
!ps         if (clo(2:2) .eq. ' ') ncy = 3
!ps         call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!ps     &                , '|F10|L|B|  |NF0|='//cso(ncx:8+ncx-2)//'|B| |NF10|L|BF18|*|NBF0|  |F0|', 0., 0 )
!ps         call plsmin ( 999., 999., 0.9*font, 0., 0 )
!ps         call pschars ( 999., 999., 0.8*font, ' '//clo(ncy:6), 0., 0 )
!ps
!ps         if (cmodel .eq. 'modbody') then
!ps           yoffset = yoffset - delta
!ps           write (cso,'(1pe9.2)') quantities(10)
!ps           cso(6:6) = 'e'
!ps           write (clo,'(f6.3)') quantities(5) / (quantities(3) + 1.0e-30)
!ps           ncy = 1
!ps           if (clo(1:1) .eq. ' ') ncy = 2
!ps           if (clo(2:2) .eq. ' ') ncy = 3
!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!ps     &                  , '|F10|N|B| |NF0|='//cso//'|B| |N|cm|S|-2|NF0|  |F0|', 0., 0 )
!ps           call plsmin ( 999., 999., 0.9*font, 0., 0 )
!ps           call pschars ( 999., 999., 0.8*font, ' '//clo(ncy:6), 0., 0 )
!ps         endif
!ps        
!ps         yoffset = yoffset - delta - 0.15
!ps         write (cso,'(f5.1)') opacity0
!ps         write (clo,'(f6.3)') opacit0relerr
!ps         ncy = 1
!ps         if (clo(1:1) .eq. ' ') ncy = 2
!ps         if (clo(2:2) .eq. ' ') ncy = 3
!ps         call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!ps     &                , '|F5|k|BF0|0|F0| |NF0|='//cso(1:5)//' |N|cm|B| |N|g|SF5|-|SF0|1|BF0|  |F0|', 0., 0 )
!ps         call plsmin ( 999., 999., 0.9*font, 0., 0 )
!ps         call pschars ( 999., 999., 0.8*font, ' '//clo(ncy:6), 0., 0 )
!ps        
!ps         yoffset = yoffset - delta
!ps         write (cso,'(f6.3)') dust2gas
!ps         write (clo,'(f6.3)') dustgasrelerr
!ps         ncy = 1
!ps         if (clo(1:1) .eq. ' ') ncy = 2
!ps         if (clo(2:2) .eq. ' ') ncy = 3
!ps         call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font, '|F5|h|BF0| |NF0|='//cso(1:6)//'  |F0|', 0., 0 )
!ps         call plsmin ( 999., 999., 0.9*font, 0., 0 )
!ps         call pschars ( 999., 999., 0.8*font, ' '//clo(ncy:6), 0., 0 )
!ps        
!!!ps         yoffset = yoffset - delta - 0.15
!!!ps         write (clo,'(1pe9.2)') quantities(2) / (quantities(1) + 1.0e-30)
!!!ps         clo(6:6) = 'e'
!!!ps         call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!!!ps     &                , '|F5|e|BF10|T|F0| |NF0|='//clo, 0., 0 )
!!!ps        
!!!ps         yoffset = yoffset - delta
!!!ps         if (distancepc .gt. 1.001d0) then
!!!ps           write (clo,'(1pe9.2)') quantities(5) / (quantities(3) + 1.0e-30)
!!!ps           clo(6:6) = 'e'
!!!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!!!ps     &                  , '|F5|e|BF10|M|F0| |NF0|='//clo, 0., 0 )
!!!ps         else
!!!ps           write (clo,'(1pe9.2)') quantities(5) / (quantities(3) + 1.0e-30)
!!!ps           clo(6:6) = 'e'
!!!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!!!ps     &                  , '|F5|e|BF10|N|F0| |NF0|='//clo, 0., 0 )
!!!ps         endif
!!!ps        
!!!ps         if (nmodparams .ge. 4) then
!!!ps           yoffset = yoffset - delta
!!!ps           write (clo,'(1pe9.2)') quantities(13) / (quantities(12) + 1.0e-30)
!!!ps           clo(6:6) = 'e'
!!!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!!!ps     &                  , '|F5|e|BF5|b|F0| |NF0|='//clo, 0., 0 )
!!!ps           if (cmodel .eq. 'modbody') then
!!!ps             yoffset = yoffset - delta
!!!ps             write (clo,'(1pe9.2)') quantities(9) / (quantities(8) + 1.0e-30)
!!!ps             clo(6:6) = 'e'
!!!ps             call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!!!ps     &                    , '|F5|e|BF5|W|F0| |NF0|='//clo, 0., 0 )
!!!ps           endif
!!!ps           yoffset = yoffset - delta
!!!ps           write (clo,'(1pe9.2)') quantities(15) / (quantities(14) + 1.0e-30)
!!!ps           clo(6:6) = 'e'
!!!ps           call pschars ( plt1x0 + xoffset, plt1y0 + yoffset, font
!!!ps     &                  , '|F5|e|BF10|D|F0| |NF0|='//clo, 0., 0 )
!!!ps         endif
!ps       endif
!ps
!ps       return
!ps       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine rfits 
       
     &            ( nx, ny, bunit, ctype1, ctype2, rp1, rp2, xr, yr, fun, dx, dy, object, ra, dec, fname, iotty, iolog, creator
     &            , beam, funmin, funmax, blank, rot1, rot2, cd11, cd12, cd21, cd22, equinox, bzero, bscale, wave, datamn, datamx
     &            , history, iverbose )
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       simple, extend, anynull, myfile
       character*(*) ctype1, ctype2, creator, fname, bunit, object, history
       character*3   cnan
       character*25  cdatetime
       character*80  comment
       integer       i, j, nx, ny, stat, blksize, group, naxis, pcount, gcount, bitpix, unit, iotty, iolog, maxdim, rw
     &             , nfound, blank, stat05, stat06, stat07, stat08, stat09, stat10, stat11, stat12, stat13, stat14, stat15, stat16
     &             , stat17, stat18, stat19, stat20, stat21, stat22, stat23, stat24, stat25, stat26, stat27, stat28, stat29, stat30
     &             , stat31, stat32, stat0, iverbose
       parameter   ( maxdim = 2 )
       integer       naxes(maxdim)
       real*8        dx, dy, ra, dec, beam, funmin, funmax, xr, yr, rp1, rp2, rot1, rot2, epoch, equinox, bzero, bscale, wave
     &             , dxdeg, dydeg, cd11, cd12, cd21, cd22, almostzero, as2deg, datamx, datamn, sign1, arg1, arg2, argm
     &             , sinrota, cosrota, pi, angle
       real*8        fun(nx,ny)
       parameter   ( pi = 3.14159265358979d0, as2deg = 1.0d0 / 3600.0d0, almostzero = 1.0d-30 )
       external      ftopen, ftgkne, ftghpr, ftgkys, ftgkyd, ftgkye, ftg2de, ftclos, printerr, ftgkyj, ftg2dd
!__________________________________________________________________________________________________________________________________
!
! Initialize variables.

       stat = 0
       unit   = 20
       group  = 1
       rw     = 0
       stat05 = 0
       stat06 = 0
       stat07 = 0
       stat08 = 0
       stat09 = 0
       stat10 = 0
       stat11 = 0
       stat12 = 0
       stat13 = 0
       stat14 = 0
       stat15 = 0
       stat16 = 0
       stat17 = 0
       stat18 = 0
       stat19 = 0
       stat20 = 0
       stat21 = 0
       stat22 = 0
       stat23 = 0
       stat24 = 0
       stat25 = 0
       stat26 = 0
       stat27 = 0
       stat28 = 0
       stat29 = 0
       stat30 = 0
       stat31 = 0
       stat32 = 0

! Open the FITS file.

       call ftopen ( unit, fname, rw, blksize, stat )

       if (stat .ne. 0) then
         if (iotty .gt. 0) write (iotty,'(/a)') ' RFITS: ERROR: Trouble opening FITS file '''//fname//'''.'
         if (iolog .gt. 0) write (iolog,'(/a)') ' RFITS: ERROR: Trouble opening FITS file '''//fname//'''.'
         stop 99
       endif

! Determine the size of the image.

       call ftgkne ( unit, 'NAXIS', 1, 2, naxes, nfound, stat )

! Check that it found both NAXIS1 and NAXIS2 keywords.

       if (nfound .ne. 2) then
         if (iotty .gt. 0) write (iotty,'(/a)') ' RFITS: ERROR: Failed to read the NAXISn keywords.'
         if (iolog .gt. 0) write (iolog,'(/a)') ' RFITS: ERROR: Failed to read the NAXISn keywords.'
         stop 99
       endif

! Get the main primary header keywords which define the array structure.

!!       call ftgphx ( unit, maxdim, simple, bitpix, naxis, naxes, pcount, gcount, extend, bscale, bzero, blank, nblank, stat )

       call ftghpr ( unit, maxdim, simple, bitpix, naxis, naxes, pcount, gcount, extend, stat )

! Read all parameters.

       call ftgkys ( unit, 'CREATOR', creator, comment, stat )
       
       if (stat .eq. 0 .and. (comment(1:21) .eq. 'Alexander Menshchikov' .or. comment(1:22) .eq. 'Alexander Men''shchikov')) then
         myfile = .true.
       else
         myfile = .false.
         stat = 0
         creator = 'FITFLUXES'
       endif

! I don't understant why blank is assigned to my array values (if blank is not zero) - I have no undefined values in arrays. Maybe
! there is a bug in the version of FITSIO I use. For now I will ignore the blank value when reading FITS files.

       blank = 0

       datamx = 0.0d0
       datamn = 0.0d0
       bunit = ' '
       ctype1 = ' '
       ctype2 = ' '
       rp1 = 0.0d0
       rp2 = 0.0d0
       rot1 = 0.0d0
       rot2 = 0.0d0
       cd11 = 0.0d0
       cd12 = 0.0d0
       cd21 = 0.0d0
       cd22 = 0.0d0
       xr = 0.0d0
       yr = 0.0d0
       dx = 0.0d0
       dy = 0.0d0
       object = ' '
       ra = 0.0d0
       dec = 0.0d0
       epoch = 0.0d0
       equinox = 0.0d0
       beam = 0.0d0
       bzero = 0.0d0
       bscale = 1.0d0
       wave = 1.0d2
       history = ' '

       call ftgkys ( unit, 'DATE', cdatetime , comment, stat05 )
       call ftgkyd ( unit, 'DATAMAX', datamx , comment, stat06 )
       call ftgkyd ( unit, 'DATAMIN', datamn , comment, stat07 )
       call ftgkyj ( unit, 'BLANK'  , blank  , comment, stat08 )
       call ftgkys ( unit, 'BUNIT'  , bunit  , comment, stat09 )
       call ftgkyd ( unit, 'BZERO'  , bzero  , comment, stat10 )
       call ftgkyd ( unit, 'BSCALE' , bscale , comment, stat11 )
       call ftgkys ( unit, 'CTYPE1' , ctype1 , comment, stat12 )
       call ftgkys ( unit, 'CTYPE2' , ctype2 , comment, stat13 )
       call ftgkyd ( unit, 'CRPIX1' , rp1    , comment, stat14 )
       call ftgkyd ( unit, 'CRPIX2' , rp2    , comment, stat15 )
       call ftgkyd ( unit, 'CROTA1' , rot1   , comment, stat16 )
       call ftgkyd ( unit, 'CROTA2' , rot2   , comment, stat17 )
       call ftgkyd ( unit, 'CD1_1'  , cd11   , comment, stat18 )
       call ftgkyd ( unit, 'CD1_2'  , cd12   , comment, stat19 )
       call ftgkyd ( unit, 'CD2_1'  , cd21   , comment, stat20 )
       call ftgkyd ( unit, 'CD2_2'  , cd22   , comment, stat21 )
       call ftgkyd ( unit, 'CRVAL1' , xr     , comment, stat22 )
       call ftgkyd ( unit, 'CRVAL2' , yr     , comment, stat23 )
       call ftgkyd ( unit, 'CDELT1' , dxdeg  , comment, stat24 )
       call ftgkyd ( unit, 'CDELT2' , dydeg  , comment, stat25 )
       call ftgkys ( unit, 'OBJECT' , object , comment, stat26 )
       call ftgkyd ( unit, 'RA'     , ra     , comment, stat27 )
       call ftgkyd ( unit, 'DEC'    , dec    , comment, stat28 )
       call ftgkyd ( unit, 'EPOCH'  , epoch  , comment, stat29 )
       call ftgkyd ( unit, 'EQUINOX', equinox, comment, stat30 )
       call ftgkyd ( unit, 'WAVE'   , wave   , comment, stat31 )
       call ftgkyd ( unit, 'BEAM'   , beam   , comment, stat )
       call ftgkys ( unit, 'HISTORY', history, history, stat32 )

       if (stat .ne. 0) then
         beam = 0.0d0
         stat = 0
       endif

       stat0 = stat05 + stat06 + stat07 + stat09 + stat10 + stat11 + stat12 + stat13 + stat14 + stat15 + stat16 + stat17 + stat18
     &       + stat19 + stat20 + stat21 + stat22 + stat23 + stat24 + stat25 + stat26 + stat27 + stat28 + stat30 + stat31

       if (iverbose .eq. 2) then
         if (stat0  .gt. 0) write (iotty,'()')
         if (stat05 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword DATE absent or error'
         if (stat06 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword DATAMAX absent or error'
         if (stat07 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword DATAMIN absent or error'
         if (stat09 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword BUNIT absent or error'
         if (stat10 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword BZERO absent or error'
         if (stat11 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword BSCALE absent or error'
         if (stat12 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CTYPE1 absent or error'
         if (stat13 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CTYPE2 absent or error'
         if (stat14 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CRPIX1 absent or error'
         if (stat15 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CRPIX2 absent or error'
         if (stat16 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CROTA1 absent or error'
         if (stat17 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CROTA2 absent or error'
         if (stat18 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CD1_1 absent or error'
         if (stat19 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CD1_2 absent or error'
         if (stat20 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CD2_1 absent or error'
         if (stat21 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CD2_2 absent or error'
         if (stat22 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CRVAL1 absent or error'
         if (stat23 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CRVAL2 absent or error'
         if (stat24 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CDELT1 absent or error'
         if (stat25 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword CDELT2 absent or error'
         if (stat26 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword OBJECT absent or error'
         if (stat27 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword RA absent or error'
         if (stat28 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword DEC absent or error'
         if (stat30 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword EQUINOX absent or error'
         if (stat31 .ne. 0) write (iotty,'(a)') '   RFITS: WARNING: Keyword WAVE absent or error'
       endif
       if (stat26 .ne. 0 .or. object(1:10) .eq. '          ') object = '<Name>'

! Read the FITS file into the 2D array.

       call ftg2dd ( unit, group, dble ( blank ), naxes(1), naxes(1), naxes(2), fun, anynull, stat )

       if (stat .ne. 0) write (iotty,'(a)') ' Error reading FITS file'

! Close the file and free the unit number.

       call ftclos ( unit, stat )

! Check for any error, and if so print out error messages.

       if (stat .gt. 0) call printerr ( stat )

       if (stat17 .ne. 0) then
         rot2 = 0.0d0
       endif
       rot1 = rot2

       if (stat29 .eq. 0 .and. stat30 .ne. 0) then
         equinox = epoch
       endif

! CDELT* is absent or error.

       if (stat24 .ne. 0 .or. stat25 .ne. 0) then

! CD1* and CD2* are defined correctly: determine CDELT* and ROTA*.

         if (stat18 .eq. 0 .and. stat19 .eq. 0 .and. stat20 .eq. 0 .and. stat21 .eq. 0) then
           sign1 = sign ( 1.0d0, cd11 * cd22 - cd12 * cd21 )
           dxdeg = sign1 * sqrt ( cd11**2 + cd21**2 )
           dydeg = sqrt ( cd12**2 + cd22**2 )
           arg1 = sign1 * cd12 / cd22
           arg2 = -sign1 * cd21 / cd11
           argm = (arg1 + arg2) / 2.0d0
           if (abs ( arg1 - arg2 ) / argm .gt. 1.0d-10) then
             write (iotty,'(/a/)') '   RFITS: WARNING: Coordinate axes are not orthogonal'
           endif
           rot2 = atan ( argm )
           sinrota = sin ( rot2 )
           cosrota = cos ( rot2 )
           rot2 = rot2 * 180.0d0 / pi
           if (cosrota .ge. 0.0d0 .and. sinrota .ge. 0.0d0) then
             rot1 = rot2
           elseif (cosrota .lt. 0.0d0 .and. sinrota .ge. 0.0d0) then
             rot1 = 180.0d0 - rot2
           elseif (cosrota .lt. 0.0d0 .and. sinrota .lt. 0.0d0) then
             rot1 = 180.0d0 + rot2
           elseif (cosrota .ge. 0.0d0 .and. sinrota .lt. 0.0d0) then
             rot1 = - rot2
           endif
           rot2 = rot1
         endif
       endif

! CD1* and CD2* are absent or error.

       if (stat18 .ne. 0 .or. stat19 .ne. 0 .or. stat20 .ne. 0 .or. stat21 .ne. 0) then

! CDELT* are defined correctly: determine CD1* and CD2*.

         if (stat24 .eq. 0 .and. stat25 .eq. 0) then
           angle = rot2 * pi / 180.0d0
           cd11 = dxdeg * cos ( angle )
           cd12 = abs ( dydeg ) * sign ( 1.0d0, dxdeg ) * sin ( angle )
           cd21 = - abs ( dxdeg ) * sign ( 1.0d0, dydeg ) * sin ( angle )
           cd22 = dydeg * cos ( angle )
           if (abs ( cd12 ) .lt. almostzero) cd12 = 0.0d0
           if (abs ( cd21 ) .lt. almostzero) cd21 = 0.0d0
         endif
       endif

       dx = abs ( dxdeg ) / as2deg
       dy = abs ( dydeg ) / as2deg

! Print out the min and max values.

       nx = naxes(1)
       ny = naxes(2)
       funmin = 1.0d30
       funmax =-1.0d30
       cnan = 'NAN'
       do j=1,ny
         do i=1,nx
           if (isnan ( fun(i,j) )) then  
             fun(i,j) = 0.0d0
           endif  
           funmin = min ( funmin, fun(i,j) )
           funmax = max ( funmax, fun(i,j) )
           if (int ( min ( abs ( fun(i,j) ), 1.0d9 ) ) .eq. blank .and. blank .ne. 0) then  
             read (cnan,*) fun(i,j)
           endif
         enddo
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine wfits 
       
     &            ( cfitsversion, nx, ny, bunit, ctype1, ctype2, rp1, rp2, xr, yr, fun, dx, dy, object, ra, dec, fname, cdate
     &            , ctime, creator, beam, blank, rot1, rot2, cd11, cd12, cd21, cd22, equinox, bzero, bscale, wave, datamn, datamx
     &            , history )
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       simple, extend
       character*(*) fname, ctype1, ctype2, object, creator, ctime, cdate, bunit, cfitsversion, history
       character*19  cdatetime
       integer       i, j, nx, ny, status, blocksze, group, naxis, pcount, gcount, bitpix, unit, naxes(2), blank
       real*8        xr, yr, dx, dy, rp1, rp2, ra, dec, beam, bzero, bscale, as2deg, rot1, rot2, equinox, datamn
     &             , datamx, dxdeg, dydeg, wave, cd11, cd12, cd21, cd22, almostzero, pi
       real*8        fun(nx,ny)
       parameter   ( pi = 3.14159265358979d0, as2deg = 1.0d0 / 3600.0d0, almostzero = 1.0d-30 )
       external      fdelete, ftgiou, ftinit, ftphpr, ftpkys, ftpcom, ftpkyd, ftpkyj, ftphis, ftp2dd, ftclos, ftfiou, printerr
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
           if (.not.isnan ( fun(i,j) )) then
             if (fun(i,j) .lt. datamn ) datamn = fun(i,j)
             if (fun(i,j) .gt. datamx ) datamx = fun(i,j)
           endif
         enddo
       enddo

! When packing numbers using the formulas below, I got values of the function different in the 5-th digit from DATAMIN and
!  DATAMAX placed in the header when writing the data. From now on I will use bzero = 0 bscale = 1 and bitpix = -32, as these
! seem to cure the problem.

cc!      bitpix = 16
cc!      bzero  = 0.5d0 * (funmin + funmax)
cc!      bscale = dabs ( funmax - bzero ) / 32767.0d0

       bitpix = -32     !<-- single precision floating-point values
       bzero  = 0.0d0
       bscale = 1.0d0
       blank = blank

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

       call ftpkys ( unit, 'CREATOR', creator, 'Alexander Menshchikov, DAp IRFU CEA Saclay', status )
       call ftpcom ( unit,' ',status)
       call ftpkys ( unit, 'DATE'  , cdatetime, 'creation date and time', status )

! Write all parameters.

       dxdeg = - dx * as2deg
       dydeg = dy * as2deg

       call ftpkyd ( unit, 'BZERO'  , bzero , 13 , 'zero point in scaling equation', status )
       call ftpkyd ( unit, 'BSCALE' , bscale, 13 , 'linear factor in scaling equation' , status )
       call ftpkyd ( unit, 'DATAMAX', datamx, 13 , 'maximum data value', status )
       call ftpkyd ( unit, 'DATAMIN', datamn, 13 , 'minimum data value', status )
!!       if  (blank .gt. 0)
!!     & call ftpkyj ( unit, 'BLANK'  , blank      , 'value used for undefined array elements', status )
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
       call ftpkyd ( unit, 'CDELT1' , dxdeg , 13 , 'coordinate increment along axis', status )
       call ftpkyd ( unit, 'CDELT2' , dydeg , 13 , 'coordinate increment along axis', status )
       call ftpkys ( unit, 'OBJECT' , object     , 'object identifier', status )
       call ftpkyd ( unit, 'RA'     , ra     , 13, 'right ascension', status )
       call ftpkyd ( unit, 'DEC'    , dec    , 13, 'declination', status )
       call ftpkyd ( unit, 'EQUINOX', equinox, 13, 'equinox', status )
       call ftpkyd ( unit, 'WAVE'   , wave   , 13, 'wavelength (microns)', status )
       call ftphis ( unit, history, status )

       if (beam .ne. 0.0d0) call ftpkyd ( unit, 'BEAM', beam, 13, 'convolved with the beam (arcsec)', status )

! Write the array to the FITS file.

       call ftp2dd ( unit, group, nx, nx, ny, fun, status )

! Close the file.

       call ftclos ( unit, status )

! Free the unit number.

       call ftfiou ( unit, status )

! Check for any error, and if so print out error messages

       if (status .gt. 0) call printerr ( status )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
