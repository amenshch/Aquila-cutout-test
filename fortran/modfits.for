
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       program modfits
!__________________________________________________________________________________________________________________________________
!
! Modify (in many ways) a single FITS image and some keywords of its header.
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extrmaction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!__________________________________________________________________________________________________________________________________
!                 
       implicit      none

       integer       inmx, iotty, nextrmax, nlevelmax
       parameter   ( iotty = 6, nextrmax = 10000000, nlevelmax = 1000000 )

       logical       inellipse, lfnama, changed, lunix, l1value, l4values, ldebug, l4conn, l8conn, lskeletons, lazimuthaverage
     &             , lnorthborder, lsouthborder, leastborder, lwestborder, l8simple, l4isolated, l8isolated, l4endpoint, l8endpoint
     &             , lfixdiags, lnoname, lmiss, lnobranches, lsknorm, lsmall, lcompact, lfilamen, lcutsources, lsources !!, lexist
     &             , lcutfilaments, lcutpeaks, lcutnegat, ltiming, ltime, lcutshapes, lshapes
                                           
       character*1   answer, cha
       character*3   dot
       character*5   cnmax
       character*6   cfitsversion, cverbose, cexpo
       character*7   clibname, cnrem, cnorm, cminpx, cx, cy
       character*8   ctime, cfpptime
       character*10  ccfact, cbsize, cbarea, cpixel, cdate, crout
       character*11  cwhat
       character*21  compda
       character*80  object, creator, ctype1, ctype2, cvalue, cvalue1, cvalue2, cvalue3, cvalue4, bunit, history
       character*128 maction, param, cline
       character*500 arg1, arg2, arg3, arg4, arg5, arg6, arg7, filename, outname, profname

       integer       firstb, lastc, lasts, i, j, k, isp1, fnlen1, ion, ndate, nx, ny, nx1, ny1, nx2, ny2, nxo, nyo, ic, jc, nconn
     &             , ia1, ia2, ia3, ia4, ia5, ia6, ia7, nbw, iact, ival, ipar, blank, icva, im1, ip1, jm1, jp1, lb, n, i1, i2, j1
     &             , ic1, ic2, ic3, icva1, icva2, icva3, icva4, nextr, irc, nrem, naddpixmax, nymin, nymax, nextrmx, j2, ii, jj
     &             , ipn, ix1, ix2, iy1, iy2, iu1, iu2, iv1, iv2, ixlength, iylength, iulength, ncon, iminp, ivlength, l, m, nlevmax
     &             , lengthmin, idirection, iwidth, nskw, numpix4, numpix8, lev, npts, ik, jk, idum, npt(10000), newrandomseq, ipt
     &             , newran, npoints, npp1, n1, n2, m1, m2, nto, kz, numa1, numb1, numa2, numa4, n248, n246, nw, ip2, jp2, im2, jm2
     &             , inrem, iter, itermax, i0, j0, ncx, ncy, ixco, jyco, npxmax, nextro, ndelta, nozero, nsrco, itgmax, iplu, iexc
     &             , nlevpos, nlevneg, nlevsmin, wchours, wcmins, wcsecs, cpuhours, cpumins, cpusecs, isl, icom, nwin, nm, ilen

       real          fitsvers
       real*8        timer, dx, dy, sfact, crpix1, crpix2, crval1, crval2, value1, value2, value3, value4, conversionfactor
     &             , rcore, ra, dec, beam, funmin, funmax, const, fJypixel, sq_arcsecs_per_sterad, fnewmin, fnewmax, value !!,thresh
     &             , sq_arcsecs_per_beam_area, pi, crota1, crota2, equinox, bzero, bscale, datamin, datamax, funfit, minpix
     &             , almostzero, deriv1a, deriv1b, deriv1c, deriv1d, deriv1e, deriv1f, deriv1g, deriv1h, deriv2a, deriv2b, deriv2c
     &             , deriv2d, wave, factMJysr2Jybeam, factMJysr2mJybeam, factJypixel2Jybeam, factJypixel2mJybeam, background, rxy
     &             , pixels_per_beam_area, factJypixel2MJysr, boupixval, newpixval, zeros, beamsize, beamarea, rad, factor, factpix
     &             , sparsity, goodtoremove, cd11, cd12, cd21, cd22, lowintensity, dlograd, tdust, dbll, dblm, clefactor !,funminlog
     &             , directioncode, levelx, xi1, xi2, yj1, yj2, tlen, cosa, sina, pt1, bbody1, bbody2, srms, ran1, gasdev, q1
     &             , totint, totrms, xpoints, rw2, rx2, ry2, fmean, variance, sigma, delta, tdc, tdmin, tdmax, dtdi, dtdj, funbmin
     &             , expo, sigm2g, sigm2m, delx, dely, rad2, argex, fitlev, hwhm, correct, radmax2, peak, dxym, radfitx
     &             , muH2, amu, speedolight, dust2gas, frequency, frequency0, beta, opacity, opacity0, levelxlog, dlpos !!, dlneg
     &             , deltpos, deltapos, elongmaxsrc, elongminfil, sparsmaxsrc, sparsminfil, funmaxlog, Rout !!, deltneg, deltaneg
     &             , sigm2g2, gauss, plaw, fitfact, numpart, numtotal, isumpct, radfit, expfit, plwfit, betax, selectk
     &             , eps, epsiln, factgeom, dgausdr, radmax, Rout2, gam1, gam2, gamma, wsize, xi, E, F, G, S, Z, Has, P, T, U, V, Q
     &             , chi

       real*8        cpu_modfits, wal_modfits, cput_modfits, walt_modfits, cputot, wctot

       integer, allocatable :: n1x(:), n2x(:), ntouching(:,:), nsx(:), nxmn(:), nxmx(:), nymn(:), nymx(:), cimask(:,:), npixfil(:)
     &             , np(:), nsourceo(:)
       real*8 , allocatable :: funa(:,:), funb(:,:), worker(:,:), nsegm(:,:), mxco(:), myco(:), deltx(:), delty(:), nonzero(:)
     &             , afwhm(:), bfwhm(:), atheta(:), equivrad(:), momxco(:), momyco(:), elongation(:), arg(:) !,func(:,:), fund(:,:)
     &             , argz(:), prof(:), profrms(:), slope(:), afoot(:), bfoot(:), xpix(:), ypix(:), pixdist(:,:), slo(:), slom(:)
     &             , origmask(:,:), xpixok(:,:), ypixok(:,:), pixd(:), nonzeroo(:), nsegmo(:,:), shapes(:,:), level(:)

       parameter   ( sq_arcsecs_per_sterad = 3282.80635d0 * 3600.0d0**2, pi = 3.14159265358979d0, almostzero = 1.0d-30 )
       parameter   ( boupixval = 0.0d0, fitlev = 1.0d0, itgmax = 99, beta = 2.0d0, muH2 = 2.8d0, speedolight = 2.99792458d10 )
       parameter   ( amu = 1.6605402D-24, frequency0 = 1.0d12, opacity0 = 10.0d0, dust2gas = 1.0d-2, dot = '•'
     &             , ldebug = .false., ltiming = .false.

! Tests with v200408:
! Note1: elongmaxsrc < 1.45 does NOT work: it suddenly drastically clips elongated sources.
! Note2: elongmaxsrc = 1.47 seems to work well, clipping stuff unrelated to clear source peaks.
! Note3: sparsmaxsrc < 1.39 does NOT work: it suddenly drastically clips all sources except their tiny peaks few pixels across.
! Note4: sparsminfil = 1.39 also seems to work well for the spiral filament of the benchmark.
! Note5: elongminfil < 1.5 cuts too much of the roundish background into filaments, not good.

     &             , elongmaxsrc = 1.47d0, sparsmaxsrc = 1.39d0  !<-- for clipping sources
     &             , elongminfil = 3.00d0, sparsminfil = 1.39d0  !<-- for clipping filaments
     &             , deltapos = 0.05d0 )                         !, deltaneg = 0.1d0

       external      firstb, lastc, lastbs, ran1, gasdev, timer, rfits, wfits, when, osystem, ftvers, lasts, getfitshead, planck
     &             , convert84, fillarea, sizemeasure, inellipse, traceskels, expandit, shrinkit, showprogress
!__________________________________________________________________________________________________________________________________
!
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

       call getarg ( 1, arg1 )
       call getarg ( 2, arg2 )
       call getarg ( 3, arg3 )
       call getarg ( 4, arg4 )
       call getarg ( 5, arg5 )
       call getarg ( 6, arg6 )
       call getarg ( 7, arg7 )

       ia1 = lastc ( arg1 )
       ia2 = lastc ( arg2 )
       ia3 = lastc ( arg3 )
       ia4 = lastc ( arg4 )
       ia5 = lastc ( arg5 )
       ia6 = lastc ( arg6 )
       ia7 = lastc ( arg7 )
       
       cverbose = '-verb2'
       if (ia7 .gt. 0 .and. arg7(1:ia7-1) .eq. '-verb') cverbose = arg7(1:ia7)
       if (ia6 .gt. 0 .and. arg6(1:ia6-1) .eq. '-verb') cverbose = arg6(1:ia6)
       if (ia5 .gt. 0 .and. arg5(1:ia5-1) .eq. '-verb') cverbose = arg5(1:ia5)
       if (ia4 .gt. 0 .and. arg4(1:ia4-1) .eq. '-verb') cverbose = arg4(1:ia4)
       if (ia3 .gt. 0 .and. arg3(1:ia3-1) .eq. '-verb') cverbose = arg3(1:ia3)
       
       if (cverbose .eq. '-verb2') then
         write (iotty,'( )')
         if (arg1(1:ia1) .eq. ':') write (iotty,'(a$)') '  '
         write (iotty,'(a$)') ' MODFITS '//dot//' Modify FITS Image or its Header '//dot//' '//compda
       endif
       if (arg1(1:ia1) .eq. ':') stop
       if (cverbose .eq. '-verb2') write (iotty,'()')
       if (cverbose .eq. '-verb2') write (iotty,'(a)') ' Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
       if (cverbose .eq. '-verb2') write (iotty,'(a)') ' Using '//clibname(lb:7)//' library version'//cfitsversion
     &                                               //' by William D Pence.'
       if (cverbose .eq. '-verb1') then
           write (iotty,'(a)') ' MODFITS: '//arg1(1:ia1)//' '//arg2(1:ia2)//' '//arg3(1:ia3)//' '//arg4(1:ia4)//' '
     &                                     //arg5(1:ia5)//' '//arg6(1:ia6)//' '//arg7(1:ia7)
       endif
       maction = ' '
       value = 0.0d0
       cvalue = ' '
       param = ' '
       icva = 0
       iact = 0
       ival = 0
       ipar = 1
       l1value = .true.
       l4values = .false.
       naddpixmax = 0
       ic2 = 0
       ic3 = 0
       correct = 0.0d0
       background = 0.0d0
       dtdi = 0.0d0
       dtdj = 0.0d0
       idirection = 0
       nextro = 0
       i1 = 0
       j1 = 0
       xi1 = 0.0d0
       xi2 = 0.0d0
       yj1 = 0.0d0
       yj2 = 0.0d0
       lsources = .false.
       lshapes = .false.
       
       if (ia4 .gt. 0 .and. (arg4(1:ia4) .ne. '-o' .and. arg3(1:ia3) .ne. '-o' .and. arg4(1:ia4-1) .ne. '-verb')) then
         maction = arg1(1:ia1)
         param = arg3(1:ia3)//' '
         filename = arg4(1:ia4)
         if (maction .ne. 'keyword' .and. maction .ne. 'key') then
! Sometimes due to mistakes on command line, the value may have a minus before a negative number, so we'll take care of that.
           if (arg2(1:2) .ne. '--') then
             read (arg2(1:ia2),*,err=90) value
           else
             read (arg2(2:ia2),*,err=90) value
           endif   
           goto 91
  90       continue         
           write (iotty,'(/a)') '   MODFITS: ERROR: Trouble reading a value from the string (1): '//arg2(1:ia2)
           stop 99
  91       continue
         endif
         cvalue = arg2(1:ia2)//' '
         icva = lastc ( cvalue ) + 1
         ipar = ia3 + 1
         ival = ia2
         iact = ia1
       elseif (ia3 .gt. 0 .and. (arg3(1:ia3) .ne. '-o' .and. arg2(1:ia2) .ne. '-o' .and. arg3(1:ia3-1) .ne. '-verb')) then
         maction = arg1(1:ia1)
         filename = arg3(1:ia3)
         param = 'n '
! Sometimes due to mistakes on command line, the value may have a minus before a negative number, so we'll take care of that.
         if (arg2(1:2) .eq. '--') arg2(1:2) = ' -'
         ic1 = index ( arg2(1:ia2), ',' )
         if (ic1 .eq. 0) then
           read (arg2(1:ia2),*,err=1) value
           l1value = .true.
           l4values = .false.
           cvalue = arg2(1:ia2)//' '
           icva = lastc ( cvalue ) + 1
           goto 3
   1       continue
             write (iotty,'(/a)') '   MODFITS: ERROR: Trouble reading a value from the string (2): '//arg2(1:ia2)
             stop 99
         else
           l1value = .false.
           l4values = .true.
           maction = arg1(1:ia1)
           iact = ia1
           cvalue1 = arg2(1:ic1-1)//' '
           if (maction(1:iact) .eq. 'proavg') then
             cvalue2 = arg2(ic1+1:ia2)//' '
           else
             ic2 = index ( arg2(ic1+1:ia2), ',' ) + ic1
             ic3 = index ( arg2(ic2+1:ia2), ',' ) + ic2
             cvalue2 = arg2(ic1+1:ic2-1)//' '
             cvalue3 = arg2(ic2+1:ic3-1)//' '
             cvalue4 = arg2(ic3+1:ia2)//' '
           endif
           icva1 = lastc ( cvalue1 ) + 1
           icva2 = lastc ( cvalue2 ) + 1
           icva3 = lastc ( cvalue3 ) + 1
           icva4 = lastc ( cvalue4 ) + 1
           read (cvalue1,*,end=2,err=2) value1
           read (cvalue2,*,end=2,err=2) value2
           if (maction(1:iact) .ne. 'proavg') then
           read (arg2(ic2+1:ic3-1),*,end=2,err=2) value3
           read (arg2(ic3+1:ia2),*,end=2,err=2) value4
           endif
           goto 3
   2       continue
             write (iotty,'(/a)') ' MODFITS: '//maction(1:iact)//' '//cvalue1(1:icva1)//cvalue2(1:icva2)
     &                                        //cvalue3(1:icva3)//cvalue4(1:icva4)//param(1:ipar)//arg3(1:ia3)
             write (iotty,'(/a)') '   MODFITS: ERROR: Trouble reading 2 or 4 values from the string: '//arg2(1:ia2)
             stop 99
         endif
   3     continue
         ipar = 2
         ival = ia2
         iact = ia1
       elseif (ia2 .gt. 0 .and. arg2(1:ia2) .ne. '-o') then  
         filename = arg2(1:ia2)
         maction = arg1(1:ia1)
         iact = ia1
         if (index ( maction(1:iact), '~' ) .gt. 0) then
           cvalue = arg2(1:ia2)//' '
           icva = lastc ( cvalue ) + 1
         endif
       elseif (ia1 .gt. 0) then
         if (index ( arg1(1:ia1), '~' ) .gt. 0) then
           maction = arg1(1:ia1)
           iact = ia1
         else
           filename = arg1(1:ia1)
         endif
       else
         write (iotty,'(a)')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USAGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'(a)')
         write (iotty,'(a)') ' modfits <maction> [<value>] [<param>] <imagein> [-o <imageout>] [-verb{0|1|2}]'
         write (iotty,'( )')   
         write (iotty,'(a)') ' This utility allows various modifications of a single FITS image; <imagein>'
         write (iotty,'(a)') ' and <imageout> are the names of the input and output images.'
         write (iotty,'(a)')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')   
         write (iotty,'(a)') ' <maction>'
         write (iotty,'( )')   
         write (iotty,'(a)') ' abs[olute] ..................... convert image to absolute values'
         write (iotty,'(a)') ' neg[ate] ....................... negate image (multiply the image by -1)'
         write (iotty,'(a)') ' nan[zero] ...................... convert all NaNs to zero value'
         write (iotty,'(a)') ' inv[ert] ....................... invert image (non-zero pixels)'
         write (iotty,'(a)') ' squ[are] ....................... convert image to squared values'
         write (iotty,'(a)') ' sqr[oot][m] .................... convert image to square root values'
         write (iotty,'(a)') ' log[10][m] ..................... convert image to logarithm values'
         write (iotty,'(a)') ' exp[onent] ..................... convert image to exponential values'
         write (iotty,'(a)') ' der[iv]1 ....................... convert image to 1st derivatives'
         write (iotty,'(a)') ' der[iv]2 ....................... convert image to 2nd derivatives'
         write (iotty,'(a)') ' edg[e] ......................... remove inner pixels of non-zero structures'
         write (iotty,'(a)') ' ann[ex] ........................ annex connected pixels by segmentation masks'
         write (iotty,'(a)') ' thi[cken] ...................... thicken structures by annexing one pixel'
         write (iotty,'(a)') ' con[tract] ..................... contract structures by removing one pixel'
         write (iotty,'(a)') ' mer[ge] ........................ merge structures separated by one pixel'
         write (iotty,'(a)') ' rec[onnect] .................... reconnect one-pixel breaks in skeletons'
         write (iotty,'(a)') ' bsr[core] ...................... background-subtract a model with rcore'
         write (iotty,'(a)') ' ren[umber] ..................... renumber segmentation image sequentially'
         write (iotty,'( )')   
         write (iotty,'(a)') ' <maction> <value>'
         write (iotty,'( )')   
         write (iotty,'(a)') ' ske[leton][1] [<width>] ........ skeletonize (Hilditch method) <width> pixels'
         write (iotty,'(a)') ' seg[ment]{4|8} [<minpix>] ...... segment image in shapes (of at least <minpix>)'
         write (iotty,'(a)') ' add[p]     <const> ............. add <const> to image or create a plane (addp)'
         write (iotty,'(a)') ' sub[tract] <const> ............. subtract <const> from image'
         write (iotty,'(a)') ' mul[tiply] <factor> ............ multiply image by <factor>'
         write (iotty,'(a)') ' div[ide]   <factor> ............ divide image by <factor>'
         write (iotty,'(a)') ' pow[er]    <expon> ............. raise image to the power <expon>'
         write (iotty,'(a)') ' expand     <number> ............ expand image at all edges by <number> pixels'
         write (iotty,'(a)') ' shrink     <number> ............ shrink image into itself by <number> pixels'
         write (iotty,'(a)') ' lay[er]    <value> ............. remove "zero" layer of pixels below <value>'
         write (iotty,'(a)') ' noi[se]    <wave> .............. add Gaussian pixel noise at <wave> (µm)'
         write (iotty,'(a)') ' surf[den]  <wave> .............. surfden: const / Bnu(<imagein>) at <wave> (µm)'
         write (iotty,'(a)') ' bbo[dy]    <wave> .............. produce Bnu(<imagein>) at <wave> (µm)'
         write (iotty,'(a)') ' pro[file|avg] <x1,y1,x2,y2> .... profile image along line or azimuth-averaged'
         write (iotty,'(a)') ' spr[ead]   <x1,y1,x2,y2> ....... spread {X|Y}-line of pixels along {Y|X}-axis'
         write (iotty,'(a)') ' fix[pix] <x1,y1,[-]beam,[-]peak> .. replace an area with a circular shape'
         write (iotty,'(a)') ' bor[der] [-]<number> ........... add or remove borders <number> pixels wide'
         write (iotty,'(a)') ' bor[der] [-]<left>,[-]<right>,[-]<lower>,[-]<upper> ... same, more flexible' 
         write (iotty,'( )')   
         write (iotty,'(a)') ' <maction> <value> <param>'
         write (iotty,'( )')   
         write (iotty,'(a)') ' max[imum] <maxval> {y|n} ........ set maximum to <maxval> (& zero pixels if "y")'
         write (iotty,'(a)') ' min[imum] <minval> {y|n} ........ set minimum to <minval> (& zero pixels if "y")'
         write (iotty,'(a)') ' maxs[tdev] <nsigm> <radpix> ..... set pixels with values > nsigm*stdev to zero'
         write (iotty,'(a)') ' key[word] <keyword> <param> ..... set FITS header <keyword> to value <param>'
         write (iotty,'(a)') ' cle[an]{4|8}[sk] <val> <par> .... remove 4- or 8-connected clusters of pixels'
         write (iotty,'(a)') ' cli[p]{4|8}{s|f} <val> <par> .... cut 4- or 8-connected sources from filaments'
         write (iotty,'(a)') ' sha[pe] <width> 0[-<R>] ......... reshape by Gaussian <width> with edge at <R>'
         write (iotty,'(a)') ' sha[pe] <width> <slope>[-<R>] ... reshape by Gaussian <width>, <slope>, [edge <R>]'
         write (iotty,'(a)') ' sha[pe] <width> <slope>+[<R>] ... reshape by Plummer <width>, <beta>, [edge <R>]'
         write (iotty,'(a)') ' sha[pe] <width> +<slope>+[<R>] .. reshape by Plummer <width>, <gamma>, [edge <R>]'
         write (iotty,'(a)') ' sha[pe] <width> <slope>+[<R>!] .. reshape image by geometry fun with <beta>'
         write (iotty,'(a)') ' sha[pe] <width> +<slope>+[<R>!] . reshape image by geometry fun with <gamma>'
         write (iotty,'(a)') ' sed[tdust] <Td> <wave> .......... scale with SED at <Td> at <wave> (µm)'
         write (iotty,'( )')   
         write (iotty,'(a)') ' Conversion of intensity units (first ~> second):'
         write (iotty,'( )')   
         write (iotty,'(a)') ' Jy/pixel~MJy/sr .............. convert units'
         write (iotty,'(a)') ' MJy/sr~Jy/pixel .............. convert units'
         write (iotty,'(a)') ' MJy/sr~Jy/beam    [-]<beam> .. convert units; <beam> is FWHM (or area if < 0)'
         write (iotty,'(a)') ' Jy/beam~MJy/sr    [-]<beam> .. convert units; <beam> is FWHM (or area if < 0)'
         write (iotty,'(a)') ' MJy/sr~mJy/beam   [-]<beam> .. convert units; <beam> is FWHM (or area if < 0)'
         write (iotty,'(a)') ' mJy/beam~MJy/sr   [-]<beam> .. convert units; <beam> is FWHM (or area if < 0)'
         write (iotty,'(a)') ' Jy/pixel~Jy/beam  [-]<beam> .. convert units; <beam> is FWHM (or area if < 0)'
         write (iotty,'(a)') ' Jy/beam~Jy/pixel  [-]<beam> .. convert units; <beam> is FWHM (or area if < 0)'
         write (iotty,'(a)') ' Jy/pixel~mJy/beam [-]<beam> .. convert units; <beam> is FWHM (or area if < 0)'
         write (iotty,'(a)') ' mJy/beam~Jy/pixel [-]<beam> .. convert units; <beam> is FWHM (or area if < 0)'
         write (iotty,'( )')   
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FITS KEYWORDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')   
         write (iotty,'(a)') ' Keywords recognized: HEADER CREATOR OBJECT CTYPE1 CTYPE2 CRPIX1 CRPIX2'
         write (iotty,'(a)') ' CRVAL1 CRVAL2 CROTA1 CROTA2 CDELT1 CDELT2 CD1_1 CD1_2 CD2_1 CD2_2 RA DEC'
         write (iotty,'(a)') ' NAXIS1 NAXIS2 BZERO BSCALE BLANK DATAMIN DATAMAX EQUINOX WAVE HISTORY'
         write (iotty,'( )')   
         write (iotty,'(a)') ' NOTE1: CDELT1 and CDELT2 must be given in arcsec.'  
         write (iotty,'(a)') ' NOTE2: beamsize must be specified in arcsec.'  
         write (iotty,'(a)') ' NOTE3: beamarea must be specified in arcsec^2.'  
         write (iotty,'(a)')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         stop 99
       endif  

! Check if the input file exists and open it.

       fnlen1 = lastc ( filename )
       inquire ( file=filename(1:fnlen1), exist=lfnama )
       if (.not.lfnama .and. filename(max(fnlen1-4,1):fnlen1) .ne. '.fits') then
         filename = filename(1:fnlen1)//'.fits'
         fnlen1 = lastc ( filename )
         inquire ( file=filename(1:fnlen1), exist=lfnama )
       endif
       isp1 = lasts ( filename ) + 1
       outname = filename(isp1:fnlen1-5)//'.m.fits'
       lnoname = .false.
       if (arg2(1:2) .eq. '-o') then
         outname = arg3
       else if (arg3(1:2) .eq. '-o') then
         outname = arg4
       else if (arg4(1:2) .eq. '-o') then
         outname = arg5
       else if (arg5(1:2) .eq. '-o') then
         outname = arg6
       else
         lnoname = .true.
       endif
       ion = lastc ( outname )
       if (outname(max(1,ion-4):ion) .ne. '.fits') then
         outname = outname(1:ion)//'.fits'
         ion = lastc ( outname )
       endif
       
!!       if (l1value) then
!!         if (cverbose .eq. '-verb1') write (iotty,'(/a)') ' MODFITS: '//maction(1:iact)//' '//cvalue(1:icva)//param(1:ipar)
!!     &                                                  //filename(1:fnlen1-5)
!!       endif
!!       if (l4values) then
!!         if (cverbose .eq. '-verb1') write (iotty,'(/a)') ' MODFITS: '//maction(1:iact)//' '//cvalue1(1:icva1)//cvalue2(1:icva2)
!!     &                                                  //cvalue3(1:icva3)//cvalue4(1:icva4)//param(1:ipar)//filename(1:fnlen1-5)
!!       endif
       if (.not.lfnama .and. index ( maction(1:iact), '~' ) .eq. 0) then
         write (iotty,'(/a)') '   MODFITS: ERROR: File '''//filename(1:fnlen1-5)//''' not found.'
         stop 1
       endif
       if (cverbose .eq. '-verb2') write (iotty,'()')
       
       dx = 1.0d-10
       dy = 1.0d-10
       rcore = 0.0d0
                                   
! Determine numbers of pixels in the 1st FITS image.
! Signal GETFITSHEAD to ignore bad values of CDELT1, CDELT2, if correcting header.

       if (maction(1:iact) .eq. 'keyword' .or. maction(1:iact) .eq. 'key') then
         dx = 1.0d100
         dy = 1.0d100
       endif

       if (maction(1:iact) .eq. 'border' .or. maction(1:iact) .eq. 'bor') then
         if (l1value) then
           naddpixmax = nint ( abs ( value ) )
         endif
         if (l4values) then
           naddpixmax = max ( abs ( nint ( value1 ) ) + abs ( nint ( value2 ) ), abs ( nint ( value3 ) ) + abs ( nint ( value4 ) ) )
         endif
       endif
       if (maction(1:iact) .eq. 'expand' .or. maction(1:iact) .eq. 'shrink') then
         naddpixmax = nint ( value )
       endif
       
       if (lfnama) then
       
         call getfitshead ( filename(1:fnlen1), nx, ny, dx, dy, bunit )

         inmx = max ( nx, ny ) + max ( 2 * naddpixmax, 0 )
       
         allocate ( funa(inmx,inmx), funb(inmx,inmx), level(nlevelmax), stat=irc )
         
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (11).'
           stop 11
         endif

         write (cx,'(i7)') nx
         write (cy,'(i7)') ny
         ncx = int ( log10 ( dble ( nx ) ) ) + 1
         ncy = int ( log10 ( dble ( ny ) ) ) + 1

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Reading ('//cx(7-ncx+1:7)//' x '//cy(7-ncy+1:7)//') '''
     &                                                 //filename(1:fnlen1)//''''

         call rfits ( inmx, nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, funa, dx, dy, object, ra, dec
     &              , filename(1:fnlen1), iotty, 0, creator, beam, funmin, funmax, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history, cverbose, rcore )
         
         if (cverbose .eq. '-verb2') write (iotty,'(a,1x,2(1pe14.7))') '   Minmax values in the image:', funmin, funmax
         changed = .false.
         nx1 = 0
         nx2 = 0
         ny1 = 0
         ny2 = 0
         sfact = 1.0d0
         const = 0.0d0
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j)
           enddo
         enddo
       endif

       if (ltiming) then
         if (maction(1:5) .eq. 'clip' .or. maction(1:3) .eq. 'cli') then
           inquire ( file='+timing.modfits', exist=ltime )
           open ( 60, file='+timing.modfits', status='unknown' )
           if (ltime) then
             do 
               read (60,'()',end=61)
             enddo
  61         backspace ( 60 )
           endif
           write (60,'(/a)') ' MODFITS started: '//cdate//' '//ctime//' Image: '//filename(1:fnlen1)
           write (60,'(a)') ' ------------------------------------------------------'
           cput_modfits = 0.0d0
           walt_modfits = 0.0d0
           cpu_modfits = timer ( 'cpu', 0.0d0 )
           wal_modfits = timer ( 'wal', 0.0d0 )
         endif
       endif
       
       fJypixel = 1.0d6 * dx * dy / sq_arcsecs_per_sterad

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'keyword' .or. maction(1:iact) .eq. 'key') then

         if (cverbose .eq. '-verb2')
     &       write (iotty,'(a)') '   Changing image header keyword '//cvalue(1:icva)//'to '//param(1:ipar)
  
         if (cvalue .eq. 'CREATOR') then
           creator = param(1:80)
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'OBJECT') then
           object = param(1:80)
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CTYPE1') then
           ctype1 = param(1:80)
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CTYPE2') then
           ctype2 = param(1:80)
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'BUNIT') then
           bunit = param(1:80)
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CRPIX1') then
           read (param,*) crpix1
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CRPIX2') then
           read (param,*) crpix2
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CRVAL1') then
           read (param,*) crval1
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CRVAL2') then
           read (param,*) crval2
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CROTA1') then
           read (param,*) crota1
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CROTA2') then
           read (param,*) crota2
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'RA') then
           read (param,*) ra
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'DEC') then
           read (param,*) dec
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'NAXIS1') then
           read (param,*) nx
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'NAXIS2') then
           read (param,*) ny
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CDELT1') then
           read (param,*) dx
           dx = abs ( dx )
!!!           if (dx .lt. 0.03d0) dx = dx * 3600.0d0
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CDELT2') then
           read (param,*) dy
!!!           if (abs ( dy ) .lt. 0.03d0) dy = dy * 3600.0d0
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CD1_1') then
           read (param,*) cd11
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CD1_2') then
           read (param,*) cd12
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CD2_1') then
           read (param,*) cd21
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'CD2_2') then
           read (param,*) cd22
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'DATAMIN') then
           read (param,*) datamin
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'DATAMAX') then
           read (param,*) datamax
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'BZERO') then
           read (param,*) bzero
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'BSCALE') then
           read (param,*) bscale
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'BLANK') then
           read (param,*) blank
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'WAVE') then
           read (param,*) wave
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'EQUINOX') then
           read (param,*) equinox
           changed = .true.
           goto 77
         endif
         if (cvalue .eq. 'HISTORY') then
           history = param(1:80)
           changed = .true.
           goto 77
         endif
 77      continue
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'absolute' .or. maction(1:iact) .eq. 'abs') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Taking absolute value of the image'
!!         thresh = value
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
!!             if (funa(i,j) .lt. thresh) funb(i,j) = - funa(i,j)
             funb(i,j) = abs ( funa(i,j) )
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.abs.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'negate' .or. maction(1:iact) .eq. 'neg') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Negating image (multiplying by -1.0)'
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             funb(i,j) = - funa(i,j)
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.neg.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
       
       if (maction(1:iact) .eq. 'nanzero' .or. maction(1:iact) .eq. 'nan') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (.not.isnan ( funa(i,j) )) then  
               if (funa(i,j) .gt. funmax) funmax = funa(i,j)
               if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             endif  
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Converting all NaNs to zero'
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) )) then  
               funb(i,j) = 0.0d0
             endif
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.nan.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'add' .or. maction(1:iact) .eq. 'addz' .or. maction(1:iact) .eq. 'addp') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (maction(1:iact) .ne. 'addp') then
           if (cverbose .eq. '-verb2') write (iotty,'(''   Adding a constant value of'',1pe11.4)') value
         else
           if (cverbose .eq. '-verb2') write (iotty,'(''   Creating a plane with liner gradient in direction'',i2)') 
     &                                               nint ( value )
         endif
         const = value
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (maction(1:iact) .ne. 'addp') then
               funb(i,j) = funa(i,j) + const
             else                            !<-- Make a surface with a linear gradient in 8 main directions.
               if (nint ( const ) .eq. 1) funb(i,j) = dble ( i - 1 ) + 1.0d0
               if (nint ( const ) .eq. 2) funb(i,j) = dble ( i - 1 ) + dble ( j - 1 ) + 1.0d0
               if (nint ( const ) .eq. 3) funb(i,j) = dble ( j - 1 ) + 1.0d0
               if (nint ( const ) .eq. 4) funb(i,j) = dble ( nx - i ) + dble ( j - 1 ) + 1.0d0
               if (nint ( const ) .eq. 5) funb(i,j) = dble ( nx - i ) + 1.0d0
               if (nint ( const ) .eq. 6) funb(i,j) = dble ( nx - i ) + dble ( ny - j ) + 1.0d0
               if (nint ( const ) .eq. 7) funb(i,j) = dble ( ny - j ) + 1.0d0
               if (nint ( const ) .eq. 8) funb(i,j) = dble ( i - 1 ) + dble ( ny - j ) + 1.0d0
             endif
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           if (maction(1:iact) .ne. 'addp') then
             outname = filename(isp1:fnlen1-5)//'.add'//cvalue(1:icva-1)//'.fits'
           else
             outname = filename(isp1:fnlen1-5)//'.plane'//cvalue(1:icva-1)//'.fits'
           endif
         endif   
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'subtract' .or. maction(1:iact) .eq. 'sub') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(''   Subtracting a constant value of '',1pe11.4)') value

         const = value
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j) - const
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.sub'//cvalue(1:icva-1)//'.fits'
         endif   
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'multiply' .or. maction(1:iact) .eq. 'mul') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (.not.isnan ( funa(i,j) )) then  
               if (funa(i,j) .gt. funmax) funmax = funa(i,j)
               if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             endif
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(''   Multiplying image by a factor of '',1pe11.4)') value

         sfact = value
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (.not.isnan ( funa(i,j) )) then  
               funb(i,j) = funa(i,j) * sfact
               if (funb(i,j) .gt. funmax) funmax = funb(i,j)
               if (funb(i,j) .lt. funmin) funmin = funb(i,j)
             else
               funb(i,j) = funa(i,j)
             endif
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.mul'//cvalue(1:icva-1)//'.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
       
       if (maction(1:iact) .eq. 'surfden' .or. maction(1:iact) .eq. 'col') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             funb(i,j) = 0.0d0
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') then
           write (iotty,'(a)') '   Computing surface density factorization map from an image at '//cvalue(1:icva)//'µm'
           write (iotty,'(a)') '   Surface densities will be produced from the factors when multiplied by a real image'
         endif
         wave = value
         funmin =  1.0d+30
         funmax = -1.0d+30

! It is useful to normalize the resulting image to have values corresponding to surface densities.
! Inu = Bnu(Td) * opacity * dust2gas * mu * mH * NH2  ~>  NH2 = Inu / (Bnu(Td) * opacity * dus2gas * mu * mH)
! Inu / Bnu(Td) = image1 / Bnu(image2) = opacity * dus2gas * mu * mH * NH2

         frequency = speedolight * 1.0d4 / wave
         opacity = opacity0 * (frequency / frequency0)**beta
         factor = 1.0d-17 / (opacity * dust2gas * muH2 * amu)
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. 3.0d0) then
               call planck ( funa(i,j), 1.0d4 / wave, bbody2 )
               funb(i,j) = factor / bbody2
             endif
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         bunit = '(H2/cm^2)/(MJy/sr)'
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.cdensfactor.'//cvalue(1:icva-1)//'um.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
       
       if (maction(1:iact) .eq. 'shape' .or. maction(1:iact) .eq. 'sha') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             funb(i,j) = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         beam = value
         iplu = 0
         icom = 0
         ilen = 0
         iexc = 0
         expo = 0.0d0
         ipar = lastc ( param )
         Rout = 0.0d0
         
         if (param(1:1) .ne. '/') then  !<-- to avoid a weird special case: reading a slash would give 0...

! Plummer function (includes power-law envelope) and optional finite outer radius.

           iplu = index ( param(1:ipar), '+' )
           if (iplu .eq. 1) then
             cwhat = 'slope=gamma'
             param = param(2:ipar)
             ipar = ipar - 1
           else
             cwhat = 'slope=beta'
           endif
           iplu = index ( param(1:ipar), '+' )
           if (iplu .gt. 0) then
             cha = '+'
             read (param(1:iplu-1),*,err=444) expo
             ilen = lastc ( param(iplu:ipar) )
             if (ilen .gt. 1) then
               iexc = index ( param(iplu+1:ipar), '!' )
               if (iexc .gt. 0) ipar = ipar - 1                 
               read (param(iplu+1:ipar),*,err=444) Rout
               if (iexc .gt. 0) ipar = ipar + 1
             endif
           else
           
! Gaussian function with optional power-law envelope fitted requiring matching derivatives and optional finite outer radius.
           
             icom = index ( param(1:ipar), '-' )
             if (icom .gt. 0) then
               cha = '-'
               read (param(1:icom-1),*,err=444) expo
               ilen = lastc ( param(icom:ipar) )
               if (ilen .gt. 1) then
                 iexc = index ( param(icom+1:ipar), '!' )
                 if (iexc .gt. 0) ipar = ipar - 1
                 read (param(icom+1:ipar),*,err=444) Rout
                 if (iexc .gt. 0) ipar = ipar + 1
               endif
             else
               cha = '-'
               read (param(1:ipar),*,err=444) expo
             endif
           endif
         endif
         goto 5555
 444     continue         
           write (iotty,'(a)') '   MODFITS: ERROR: Trouble reading EXPO from '//param(1:ipar)
           stop 5
 5555    continue
         if (expo .lt. -almostzero) then
           write (iotty,'(a)') '   MODFITS: ERROR: Shape slope is negative, but must be positive.'
           stop 6
         endif
         write (cexpo,'(f6.2)',err=555) expo
         if (cexpo(1:1) .eq. ' ') cexpo = cexpo(2:)
         if (cexpo(1:1) .eq. ' ') cexpo = cexpo(2:)
         if (cexpo(1:1) .eq. ' ') cexpo = cexpo(2:)
         isl = lastc ( cexpo )
         goto 6666
 555     continue         
           isl = lastc ( cexpo )
           write (iotty,'(a)') '   MODFITS: ERROR: Trouble writing SLOPE to '//cexpo(1:isl)
           stop 7
 6666    continue
 
         if (cverbose .eq. '-verb2') then
           if (iexc .eq. 0) then
             if (iplu .gt. 0) then
               if (cwhat .eq. 'slope=beta') then
                 if (Rout .lt. almostzero) then
                   write (iotty,'(a)') '   Reshaping image by Plummer (h='//cvalue(1:icva-1)//' as, beta='//cexpo(1:isl)//')'
                 else
                   write (iotty,'(a)') '   Reshaping image by Plummer (h='//cvalue(1:icva-1)//' as, beta='//cexpo(1:isl)
     &                               //', Rout='//param(iplu+1:ipar)//' as)'
                 endif
               else
                 if (Rout .lt. almostzero) then
                   write (iotty,'(a)') '   Reshaping image by Plummer (w='//cvalue(1:icva-1)//' as, gamma='//cexpo(1:isl)//')'
                 else
                   write (iotty,'(a)') '   Reshaping image by Plummer (w='//cvalue(1:icva-1)//' as, gamma='//cexpo(1:isl)
     &                               //', Rout='//param(iplu+1:ipar)//' as)'
                 endif
               endif
             endif
             if (iplu .eq. 0 .and. expo .gt. almostzero) then
               if (cwhat .eq. 'slope=beta') then
                 if (Rout .lt. almostzero) then
                   write (iotty,'(a)') '   Reshaping image by Gausspl ('//cvalue(1:icva-1)//' as, beta='//cexpo(1:isl)//')'
                 else
                   write (iotty,'(a)') '   Reshaping image by Gausspl ('//cvalue(1:icva-1)//' as, beta='//cexpo(1:isl)
     &                               //', Rout='//param(icom+1:ipar)//' as)'
                 endif
               else
                 if (Rout .lt. almostzero) then
                   write (iotty,'(a)') '   Reshaping image by Gausspl ('//cvalue(1:icva-1)//' as, gamma='//cexpo(1:isl)//')'
                 else
                   write (iotty,'(a)') '   Reshaping image by Gausspl ('//cvalue(1:icva-1)//' as, gamma='//cexpo(1:isl)
     &                               //', Rout='//param(icom+1:ipar)//' as)'
                 endif
               endif
             endif
             if (iplu .eq. 0 .and. expo .le. almostzero) then
               if (Rout .lt. almostzero) then
                 write (iotty,'(a)') '   Reshaping image by Gauss ('//cvalue(1:icva-1)//' as)'
               else
                 write (iotty,'(a)') '   Reshaping image by Gauss ('//cvalue(1:icva-1)
     &                             //', Rout='//param(icom+1:ipar)//' as)'
               endif
             endif
           else
             if (iplu .gt. 0) then
               if (cwhat .eq. 'slope=beta') then
                 write (iotty,'(a)') '   Reshaping image by finite geometry function (beta='//cexpo(1:isl)
     &                             //', Rout='//param(iplu+1:ipar)//' as)'
               else
                 write (iotty,'(a)') '   Reshaping image by finite geometry function (gamma='//cexpo(1:isl)
     &                             //', Rout='//param(iplu+1:ipar)//' as)'
               endif
             elseif (icom .gt. 0) then
               write (iotty,'(a)') '   Reshaping image by finite geometry function (slope='//cexpo(1:isl)
     &                           //', Rout='//param(icom+1:ipar)//' as)'
             else
               write (iotty,'(a)') '   Reshaping image by finite geometry function (Rout='//param(icom+1:ipar)//' as)'
             endif
           endif
         endif

         if (Rout .gt. almostzero) then

           allocate ( worker(nx,ny), stat=irc )
  
           if (irc .ne. 0) then
             write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (12).'
             stop 112
           endif
  
           do j=1,ny
             do i=1,nx
               worker(i,j) = funa(i,j)
             enddo
           enddo
         endif
         
         hwhm = beam / 2.0d0
         sigm2m = hwhm**2
         sigm2g = hwhm**2 / log ( 4.0d0 )
         sigm2g2 = 2.0d0 * sigm2g
         Rout2 = Rout**2

! Fit parameters for the Gaussian and power-law functions at radius where they derivatives match.

         if (expo .gt. almostzero) then
           radfit = sqrt ( expo * sigm2g )
           radfitx = expo * sigm2g
           expfit = exp ( -min ( radfit**2 / sigm2g2, 50.0d0 ))
           plwfit = 1.0d0 / max ( radfit, dx )**expo
           fitfact = expfit / plwfit
!!           dplawdr = expo
           gamma = dble ( expo )
           wsize = beam

           if (Rout .gt. almostzero) then

             if (cwhat .eq. 'slope=beta') then

! Using bisection method to solve the function gamma(beta,gamma) for gamma.

               betax = dble ( expo )
               xi = Rout / beam
               gam1 = 0.01d0
               gam2 = 20.0d0

               do k=1,100
                 gamma = (gam1 + gam2) / 2.0d0
                 eps = betax - gamma - 1.529d0 / (1.0d0 + xi**(-0.199d0) * exp(-2.725d0 * (gamma * (xi**(0.03d0 * gamma**(-0.7d0))
     &               + xi**(0.26d0 * xi**(0.03d0) * gamma**(-0.1d0))) * 0.5d0 - 0.319d0))) + 0.541d0

                 if (gam1 .gt. 0.0d0 .and. eps .gt. 0.0d0) then
                   gam1 = gamma
                 else
                   gam2 = gamma
                 endif
!!                 write (*,*) gam1, gamma, gam2, eps
                 if (abs ( eps ) .lt. almostzero) exit
                 if (abs ((abs ( gam1 ) - abs ( gam2 )) * 2.0d0 / (abs ( gam1 ) + abs ( gam2 ))) .lt. 1.0d-7) exit
               enddo

               epsiln = 40.115d0 / (1.0d0 + exp(-0.3782d0 * (betax + 2.542d0)) * xi**(0.012d0))
     &                + 8.21d0 * exp(-0.7497d0 * betax * xi**(0.0475d0)) * xi**(-0.015d0) - 34.104d0
                 
               E = 0.77149d0 * (1.0d0 - exp ( -((xi - 0.52709d0) / 0.7156d0)**(-0.9095d0) )) + 0.0026857d0
               F = -0.31586d0 * xi**(-1.9388d0) + 0.57344d0 * xi**(-0.96778d0) + 6.5472d0 * xi**(0.14471d0) - 5.1551d0
               G = 1.0811d0 * xi**(-0.62466d0) - 1.4375d0 * xi**(-1.0203d0) + 4.0626d0 * xi**(0.00955d0) - 3.1165d0
               S = 0.034613d0 * exp ( -0.014394d0 * xi ) + 0.036328d0 * exp ( -0.27696d0 * xi ) + 2.1537d0 * exp ( -5.104d0 * xi )
     &           + 0.94275d0
               Z = 0.26355d0 * exp ( -0.059635d0 * xi ) + 8.8497d0 * exp ( -4.824d0 * xi ) - 645.57d0 * exp ( -12.634d0 * xi )
     &           + 0.24014d0
              
               Has = beam / (S + (E - S) * ( 1.0d0 + (betax / G)**F )**(-Z))
              
               wsize = Has * (235.7d0 * exp(-20.0d0 * betax * xi**(-0.2d0)) + 5.0d-5 * xi**(0.5d0) * beta**(-6.0d0 * xi**(0.21d0))
     &               + 2.878d0 * exp(-1.069d0 * betax * xi**(0.22d0 * betax)) + 113.0d0 * exp(-10.88d0 * betax * xi**(-0.2d0))
     &               + 1.022d0 * xi**(-0.0077d0))
               
               P = 0.93537d0 * (1.0d0 - exp ( -((xi - 1.0d0) / 161.325d0)**(0.64632d0) )) + 0.59532d0
               Q = 10.843d0 * (1.0d0 - exp ( -((xi - 0.51536d0) / 86.911d0)**(0.415d0) )) + 0.46752d0
               V = -4.2115d-2 * exp ( -9.8242d-2 * xi ) - 6.915d-1 * exp ( -7.56d-1 * xi ) - 2.8652d0 * exp ( -2.5807d0 * xi ) + 1.1694d0
               T = 1.1112d0 * (xi + 27.889d0)**(0.97013d0) - 28.0154d0
               U = 0.51224d0 * xi**(-2.2126d0) - 1.3366d0 * xi**(-1.1266d0) + 1.245d0 * xi**(-0.5599d0) + 0.13369d0
     
               chi = T + (P - T) * ( 1.0d0 + (betax / V)**Q )**(-U)
              
               write (*,'(a)') '   '//cwhat
               write (*,'(a,1pe12.5)') '      beta:', betax
               write (*,'(a,1pe12.5)') '     gamma:', gamma
               write (*,'(a,1pe12.5)') '    epsiln:', epsiln
               write (*,'(a,1pe12.5)') '   correct:', correct
               write (*,'(a,1pe12.5)') '         h:', beam
               write (*,'(a,1pe12.5)') '         w:', wsize
               write (*,'(a,1pe12.5)') '         H:', Has
               write (*,'(a,1pe12.5)') '       w/h:', wsize / beam
               write (*,'(a,1pe12.5)') '       H/h:', Has / beam
               write (*,'(a,1pe12.5)') '    xi=R/h:', xi
               write (*,'(a,1pe12.5)') '   chi=R/H:', chi
               write (*,'(a,1pe12.5)') '       R/H:', Rout / Has
               write (*,'(a,1pe12.5)') '       R/w:', Rout / wsize
               write (*,'(a,1pe12.5)') '       w/H:', wsize / Has
             else
               betax = 0.0d0
               epsiln = 0.0d0
               Has = 0.0d0
               xi = 0.0d0
               chi = 0.0d0
!!               betax = gamma + 1.529d0 / (1.0d0 + xi**(-0.199d0) * exp(-2.725d0 * (gamma * (xi**(0.03d0 * gamma**(-0.7d0))
!!     &               + xi**(0.26d0 * xi**(0.03d0) * gamma**(-0.1d0))) * 0.5d0 - 0.319d0))) - 0.541d0
             endif
           endif
           correct = (2.0d0**(2.0d0 / gamma) - 1.0d0)
         endif

         isumpct = 0.0d0
         numtotal = dble ( ny )
         numpart = 0.0d0

         if (iexc .eq. 0) then

           if (cverbose .eq. '-verb2') then
             if (iotty .gt. 0) then
               write (iotty,'(a)', advance='no') '   Progress: '
               if (iotty .gt. 0) endfile   ( iotty, err=777 )
 777           continue 
               if (iotty .gt. 0) backspace ( iotty )
             endif       
           endif

           if (Rout .gt. almostzero) then
             radmax = 1.1d0 * Rout
             nw = nint ( radmax / dx )
             nw = min ( nw, max ( nx, ny ) )
           else
             nw = max ( nx, ny )
           endif

           do jj=1,ny
             if (cverbose .eq. '-verb2') then
               numpart = dble ( jj )
               call showprogress ( iotty, 0, lunix, numpart, numtotal, isumpct )
             endif

             do ii=1,nx
               if (funa(ii,jj) .gt. almostzero) then
                 do j=max(jj-nw,1),min(jj+nw,ny)
!!                 do j=1,ny
                   dely = (dble ( j - jj ) * dy)**2
                   do i=max(ii-nw,1),min(ii+nw,nx)
!!                   do i=1,nx
                     delx = (dble ( i - ii ) * dx)**2
                     rad2 = delx + dely
                     rad = sqrt ( rad2 )
                     if (iplu .gt. 0) then
                       funb(i,j) = max ( funb(i,j), funa(ii,jj) * (1.0d0 + correct * (2.0d0 * rad / wsize)**2)**(-gamma / 2.0d0) )
                     else
                       argex = min ( rad2 / sigm2g2, 50.0d0 )
                       gauss = exp ( -argex )
!!                       dgausdr = rad2 / sigm2g
                       if (expo .gt. almostzero) then
                         plaw = 1.0d0 / sqrt ( max ( rad2, dx ) )**expo
                         if (rad2 .lt. radfitx) then
                           funb(i,j) = max ( funb(i,j), funa(ii,jj) * gauss )
                         else
                           funb(i,j) = max ( funb(i,j), funa(ii,jj) * plaw * fitfact )
                         endif
                       else
                         funb(i,j) = max ( funb(i,j), funa(ii,jj) * gauss )
                       endif
                     endif
                   enddo
                 enddo
               endif
             enddo
           enddo
         endif

         if (Rout .gt. almostzero .and. cwhat .eq. 'slope=beta') then

           isumpct = 0.0d0
           numtotal = dble ( ny )
           if (cverbose .eq. '-verb2') then
             if (iotty .gt. 0) then
               write (iotty,'(a)', advance='no') '   Progress: '
               if (iotty .gt. 0) endfile   ( iotty, err=888 )
 888           continue 
               if (iotty .gt. 0) backspace ( iotty )
             endif       
           endif

           radmax = 1.1d0 * Rout
           nw = nint ( radmax / dx )
           nw = min ( nw, max ( nx, ny ) )
         
           do jj=1,ny
             if (cverbose .eq. '-verb2') then
               numpart = dble ( jj )
               call showprogress ( iotty, 0, lunix, numpart, numtotal, isumpct )
             endif

              do ii=1,nx
               if (funa(ii,jj) .gt. almostzero) then
                 do j=max(jj-nw,1),min(jj+nw,ny)
                   dely = (dble ( j - jj ) * dy)**2
                   do i=max(ii-nw,1),min(ii+nw,nx)
                     delx = (dble ( i - ii ) * dx)**2
                     rad2 = delx + dely
                     if (expo .gt. almostzero) then
                       factgeom = sqrt ( 1.0d0 - (sqrt ( min ( rad2, Rout2 ) ) / Rout )**(epsiln) )
                       worker(i,j) = max ( worker(i,j), funa(ii,jj) * factgeom )
                     else
                       dgausdr = rad2 / sigm2g
                       factgeom = sqrt ( 1.0d0 - (sqrt ( min ( rad2, Rout2 ) ) / Rout )**(2) )    !!dgausdr
                       worker(i,j) = max ( worker(i,j), funa(ii,jj) * factgeom )
                     endif
                   enddo
                 enddo
               endif
             enddo
           enddo
         endif

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (Rout .gt. almostzero .and. cwhat .eq. 'slope=beta') then
               if (iexc .eq. 0) then
                 funb(i,j) = funb(i,j) * worker(i,j)
               else
                 funb(i,j) = worker(i,j)
               endif
             endif
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax

         if (Rout .gt. almostzero) then
           deallocate ( worker )
         endif

         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           ipt = index ( cexpo(1:isl), '.' )
           if (ipt .gt. 0) cexpo(ipt:ipt) = 'p'
           ipt = index ( cexpo(ipt+1:isl), '.' ) + ipt
           if (ipt .gt. 0) cexpo(ipt:ipt) = 'p'

           if (Rout .gt. almostzero) then
             if (iplu .gt. 0) crout = param(iplu+1:ipar)
             if (icom .gt. 0) crout = param(icom+1:ipar)
             ilen = lastc ( crout )
             ipt = index ( crout(1:ilen), '.' )
             if (ipt .gt. 0) crout(ipt:ipt) = 'p'
           endif

           if (Rout .lt. almostzero) then
             if (expo .gt. almostzero) then
               outname = filename(isp1:fnlen1-5)//'.sha'//cvalue(1:icva-1)//'x'//cexpo(1:isl)//cha//'.fits'
             else
               outname = filename(isp1:fnlen1-5)//'.sha'//cvalue(1:icva-1)//'.fits'
             endif
           else
             if (expo .gt. almostzero) then
               if (cwhat .eq. 'slope=beta') then
                 outname = filename(isp1:fnlen1-5)//'.sha'//cvalue(1:icva-1)//'x'//cexpo(1:isl)//cha//crout(1:ilen)//'.fits'
                 ilen = lastc ( outname )
                 open (44, file=outname(1:ilen-5)//'.txt', status='unknown')
                 write (44,'(a)') ' '//cwhat
                 write (44,'(a)') '  betax       gamma       epsiln      correct     h           w           H           w/h'
     &                          //'         H/h         xi=R/h      chi=R/H     R/H         R/w         w/H'
                 write (44,'(20(1pe12.5))') betax, gamma, epsiln, correct, beam, wsize, Has, wsize / beam, Has / beam, xi
     &                                    , chi, Rout / Has, Rout / wsize, wsize / Has
                 close (44)
               else
                 outname = filename(isp1:fnlen1-5)//'.sha'//cvalue(1:icva-1)//'x+'//cexpo(1:isl)//cha//crout(1:ilen)//'.fits'
               endif
             else
               outname = filename(isp1:fnlen1-5)//'.sha'//cvalue(1:icva-1)//cha//crout(1:ilen)//'.fits'
             endif
           endif
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

       !!             eps = (beam / Rout) / (0.9d0 + 1.0d0 / dble ( expo )) / (dble ( expo ) + 1.0d0)
       !!             expoeps = dble ( expo ) + eps

            !!OLD               gam1 = 0.001d0
            !!OLD               gam2 = betax
            !!OLD               gam1 = betax - gam1 - 1.0d0 + 0.7485d0 * exp ( -5.19d0 * gam1 )
            !!OLD               gam2 = betax - gam2 - 1.0d0 + 0.7485d0 * exp ( -5.19d0 * gam2 )
            !!OLD               
            !!OLD               do k=1,100
            !!OLD                 gamma = (gam1 + gam2) / 2.0d0
            !!OLD                 eps = betax - gamma - 1.0d0 + 0.7485d0 * exp ( -5.19d0 * gamma )
            !!OLD                 if (gam1 .gt. 0.0d0 .and. eps .lt. 0.0d0) then
            !!OLD                   gam1 = gamma
            !!OLD                 else
            !!OLD                   gam2 = gamma
            !!OLD                 endif
            !!OLD                 if (abs ( eps ) .lt. almostzero) exit
            !!OLD                 if (abs ((abs ( gam1 ) - abs ( gam2 )) * 2.0d0 / (abs ( gam1 ) + abs ( gam2 ))) .lt. 1.0d-7) exit
            !!OLD               enddo
       !!OLD               betax = gamma + 1.0d0 - 0.7485d0 * exp ( -5.19d0 * gamma )
       !!OLD             epsiln = 1.54d0 * (-betax + 0.855d0 * (xi)**(0.025d0))**2
       !!OLD     &              / (1.1d0 * (betax - 0.3d0)**(1.17d0) + 0.035d0 * betax**2 + 0.5d0) + 0.33d0
       !!OLD     &              + 1.6d0 * exp ( -(xi**(1.02d0) - 0.33d0)**(0.18d0) )**(0.86d0) - 0.265d0
       !!             wsize = Has * (235.7d0 * exp(-20.0d0 * betax * xi**(-0.2d0)) + 2.878d0 * exp(-1.069d0 * betax * xi**(0.22d0 * betax))
       !!     &             + 113.0d0 * exp(-10.88d0 * betax * xi**(-0.2d0)) + 1.022d0 * xi**(-0.0077d0))
             
       !!OLD             wsize = Has * (0.6975d0 * (betax**(-14) + 0.6d0 * betax**(-13) + 1.3d0 * betax**(-12) + betax**(-11) + betax**(-10)
       !!OLD     &           + betax**(-9) + betax**(-8) + 2.0d0 * betax**(-7) + 2.0d0 * betax**(-6) + 3.0d0 * betax**(-5) + 1.3d0 * betax**(-4)
       !!OLD     &           + 0.6d0 * xi**(-1) * (betax**(-3) + betax**(-2))) * (xi**(-1.8d0) + 0.25d0) + 1.0d0)

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
       
       if (maction(1:iact) .eq. 'bbody' .or. maction(1:iact) .eq. 'bbo') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             funb(i,j) = 0.0d0
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Converting Td image to Bnu(Td) at '//cvalue(1:icva)//'µm'
         
         wave = value
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             call planck ( min ( max ( funa(i,j), 3.0d0 ), 500.0d0 ), 1.0d4 / wave, funb(i,j) )
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.bbody.'//cvalue(1:icva-1)//'um.fits'
         endif 
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
       
       if (maction(1:iact) .eq. 'sedtdust' .or. maction(1:iact) .eq. 'sed') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') 
     &       write (iotty,'(a)') '   Computing pixel SED for uniform Tdust of '//cvalue(1:icva)//'K at '//param(1:ipar)//'µm' 

         tdust = value
         read (param(1:ipar),*,err=434) wave
         goto 5454
 434     continue         
           write (iotty,'(/a)') '   MODFITS: ERROR: Trouble reading WAVE from operand: '//param(1:ipar)
           stop 5
 5454    continue
         funmin =  1.0d+30
         funmax = -1.0d+30
         if (tdust .gt. 0.0d0) then
           call planck ( tdust, 1.0d4 / wave, bbody1 )
           call planck ( tdust, 1.0d4 / 100.0d0, bbody2 )
         else
           tdmin = 15.0d0
           tdmax = 20.0d0
           dtdi = (tdmax - tdmin) * 0.5d0 / dble ( nx - 1 )
           dtdj = (tdmax - tdmin) * 0.5d0 / dble ( ny - 1 )                 
         endif                                                      
                                                                    
! Scaling with frequency squared is to account for dust emissivity (lambda^-2 at long waves).
                                                                    
         do j=1,ny                                                  
           do i=1,nx                                                
             if (tdust .gt. 0.0d0) then                             
               funb(i,j) = funa(i,j) * (bbody1 / bbody2) * (100.0d0 / wave)**2
             else                                                   
               tdc = tdmin + dble ( nx - i ) * dtdi + dble ( j - 1 ) * dtdj
               call planck ( tdc, 1.0d4 / wave, bbody1 )
               call planck ( tdc, 1.0d4 / 100.0d0, bbody2 )
               funb(i,j) = funa(i,j) * (bbody1 / bbody2) * (100.0d0 / wave)**2
             endif
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.sed.Td'//cvalue(1:icva-1)//'K.'//param(1:ipar-1)//'um.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (index ( maction(1:iact), '~' ) .gt. 0) then
         
         if (lfnama) then
           funmin =  1.0d+30
           funmax = -1.0d+30
           do j=1,ny
             do i=1,nx
               if (funa(i,j) .gt. funmax) funmax = funa(i,j)
               if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             enddo
           enddo
           if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
           if (cverbose .eq. '-verb2') write (iotty,*) funmax
           if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
           if (cverbose .eq. '-verb2') write (iotty,*) funmin
         endif

         if (index ( maction(1:iact), 'beam' ) .gt. 0) then
           if (ia2 .eq. 0) then
             write (iotty,'(''   Enter beam size (arcsec): '',$)')
             read (iotty,*) value
           else
             if (ia2 .gt. 0) read (arg2(1:ia2),*,err=46) value
             goto 47
  46         continue
             if (cverbose .eq. '-verb2') write (iotty,'(/a)') '   MODFITS: ERROR: Trouble reading beam size from command line.'
             stop 99
  47         continue
           endif
           if (cverbose .eq. '-verb2')
     &         write (iotty,'(''   Convert image physical units: '',a)') maction(1:iact)//' (beam='//cvalue(1:icva-1)//')'
         else
           if (cverbose .eq. '-verb2') write (iotty,'(''   Converting image physical units: '',a)') maction(1:iact)
         endif
         beamsize = 0.0d0
         beamarea = 0.0d0
         if (value .gt. 0.0d0) then
           beamsize = value
         else
           beamarea = abs ( value )
         endif
         if (beamsize .gt. almostzero) then
           sq_arcsecs_per_beam_area = pi * (beamsize / 2.0)**2 / log ( 2.0 )
           pixels_per_beam_area = pi * (beamsize / 2.0)**2 / (dx * dy) / log ( 2.0 )
         else
           sq_arcsecs_per_beam_area = beamarea
           pixels_per_beam_area = beamarea
         endif

         factMJysr2Jybeam = 1.0d6 * sq_arcsecs_per_beam_area / sq_arcsecs_per_sterad

         if (maction(1:iact) .eq. 'MJy/sr~Jy/beam') conversionfactor = factMJysr2Jybeam
         if (maction(1:iact) .eq. 'Jy/beam~MJy/sr') conversionfactor = 1.0d0 / factMJysr2Jybeam

         factMJysr2mJybeam = 1.0d9 * sq_arcsecs_per_beam_area / sq_arcsecs_per_sterad

         if (maction(1:iact) .eq. 'MJy/sr~mJy/beam') conversionfactor = factMJysr2mJybeam
         if (maction(1:iact) .eq. 'mJy/beam~MJy/sr') conversionfactor = 1.0d0 / factMJysr2mJybeam

         factJypixel2Jybeam = pixels_per_beam_area

         if (maction(1:iact) .eq. 'Jy/pixel~Jy/beam') conversionfactor = factJypixel2Jybeam
         if (maction(1:iact) .eq. 'Jy/beam~Jy/pixel') conversionfactor = 1.0d0 / factJypixel2Jybeam

         factJypixel2mJybeam = 1.0d3 * pixels_per_beam_area

         if (maction(1:iact) .eq. 'Jy/pixel~mJy/beam') conversionfactor = factJypixel2mJybeam
         if (maction(1:iact) .eq. 'mJy/beam~Jy/pixel') conversionfactor = 1.0d0 / factJypixel2mJybeam
                                                    
         factJypixel2MJysr = sq_arcsecs_per_sterad / (dx * dy) / 1.0d6

         if (maction(1:iact) .eq. 'Jy/pixel~MJy/sr') conversionfactor = factJypixel2MJysr
         if (maction(1:iact) .eq. 'MJy/sr~Jy/pixel') conversionfactor = 1.0d0 / factJypixel2MJysr
                                                    
         if (cverbose .eq. '-verb1') write (iotty,*)
         if (cverbose .ne. '-verb0') write (iotty,*) '  conversion factor:', conversionfactor

         write (ccfact,'(1pe10.3)') conversionfactor

         if (beamsize .gt. almostzero) then
           write (cbsize,'(1pe10.3)') beamsize
           history = maction(1:iact)//': beam size:'//cbsize//'; conversion factor:'//ccfact
         elseif (beamarea .gt. almostzero) then
           write (cbarea,'(1pe10.3)') beamarea
           history = maction(1:iact)//': beam area:'//cbarea//'; conversion factor:'//ccfact
         else
           write (cpixel,'(1pe10.3)') sqrt ( dx * dy )
           history = maction(1:iact)//': pixel size:'//cpixel//'; conversion factor:'//ccfact
         endif
         
         if (.not.lfnama .and. index ( maction(1:iact), '~' ) .gt. 0) then
           write (iotty,'()')
!!           write (iotty,'(a)') '   MODFITS: ERROR: File '''//filename(1:fnlen1-5)//''' not found.'
           stop
         endif
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j) * conversionfactor
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'divide' .or. maction(1:iact) .eq. 'div') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(''   Dividing image by a factor of '',1pe11.4)') value
         sfact = value
         if (sfact .eq. 0.0d0) then
           write (iotty,'(/a)') '   MODFITS: ERROR: Division by zero not allowed.'
           stop 99
         endif
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j) / sfact
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.div'//cvalue(1:icva-1)//'.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'square' .or. maction(1:iact) .eq. 'squ') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing square of the image'
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j)**2
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.squ.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'power' .or. maction(1:iact) .eq. 'pow') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(''   Computing power of the image:'',i3)') nint ( value )
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j)**nint ( value )
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.pow'//cvalue(1:icva-1)//'.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'thicken' .or. maction(1:iact) .eq. 'thi') then
             
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Thickening structures by annexing one boundary pixel'

         do j=1,ny
           jm1 = max (j - 1, 1)
           jp1 = min (j + 1, ny)
           do i=1,nx
             im1 = max (i - 1, 1)
             ip1 = min (i + 1, nx)
             if (funa(i,j) .lt. almostzero) then
               if (funa(im1,j  ) .gt. almostzero) funb(i,j) = funa(im1,j  )
               if (funa(ip1,j  ) .gt. almostzero) funb(i,j) = funa(ip1,j  )
               if (funa(i  ,jm1) .gt. almostzero) funb(i,j) = funa(i  ,jm1)
               if (funa(i  ,jp1) .gt. almostzero) funb(i,j) = funa(i  ,jp1)
               if (funa(im1,jm1) .gt. almostzero) funb(i,j) = funa(im1,jm1)
               if (funa(ip1,jp1) .gt. almostzero) funb(i,j) = funa(ip1,jp1)
               if (funa(im1,jp1) .gt. almostzero) funb(i,j) = funa(im1,jp1)
               if (funa(ip1,jp1) .gt. almostzero) funb(i,j) = funa(ip1,jp1)
             endif
           enddo
         enddo
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.thi.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'contract' .or. maction(1:iact) .eq. 'con') then
             
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Contracting structures by removing one boundary pixel'

         do j=1,ny
           jm1 = max (j - 1, 1)
           jp1 = min (j + 1, ny)
           do i=1,nx
             im1 = max (i - 1, 1)
             ip1 = min (i + 1, nx)
             if (funa(i,j) .gt. almostzero) then
               if (funa(im1,j  ) .lt. almostzero) funb(i,j) = funa(im1,j  )
               if (funa(ip1,j  ) .lt. almostzero) funb(i,j) = funa(ip1,j  )
               if (funa(i  ,jm1) .lt. almostzero) funb(i,j) = funa(i  ,jm1)
               if (funa(i  ,jp1) .lt. almostzero) funb(i,j) = funa(i  ,jp1)
               if (funa(im1,jm1) .lt. almostzero) funb(i,j) = funa(im1,jm1)
               if (funa(ip1,jp1) .lt. almostzero) funb(i,j) = funa(ip1,jp1)
               if (funa(im1,jp1) .lt. almostzero) funb(i,j) = funa(im1,jp1)
               if (funa(ip1,jp1) .lt. almostzero) funb(i,j) = funa(ip1,jp1)
             endif
           enddo
         enddo
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.con.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'merge' .or. maction(1:iact) .eq. 'mer') then
       
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Merging structures separated by 1 pixel'

         do j=1,ny
           jm1 = max (j - 1, 1)
           jp1 = min (j + 1, ny)
           do i=1,nx
             im1 = max (i - 1, 1)
             ip1 = min (i + 1, nx)
             if (funa(i,j) .lt. almostzero) then
               if (funa(im1,j) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero) funb(i,j) = (funa(im1,j) + funa(ip1,j)) / 2.0d0
               if (funa(i,jm1) .gt. almostzero .and. 
     &             funa(i,jp1) .gt. almostzero) funb(i,j) = (funa(i,jm1) + funa(i,jp1)) / 2.0d0
             endif
           enddo
         enddo
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.mer.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:6) .eq. 'sqroot' .or. maction(1:3) .eq. 'sqr') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (maction(1:iact) .eq. 'sqrootm' .or. maction(1:iact) .eq. 'sqrm') then 
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing square root (modified) of the image'
         else
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing square root of the image'
         endif

         if (maction(1:7) .eq. 'sqrootm' .or. maction(1:4) .eq. 'sqrm') then 
!!           if (maction(1:iact) .eq. 'sqrootmp' .or. maction(1:iact) .eq. 'sqrmp') then 
                    
! The critical value of 0.32768 is found by requiring that the first derivative of the power function y=x**0.8 
! is equal to the first derivative of y=x, i.e. that the resulting first derivative is smooth.
! This is easily calculated as x = (1/0.8)**(1/(0.8-1)) = 0.8**5 = 0.32768. In order to have the function
! itself smooth, one needs to shift the power function down by 0.08192 = 0.32768**0.8 - 0.32768 = 0.8**4 - 0.8**5
!
! For function y=x**0.5 it's calculated as x = (1/0.5)**(1/(0.5-1)) = 0.5**2 = 0.25. In order to have the function
! itself smooth, one needs to shift the power function down by delta = 0.25**0.5 - 0.5**2
                    
!!              funfit = 0.8d0**5    !< = 0.32768
             funfit = 0.5d0**2    !< = 0.25
             
             do j=1,ny
               do i=1,nx
                 if (funa(i,j) .ge. funfit) then
                   funb(i,j) = funa(i,j)**(0.5d0) - (0.25d0**(0.5d0) - funfit)
                 else
                   if (funa(i,j) .le. -funfit) then
                     funb(i,j) = -abs ( funa(i,j) )**(0.5d0) + (0.25d0**(0.5d0) - abs ( funfit ))
                   else
                     funb(i,j) = funa(i,j)
                   endif
                 endif
               enddo
             enddo
         else
           do j=1,ny
             do i=1,nx
               if (funa(i,j) .ge. 0.0d0) then
                 funb(i,j) = sqrt ( funa(i,j) )
               else
                 funb(i,j) = 0.0d0
               endif
             enddo
           enddo
         endif
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.sqr.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:5) .eq. 'log10' .or. maction(1:3) .eq. 'log') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (maction(1:iact) .eq. 'log10m' .or. maction(1:iact) .eq. 'logm') then 
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing log10 (modified) of the image'
         else
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing log10 of the image'
         endif

         if (maction(1:iact) .eq. 'log10m' .or. maction(1:iact) .eq. 'logm') then 
           do j=1,ny
             do i=1,nx
               if (funa(i,j) .ge. 1.0d0) then
                 funb(i,j) = log10 ( funa(i,j) ) + 1.0d0
               else
                 if (funa(i,j) .le. -1.0d0) then
                   funb(i,j) = -log10 ( -funa(i,j) ) - 1.0d0
                 else
                   funb(i,j) = funa(i,j)
                 endif
               endif
             enddo
           enddo
         else
           if (funmin .le. 1.0d-10) then 
             do j=1,ny
               do i=1,nx
                 if (funa(i,j) .ge. 1.0d0) then
                   funb(i,j) = max ( log10 ( funa(i,j) ), 0.0d0 )
                 else
                   if (funa(i,j) .ge. 0.0d0) then
                     funb(i,j) = 0.0d0
                   endif
                 endif
               enddo
             enddo
             do j=1,ny
               do i=1,nx
                 if (funa(i,j) .le. -1.0d0) then
                   funb(i,j) = min ( -log10 ( -funa(i,j) ), 0.0d0 )
                 else
                   if (funa(i,j) .le. 0.0d0) then
                     funb(i,j) = 0.0d0
                   endif
                 endif
               enddo
             enddo
           else
             do j=1,ny
               do i=1,nx
                 funb(i,j) = log10 ( funa(i,j) )
               enddo
             enddo
           endif
         endif
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamax = funmax
         datamin = funmin
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.log.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'exponent' .or. maction(1:iact) .eq. 'exp') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing exponent of the image'

         do j=1,ny
           do i=1,nx
             if (funa(i,j) .lt. 88.0d0) then
               funb(i,j) = exp ( funa(i,j) )
             else
               funb(i,j) = 0.0d0  
             endif
           enddo
         enddo
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamax = funmax
         datamin = funmin
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.exp.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:5) .eq. 'layer' .or. maction(1:3) .eq. 'lay') then
 
         sfact = value
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(''   Removing ''''zero layer'''' of +/-'',1pe11.4)') value

         do j=1,ny
           do i=1,nx
             if (funa(i,j) .ge. sfact .or. funa(i,j) .le. -sfact) then 
               funb(i,j) = funa(i,j)
             else
               funb(i,j) = 0.0d0
!!               funb(i,j) = sfact
             endif
           enddo
         enddo
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamax = funmax
         datamin = funmin
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.lay'//cvalue(1:icva-1)//'.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:5) .eq. 'noise' .or. maction(1:3) .eq. 'noi') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(''   Adding Gaussian noise at wavelength'',1pe11.4)') value

! Define random noise at each pixel using a random number generator 'gasdev' with a normal (Gaussian) distribution.
! Re-initialize the random sequence first (idum = -1).

         write (cline,*) nint ( value )
         idum = -1
         read (cline,*) newrandomseq
         srms = 0.0d0

! Inner loop is to make a new independent random sequence for every wavelength, to avoid correlation between the quasi-random
! noise at any wavelength. The initialization consists of multiplying the integer number equal to the wavelength by a uniform 
! random number in the range [0,1]. So many cycles will be skept out of the Gaussian random sequence between the rows of pixels.
         
         do j=1,ny
           newran = int ( newrandomseq * ran1 ( idum ) )
           do i=1,newran
             q1 = gasdev ( idum )
           enddo
           do i=1,nx
             q1 = gasdev ( idum )
             funb(i,j) = funb(i,j) + q1
             srms = srms + q1**2
           enddo
         enddo
         srms = sqrt ( srms / real ( nx * ny ) )
         if (cverbose .eq. '-verb2') write (iotty,'(''     Noise standard deviation: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) srms
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamax = funmax
         datamin = funmin
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.noi.'//cvalue(1:icva-1)//'um.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'getfilaments' .or. maction(1:iact) .eq. 'get') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             funb(i,j) = 0.0d0
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         
         nskw = max ( nint ( value ), 1 )
         if (mod ( nskw, 2 ) .eq. 0) nskw = nskw + 1
         if (cverbose .eq. '-verb2') write (iotty,'(''   Skeletonizing image using my algorithm: '',i1,'' pixels wide'')') nskw

         do j=1,ny
           jm1 = max ( j - 1, 1 )
           jp1 = min ( j + 1, ny )
           do i=1,nx
             if (funa(i,j) .gt. almostzero) then
               im1 = max ( i - 1, 1 )
               ip1 = min ( i + 1, nx )
               numpix4 = 0
               numpix8 = 0
               if (funa(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
               if (funa(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
               if (funa(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
               if (funa(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
               if (funa(im1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
               if (funa(ip1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
               if (funa(im1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
               if (funa(ip1,jm1) .gt. almostzero) numpix8 = numpix8 + 1

               if (numpix4 .eq. 4 .or. numpix8 .ge. 3) then

                 lowintensity = 0.5d0 * funa(i,j)
                 ix1 = 1
                 do ii=i,1,-1
                   if (funa(ii,j) .lt. lowintensity) then
                     ix1 = ii + 1
                     exit
                   endif
                 enddo
                 ix2 = nx
                 do ii=i,nx
                   if (funa(ii,j) .lt. lowintensity) then
                     ix2 = ii - 1
                     exit
                   endif
                 enddo
                 ixlength = ix2 - ix1 + 1
                 iy1 = 1
                 do jj=j,1,-1
                   if (funa(i,jj) .lt. lowintensity) then
                     iy1 = jj + 1
                     exit
                   endif
                 enddo
                 iy2 = ny
                 do jj=j,ny
                   if (funa(i,jj) .lt. lowintensity) then
                     iy2 = jj - 1
                     exit
                   endif
                 enddo
                 iylength = iy2 - iy1 + 1
                 iu1 = 1
                 do m=1,min(i-1,j-1)
                   if (funa(i-m,j-m) .lt. lowintensity) then
                     iu1 = m
                     exit
                   endif
                 enddo
                 iu2 = 1
                 do n=1,min(nx-i,ny-j)
                   if (funa(i+n,j+n) .lt. lowintensity) then
                     iu2 = n
                     exit
                   endif
                 enddo
                 iulength = iu2 + iu1 + 1
                 iv1 = 1
                 do m=1,min(i-1,ny-j)
                   if (funa(i-m,j+m) .lt. lowintensity) then
                     iv1 = m
                     exit
                   endif
                 enddo
                 iv2 = 1
                 do n=1,min(nx-i,j-1)
                   if (funa(i+n,j-n) .lt. lowintensity) then
                     iv2 = n
                     exit
                   endif
                 enddo
                 ivlength = iv2 + iv1 + 1

! Sort directions of the filament profiles in a sequence of increasing widths.

                 lengthmin = min ( ixlength, iylength, iulength, ivlength )
                 if (lengthmin .eq. 0) then
                   write (iotty,'(a)') '   MODFITS: SKELETON: ERROR: lengthmin = 0'
                   stop 99
                 endif
                 if (lengthmin .eq. ixlength) then
                   idirection = 1
                   iwidth = ixlength
                 else if (lengthmin .eq. iylength) then
                   idirection = 2
                   iwidth = iylength
                 else if (lengthmin .eq. iulength) then
                   idirection = 3
                   iwidth = iulength
                 else if (lengthmin .eq. ivlength) then
                   idirection = 4
                   iwidth = ivlength
                 endif

! The non-zero values of the skeleton image will contain the direction, where the corresponding filament has smallest width. 
! The coded information is in the form (1.0d0 + directioncode), where directioncode = direction / 1000.0d0, where 
! direction = {1|2|3|4} is one of the four main directions. The coded direction value is transferred to the skeleton image,
! which is expected to always be normalized to 1.0 (plus the direction code). **NOT USED ANYMORE**

                 newpixval = 1.0d0
!!!                 directioncode = dble ( idirection ) / 1000.0d0
                 directioncode = 0.0d0
                   
                 if (idirection .eq. 1) then
                   if (funa(im1,j) .lt. funa(i,j) .and. funa(ip1,j) .lt. funa(i,j)) then
!!!                      funb(i,j) = newpixval
                     do l=max(j-nskw/2,1),min(j+nskw/2,ny)
                       do k=max(i-nskw/2,1),min(i+nskw/2,nx)
                         if (funb(k,l) .lt. almostzero) funb(k,l) = newpixval + directioncode
                       enddo
                     enddo
                     goto 10
                   endif
                 endif
                 if (idirection .eq. 2) then
                   if (funa(i,jm1) .lt. funa(i,j) .and. funa(i,jp1) .lt. funa(i,j)) then
!!!                      funb(i,j) = newpixval
                     do l=max(j-nskw/2,1),min(j+nskw/2,ny)
                       do k=max(i-nskw/2,1),min(i+nskw/2,nx)
                         if (funb(k,l) .lt. almostzero) funb(k,l) = newpixval + directioncode
                       enddo
                     enddo
                     goto 10
                   endif
                 endif
                 if (idirection .eq. 3) then
                   if (funa(im1,jm1) .lt. funa(i,j) .and. funa(ip1,jp1) .lt. funa(i,j) .and. 
     &                  funa(im1,jm1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero) then
!!!                      funb(i,j) = newpixval
                     do l=max(j-nskw/2,1),min(j+nskw/2,ny)
                       do k=max(i-nskw/2,1),min(i+nskw/2,nx)
                         if (funb(k,l) .lt. almostzero) funb(k,l) = newpixval + directioncode
                       enddo
                     enddo
                     goto 10
                   endif
                 endif
                 if (idirection .eq. 4) then
                   if (funa(im1,jp1) .lt. funa(i,j) .and. funa(ip1,jm1) .lt. funa(i,j) .and.
     &                  funa(im1,jp1) .gt. almostzero .and. funa(ip1,jm1) .gt. almostzero) then
!!!                      funb(i,j) = newpixval
                     do l=max(j-nskw/2,1),min(j+nskw/2,ny)
                       do k=max(i-nskw/2,1),min(i+nskw/2,nx)
                         if (funb(k,l) .lt. almostzero) funb(k,l) = newpixval + directioncode
                       enddo
                     enddo
                     goto 10
                   endif
                 endif
  10             continue
               endif
             endif
           enddo
         enddo
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.get.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:8) .eq. 'skeleton' .or. maction(1:3) .eq. 'ske') then

         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         itermax = 2000

         lsknorm = .false.
         if (maction(9:9) .eq. '1' .or. maction(4:4) .eq. '1') lsknorm = .true.
         maction = 'skeleton'
         iact = lastc ( maction )
         if (lsknorm) then
           cnorm = ' (norm)'
         else
           cnorm = ''
         endif
         nskw = max ( nint ( value ), 1 )
         if (mod ( nskw, 2 ) .eq. 0) nskw = nskw + 1
         inrem = int ( log10 ( dble ( max ( nskw, 1 ) ) ) ) + 1
         write (cnrem,'(i7)') nskw
         if (cverbose .eq. '-verb2')
     &       write (iotty,'(a)') '   Skeletonizing image using Hilditch algorithm: '//cnrem(7-inrem+1:7)//' pixels wide'//cnorm

! My implementation of the Hilditch's algorithm for making skeletons by thinning structures (iteratively).

         do iter=1,itermax
           nrem = 0

           do j=1,ny
             jm1 = max ( j - 1, 1 )
             jm2 = max ( j - 2, 1 )
             jp1 = min ( j + 1, ny )
             jp2 = min ( j + 2, ny )
             do i=1,nx
               im1 = max ( i - 1, 1 )
               im2 = max ( i - 2, 1 )
               ip1 = min ( i + 1, nx )
               ip2 = min ( i + 2, nx )

               if (funa(i,j) .gt. almostzero) then
                 numpix4 = 0
                 if (funa(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                 if (funa(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                 if (funa(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                 if (funa(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                 numpix8 = 0
                 if (funa(im1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
                 if (funa(ip1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
                 if (funa(im1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
                 if (funa(ip1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
                 numa1 = 0
                 if (funa(im1,j  ) .lt. almostzero .and. funa(im1,jp1) .gt. almostzero) numa1 = numa1 + 1
                 if (funa(im1,jp1) .lt. almostzero .and. funa(i  ,jp1) .gt. almostzero) numa1 = numa1 + 1
                 if (funa(i  ,jp1) .lt. almostzero .and. funa(ip1,jp1) .gt. almostzero) numa1 = numa1 + 1
                 if (funa(ip1,jp1) .lt. almostzero .and. funa(ip1,j  ) .gt. almostzero) numa1 = numa1 + 1
                 if (funa(ip1,j  ) .lt. almostzero .and. funa(ip1,jm1) .gt. almostzero) numa1 = numa1 + 1
                 if (funa(ip1,jm1) .lt. almostzero .and. funa(i  ,jm1) .gt. almostzero) numa1 = numa1 + 1
                 if (funa(i  ,jm1) .lt. almostzero .and. funa(im1,jm1) .gt. almostzero) numa1 = numa1 + 1
                 if (funa(im1,jm1) .lt. almostzero .and. funa(im1,j  ) .gt. almostzero) numa1 = numa1 + 1
                 numa2 = 0
                 if (funa(im1,jp1) .lt. almostzero .and. funa(im1,jp2) .gt. almostzero) numa2 = numa2 + 1
                 if (funa(im1,jp2) .lt. almostzero .and. funa(i  ,jp2) .gt. almostzero) numa2 = numa2 + 1
                 if (funa(i  ,jp2) .lt. almostzero .and. funa(ip1,jp2) .gt. almostzero) numa2 = numa2 + 1
                 if (funa(ip1,jp2) .lt. almostzero .and. funa(ip1,jp1) .gt. almostzero) numa2 = numa2 + 1
                 if (funa(ip1,jp1) .lt. almostzero .and. funa(ip1,j  ) .gt. almostzero) numa2 = numa2 + 1
                 if (funa(ip1,j  ) .lt. almostzero .and. funa(i  ,j  ) .gt. almostzero) numa2 = numa2 + 1
                 if (funa(i  ,j  ) .lt. almostzero .and. funa(im1,j  ) .gt. almostzero) numa2 = numa2 + 1
                 if (funa(im1,j  ) .lt. almostzero .and. funa(im1,jp1) .gt. almostzero) numa2 = numa2 + 1
                 numa4 = 0
                 if (funa(i  ,j  ) .lt. almostzero .and. funa(i  ,jp1) .gt. almostzero) numa4 = numa4 + 1
                 if (funa(i  ,jp1) .lt. almostzero .and. funa(ip1,jp1) .gt. almostzero) numa4 = numa4 + 1
                 if (funa(ip1,jp1) .lt. almostzero .and. funa(ip2,jp1) .gt. almostzero) numa4 = numa4 + 1
                 if (funa(ip2,jp1) .lt. almostzero .and. funa(ip2,j  ) .gt. almostzero) numa4 = numa4 + 1
                 if (funa(ip2,j  ) .lt. almostzero .and. funa(ip2,jm1) .gt. almostzero) numa4 = numa4 + 1
                 if (funa(ip2,jm1) .lt. almostzero .and. funa(ip1,jm1) .gt. almostzero) numa4 = numa4 + 1
                 if (funa(ip1,jm1) .lt. almostzero .and. funa(i  ,jm1) .gt. almostzero) numa4 = numa4 + 1
                 if (funa(i  ,jm1) .lt. almostzero .and. funa(i  ,j  ) .gt. almostzero) numa4 = numa4 + 1
                 n246 = 1
                 if (funa(i,jp1) .lt. almostzero .or. funa(ip1,j) .lt. almostzero .or. funa(i,jm1) .lt. almostzero) n246 = 0
                 n248 = 1
                 if (funa(i,jp1) .lt. almostzero .or. funa(ip1,j) .lt. almostzero .or. funa(im1,j) .lt. almostzero) n248 = 0
                 numb1 = numpix4 + numpix8

! The original Hilditch's algorithm does not work for skeletons along the main diagonals, eating them up pixel by pixel
! during iterative passes. The next statements defining LFIXDIAGS fixes that deficiency, preserving the diagonal skeletons.
! 2019-02-05: Fixed mistake in the 5-th condition below.

                 lfixdiags = .false.
             
                 if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.   !  1 1 0 0
     &               funa(i,jp1) .lt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.   !  * 1 + 0
     &               funa(ip1,j) .lt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.   !  * 0 0 0
     &               funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero .and.   !  * * * *
     &               funa(im2,jp1) .gt. almostzero) lfixdiags = .true.
                 
                 if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.   !  1 * * *
     &               funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.   !  1 1 0 *
     &               funa(ip1,j) .lt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.   !  0 + 0 *
     &               funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero .and.   !  0 0 0 *
     &               funa(im1,jp2) .gt. almostzero) lfixdiags = .true.
             
                 if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.   !  0 0 1 1
     &               funa(i,jp1) .lt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.   !  0 + 1 *
     &               funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.   !  0 0 0 *
     &               funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero .and.   !  * * * *
     &               funa(ip2,jp1) .gt. almostzero) lfixdiags = .true.
                 
                 if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.   !  * * * 1
     &               funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.   !  * 0 1 1
     &               funa(ip1,j) .lt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.   !  * 0 + 0
     &               funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero .and.   !  * 0 0 0
     &               funa(ip1,jp2) .gt. almostzero) lfixdiags = .true.
             
                 if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.   !  * * * *
     &               funa(i,jp1) .lt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.   !  * 0 0 0
     &               funa(ip1,j) .lt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.   !  * 1 + 0
     &               funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .gt. almostzero .and.   !  1 1 0 0
     &               funa(im2,jm1) .gt. almostzero) lfixdiags = .true.
                 
                 if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.   !  0 0 0 *
     &               funa(i,jp1) .lt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.   !  0 + 0 *
     &               funa(ip1,j) .lt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.   !  1 1 0 *
     &               funa(i,jm1) .gt. almostzero .and. funa(im1,jm1) .gt. almostzero .and.   !  1 * * *
     &               funa(im1,jm2) .gt. almostzero) lfixdiags = .true.
             
                 if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.   !  * * * *
     &               funa(i,jp1) .lt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.   !  0 0 0 *
     &               funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .gt. almostzero .and.   !  0 + 1 *
     &               funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero .and.   !  0 0 1 1
     &               funa(ip2,jm1) .gt. almostzero) lfixdiags = .true.
                 
                 if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.   !  * 0 0 0
     &               funa(i,jp1) .lt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.   !  * 0 + 0
     &               funa(ip1,j) .lt. almostzero .and. funa(ip1,jm1) .gt. almostzero .and.   !  * 0 1 1
     &               funa(i,jm1) .gt. almostzero .and. funa(im1,jm1) .lt. almostzero .and.   !  * * * 1
     &               funa(ip1,jm2) .gt. almostzero) lfixdiags = .true.

                 if (numb1 .ge. 2 .and. numb1 .le. 6 .and. numa1 .eq. 1 .and.
     &               (n248 .eq. 0 .or. numa2 .ne. 1) .and. (n246 .eq. 0 .or. numa4 .ne. 1)) then
                   if (.not.lfixdiags) then
                     nrem = nrem + 1
                     funb(i,j) = 0.0d0
                   endif
                 endif
               endif
             enddo
           enddo
           
           if (cverbose .eq. '-verb2') then
             inrem = int ( log10 ( dble ( max ( nrem, 1 ) ) ) ) + 1
             write (cnrem,'(i7)') nrem
             write (iotty,'(a)') '   Removed '//cnrem(7-inrem+1:7)//' pixels'
           endif
           do j=1,ny
             do i=1,nx
               funa(i,j) = funb(i,j)
             enddo
           enddo
           if (iter .eq. itermax) then
             write (iotty,'(/a)') '   MODFITS: ERROR: No convergence in the Hilditch algorithm'
             stop 99
           endif
           if (nrem .eq. 0) exit
         enddo
         
! Break loops and branches from the derived skeletons, based on the surrounding pixel values;
! a pixel with a minimum value serves as a breaking point.
       
         nrem = 0

         do j=1,ny
           jm1 = max ( j - 1, 1 )
           jm2 = max ( j - 2, 1 )
           jp1 = min ( j + 1, ny )
           jp2 = min ( j + 2, ny )
           do i=1,nx
             im1 = max ( i - 1, 1 )
             im2 = max ( i - 2, 1 )
             ip1 = min ( i + 1, nx )
             ip2 = min ( i + 2, nx )

             if (funb(i,j) .gt. almostzero) then
! case 1:
               if (funb(im1,j) .gt. almostzero .and. funb(im1,jm1) .gt. almostzero .and. funb(im1,jp1) .gt. almostzero .and.
     &             funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(im1,j), funb(im1,jm1), funb(im1,jp1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(im1,j) - funbmin ) .lt. almostzero) then
                   funb(im1,j) = 0.0d0
                 else if (abs ( funb(im1,jm1) - funbmin ) .lt. almostzero) then
                   funb(im1,jm1) = 0.0d0
                 else if (abs ( funb(im1,jp1) - funbmin ) .lt. almostzero) then
                   funb(im1,jp1) = 0.0d0
                 endif
               endif
               if (funb(ip1,j) .gt. almostzero .and. funb(ip1,jm1) .gt. almostzero .and. funb(ip1,jp1) .gt. almostzero .and.
     &             funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(ip1,j), funb(ip1,jm1), funb(ip1,jp1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(ip1,j) - funbmin ) .lt. almostzero) then
                   funb(ip1,j) = 0.0d0
                 else if (abs ( funb(ip1,jm1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jm1) = 0.0d0
                 else if (abs ( funb(ip1,jp1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jp1) = 0.0d0
                 endif
               endif
               if (funb(i,jm1) .gt. almostzero .and. funb(ip1,jm1) .gt. almostzero .and. funb(im1,jm1) .gt. almostzero .and.
     &             funb(im1,j) .lt. almostzero .and. funb(ip1,j) .lt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(i,jm1), funb(ip1,jm1), funb(im1,jm1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(i,jm1) - funbmin ) .lt. almostzero) then
                   funb(i,jm1) = 0.0d0
                 else if (abs ( funb(ip1,jm1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jm1) = 0.0d0
                 else if (abs ( funb(im1,jm1) - funbmin ) .lt. almostzero) then
                   funb(im1,jm1) = 0.0d0
                 endif
               endif
               if (funb(i,jp1) .gt. almostzero .and. funb(ip1,jp1) .gt. almostzero .and. funb(im1,jp1) .gt. almostzero .and.
     &             funb(im1,j) .lt. almostzero .and. funb(ip1,j) .lt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(i,jp1), funb(ip1,jp1), funb(im1,jp1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(i,jp1) - funbmin ) .lt. almostzero) then
                   funb(i,jp1) = 0.0d0
                 else if (abs ( funb(ip1,jp1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jp1) = 0.0d0
                 else if (abs ( funb(im1,jp1) - funbmin ) .lt. almostzero) then
                   funb(im1,jp1) = 0.0d0
                 endif
               endif
! case 2:
               if (funb(im1,jm1) .gt. almostzero .and. funb(i,jp1) .gt. almostzero .and. funb(ip1,j) .gt. almostzero .and.
     &             funb(im1,jp1) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and. funb(ip1,jp1) .lt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(im1,jm1), funb(i,jp1), funb(ip1,j) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(im1,jm1) - funbmin ) .lt. almostzero) then
                   funb(im1,jm1) = 0.0d0
                 else if (abs ( funb(i,jp1) - funbmin ) .lt. almostzero) then
                   funb(i,jp1) = 0.0d0
                 else if (abs ( funb(ip1,j) - funbmin ) .lt. almostzero) then
                   funb(ip1,j) = 0.0d0
                 endif
               endif
               if (funb(im1,jp1) .gt. almostzero .and. funb(i,jm1) .gt. almostzero .and. funb(ip1,j) .gt. almostzero .and.
     &             funb(im1,jm1) .lt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(im1,jp1), funb(i,jm1), funb(ip1,j) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(im1,jp1) - funbmin ) .lt. almostzero) then
                   funb(im1,jp1) = 0.0d0
                 else if (abs ( funb(i,jm1) - funbmin ) .lt. almostzero) then
                   funb(i,jm1) = 0.0d0
                 else if (abs ( funb(ip1,j) - funbmin ) .lt. almostzero) then
                   funb(ip1,j) = 0.0d0
                 endif
               endif
               if (funb(ip1,jp1) .gt. almostzero .and. funb(im1,j) .gt. almostzero .and. funb(i,jm1) .gt. almostzero .and.
     &             funb(im1,jp1) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and. funb(im1,jm1) .lt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(ip1,jp1), funb(im1,j), funb(i,jm1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(ip1,jp1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jp1) = 0.0d0
                 else if (abs ( funb(im1,j) - funbmin ) .lt. almostzero) then
                   funb(im1,j) = 0.0d0
                 else if (abs ( funb(i,jm1) - funbmin ) .lt. almostzero) then
                   funb(i,jm1) = 0.0d0
                 endif
               endif
               if (funb(ip1,jm1) .gt. almostzero .and. funb(im1,j) .gt. almostzero .and. funb(i,jp1) .gt. almostzero .and.
     &             funb(im1,jm1) .lt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and. funb(im1,jp1) .lt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(ip1,jm1), funb(im1,j), funb(i,jp1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(ip1,jm1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jm1) = 0.0d0
                 else if (abs ( funb(im1,j) - funbmin ) .lt. almostzero) then
                   funb(im1,j) = 0.0d0
                 else if (abs ( funb(i,jp1) - funbmin ) .lt. almostzero) then
                   funb(i,jp1) = 0.0d0
                 endif
               endif
! case 3:
               if (funb(im1,jm1) .gt. almostzero .and. funb(ip1,jm1) .gt. almostzero .and. funb(i,jp1) .gt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(im1,jm1), funb(ip1,jm1), funb(i,jp1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(im1,jm1) - funbmin ) .lt. almostzero) then
                   funb(im1,jm1) = 0.0d0
                 else if (abs ( funb(ip1,jm1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jm1) = 0.0d0
                 else if (abs ( funb(i,jp1) - funbmin ) .lt. almostzero) then
                   funb(i,jp1) = 0.0d0
                 endif
               endif
               if (funb(im1,jp1) .gt. almostzero .and. funb(ip1,jp1) .gt. almostzero .and. funb(i,jm1) .gt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(im1,jp1), funb(ip1,jp1), funb(i,jm1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(im1,jp1) - funbmin ) .lt. almostzero) then
                   funb(im1,jp1) = 0.0d0
                 else if (abs ( funb(ip1,jp1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jp1) = 0.0d0
                 else if (abs ( funb(i,jm1) - funbmin ) .lt. almostzero) then
                   funb(i,jm1) = 0.0d0
                 endif
               endif
               if (funb(ip1,jm1) .gt. almostzero .and. funb(ip1,jp1) .gt. almostzero .and. funb(im1,j) .gt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(ip1,jm1), funb(ip1,jp1), funb(im1,j) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(ip1,jm1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jm1) = 0.0d0
                 else if (abs ( funb(ip1,jp1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jp1) = 0.0d0
                 else if (abs ( funb(im1,j) - funbmin ) .lt. almostzero) then
                   funb(im1,j) = 0.0d0
                 endif
               endif
               if (funb(im1,jm1) .gt. almostzero .and. funb(im1,jp1) .gt. almostzero .and. funb(ip1,j) .gt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(im1,jm1), funb(im1,jp1), funb(ip1,j) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(im1,jm1) - funbmin ) .lt. almostzero) then
                   funb(im1,jm1) = 0.0d0
                 else if (abs ( funb(im1,jp1) - funbmin ) .lt. almostzero) then
                   funb(im1,jp1) = 0.0d0
                 else if (abs ( funb(ip1,j) - funbmin ) .lt. almostzero) then
                   funb(ip1,j) = 0.0d0
                 endif
               endif
! case 4:
               if (funb(im1,jm1) .gt. almostzero .and. funb(im1,jp1) .gt. almostzero .and. funb(ip1,jp1) .gt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(im1,jm1), funb(im1,jp1), funb(ip1,jp1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(im1,jm1) - funbmin ) .lt. almostzero) then
                   funb(im1,jm1) = 0.0d0
                 else if (abs ( funb(im1,jp1) - funbmin ) .lt. almostzero) then
                   funb(im1,jp1) = 0.0d0
                 else if (abs ( funb(ip1,jp1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jp1) = 0.0d0
                 endif
               endif
               if (funb(im1,jp1) .gt. almostzero .and. funb(ip1,jp1) .gt. almostzero .and. funb(ip1,jm1) .gt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(im1,jp1), funb(ip1,jp1), funb(ip1,jm1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(im1,jp1) - funbmin ) .lt. almostzero) then
                   funb(im1,jp1) = 0.0d0
                 else if (abs ( funb(ip1,jp1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jp1) = 0.0d0
                 else if (abs ( funb(ip1,jm1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jm1) = 0.0d0
                 endif
               endif
               if (funb(im1,jm1) .gt. almostzero .and. funb(im1,jp1) .gt. almostzero .and. funb(ip1,jm1) .gt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(im1,jm1), funb(im1,jp1), funb(ip1,jm1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(im1,jm1) - funbmin ) .lt. almostzero) then
                   funb(im1,jm1) = 0.0d0
                 else if (abs ( funb(im1,jp1) - funbmin ) .lt. almostzero) then
                   funb(im1,jp1) = 0.0d0
                 else if (abs ( funb(ip1,jm1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jm1) = 0.0d0
                 endif
               endif
               if (funb(im1,jm1) .gt. almostzero .and. funb(ip1,jp1) .gt. almostzero .and. funb(ip1,jm1) .gt. almostzero) then
                 nrem = nrem + 1
                 funbmin = min ( funb(i,j), funb(im1,jm1), funb(ip1,jp1), funb(ip1,jm1) )
                 if (abs ( funb(i,j) - funbmin ) .lt. almostzero) then
                   funb(i,j) = 0.0d0
                 else if (abs ( funb(im1,jm1) - funbmin ) .lt. almostzero) then
                   funb(im1,jm1) = 0.0d0
                 else if (abs ( funb(ip1,jp1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jp1) = 0.0d0
                 else if (abs ( funb(ip1,jm1) - funbmin ) .lt. almostzero) then
                   funb(ip1,jm1) = 0.0d0
                 endif
               endif
             endif
           enddo
         enddo

         if (cverbose .eq. '-verb2') then
           inrem = int ( log10 ( dble ( max ( nrem, 1 ) ) ) ) + 1
           write (cnrem,'(i7)') nrem
           write (iotty,'(a)') '   Fixed '//cnrem(7-inrem+1:7)//' pixels'
         endif

         do j=1,ny
           do i=1,nx
             funa(i,j) = funb(i,j)
           enddo
         enddo
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. almostzero) then
               do l=max(j-nskw/2,1),min(j+nskw/2,ny)
                 do k=max(i-nskw/2,1),min(i+nskw/2,nx)
                   if (lsknorm) then
                     funb(k,l) = 1.0d0
                   endif
                 enddo
               enddo
             endif
           enddo
         enddo
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.ske'//cvalue(1:icva-1)//'.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif
       
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'rosenfeld' .or. maction(1:iact) .eq. 'ros') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Skeletonizing image using Rosenfeld algorithm'

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. almostzero) then
               im1 = max ( i - 1, 1 )
               jm1 = max ( j - 1, 1 )
               ip1 = min ( i + 1, nx )
               jp1 = min ( j + 1, ny )

               numpix4 = 0
               if (funa(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
               if (funa(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
               if (funa(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
               if (funa(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
               numpix8 = 0
               if (funa(im1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
               if (funa(ip1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
               if (funa(im1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
               if (funa(ip1,jm1) .gt. almostzero) numpix8 = numpix8 + 1

               lnorthborder = funa(i,jp1) .lt. almostzero
               lsouthborder = funa(i,jm1) .lt. almostzero
               leastborder  = funa(ip1,j) .lt. almostzero
               lwestborder  = funa(im1,j) .lt. almostzero
               l4endpoint   = numpix4 .eq. 1
               l8endpoint   = numpix8 .eq. 1
               l4isolated   = numpix4 .eq. 0
               l8isolated   = numpix8 .eq. 0

!!               l8simple = .false.
!!               if (numpix8 .le. 1) l8simple = .true.

               if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero) l8simple = .true.

               if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero) l8simple = .true.

               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero) l8simple = .true.

               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero) l8simple = .true.

               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .gt. almostzero) l8simple = .true.

               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .lt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .gt. almostzero) l8simple = .true.

               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .gt. almostzero) l8simple = .true.

               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .lt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .gt. almostzero) l8simple = .true.


               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(i,jp1) .lt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .gt. almostzero .and. funa(im1,jm1) .lt. almostzero) l8simple = .true.

               if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .gt. almostzero .and. funa(im1,jm1) .gt. almostzero) l8simple = .true.


               if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.
     &             funa(i,jp1) .lt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero) l8simple = .true.

               if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .lt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero) l8simple = .true.

               if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .lt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero) l8simple = .true.

               if (funa(im1,j) .lt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero) l8simple = .true.

               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .lt. almostzero) l8simple = .true.

               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .lt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .gt. almostzero) l8simple = .true.

               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .lt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .gt. almostzero) l8simple = .true.

               if (funa(im1,j) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and.
     &             funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .gt. almostzero .and.
     &             funa(i,jm1) .lt. almostzero .and. funa(im1,jm1) .gt. almostzero) l8simple = .true.

               if (lnorthborder .and. l8simple .and. .not.l8isolated .and. .not.l8endpoint) then
                 funb(i,j) = 0.0d0
               endif
               if (lsouthborder .and. l8simple .and. .not.l8isolated .and. .not.l8endpoint) then
                 funb(i,j) = 0.0d0
               endif
               if (leastborder .and. l8simple .and. .not.l8isolated .and. .not.l8endpoint) then
                 funb(i,j) = 0.0d0
               endif
               if (lwestborder .and. l8simple .and. .not.l8isolated .and. .not.l8endpoint) then
                 funb(i,j) = 0.0d0
               endif
               if (funb(i,j) .gt. funmax) funmax = funb(i,j)
               if (funb(i,j) .lt. funmin) funmin = funb(i,j)
             endif
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.ske.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif
       
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:5) .eq. 'clean' .or. maction(1:3) .eq. 'cle') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         l4conn = .false.
         l8conn = .false.
         lskeletons = .false.
         if (maction(6:6) .eq. '4' .or. maction(4:4) .eq. '4') l4conn = .true.
         if (maction(6:6) .eq. '8' .or. maction(4:4) .eq. '8') l8conn = .true.
         if (maction(7:8) .eq. 'sk' .or. maction(5:6) .eq. 'sk') lskeletons = .true.
         if (.not.l4conn .and. .not.l8conn) then
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   MODFITS: ERROR: Invalid action: "'//maction(1:iact)//'".'
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '            Known actions: "clean4", "cle4", "clean8", "cle8",'
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '            "clean4sk", "cle4sk", "clean8sk", "cle8sk"'
           stop 99
         endif         
         minpix = dble ( nint ( value ) )
         read (param(1:ipar),*,err=48) clefactor
         goto 49
  48     continue
         if (cverbose .eq. '-verb2') write (iotty,'(/a)') '   MODFITS: ERROR: Trouble reading CLEFACTOR from command line.'
         stop 99
  49     continue

         if (cverbose .eq. '-verb2') then
           if (l4conn) then
             write (iotty,'(a,1pe10.3,a)') '   Removing 4-conn clusters smaller than', minpix, ' pixels ('//param(1:ipar-1)//')' 
           endif
           if (l8conn) then
             write (iotty,'(a,1pe10.3,a)') '   Removing 8-conn clusters smaller than', minpix, ' pixels ('//param(1:ipar-1)//')' 
           endif
         endif  

         allocate ( worker(nx,ny), nsegm(nx,ny), origmask(nx,ny), n1x(ny), n2x(ny), mxco(1), myco(1), deltx(1), delty(1)
     &            , ntouching(1,1), nsx(1), nxmn(nextrmax), nxmx(nextrmax), nymn(nextrmax), nymx(nextrmax), stat=irc )

         if (irc .ne. 0) then
           write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (12).'
           stop 12
         endif
         changed = .false.

         if (funmax .gt. almostzero) then

           do j=1,ny
             n1x(j) = 1
             n2x(j) = nx
           enddo
           nymin = 1
           nymax = ny
           do j=1,ny
             do i=1,nx
               if (funa(i,j) .gt. almostzero) then
                 origmask(i,j) = 1.0d0
               else
                 origmask(i,j) = 0.0d0
               endif
               worker(i,j) = funa(i,j)
               funb(i,j) = worker(i,j)
               nsegm(i,j) = 0.0d0
             enddo
           enddo
         
! Initially remove isolated pixels.

           if (l8conn) then

             do j=1,ny
               jm1 = max ( j - 1, 1 )
               jp1 = min ( j + 1, ny )
               do i=1,nx
                 im1 = max ( i - 1, 1 )
                 ip1 = min ( i + 1, nx )
                 if (funa(i,j) .gt. almostzero) then
                   nozero = 0
                   if (funa(im1,jm1) .gt. almostzero) nozero = nozero + 1
                   if (funa(im1,j  ) .gt. almostzero) nozero = nozero + 1
                   if (funa(im1,jp1) .gt. almostzero) nozero = nozero + 1
                   if (funa(i  ,jp1) .gt. almostzero) nozero = nozero + 1
                   if (funa(ip1,jp1) .gt. almostzero) nozero = nozero + 1
                   if (funa(ip1,j  ) .gt. almostzero) nozero = nozero + 1
                   if (funa(ip1,jm1) .gt. almostzero) nozero = nozero + 1
                   if (funa(i  ,jm1) .gt. almostzero) nozero = nozero + 1
                   if (nozero .eq. 0) then
                     funa(i,j) = 0.0d0
                     funb(i,j) = 0.0d0
                     worker(i,j) = 0.0d0
                     if (cverbose .eq. '-verb2') then
                       write (iotty,'(a,2i5,1x,i2,a)') '   Removing isolated pixel', i, j, nozero
                     endif
                   endif
                 endif
               enddo
             enddo
           
! To clean clusters, I use TintFill which works for 4-conn clusters only (in my implementation).
! The following call converts all 8-conn pixels to 4-conn before calling TintFill.
          
             if (lskeletons) then
               call convert84 ( inmx, nx, ny, funb, worker, filename(1:fnlen1), cverbose, iotty, almostzero, lskeletons )
             else
               do j=1,ny
                 jm1 = max ( j - 1, 1 )
                 jp1 = min ( j + 1, ny )
                 do i=1,nx
                   if (worker(i,j) .gt. almostzero) then
                     im1 = max ( i - 1, 1 )
                     ip1 = min ( i + 1, nx )
                     if (worker(im1,jm1) .gt. almostzero .and.
     &                   worker(im1,j) .lt. almostzero .and. worker(i,jm1) .lt. almostzero ) then
                       worker(i,jm1) = worker(i,j)
                       funb(i,jm1) = worker(i,j)
                     endif
                     if (worker(im1,jp1) .gt. almostzero .and.
     &                   worker(im1,j) .lt. almostzero .and. worker(i,jp1) .lt. almostzero ) then
                       worker(im1,j) = worker(i,j)
                       funb(im1,j) = worker(i,j)
                     endif
                     if (worker(ip1,jp1) .gt. almostzero .and.
     &                   worker(ip1,j) .lt. almostzero .and. worker(i,jp1) .lt. almostzero ) then
                       worker(i,jp1) = worker(i,j)
                       funb(i,jp1) = worker(i,j)
                     endif
                     if (worker(ip1,jm1) .gt. almostzero .and.
     &                   worker(ip1,j) .lt. almostzero .and. worker(i,jm1) .lt. almostzero ) then
                       worker(ip1,j) = worker(i,j)
                       funb(ip1,j) = worker(i,j)
                     endif
                   endif
                 enddo
               enddo
             endif
           endif

!!           if (l4conn) then
!!             do j=1,ny
!!               jm1 = max ( j - 1, 1 )
!!               jp1 = min ( j + 1, ny )
!!               do i=1,nx
!!                 if (worker(i,j) .gt. almostzero) then
!!                   im1 = max ( i - 1, 1 )
!!                   ip1 = min ( i + 1, nx )
!!                   numpix = 0
!!                   if (worker(im1,j  ) .gt. almostzero) numpix = numpix + 1
!!                   if (worker(ip1,j  ) .gt. almostzero) numpix = numpix + 1
!!                   if (worker(i  ,jm1) .gt. almostzero) numpix = numpix + 1
!!                   if (worker(i  ,jp1) .gt. almostzero) numpix = numpix + 1
!!                   if (worker(im1,jm1) .gt. almostzero) numpix = numpix + 1
!!                   if (worker(ip1,jp1) .gt. almostzero) numpix = numpix + 1
!!                   if (worker(im1,jp1) .gt. almostzero) numpix = numpix + 1
!!                   if (worker(ip1,jm1) .gt. almostzero) numpix = numpix + 1
!!                   if (numpix .le. 3) then
!!                     worker(i,j) = 0.0d0
!!                   endif
!!                 endif
!!               enddo
!!             enddo
!!           endif
           
           do j=1,ny
             do i=1,nx
               funb(i,j) = worker(i,j)
               worker(i,j) = funb(i,j) / funmax / 1.1d0
               funa(i,j) = funb(i,j)
             enddo
           enddo
           
! Use TintFill to find connected clusters and numbers of connected pixels in them.

           if (cverbose .eq. '-verb2') write (iotty,'(a$)') '   Applying 4-conn TintFill algorithm '
           nextrmx = 1
           zeros = 1.0d-20
           mxco(1) = 1.0d0
           myco(1) = 1.0d0
           deltx(1) = 0.0d0
           delty(1) = 0.0d0
           ntouching(1,1) = 0
           nsx(1) = 0
           nextr = 0

           do j=1,ny
             do i=1,nx
               if (worker(i,j) .gt. almostzero .and. nsegm(i,j) .le. almostzero) then
                 nextr = nextr + 1
                 if (nextr .gt. nextrmax) then
                   write (iotty,'(/a,i8)') '   MODFITS: ERROR: Number of isolated clusters reached the maximum of ', nextrmax
                   write (iotty,'(a    )') '            Increase the value and recompile MODFITS if you wish to continue.'
                   stop 89
                 endif
                 newpixval = dble ( nextr )
           
                 call fillarea ( .false., i, j, nextrmx, nextr, newpixval, 'positive', zeros, boupixval, nx, ny, n1x, n2x, nymin
     &                         , nymax, nxmn(nextr), nxmx(nextr), nymn(nextr), nymx(nextr), worker, nsegm, nsegm, mxco, myco
     &                         , deltx, delty, ntouching, nsx )
               endif
             enddo
           enddo                   

           if (cverbose .eq. '-verb2') then
             inrem = int ( log10 ( dble ( max ( nextr, 1 ) ) ) ) + 1
             write (cnrem,'(i7)') nextr
             write (iotty,'(a)') cnrem(7-inrem+1:7)//' clusters'
           endif
           
           allocate ( momxco(nextr), momyco(nextr), afwhm(nextr), bfwhm(nextr), atheta(nextr), equivrad(nextr), elongation(nextr)
     &              , nonzero(nextr), afoot(nextr), bfoot(nextr), stat=irc )
           
           if (irc .ne. 0) then
             write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (13).'
             stop 13
           endif
           
           do k=1,nextr
             nonzero(k) = 0.0d0
             do j=nymn(k),nymx(k)
               do i=nxmn(k),nxmx(k)
                 if (nint ( nsegm(i,j) ) .eq. k) then
                   nonzero(k) = nonzero(k) + 1.0d0
                 endif
               enddo
             enddo                   
           enddo

           call sizemeasure ( iotty, nextr, nxmn, nxmx, nymn, nymx, inmx, nx, ny, funb, nsegm, dx, dy, momxco, momyco, afwhm, bfwhm
     &                      , atheta, afoot, bfoot, equivrad, elongation, nonzero, almostzero )
           
           if (ldebug) then
             open ( 15, file=filename(1:fnlen1)//'.cleaning', status='unknown' )
             write (15,'(a)') '#'
             write (15,'(a)') '# Results of cleaning of connected clusters of pixels by MODFITS.'
             write (15,'(a)') '# The value GOOD = 1.00 has a meaning of GOOD FOR REMOVAL.'
             write (15,'(a)') '# The values of AFHM01 and BFWHM01 were multiplied by 2.'
             write (15,'(a)') '#'
             write (15,'(a,2(1pe11.3))') '# Input parameters (minpix, clefactor):', minpix, clefactor
             write (15,'(a)') '#'
             write (15,'(a)') '#      N    XCO_P    YCO_P  GOOD  AFWH01x2   BFWH01x2   THEP01     ELONGATIO'
     &                      //'     NONZERO  FXP_BEST01 FXP_ERRO01 FXT_BEST01 FXT_ERRO01 SIG_MONO01 SN_RATIO01'
             write (15,'(a)') '#'
           endif
           
           if (cverbose .eq. '-verb2') then
             write (iotty,'(a)') '     k  nonzero  elongatio  sparsity  nxmn nxmx toremove'
           endif

! Removal of structures whose area is smaller than an input value of MINPIX. If a structure (at some level) is not very elongated
! and it is also not very sparse, it belongs to the compact sources component and must be completely cut off, even if its area is
! larger than MINPIX. 
           
           nrem = 0
           do k=1,nextr
             sparsity = afwhm(k) * bfwhm(k) / equivrad(k)**2
             goodtoremove = 0.0d0

! Only small and compact shapes are removed when cutting sources.
! If a shape is compact, it will also be considered small, even if its area is larger than MINPIX.

             lcompact = elongation(k) .lt. clefactor * elongmaxsrc .and. sparsity .lt. clefactor * sparsmaxsrc
             lsmall = nonzero(k) .le. minpix .or. lcompact

             if (lsmall) then
               nrem = nrem + 1
               goodtoremove = 1.0d0
               do j=nymn(k),nymx(k)
                 do i=nxmn(k),nxmx(k)
                   if (nint ( nsegm(i,j) ) .eq. k) then
                     funa(i,j) = 0.0d0
                   endif
                 enddo
               enddo
             endif

             if (cverbose .eq. '-verb2') then
               write (iotty,'(i6,3(1pe10.3),3i5,2l)') k, nonzero(k), elongation(k), sparsity
     &               , nxmn(k), nxmx(k), nint ( goodtoremove ), lsmall, lcompact
             endif

! To debug, draw ellipses around clusters of connected pixels.

             if (ldebug) then
               dxym = 2.0d0 * sqrt ( dx * dy )
               do j=nymn(k),nymx(k)
                 dblm = dble ( j )
                 do i=nxmn(k),nxmx(k)
                   dbll = dble ( i )
                   if (inellipse ( dbll, dblm, dx, dy, momxco(k), momyco(k), afoot(k), bfoot(k), atheta(k) )) then
                     if (.not.inellipse ( dbll, dblm, dx, dy, momxco(k), momyco(k), afoot(k)-dxym, bfoot(k)-dxym, atheta(k) )) then
                       funa(i,j) = 10.0d0
                     endif
                   endif
                 enddo
               enddo                   
               write (15,'(i8,2f9.1,f6.2,13(1pe11.3))') 
     &                   k, momxco(k), momyco(k), goodtoremove, afoot(k), bfoot(k), atheta(k), elongation(k)
     &                 , nonzero(k), 1.0d0, 0.0d0, 1.0d0, 0.0d0, 100.0d0, 100.0d0
             endif
           enddo
           
           if (ldebug) then
             close ( 15 )
           endif
           
           deallocate ( momxco, momyco, afwhm, bfwhm, atheta, afoot, bfoot, equivrad, elongation, nonzero )

!!           if (lskeletons) then
!!             do j=1,ny
!!               jm1 = max ( j - 1, 1 )
!!               jp1 = min ( j + 1, ny )
!!               do i=1,nx
!!                 im1 = max ( i - 1, 1 )
!!                 ip1 = min ( i + 1, nx )
!!                 ip2 = min ( i + 2, nx )
!!                 if (funb(i,j) .gt. almostzero) then
!!                   numpix4 = 0
!!                   numpix8 = 0
!!                   if (funb(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
!!                   if (funb(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
!!                   if (funb(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
!!                   if (funb(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
!!                   if (funb(im1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
!!                   if (funb(ip1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
!!                   if (funb(im1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
!!                   if (funb(ip1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
!!
!!                   if (numpix4 .eq. 2 .and. numpix8 .eq. 1) then
!!                     funa(i,j) = 0.0d0
!!                     if (cverbose .eq. '-verb2') then
!!                       write (iotty,'(a,2i5,1x,4i2)') '   Additionally disconnected branch: ', i, j, numpix4, numpix8
!!                     endif
!!                   endif
!!                 endif
!!               enddo
!!             enddo
!!           endif
           
           if (cverbose .eq. '-verb2') then
             inrem = int ( log10 ( dble ( max ( nrem, 1 ) ) ) ) + 1
             write (cnrem,'(i7)') nrem
             write (iotty,'(a)') '   Removed '//cnrem(7-inrem+1:7)//' clusters of connected pixels'
           endif

! Multiplication by the original input image is necessary to restore the original pixels connections, because theyr were
! damaged above by adding 4-conn pixels when replacing 8-conn clusters with 4-conn to make TittFill work.
! Just a multiplication must be sufficient to restore the original connections.
           
           funmin =  1.0d+30
           funmax = -1.0d+30
           do j=1,ny
             do i=1,nx
               funb(i,j) = funa(i,j)
               if (lskeletons) then
                 funb(i,j) = funb(i,j)    !!* origmask(i,j)
               endif
               if (funb(i,j) .gt. funmax) funmax = funb(i,j)
               if (funb(i,j) .lt. funmin) funmin = funb(i,j)
             enddo
           enddo
           datamax = funmax
           datamin = funmin

           deallocate ( worker, nsegm, origmask, n1x, n2x, mxco, myco, deltx, delty, ntouching, nsx, nxmn, nxmx, nymn, nymx )

           if (lnoname) then
             ic1 = index ( cvalue(1:icva-1), '.' )
             if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
             outname = filename(isp1:fnlen1-5)//'.cle'//cvalue(1:icva-1)//'.fits'
           endif
           if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
           if (cverbose .eq. '-verb2') write (iotty,*) datamax
           if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
           if (cverbose .eq. '-verb2') write (iotty,*) datamin
           changed = .true.
         endif
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:5) .eq. 'clip' .or. maction(1:3) .eq. 'cli') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         l4conn = .false.
         l8conn = .false.
         lcutsources = .false.
         lcutfilaments = .false.
         lcutshapes = .false.
         if (maction(5:5) .eq. '4' .or. maction(3:3) .eq. '4') l4conn = .true.
         if (maction(5:5) .eq. '8' .or. maction(3:3) .eq. '8') l8conn = .true.
         if (maction(6:6) .eq. 's' .or. maction(4:4) .eq. 's') lcutsources = .true.
         if (maction(6:6) .eq. 'f' .or. maction(4:4) .eq. 'f') lcutfilaments = .true.
         if (maction(6:6) .eq. 'x' .or. maction(4:4) .eq. 'x') lcutshapes = .true.
         if ((.not.l4conn .and. .not.l8conn) .or. (.not.lcutsources .and. .not.lcutfilaments .and. .not.lcutshapes)) then
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   MODFITS: ERROR: Invalid action: "'//maction(1:iact)//'".'
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '            Known actions: "clip4s", "cli4s", "clip8s", "cli8s",'
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '                           "clip4f", "cli4f", "clip8f", "cli8f",'
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '                           "clip4x", "cli4x", "clip8x", "cli8x".'
           stop 99
         endif         
         minpix = dble ( nint ( value ) )
         read (param(1:ipar),*,err=58) levelx
         if (param(1:1) .eq. '+') then
           lcutpeaks = .true.
           lcutnegat = .false.
         else
           if (param(1:1) .eq. '-') then
             levelx = abs ( levelx )
             lcutpeaks = .false.
             lcutnegat = .true.
           else
             lcutpeaks = .false.
             lcutnegat = .false.
           endif
         endif
         goto 59
  58     continue
         if (cverbose .eq. '-verb2') write (iotty,'(/a)') '   MODFITS: ERROR: Trouble reading LEVELX from command line.'
         stop 99
  59     continue
         if (levelx .lt. almostzero) then
           if (cverbose .eq. '-verb2') write (iotty,'(/a)') '   MODFITS: ERROR: Incorrect value of LEVELX = 0 (| LEVELX | > 0).'
           stop 99
         endif
         if (cverbose .eq. '-verb2') then
           if (lcutpeaks) then
             if (l4conn) then
               write (iotty,'(a,1pe10.3,a)') '   Cutting off 4-conn peaks:', minpix, ' pixels ('//param(1:ipar-1)//')' 
             endif
             if (l8conn) then
               write (iotty,'(a,1pe10.3,a)') '   Cutting off 8-conn peaks:', minpix, ' pixels ('//param(1:ipar-1)//')' 
             endif
           else
             if (l4conn) then
               write (iotty,'(a,1pe10.3,a)') '   Cutting 4-conn peaks from filaments:', minpix, ' pixels ('//param(1:ipar-1)//')' 
             endif
             if (l8conn) then
               write (iotty,'(a,1pe10.3,a)') '   Cutting 8-conn peaks from filaments:', minpix, ' pixels ('//param(1:ipar-1)//')' 
             endif
           endif
         endif  
         
         allocate ( worker(nx,ny), nsegm(nx,ny), shapes(nx,ny), n1x(ny), n2x(ny), mxco(1), myco(1), deltx(1), delty(1)
     &            , ntouching(1,1), nsx(1), nxmn(nextrmax), nxmx(nextrmax), nymn(nextrmax), nymx(nextrmax), nonzero(nextrmax) 
     &            , nonzeroo(nextrmax), nsegmo(nx,ny), nsourceo(nextrmax), stat=irc )
     
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (12a).'
           stop 12
         endif

         do j=1,ny
           n1x(j) = 1
           n2x(j) = nx
         enddo
         do j=1,ny
           do i=1,nx
             nsegm(i,j) = 0.0d0
             if (lcutshapes) then
               shapes(i,j) = funa(i,j)
             else
               shapes(i,j) = 0.0d0
             endif
           enddo
         enddo
         nymin = 1
         nymax = ny
         nextrmx = 1
         zeros = 1.0d-20
         mxco(1) = 1.0d0
         myco(1) = 1.0d0
         deltx(1) = 0.0d0
         delty(1) = 0.0d0
         ntouching(1,1) = 0
         nsx(1) = 0
         
! Levels are logarithmically equidistant above levelx and below -levelx.

         nlevsmin = 2
         levelxlog = log ( levelx )
         funmaxlog = log ( funmax )
         dlpos = max ( funmaxlog - levelxlog, 0.0d0 )
         deltpos = deltapos
         if (lcutpeaks) then
           deltpos = 0.01d0
         endif
         nlevpos = nint ( dlpos / deltpos )
         if (nlevpos .lt. nlevsmin) then
           nlevpos = nlevsmin
           deltpos = dlpos / dble ( nlevpos )
         endif
         if (cverbose .eq. '-verb2') 
     &      write (iotty,'(a,i5,a,1pe10.3)') '   Positive levels:', nlevpos, ' delta:', deltpos         !! exp ( levelxlog ) * 

         do lev=1,nlevpos
           level(lev) = exp ( levelxlog + dble ( nlevpos - lev ) * deltpos )
         enddo
         nlevmax = nlevpos
         nlevneg = 0

! UNTRUE: No sense to clip filaments at zero and negative levels: shapes at those levels are interconnected.
! UPDATE 200305: It DOES make sense to clip negative levels even for filaments, because clipping just the positives
! does not dig deep enough, losing the faint outskirts of the filaments and making their profiles steeper than they are,
! which effectively means the filament background becomes overestimated.

!!         if (.not.lcutpeaks .and. .not.lcutshapes) then   !UPD200305 !! .and. .not.lcutfilaments
!!!!            nlevmax = nlevmax + 1
!!!!            level(nlevmax) = 0.0d0
!!           funminlog = log ( abs ( funmin ) )
!!           dlneg = max ( funminlog - levelxlog, 0.0d0 )
!!           deltneg = deltaneg           
!!           nlevneg = nint ( dlneg / deltaneg )
!!           if (nlevneg .lt. nlevsmin) then
!!             nlevneg = nlevsmin
!!             deltneg = dlneg / dble ( nlevneg )
!!           endif
!!           if (cverbose .eq. '-verb2') 
!!     &        write (iotty,'(a,i5,a,1pe10.3)') '   Negative levels:', nlevneg, ' delta:', deltneg      !! exp ( levelxlog ) * 
!!        
!!           do lev=1,nlevneg
!!             level(lev+nlevmax) = - exp ( levelxlog + dble ( lev - 1 ) * deltneg )
!!           enddo
!!           nlevmax = nlevmax + nlevneg
!!         endif

! Work down over all levels.

         nextr = 0
         do lev=1,nlevmax

           do j=1,ny
             do i=1,nx
               if (funa(i,j) .ge. level(lev)) then
                 worker(i,j) = (funa(i,j) - funmin) / (funmax - funmin) / 1.01d0
               else
                 worker(i,j) = 0.0d0
               endif
               funb(i,j) = worker(i,j)  !! funb(i,j) = min ( worker(i,j), level(lev) ) << 2019-12-23 old error from Feb 2019 ?? 
             enddo
           enddo

! To clean clusters, I use TintFill for 4-conn clusters (in my implementation).
! The following fragment converts all 8-conn pixels to 4-conn before calling TintFill.

           if (l8conn) then
             call convert84 ( inmx, nx, ny, funb, worker, filename(1:fnlen1), cverbose, iotty, almostzero, .false. )
           endif
           
           do j=1,ny
             do i=1,nx
               funb(i,j) = worker(i,j)  !<-- Remember worker to restore it before size measurements
               nsegmo(i,j) = nsegm(i,j)
               nsegm(i,j) = 0.0d0
             enddo
           enddo
           nextro = nextr
           nextr = 0

           do j=1,ny
             do i=1,nx
               if (worker(i,j) .gt. almostzero .and. nsegm(i,j) .le. almostzero .and. nsegmo(i,j) .le. almostzero) then
                 nextr = nextr + 1
                 if (nextr .gt. nextrmax) then
                   write (iotty,'(/a,i8)') '   MODFITS: ERROR: Number of isolated clusters reached the maximum of ', nextrmax
                   write (iotty,'(a)'    ) '      HINT: Increase NEXTRMAX and recompile MODFITS to continue.'
                   stop 89
                 endif
                 newpixval = dble ( nextr )
    
                 call fillarea ( .false., i, j, nextrmx, nextr, newpixval, 'positive', zeros, boupixval, nx, ny, n1x, n2x
     &                         , nymin, nymax, nxmn(nextr), nxmx(nextr), nymn(nextr), nymx(nextr), worker, nsegm, nsegm
     &                         , mxco, myco, deltx, delty, ntouching, nsx )
               endif
             enddo
           enddo
           
           allocate ( momxco(nextr), momyco(nextr), afwhm(nextr), bfwhm(nextr), atheta(nextr), afoot(nextr), bfoot(nextr)
     &              , equivrad(nextr), elongation(nextr), stat=irc )
     
           if (irc .ne. 0) then
             write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (13).'
             stop 13
           endif
    
           do k=1,nextr
             nonzero(k) = 0.0d0
             do j=nymn(k),nymx(k)
               do i=nxmn(k),nxmx(k)
                 if (nint ( nsegm(i,j) ) .eq. k) then
                   nonzero(k) = nonzero(k) + 1.0d0
                 endif
               enddo
             enddo                   
           enddo

           if (nextr .gt. 0) then
             call sizemeasure ( iotty, nextr, nxmn, nxmx, nymn, nymx, inmx, nx, ny, funb, nsegm, dx, dy, momxco, momyco
     &                        , afwhm, bfwhm, atheta, afoot, bfoot, equivrad, elongation, nonzero, almostzero )
     
             if (cverbose .eq. '-verb2') then
               write (iotty,'(a)') '   lev       k  nonzero    level   elongatio  sparsity  factpix nsrco ixco jyco'
     &                           //' nxmn nxmx toremove'
             endif
           endif

           do k=1,nextr
             ixco = nint ( momxco(k) )
             jyco = nint ( momyco(k) )
             nsrco = 0
             goodtoremove = 0.0d0
             factpix = 1.0d0
             
             sparsity = afwhm(k) * bfwhm(k) / equivrad(k)**2
             lsmall = nonzero(k) .le. minpix

             if (lcutpeaks) then
               lsources = lsmall
             endif

             if (lcutsources) then
               lfilamen = .false.
               lcompact = elongation(k) .le. elongmaxsrc .and. sparsity .le. sparsmaxsrc
               if (minpix .lt. almostzero) then
                 lsources = lcompact
               else
                 lsources = lcompact .and. lsmall
               endif
             endif
             
             if (lcutfilaments) then
               lcompact = .false.
               lsmall = .false.
               lfilamen = elongation(k) .gt. elongminfil .or. sparsity .gt. value * sparsminfil
             endif

             if (lcutshapes) then
               lcompact = .false.
               lshapes = elongation(k) .le. elongminfil .and. sparsity .le. sparsminfil
               lshapes = lshapes .or. lsmall
             endif

! Only small and compact shapes are removed when cutting sources; and filamentary shapes when cutting filaments.

             if ((lcutsources .and. lsources) .or. 
     &           (lcutfilaments .and. lfilamen) .or.
     &           (lcutshapes .and. lshapes)) then

               if (lcutshapes) then
                 do j=nymn(k),nymx(k)
                   do i=nxmn(k),nxmx(k)
                     if (nint ( nsegm(i,j) ) .eq. k) then
                       if (lev .eq. nlevpos) then
                         funa(i,j) = 0.0d0
                       else
                         funa(i,j) = level(lev)
                       endif
                     endif
                   enddo
                 enddo
               else
                 do j=nymn(k),nymx(k)
                   do i=nxmn(k),nxmx(k)
                     if (nint ( nsegm(i,j) ) .eq. k) then
                       if (lev .eq. 1) then
                         shapes(i,j) = funa(i,j) - level(lev)
                       else
                         if (lev .eq. nlevpos) then
                           shapes(i,j) = shapes(i,j) + min ( max ( funa(i,j), 0.0d0 ), level(lev-1) )
                         else
                           shapes(i,j) = shapes(i,j) + min ( funa(i,j), level(lev-1) ) - level(lev)
                         endif
                       endif
                     endif
                   enddo
                 enddo
               endif
               goodtoremove = 1.0d0
             endif
             
             if (cverbose .eq. '-verb2') then
               write (iotty,'(i6,i8,5(1pe10.3),6i5,5l)') lev, k, nonzero(k), level(lev), elongation(k), sparsity, factpix, nsrco
     &               , ixco, jyco, nxmn(k), nxmx(k), nint ( goodtoremove ), lcutsources, lcutfilaments, lsmall, lcompact, lfilamen
             endif

! To debug, draw ellipses around clusters of connected pixels.

!!             if (ldebug .and. level(max(lev-1,1)) .ge. levelx .and. level(lev) .lt. levelx) then
!!               dxym = 4.0d0 * sqrt ( dx * dy )
!!               do j=nymn(k),nymx(k)
!!                 dblm = dble ( j )
!!                 do i=nxmn(k),nxmx(k)
!!                   dbll = dble ( i )
!!                   if (inellipse ( dbll, dblm, dx, dy, momxco(k), momyco(k), afoot(k), bfoot(k), atheta(k) )) then
!!                     funa(i,j) = max ( funa(i,j), level(lev) ) + funmax / 3.0d0
!!                   endif
!!                 enddo
!!               enddo
!!             endif
           enddo

           deallocate ( momxco, momyco, afwhm, bfwhm, afoot, bfoot, atheta, equivrad, elongation )
    
!!           if (ldebug .and. level(max(lev-1,1)) .ge. levelx .and. level(lev) .lt. levelx) then
!!             exit
!!           endif
         enddo

         if (lcutshapes) then
           do j=1,ny
             do i=1,nx
               if (funa(i,j) .gt. almostzero) then
                 funb(i,j) = shapes(i,j)
               else
                 funb(i,j) = 0.0d0
               endif
             enddo
           enddo
         else
           do j=1,ny
             do i=1,nx
               funb(i,j) = funa(i,j) - shapes(i,j)
             enddo
           enddo
         endif
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamax = funmax
         datamin = funmin

         deallocate ( worker, nsegm, shapes, n1x, n2x, mxco, myco, deltx, delty, ntouching, nsx, nxmn, nxmx, nymn, nymx, nonzero
     &              , nonzeroo, nsegmo, nsourceo )

         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           ic1 = index ( param(1:ipar-1), '.' )
           if (ic1 .gt. 0) param(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.cli'//cvalue(1:icva-1)//'l'//param(1:ipar-1)//'.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:7) .eq. 'segment' .or. maction(1:3) .eq. 'seg') then !! .or.
!!     &     maction(1:iact) .eq. 'simfil' .or. maction(1:iact) .eq. 'sim') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         l4conn = .false.
         l8conn = .false.
         lskeletons = .false.
         if (maction(8:8) .eq. '4' .or. maction(4:4) .eq. '4') l4conn = .true.
         if (maction(8:8) .eq. '8' .or. maction(4:4) .eq. '8') l8conn = .true.
         if (maction(9:10) .eq. 'sk') lskeletons = .true.
         minpix = dble ( nint ( value ) )
!!         if (maction(1:iact) .eq. 'simfil' .or. maction(1:iact) .eq. 'sim') then
!!           l8conn = .true.
!!           minpix = 10
!!         endif
         if (.not.l4conn .and. .not.l8conn) then
           write (iotty,'(/a)') '   MODFITS: ERROR: Invalid action: "'//maction(1:iact)//'".'
           write (iotty,'(a)' ) '      HINT: Known actions: "segment4", "segment8", "segment8sk", or "seg4", "seg8", "seg8sk".'
           stop 99
         endif
         if (lskeletons .and. abs ( minpix ) .lt. almostzero) then
           write (iotty,'(/a)') '   MODFITS: ERROR: Invalid command-line value of MINPIX for: "'//maction(1:iact)//'".'
           write (iotty,'(a)' ) '      HINT: For segmentation of skeletons, MINPIX must not be zero (e.g., 10).'
           stop 99
         endif
         if (cverbose .eq. '-verb2') then
           iminp = int ( log10 ( dble ( max ( nint ( minpix ), 1 ) ) ) ) + 1
           write (cminpx,'(i7)') nint ( minpix )
           if (lskeletons) then
             write (iotty,'(a)') '   Segment image in skeletons connecting at least '//cminpx(7-iminp+1:7)//' pixels'
           else
             write (iotty,'(a)') '   Segment image in shapes connecting at least '//cminpx(7-iminp+1:7)//' pixels'
           endif
         endif  
         changed = .false.

         if (lskeletons) then
           nextro = 0
           itermax = 40
         else
           itermax = 1
         endif

         if (abs ( funmax ) .gt. almostzero) then

           allocate ( worker(nx,ny), nsegm(nx,ny), n1x(ny), n2x(ny), mxco(1), myco(1), deltx(1), delty(1)
     &              , nxmn(nextrmax), nxmx(nextrmax), nymn(nextrmax), nymx(nextrmax), ntouching(1,1), nsx(1), stat=irc )

           if (irc .ne. 0) then
             write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (12).'
             stop 12
           endif
          
           do j=1,ny
             do i=1,nx
               worker(i,j) = funb(i,j)
             enddo
           enddo
           
! Repeat segmentation several times, until the number of segments does not change.

           do iter=1,itermax

! To segment images, I use TintFill which works for 4-conn clusters only (in my implementation).
! The following call converts all 8-conn pixels to 4-conn before calling TintFill.

             if (l8conn) then
               call convert84 ( inmx, nx, ny, funb, worker, filename(1:fnlen1), cverbose, iotty, almostzero, lskeletons )
             endif
             if (cverbose .eq. '-verb2') then
               write (iotty,'(a,i2,a$)') '   Applying 4-conn TintFill algorithm (iter:', iter, ') '
             endif
             nymin = 1
             nymax = ny
             do j=1,ny
               n1x(j) = 1
               n2x(j) = nx
               do i=1,nx
                 funb(i,j) = worker(i,j)
                 worker(i,j) = funb(i,j) / funmax / 1.1d0
                 funa(i,j) = funb(i,j)
                 nsegm(i,j) = 0.0d0
               enddo
             enddo
             
             nextrmx = 1
             zeros = 1.0d-20
             mxco(1) = 1.0d0
             myco(1) = 1.0d0
             deltx(1) = 0.0d0
             delty(1) = 0.0d0
             ntouching(1,1) = 0
             nsx(1) = 0
             nextr = 0
             do j=1,ny
               do i=1,nx
                 if (worker(i,j) .gt. almostzero .and. nsegm(i,j) .le. almostzero) then
                   nextr = nextr + 1
                   newpixval = dble ( nextr )
            
                   call fillarea ( .false., i, j, nextrmx, nextr, newpixval, 'positive', zeros, boupixval, nx, ny, n1x, n2x, nymin
     &                           , nymax, nxmn(nextr), nxmx(nextr), nymn(nextr), nymx(nextr), worker, nsegm, nsegm, mxco, myco
     &                           , deltx, delty, ntouching, nsx )
                 endif
               enddo
             enddo                   
            
             if (cverbose .eq. '-verb2') then
               inrem = int ( log10 ( dble ( max ( nextr, 1 ) ) ) ) + 1
               write (cnrem,'(i7)') nextr
               write (iotty,'(a)') cnrem(7-inrem+1:7)//' clusters'
             endif

             allocate ( nonzero(nextr), np(nextr), stat=irc )
            
             if (irc .ne. 0) then
               write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (39).'
               stop 39
             endif
            
             npxmax = 0
             do k=1,nextr
               nonzero(k) = 0.0d0
               np(k) = 0
               do j=nymn(k),nymx(k)
                 do i=nxmn(k),nxmx(k)
                   if (nint ( nsegm(i,j) ) .eq. k) then
                     nonzero(k) = nonzero(k) + 1.0d0
                     np(k) = np(k) + 1
                   endif
                 enddo
               enddo                   
               if (np(k) .gt. npxmax) npxmax = np(k)
             enddo
             write (cnmax,'(i5)') npxmax

             if (lskeletons) then
              
               allocate ( npixfil(nextr), xpix(npxmax), ypix(npxmax), pixd(npxmax), xpixok(nextr,npxmax), ypixok(nextr,npxmax)
     &                  , pixdist(nx,ny), stat=irc )
            
               if (irc .ne. 0) then
                 write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (40).'
                 stop 40
               endif
            
! Trace all pixels of skeletons and separate their branches.
            
               call traceskels ( nx, ny, nsegm, filename(1:fnlen1), iotty, nextr, npixfil, np, xpix, ypix, xpixok, ypixok, pixd
     &                         , pixdist, cverbose, almostzero, lnobranches )
            
! The above call to TRACESKELS may remove single pixels off skeletons, making them again 8-conn.
! Next removes those 8-conn pixels and separates branches completely, making skeletons 4-conn.

               do j=1,ny
                 jm1 = max ( j - 1, 1 )
                 jp1 = min ( j + 1, ny )
                 do i=1,nx
                   im1 = max ( i - 1, 1 )
                   ip1 = min ( i + 1, nx )
                   if (nsegm(i,j) .gt. almostzero) then
                     numpix4 = 0
                     numpix8 = 0
                     if (nsegm(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                     if (nsegm(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                     if (nsegm(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                     if (nsegm(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                     if (nsegm(im1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
                     if (nsegm(ip1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
                     if (nsegm(im1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
                     if (nsegm(ip1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
              
                     if ((numpix4 .eq. 1 .and. numpix8 .eq. 2) .or. (numpix4 .eq. 0 .and. numpix8 .eq. 1) .or.
     &                   (nsegm(i,jm1) .gt. almostzero .and. nsegm(im1,jp1) .gt. almostzero .and. nsegm(ip1,jp1) + 
     &                    nsegm(ip1,j) + nsegm(ip1,jm1) + nsegm(im1,j) + nsegm(im1,jm1) + nsegm(i,jp1) .lt. almostzero) .or.
     &                   (nsegm(i,jm1) .gt. almostzero .and. nsegm(ip1,jp1) .gt. almostzero .and. nsegm(im1,jp1) + 
     &                    nsegm(im1,j) + nsegm(im1,jm1) + nsegm(ip1,j) + nsegm(ip1,jm1) + nsegm(i,jp1) .lt. almostzero) .or.
     &                   (nsegm(i,jp1) .gt. almostzero .and. nsegm(im1,jm1) .gt. almostzero .and. nsegm(ip1,jp1) + 
     &                    nsegm(ip1,j) + nsegm(ip1,jm1) + nsegm(im1,j) + nsegm(im1,jp1) + nsegm(i,jm1) .lt. almostzero) .or.
     &                   (nsegm(i,jp1) .gt. almostzero .and. nsegm(ip1,jm1) .gt. almostzero .and. nsegm(im1,jp1) + 
     &                    nsegm(im1,j) + nsegm(im1,jm1) + nsegm(ip1,j) + nsegm(ip1,jp1) + nsegm(i,jm1) .lt. almostzero) .or.
     &                   (nsegm(im1,j) .gt. almostzero .and. nsegm(ip1,jp1) .gt. almostzero .and. nsegm(im1,jm1) + 
     &                    nsegm(i,jm1) + nsegm(ip1,jm1) + nsegm(im1,jp1) + nsegm(i,jp1) + nsegm(ip1,j) .lt. almostzero) .or.
     &                   (nsegm(im1,j) .gt. almostzero .and. nsegm(ip1,jm1) .gt. almostzero .and. nsegm(im1,jp1) + 
     &                    nsegm(i,jp1) + nsegm(ip1,jp1) + nsegm(im1,jm1) + nsegm(i,jm1) + nsegm(ip1,j) .lt. almostzero) .or.
     &                   (nsegm(ip1,j) .gt. almostzero .and. nsegm(im1,jm1) .gt. almostzero .and. nsegm(im1,jp1) + 
     &                    nsegm(i,jp1) + nsegm(ip1,jp1) + nsegm(i,jm1) + nsegm(ip1,jm1) + nsegm(im1,j) .lt. almostzero) .or.
     &                   (nsegm(ip1,j) .gt. almostzero .and. nsegm(im1,jp1) .gt. almostzero .and. nsegm(im1,jm1) + 
     &                    nsegm(i,jm1) + nsegm(ip1,jm1) + nsegm(i,jp1) + nsegm(ip1,jp1) + nsegm(im1,j) .lt. almostzero)) then
              
                       if (cverbose .eq. '-verb2' .and. iotty .gt. 0) then
                         write (iotty,'(a,2(i4,a))') '   Zeroing 8-conn pixel (', i, ',', j,') induced by branch separation'
                         lnobranches = .false.
                       endif
                       nsegm(i,j) = 0.0d0
                     endif
                   endif
                 enddo
               enddo
             endif

! The above call to TRACESKELS may separate short segments of skeletons; next removes such short skeletons.

             nrem = 0
             do k=1,nextr
               if (nonzero(k) .lt. minpix) then
                 nrem = nrem + 1
                 do j=nymn(k),nymx(k)
                   do i=nxmn(k),nxmx(k)
                     if (nint ( nsegm(i,j) ) .eq. k) nsegm(i,j) = 0.0d0
                   enddo
                 enddo                   
               endif
             enddo

             if (cverbose .eq. '-verb2') then
               inrem = int ( log10 ( dble ( max ( nrem, 1 ) ) ) ) + 1
               iminp = int ( log10 ( dble ( max ( nint ( minpix ), 1 ) ) ) ) + 1
               write (cnrem,'(i7)') nrem
               write (cminpx,'(i7)') nint ( minpix )
               write (iotty,'(a,i3,a)') '   Removed '//cnrem(7-inrem+1:7)//' clusters connecting fewer than '//cminpx(7-iminp+1:7)
     &                                //' pixels'
             endif

             funmin =  1.0d+30
             funmax = -1.0d+30
             do j=1,ny
               do i=1,nx
                 worker(i,j) = nsegm(i,j)
                 funb(i,j) = nsegm(i,j)
                 if (funa(i,j) .gt. funmax) funmax = worker(i,j)
                 if (funa(i,j) .lt. funmin) funmin = worker(i,j)
               enddo
             enddo

             deallocate ( nonzero, np )
             if (lskeletons) then
               deallocate ( npixfil, xpix, ypix, pixd, xpixok, ypixok, pixdist )
               if (nextr .eq. nextro .and. lnobranches) exit
               nextro = nextr
             endif               
           enddo

           if (lskeletons) then
             if (nextr .ne. nextro .or. .not.lnobranches) then
               write (iotty,'(/a)') '   MODFITS: ERROR: No convergence of segmentation iterations.'
               write (iotty,'(a)' ) '   MODFITS: Image: '//filename(1:fnlen1)
               stop 99
             endif
           endif

           do j=1,ny
             do i=1,nx
               funb(i,j) = nsegm(i,j)
             enddo
           enddo
           datamax = funmax
           datamin = funmin
           
           deallocate ( worker, nsegm, n1x, n2x, mxco, myco, deltx, delty, nxmn, nxmx, nymn, nymx, ntouching, nsx )
          
           if (lnoname) then
             ic1 = index ( cvalue(1:icva-1), '.' )
             if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
             outname = filename(isp1:fnlen1-5)//'.seg'//cvalue(1:icva-1)//'.fits'
           endif
           if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
           if (cverbose .eq. '-verb2') write (iotty,*) datamax
           if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
           if (cverbose .eq. '-verb2') write (iotty,*) datamin
           changed = .true.
         endif
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:5) .eq. 'erase' .or. maction(1:3) .eq. 'era') then
       
         nconn = 0
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Erasing some pixels in skeletons'

! Erase some pixels in skeleton.

!!         do j=1,ny
!!           jm1 = max ( j - 1, 1 )
!!           jp1 = min ( j + 1, ny )
!!           do i=1,nx
!!             if (funa(i,j) .gt. almostzero) then
!!               im1 = max ( i - 1, 1 )
!!               ip1 = min ( i + 1, nx )
!!               numpix4 = 0
!!               numpix8 = 0
!!               if (funa(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
!!               if (funa(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
!!               if (funa(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
!!               if (funa(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
!!               if (funa(im1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
!!               if (funa(ip1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
!!               if (funa(im1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
!!               if (funa(ip1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
!!         
!!               if (.not.(numpix4 .le. 2 .and. numpix8 .le. 2)) then
!!!!               if (.not.(numpix4 .le. 3 .and. numpix8 .le. 2)) then
!!                 funa(i,j) = 0.0d0
!!                 funb(i,j) = 0.0d0
!!                 nconn = nconn + 1
!!               endif
!!             endif
!!           enddo
!!         enddo

         do j=1,ny
           jm1 = max ( j - 1, 1 )
           jp1 = min ( j + 1, ny )
           do i=1,nx
             im1 = max ( i - 1, 1 )
             ip1 = min ( i + 1, nx )
             if (funa(i,j) .gt. almostzero) then
               if (funa(im1,j) .gt. almostzero .and. funa(im1,jm1) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(i,j) = 0.0d0
               else if (funa(ip1,j) .gt. almostzero .and. funa(ip1,jm1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(i,j) = 0.0d0
               else if (funa(i,jm1) .gt. almostzero .and. funa(ip1,jm1) .gt. almostzero .and. funa(im1,jm1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(i,j) = 0.0d0
               else if (funa(i,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(i,j) = 0.0d0
               endif
             endif
           enddo
         enddo
         do j=1,ny
           jm1 = max ( j - 1, 1 )
           jp1 = min ( j + 1, ny )
           do i=1,nx
             im1 = max ( i - 1, 1 )
             ip1 = min ( i + 1, nx )
             if (funa(i,j) .gt. almostzero) then
               if (funa(im1,jm1) .gt. almostzero .and. funa(i,jp1) .gt. almostzero .and. funa(ip1,j) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(ip1,jp1) = funa(i,j)
                 funa(i,j) = 0.0d0
               else if (funa(im1,jp1) .gt. almostzero .and. funa(i,jm1) .gt. almostzero .and. funa(ip1,j) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(ip1,jm1) = funa(i,j)
                 funa(i,j) = 0.0d0
               else if (funa(ip1,jp1) .gt. almostzero .and. funa(im1,j) .gt. almostzero .and. funa(i,jm1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(im1,jm1) = funa(i,j)
                 funa(i,j) = 0.0d0
               else if (funa(ip1,jm1) .gt. almostzero .and. funa(im1,j) .gt. almostzero .and. funa(i,jp1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(im1,jp1) = funa(i,j)
                 funa(i,j) = 0.0d0
               endif
             endif
           enddo
         enddo
         do j=1,ny
           jm1 = max ( j - 1, 1 )
           jp1 = min ( j + 1, ny )
           do i=1,nx
             im1 = max ( i - 1, 1 )
             ip1 = min ( i + 1, nx )
             if (funa(i,j) .gt. almostzero) then
               if (funa(im1,jm1) .gt. almostzero .and. funa(ip1,jm1) .gt. almostzero .and. funa(i,jp1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(i,jm1) = funa(i,j)
                 funa(i,j) = 0.0d0
               else if (funa(im1,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and. funa(i,jm1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(i,jp1) = funa(i,j)
                 funa(i,j) = 0.0d0
               else if (funa(ip1,jm1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and. funa(im1,j) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(ip1,j) = funa(i,j)
                 funa(i,j) = 0.0d0
               else if (funa(im1,jm1) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and. funa(ip1,j) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(im1,j) = funa(i,j)
                 funa(i,j) = 0.0d0
               endif
             endif
           enddo
         enddo
         do j=1,ny
           jm1 = max ( j - 1, 1 )
           jp1 = min ( j + 1, ny )
           do i=1,nx
             im1 = max ( i - 1, 1 )
             ip1 = min ( i + 1, nx )
             if (funa(i,j) .gt. almostzero) then
               if (funa(im1,jm1) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(im1,jp1) = 0.0d0
               else if (funa(im1,jp1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and. funa(ip1,jm1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(ip1,jp1) = 0.0d0
               else if (funa(im1,jm1) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and. funa(ip1,jm1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(im1,jm1) = 0.0d0
               else if (funa(im1,jm1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero .and. funa(ip1,jm1) .gt. almostzero) then
                 nconn = nconn + 1
                 funa(ip1,jm1) = 0.0d0
               endif
             endif
           enddo
         enddo
         if (cverbose .eq. '-verb2') then
           inrem = int ( log10 ( dble ( max ( nconn, 1 ) ) ) ) + 1
           write (cnrem,'(i7)') nconn
           write (iotty,'(a)') '   Erased '//cnrem(7-inrem+1:7)//' pixels off skeletons'
         endif
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j)
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamax = funmax
         datamin = funmin
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.era.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:7) .eq. 'reconnect' .or. maction(1:3) .eq. 'rec') then
       
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Reconnecting skeletons separated by one-pixel breaks'

         ncon = max ( nint ( value ), 1 )
         nconn = 0

! Reconnect some breaks in a skeleton, if the breaks are just one pixel wide.

         do j=1,ny
           jm1 = max ( j - 1, 1 )
           jp1 = min ( j + 1, ny )
           jm2 = max ( j - 2, 1 )
           jp2 = min ( j + 2, ny )
           do i=1,nx
             if (funb(i,j) .lt. almostzero) then
               im1 = max ( i - 1, 1 )
               ip1 = min ( i + 1, nx )
               im2 = max ( i - 2, 1 )
               ip2 = min ( i + 2, nx )
               numpix4 = 0
               numpix8 = 0
               if (funb(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
               if (funb(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
               if (funb(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
               if (funb(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
               if (funb(im1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
               if (funb(ip1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
               if (funb(im1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
               if (funb(ip1,jm1) .gt. almostzero) numpix8 = numpix8 + 1

! Reconnect only *two* disconnected pixels in the 8-pixel environment of the pixel (i,j) being considered.

               if (numpix4 + numpix8 .eq. 2) then

! 4-connect pixels.

                 if (funb(im1,j) .gt. almostzero .and. funb(im1,jp1) .lt. almostzero .and. funb(im1,jm1) .lt. almostzero .and.
     &               funb(ip1,j) .gt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and.
     &               funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then 
                   nconn = nconn + 1
                   funa(i,j) = 1.0d0
                 endif
                 if (funb(i,jm1) .gt. almostzero .and. funb(im1,jp1) .lt. almostzero .and. funb(im1,jm1) .lt. almostzero .and.
     &               funb(i,jp1) .gt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and.
     &               funb(im1,j) .lt. almostzero .and. funb(ip1,j) .lt. almostzero) then 
                   nconn = nconn + 1
                   funa(i,j) = 1.0d0
                 endif
          
! 8-connect pixels.

                 if (funb(im1,jp1) .gt. almostzero .and. funb(im1,j) .lt. almostzero .and. funb(im1,jm1) .lt. almostzero .and.
     &               funb(ip1,jm1) .gt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and. funb(ip1,j) .lt. almostzero .and.
     &               funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then 
                   nconn = nconn + 1
                   funa(i,j) = 1.0d0
                 endif
                 if (funb(im1,jp1) .gt. almostzero .and. funb(im1,j) .lt. almostzero .and. funb(im1,jm1) .lt. almostzero .and.
     &               funb(ip1,jp1) .gt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and. funb(ip1,j) .lt. almostzero .and.
     &               funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then 
                   nconn = nconn + 1
                   funa(i,j) = 1.0d0
                 endif
                 if (funb(im1,jm1) .gt. almostzero .and. funb(im1,j) .lt. almostzero .and. funb(im1,jp1) .lt. almostzero .and.
     &               funb(ip1,jp1) .gt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and. funb(ip1,j) .lt. almostzero .and.
     &               funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then 
                   nconn = nconn + 1
                   funa(i,j) = 1.0d0
                 endif
                 if (funb(im1,jm1) .gt. almostzero .and. funb(im1,j) .lt. almostzero .and. funb(im1,jp1) .lt. almostzero .and.
     &               funb(ip1,jm1) .gt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and. funb(ip1,j) .lt. almostzero .and.
     &               funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then 
                   nconn = nconn + 1
                   funa(i,j) = 1.0d0
                 endif

! 4- and 8-connect mixed pixels.

                 if (funb(im1,j) .gt. almostzero .and. funb(im1,jp1) .lt. almostzero .and. funb(im1,jm1) .lt. almostzero .and.
     &               funb(ip1,j) .gt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and.
     &               funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then 
                   numpix4 = 0
                   numpix8 = 0
                   if (funb(i  ,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip2,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,jm1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip2,jp1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,jp1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip2,jm1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 + numpix8 .eq. 1) then
                     nconn = nconn + 1
                     funa(i,j) = 1.0d0
                   endif
                 endif
                 if (funb(im1,j) .gt. almostzero .and. funb(im1,jp1) .lt. almostzero .and. funb(im1,jm1) .lt. almostzero .and.
     &               funb(ip1,jm1) .gt. almostzero .and. funb(ip1,j) .lt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and.
     &               funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then 
                   numpix4 = 0
                   numpix8 = 0
                   if (funb(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip2,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,jm2) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,jm2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip2,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip2,jm2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 + numpix8 .eq. 1) then
                     nconn = nconn + 1
                     funa(i,j) = 1.0d0
                   endif
                 endif
                 if (funb(im1,j) .gt. almostzero .and. funb(im1,jp1) .lt. almostzero .and. funb(im1,jm1) .lt. almostzero .and.
     &               funb(ip1,jp1) .gt. almostzero .and. funb(ip1,j) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and.
     &               funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then 
                   numpix4 = 0
                   numpix8 = 0
                   if (funb(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip2,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,jp2) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip2,jp2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,jp2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip2,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 + numpix8 .eq. 1) then
                     nconn = nconn + 1
                     funa(i,j) = 1.0d0
                   endif
                 endif
                 if (funb(ip1,j) .gt. almostzero .and. funb(im1,j) .lt. almostzero .and. funb(im1,jp1) .lt. almostzero .and.
     &               funb(im1,jm1) .gt. almostzero .and. funb(im1,j) .lt. almostzero .and. funb(im1,jp1) .lt. almostzero .and.
     &               funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then 
                   numpix4 = 0
                   numpix8 = 0
                   if (funb(im2,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im1,jm2) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im2,jm2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(im2,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,jm2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 + numpix8 .eq. 1) then
                     nconn = nconn + 1
                     funa(i,j) = 1.0d0
                   endif
                 endif
                 if (funb(ip1,j) .gt. almostzero .and. funb(im1,j) .lt. almostzero .and. funb(im1,jp1) .lt. almostzero .and.
     &               funb(im1,jp1) .gt. almostzero .and. funb(im1,j) .lt. almostzero .and. funb(im1,jm1) .lt. almostzero .and.
     &               funb(i,jm1) .lt. almostzero .and. funb(i,jp1) .lt. almostzero) then 
                   numpix4 = 0
                   numpix8 = 0
                   if (funb(im2,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im1,jp2) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im2,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,jp2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(im2,jp2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 + numpix8 .eq. 1) then
                     nconn = nconn + 1
                     funa(i,j) = 1.0d0
                   endif
                 endif

                 if (funb(i,jm1) .gt. almostzero .and. funb(im1,jm1) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and.
     &               funb(i,jp1) .gt. almostzero .and. funb(im1,jp1) .lt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and. 
     &               funb(im1,j) .lt. almostzero .and. funb(ip1,j) .lt. almostzero) then 
                   numpix4 = 0
                   numpix8 = 0
                   if (funb(im1,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,jp2) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im1,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip1,jp2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(im1,jp2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip1,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 + numpix8 .eq. 1) then
                     nconn = nconn + 1
                     funa(i,j) = 1.0d0
                   endif
                 endif
                 if (funb(i,jm1) .gt. almostzero .and. funb(im1,jm1) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and.
     &               funb(ip1,jp1) .gt. almostzero .and. funb(i,jp1) .lt. almostzero .and. funb(im1,jp1) .lt. almostzero .and.
     &               funb(im1,j) .lt. almostzero .and. funb(ip1,j) .lt. almostzero) then 
                   numpix4 = 0
                   numpix8 = 0
                   if (funb(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip2,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,jp2) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip2,jp2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,jp2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip2,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 + numpix8 .eq. 1) then
                     nconn = nconn + 1
                     funa(i,j) = 1.0d0
                   endif
                 endif
                 if (funb(i,jm1) .gt. almostzero .and. funb(im1,jm1) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and.
     &               funb(im1,jp1) .gt. almostzero .and. funb(i,jp1) .lt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and.
     &               funb(im1,j) .lt. almostzero .and. funb(ip1,j) .lt. almostzero) then 
                   numpix4 = 0
                   numpix8 = 0
                   if (funb(im2,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im1,jp2) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im2,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,jp2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(im2,jp2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 + numpix8 .eq. 1) then
                     nconn = nconn + 1
                     funa(i,j) = 1.0d0
                   endif
                 endif
                 if (funb(i,jp1) .gt. almostzero .and. funb(im1,jp1) .lt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and.
     &               funb(ip1,jm1) .gt. almostzero .and. funb(i,jm1) .lt. almostzero .and. funb(im1,jm1) .lt. almostzero .and. 
     &               funb(im1,j) .lt. almostzero .and. funb(ip1,j) .lt. almostzero) then 
                   numpix4 = 0
                   numpix8 = 0
                   if (funb(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip2,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,jm2) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,jm2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip2,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(ip2,jm2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 + numpix8 .eq. 1) then
                     nconn = nconn + 1
                     funa(i,j) = 1.0d0
                   endif
                 endif
                 if (funb(i,jp1) .gt. almostzero .and. funb(im1,jp1) .lt. almostzero .and. funb(ip1,jp1) .lt. almostzero .and.
     &               funb(im1,jm1) .gt. almostzero .and. funb(i,jm1) .lt. almostzero .and. funb(ip1,jm1) .lt. almostzero .and.
     &               funb(im1,j) .lt. almostzero .and. funb(ip1,j) .lt. almostzero) then 
                   numpix4 = 0
                   numpix8 = 0
                   if (funb(im2,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(i  ,jm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im1,jm2) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (funb(im2,jm2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(im2,j  ) .gt. almostzero) numpix8 = numpix8 + 1
                   if (funb(i  ,jm2) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 + numpix8 .eq. 1) then
                     nconn = nconn + 1
                     funa(i,j) = 1.0d0
                   endif
                 endif

               endif
             endif
           enddo
         enddo
         if (cverbose .eq. '-verb2') then
           inrem = int ( log10 ( dble ( max ( nconn, 1 ) ) ) ) + 1
           write (cnrem,'(i7)') nconn
           write (iotty,'(a)') '   Inserted '//cnrem(7-inrem+1:7)//' pixels reconnecting skeletons'
         endif
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. almostzero) then
               funb(i,j) = funa(i,j) - dble ( nint ( funa(i,j) ) ) + 1.0d0
             else
               funb(i,j) = 0.0d0
             endif
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamax = funmax
         datamin = funmin
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.rec.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:5) .eq. 'annex' .or. maction(1:3) .eq. 'ann') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Annexing connected pixels by segmentation masks'
         nextr = nint ( funmax )

         allocate ( worker(nx,ny), nsegm(nx,ny), n1x(ny), n2x(ny), mxco(nextr), myco(nextr), deltx(nextr), delty(nextr)
     &            , ntouching(nextr,nextr), nsx(nextr), nxmn(nextr), nxmx(nextr), nymn(nextr), nymx(nextr), stat=irc )

         if (irc .ne. 0) then
           write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (12).'
           stop 12
         endif

         do k=1,nextr
           do j=1,ny
             do i=1,nx
               if (abs ( funa(i,j) - dble ( k ) ) .lt. almostzero) then
                 mxco(k) = dble ( i )
                 myco(k) = dble ( j )
                 goto 111
               endif
             enddo
           enddo
 111       continue
         enddo
         do k=1,nextr
           deltx(k) = 0.0d0
           delty(k) = 0.0d0
           nsx(k) = 0
           do n=1,nextr
             ntouching(k,n) = 0
           enddo
         enddo
         do j=1,ny
           n1x(j) = 1
           n2x(j) = nx
           do i=1,nx
             funb(i,j) = 0.0d0
             worker(i,j) = min ( funa(i,j), 0.999d0 )
             nsegm(i,j) = 0.0d0
           enddo
         enddo
         nymin = 1
         nymax = ny
         zeros = 1.0d-20
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         do k=1,nextr
           newpixval = dble ( k )
           i = nint ( mxco(k) )
           j = nint ( myco(k) )
         
           call fillarea ( .false., i, j, nextr, nextr, newpixval, 'positive', zeros, boupixval, nx, ny, n1x, n2x, nymin, nymax
     &                   , nxmn(k), nxmx(k), nymn(k), nymx(k), worker, nsegm, nsegm, mxco, myco, deltx, delty, ntouching, nsx )
         enddo
         do j=1,ny
           do i=1,nx
             funb(i,j) = nsegm(i,j)
           enddo
         enddo

         deallocate ( worker, nsegm, n1x, n2x, mxco, myco, deltx, delty, ntouching, nsx, nxmn, nxmx, nymn, nymx )

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamax = funmax
         datamin = funmin
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.ann.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:8) .eq. 'renumber' .or. maction(1:3) .eq. 'ren') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             funb(i,j) = 0.0d0
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Renumbering segmentation image sequentially'
         nextr = nint ( funmax )

         ndelta = 0
         do k=1,nextr
           lmiss = .true.
           do j=1,ny
             do i=1,nx
               if (abs ( funa(i,j) - dble ( k ) ) .lt. almostzero) then
                 funb(i,j) = dble ( k - ndelta )
                 lmiss = .false.
               endif
             enddo
           enddo
           if (lmiss) ndelta = ndelta + 1
         enddo

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamax = funmax
         datamin = funmin
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:5) .eq. 'deriv' .or. maction(1:3) .eq. 'der') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             funb(i,j) = 0.0d0
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (maction(1:6) .eq. 'deriv1' .or. maction(1:4) .eq. 'der1') then
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing 1st derivative of the image'
         endif
         if (maction(1:6) .eq. 'deriv2' .or. maction(1:4) .eq. 'der2') then
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing 2nd derivative of the image'
         endif

         deriv1a = 0.0d0
         deriv1b = 0.0d0
         deriv1c = 0.0d0
         deriv1d = 0.0d0
         deriv1e = 0.0d0
         deriv1f = 0.0d0
         deriv1g = 0.0d0
         deriv1h = 0.0d0
         do j=1,ny
           jm1 = max (  1, j - 1 )
           jp1 = min ( ny, j + 1 )
           do i=1,nx
             im1 = max (  1, i - 1 )
             ip1 = min ( nx, i + 1 )
             deriv1a = (funa(i,j) - funa(im1,j)) / dx
             deriv1b = (funa(ip1,j) - funa(i,j)) / dx
             deriv1c = (funa(i,j) - funa(i,jm1)) / dy
             deriv1d = (funa(i,jp1) - funa(i,j)) / dy
             deriv1e = (funa(i,j) - funa(im1,jm1)) / sqrt ( dx**2 + dy**2 )
             deriv1f = (funa(ip1,jp1) - funa(i,j)) / sqrt ( dx**2 + dy**2 )
             deriv1g = (funa(i,j) - funa(im1,jp1)) / sqrt ( dx**2 + dy**2 )
             deriv1h = (funa(ip1,jm1) - funa(i,j)) / sqrt ( dx**2 + dy**2 )

             if (maction(1:6) .eq. 'deriv1' .or. maction(1:4) .eq. 'der1') then
               funb(i,j) = (abs ( deriv1a ) + abs ( deriv1b ) + abs ( deriv1c ) + abs ( deriv1d )
     &                    + abs ( deriv1e ) + abs ( deriv1f ) + abs ( deriv1g ) + abs ( deriv1h )) / 8.0d0

!!               deriv1a = max ( deriv1a, 0.0d0 )
!!               deriv1b = max ( deriv1b, 0.0d0 )
!!               deriv1c = max ( deriv1c, 0.0d0 )
!!               deriv1d = max ( deriv1d, 0.0d0 )
!!               deriv1e = max ( deriv1e, 0.0d0 )
!!               deriv1f = max ( deriv1f, 0.0d0 )
!!               deriv1g = max ( deriv1g, 0.0d0 )
!!               deriv1h = max ( deriv1h, 0.0d0 )
!!               
!!               funb(i,j) = min ( deriv1a, deriv1b, deriv1c, deriv1d, deriv1e, deriv1f, deriv1g, deriv1h )
             endif

             if (maction(1:6) .eq. 'deriv2' .or. maction(1:4) .eq. 'der2') then
               deriv2a = (deriv1b - deriv1a) / dx**2
               deriv2b = (deriv1d - deriv1c) / dy**2
               deriv2c = (deriv1f - deriv1e) / sqrt ( dx**2 + dy**2 )
               deriv2d = (deriv1h - deriv1g) / sqrt ( dx**2 + dy**2 )

               funb(i,j) = -(deriv2a + deriv2b + deriv2c + deriv2d) / 4.0d0

!!               funb(i,j) = min ( deriv2a, deriv2b, deriv2c, deriv2d )
               
!!               if (funa(i,j) .gt. 1.0d-10) funb(i,j) = funb(i,j) / funa(i,j)
             endif
           enddo
         enddo
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamax = funmax
         datamin = funmin
         if (lnoname) then
           if (maction(1:6) .eq. 'deriv1' .or. maction(1:4) .eq. 'der1') then
             outname = filename(isp1:fnlen1-5)//'.der1.fits'
           else
             outname = filename(isp1:fnlen1-5)//'.der2.fits'
           endif
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) datamin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'invert' .or. maction(1:iact) .eq. 'inv') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Inverting the image'

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (abs ( funa(i,j) ) .gt. almostzero) then
               funb(i,j) = 1.0d0 / funa(i,j)
               if (funb(i,j) .gt. funmax) funmax = funb(i,j)
               if (funb(i,j) .lt. funmin) funmin = funb(i,j)
             endif
           enddo
         enddo
         do j=1,ny
           do i=1,nx
             if (abs ( funa(i,j) ) .le. almostzero) then
               funb(i,j) = funmax
             endif
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.inv.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'edge' .or. maction(1:iact) .eq. 'edg') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Removing inner pixels of all structures'

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j)
             if (funa(i,j) .gt. almostzero) then
               im1 = max ( i - 1, 1 )
               jm1 = max ( j - 1, 1 )
               ip1 = min ( i + 1, nx )
               jp1 = min ( j + 1, ny )
               if (funa(im1,j  ) .gt. almostzero .and. funa(ip1,j  ) .gt. almostzero .and.
     &             funa(i  ,jm1) .gt. almostzero .and. funa(i  ,jp1) .gt. almostzero .and.
     &             funa(im1,jm1) .gt. almostzero .and. funa(im1,jp1) .gt. almostzero .and.
     &             funa(ip1,jm1) .gt. almostzero .and. funa(ip1,jp1) .gt. almostzero) then
                 funb(i,j) = 0.0d0
               endif
               if (funb(i,j) .gt. funmax) funmax = funb(i,j)
               if (funb(i,j) .lt. funmin) funmin = funb(i,j)
             endif
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.edg.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif
       
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'bsrcore' .or. maction(1:iact) .eq. 'bsrc') then
 
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         if (rcore .lt. almostzero) then
           write (iotty,'(/a)') '   MODFITS: ERROR: Core radius is zero: check the FITS header.'
           stop 99
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Background-subtracting model object'

         ic = nx / 2 + 1
         jc = ny / 2 + 1
         do i=ic,nx
           rad = dble ( i - ic ) * sqrt ( dx * dy )
           if (rad .gt. rcore) then
             background = funa(i,jc)
             goto 50
           endif
         enddo
 50      continue
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             rad = sqrt ( (dble ( i - ic ) * dx)**2 + (dble ( j - jc ) * dy)**2 )
             if (rad .le. rcore) then
               funb(i,j) = funa(i,j) - background
             else
               funb(i,j) = 0.0d0
             endif
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.bsr.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'minimum' .or. maction(1:iact) .eq. 'min' .or.
     &     maction(1:iact) .eq. 'maximum' .or. maction(1:iact) .eq. 'max' .or.
     &     maction(1:iact) .eq. 'maxstdev' .or. maction(1:iact) .eq. 'maxs') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
           enddo
         enddo
         fnewmax = funmax
         fnewmin = funmin
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin

         if (maction(1:iact) .eq. 'minimum' .or. maction(1:iact) .eq. 'min' .or. 
     &       maction(1:iact) .eq. 'maximum' .or. maction(1:iact) .eq. 'max') then
           if (cverbose .eq. '-verb2')
     &         write (iotty,'(a)') '   Setting new minimum and/or maximum values (param='//param(1:ipar-1)//')'
         endif
         if (maction(1:iact) .eq. 'maxstdev' .or. maction(1:iact) .eq. 'maxs') then
           if (cverbose .eq. '-verb2')
     &         write (iotty,'(a)') '   Setting new value according to maxstdev (param='//param(1:ipar-1)//')'
         endif

         if (maction(1:iact) .eq. 'minimum' .or. maction(1:iact) .eq. 'min') then

           fnewmin = value           
           do i=1,80
             cline(i:i) = ' '
           enddo
           answer = 'Y'
           read (param(1:ipar),*) answer
           funmin =  1.0d+30
           funmax = -1.0d+30
           do j=1,ny
             do i=1,nx
               if (funa(i,j) .gt. funmax) funmax = funa(i,j)
               if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             enddo
           enddo
           if (fnewmin .lt. funmin .and. abs ( fnewmin - funmin ) .lt. 1.0d-5 * abs ( funmin )) fnewmin = funmin * 1.00001
           funmin =  1.0d+30
           funmax = -1.0d+30
           do j=1,ny
             do i=1,nx
               funb(i,j) = max ( funa(i,j), fnewmin )
               if ((answer .eq. 'Y' .or. answer .eq. 'y') .and. funb(i,j) .le. fnewmin) funb(i,j) = 0.0d0
               if (funb(i,j) .gt. funmax) funmax = funb(i,j)
               if (funb(i,j) .lt. funmin) funmin = funb(i,j)
             enddo
           enddo
           datamin = funmin
           datamax = funmax

         elseif (maction(1:iact) .eq. 'maximum' .or. maction(1:iact) .eq. 'max') then

           fnewmax = value
           do i=1,80
             cline(i:i) = ' '
           enddo
           answer = 'Y'
           read (param(1:ipar),*) answer
           funmin =  1.0d+30
           funmax = -1.0d+30
           do j=1,ny
             do i=1,nx
               funb(i,j) = min ( funa(i,j), fnewmax )
               if ((answer .eq. 'Y' .or. answer .eq. 'y') .and. funb(i,j) .ge. fnewmax) funb(i,j) = 0.0d0
               if (funb(i,j) .gt. funmax) funmax = funb(i,j)
               if (funb(i,j) .lt. funmin) funmin = funb(i,j)
             enddo
           enddo
           datamin = funmin
           datamax = funmax
         
         elseif (maction(1:iact) .eq. 'maxstdev' .or. maction(1:iact) .eq. 'maxs') then

           fnewmax = value
           do i=1,80
             cline(i:i) = ' '
           enddo
           answer = 'Y'
           read (param(1:ipar),*) xpoints
           
           npoints = nint ( xpoints )
           if (xpoints .le. 2.0d0) then
             rw2 = (xpoints + 0.5d0)**2
           else
             rw2 = (xpoints)**2
           endif
           npp1 = npoints + 1
  
           allocate ( cimask(-npoints:npoints,-npoints:npoints), stat=irc )
           if (irc .ne. 0) then
             write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (160).'
             stop 160
           endif

! The circular mask.

!!           if (cverbose .eq. '-verb2') write (iotty,'()')
           do l=-npoints,npoints
             ry2 = (dble ( l ))**2
             do k=-npoints,npoints
               rx2 = (dble ( k ))**2
               if (rx2 + ry2 .le. rw2) then
                 cimask(k,l) = 1 
               else
                 cimask(k,l) = 0
               endif
             enddo
!!             if (cverbose .eq. '-verb2') write (iotty,'(1000i2)') (cimask(k,l),k=-npoints,min(-npoints+85,npoints))
           enddo
!!           if (cverbose .eq. '-verb2') write (iotty,'()')

           funmin =  1.0d+30
           funmax = -1.0d+30
           do j=1,ny
             n1 = max ( j - npoints, 1 )
             n2 = min ( j + npoints, ny )
             do i=1,nx
               m1 = max ( i - npoints, 1 )
               m2 = min ( i + npoints, nx )
             
! First compute the mean value.

               nto = 0
               fmean = 0.0d0
               do m=m1,m2
                 do n=n1,n2
                   if (cimask(m-i,n-j) .ge. 1) then
                     nto = nto + 1
                     fmean = fmean + funa(m,n)
                   endif
                 enddo
               enddo
               if (nto .ge. 1) then
                 fmean = fmean / dble ( nto )
               endif

! Next compute the standard deviation.

               variance = 0.0d0
               do m=m1,m2
                 do n=n1,n2
                   if (cimask(m-i,n-j) .ge. 1) then
                     variance = variance + (funa(m,n) - fmean)**2
                   endif
                 enddo
               enddo
               if (nto .ge. 2) then
                 sigma = sqrt ( variance / dble ( nto - 1 ) )
               else
                 sigma = 0.0d0
               endif
               
               funb(i,j) = funa(i,j)
               if ((answer .eq. 'Y' .or. answer .eq. 'y') .and. funa(i,j) .gt. fnewmax * sigma) funb(i,j) = 0.0d0
               if (funb(i,j) .gt. funmax) funmax = funb(i,j)
               if (funb(i,j) .lt. funmin) funmin = funb(i,j)
             enddo
           enddo
           datamin = funmin
           datamax = funmax

           deallocate ( cimask )
         endif

         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.'//maction(1:iact)//'.fits'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'( ''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'( ''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'border' .or. maction(1:iact) .eq. 'bor') then

         if (l1value) then
           inrem = int ( log10 ( dble ( max ( nint ( abs ( value ) ), 1 ) ) ) ) + 1
           write (cnrem,'(i7)') nint ( abs ( value ) )
           if (nint(value) .ge. 0.0d0) then
             if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Adding border of '//cnrem(7-inrem+1:7)//' pixels around image'
           else
             if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Removing border of '//cnrem(7-inrem+1:7)//' pixels around image'
           endif
         endif
         if (l4values) then
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Changing image borders by '//arg2(1:ia2)//' pixels'
         endif

         nxo = nx
         nyo = ny
         if (l1value) then
           nx1 = nint ( value )
           nx2 = nx1
           ny1 = nx1
           ny2 = nx1
         endif
         if (l4values) then
           nx1 = nint ( value1 )
           nx2 = nint ( value2 )
           ny1 = nint ( value3 )
           ny2 = nint ( value4 )
         endif
         do j=1,inmx
           do i=1,inmx
             funa(i,j) = 0.0d0
           enddo
         enddo
         nx = nx + nx1 + nx2
         ny = ny + ny1 + ny2
         if (nx1 .ge. 0 .and. ny1 .ge. 0) then
           do j=1,nyo
             do i=1,nxo
               funa(i+nx1,j+ny1) = funb(i,j)
             enddo
           enddo
         endif
         if (nx1 .lt. 0 .and. ny1 .ge. 0) then
           do j=1,nyo
             do i=1,nxo+nx1
               funa(i,j+ny1) = funb(i-nx1,j)
             enddo
           enddo
         endif
         if (nx1 .ge. 0 .and. ny1 .lt. 0) then
           do j=1,nyo+ny1
             do i=1,nxo
               funa(i+nx1,j) = funb(i,j-ny1)
             enddo
           enddo
         endif
         if (nx1 .lt. 0 .and. ny1 .lt. 0) then
           do j=1,nyo+ny1
             do i=1,nxo+nx1
               funa(i,j) = funb(i-nx1,j-ny1)
             enddo
           enddo
         endif
         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j)
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax

         crpix1 = crpix1 + dble ( nx1 )
         crpix2 = crpix2 + dble ( ny1 )
         if (cverbose .eq. '-verb2') then
           write (cx,'(i7)') nxo
           write (cy,'(i7)') nyo
           ncx = int ( log10 ( dble ( nxo ) ) ) + 1
           ncy = int ( log10 ( dble ( nyo ) ) ) + 1
           write (iotty,'(a)') '   Old image dimensions: '//cx(7-ncx+1:7)//' '//cy(7-ncy+1:7)//' pixels'
         endif
         if (cverbose .eq. '-verb2') then
           write (cx,'(i7)') nx
           write (cy,'(i7)') ny
           ncx = int ( log10 ( dble ( nx ) ) ) + 1
           ncy = int ( log10 ( dble ( ny ) ) ) + 1
           write (iotty,'(a)') '   New image dimensions: '//cx(7-ncx+1:7)//' '//cy(7-ncy+1:7)//' pixels'
           write (iotty,'(a,1x,2(1pe14.7))') '   Minmax values in the image:', funmin, funmax
         endif

         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.bor'//cvalue(1:icva-1)//'.fits'
         endif
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'profile' .or. maction(1:iact) .eq. 'pro' .or. maction(1:iact) .eq. 'proavg') then

         if (maction(1:iact) .eq. 'proavg') then
           lazimuthaverage = .true.
           if (cverbose .eq. '-verb2')
     &         write (iotty,'(a)') '   Creating azimutally-averaged profile about (x1,y1): '//arg2(1:ia2)
         else
           lazimuthaverage = .false.
           if (cverbose .eq. '-verb2')
     &         write (iotty,'(a)') '   Creating radial profile along line (x1,y1,x2,y2): '//arg2(1:ia2)
         endif
         i2 = 1
         j2 = 1
         if (l4values) then
           xi1 = value1
           yj1 = value2
           i1 = nint ( xi1 )
           j1 = nint ( yj1 )
           if (.not.lazimuthaverage) then
             xi2 = value3
             yj2 = value4
             i2 = nint ( xi2 )
             j2 = nint ( yj2 )
           endif           
         endif
         if (i1 .lt. 1 .or. i1 .gt. nx .or. j1 .lt. 1 .or. j1 .gt. ny .or. 
     &       i2 .lt. 1 .or. i2 .gt. nx .or. j2 .lt. 1 .or. j2 .gt. ny) then
           write (iotty,'(/a)') '   MODFITS: ERROR: Line coordinates outside of image dimensions.'
           stop 99
         endif
         if (lazimuthaverage) then
           npts = nint ( max ( xi1 - 1, yj1 - 1, nx - xi1, ny - yj1 ) ) + 1
         else
           npts = nint ( max ( abs ( xi2 - xi1 ), abs ( yj2 - yj1 ) ) ) + 1
         endif

         allocate ( arg(npts), argz(npts), prof(npts), profrms(npts), slope(npts), slo(npts), slom(npts), stat=irc )
         
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   MODFITS: ERROR: Trouble allocating memory (15).'
           stop 15
         endif

         if (lazimuthaverage) then
           pt1 = 0.0d0
           do k=1,npts
             arg(k) = pt1 + sqrt ( dx * dy ) * dble ( k - 1 )
             prof(k) = 0.0d0
             profrms(k) = 0.0d0
             npt(k) = 0
           enddo
           do k=1,npts
             do m=max(j1-k-1,1),min(j1+k+1,ny)
               do l=max(i1-k-1,1),min(i1+k+1,nx)
                 
! 2024/01/10 ABM: Next IF avoids averaging the horizontal and vertical spikes when working on the Fourier amplitude image.
                 
                 if ((l .lt. i1 - 1 .or. l .gt. j1 + 1) .and. (m .lt. j1 - 1 .or. m .gt. j1 + 1)) then
                   rxy = sqrt ( dx * dy ) * sqrt ( (dble ( l ) - xi1)**2 + (dble ( m ) - yj1)**2 )
                   if (rxy .ge. arg(k) .and. rxy .lt. arg(k+1)) then 
                     npt(k) = npt(k) + 1
                     prof(k) = prof(k) + funa(l,m)
                   endif
                 endif
               enddo
             enddo
           enddo
           prof(1) = funa(i1,j1)
           prof(2) = 0.25d0 * (funa(i1,j1+1) + funa(i1-1,j1) + funa(i1-1,j1-1) + funa(i1+1,j1))
           do k=3,npts
             if (npt(k) .gt. 0) prof(k) = prof(k) / npt(k)
           enddo
           totint = 0.0d0
           do k=1,npts
             totint = totint + prof(k) * 2.0d0 * pi * dble ( k )
           enddo
           do k=1,npts
             npt(k) = 0
             do m=max(j1-k-1,1),min(j1+k+1,ny)
               do l=max(i1-k-1,1),min(i1+k+1,nx)
                 
! 2024/01/10 ABM: Next IF avoids averaging the horizontal and vertical spikes when working on the Fourier amplitude image.

                 if ((l .lt. i1 - 1 .or. l .gt. j1 + 1) .and. (m .lt. j1 - 1 .or. m .gt. j1 + 1)) then
                   rxy = sqrt ( dx * dy ) * sqrt ( (dble ( l ) - xi1)**2 + (dble ( m ) - yj1)**2 )
                   if (rxy .ge. arg(k) .and. rxy .lt. arg(k+1)) then 
                     npt(k) = npt(k) + 1
                     profrms(k) = profrms(k) + (funa(l,m) - prof(k))**2
                   endif
                 endif
               enddo
             enddo
           enddo
           totrms = 0.0d0
           do k=1,npts
             if (npt(k) .gt. 1) profrms(k) = sqrt ( profrms(k) / dble ( npt(k) - 1 ) )
             if (npt(k) .gt. 1) totrms = totrms + max ( profrms(k)**2, 1.0d-30 ) * 2.0d0 * pi * dble ( k )
             if (abs ( prof(k) ) .lt. almostzero) prof(k) = 1.0d-33
             if (abs ( profrms(k) ) .lt. almostzero) profrms(k) = 1.0d-33
           enddo
           totrms = sqrt ( totrms )
           funb(i1,j1) = funmax * 1.3d0
         else
           tlen = sqrt ( (xi2 - xi1)**2 + (yj2 - yj1)**2 )
           if (tlen .gt. almostzero) then
             cosa = (xi2 - xi1) / tlen
             sina = (yj2 - yj1) / tlen
           else
             cosa = 1.0d0
             sina = 0.0d0
           endif
           delta = sqrt ( dx * dy ) / max ( abs ( cosa ), abs ( sina ) )
!!           write (*,*) tlen, npts, cosa, sina, dx, delta
           
           pt1 = 0.01d0 * dx
           do k=1,npts
             ik = i1 + nint ( cosa * (delta / dx) * dble ( k - 1 ) )
             jk = j1 + nint ( sina * (delta / dy) * dble ( k - 1 ) )
             arg(k) = max ( delta * dble ( k - 1 ), pt1 )
             prof(k) = funa(ik,jk)
             if (abs ( prof(k) ) .lt. almostzero) prof(k) = 1.0d-33
             if (abs ( funmin ) .lt. almostzero .and. abs ( funmax ) .lt. almostzero) then
               funb(ik,jk) = 1.0d0
             else
               funb(ik,jk) = funmax * 1.3d0
             endif
           enddo
           kz = npts / 2 + 1
           do k=1,npts
             argz(k) = (k - kz) * delta
           enddo
         endif

! Compute logarithmic slopes of the profiles.

         slope(1) = 1.0d-31
         if (npts .gt. 1) then
           do m=2,npts-1
             dlograd = log10 ( arg(m) ) - log10 ( arg(m-1) )
             if (prof(m-1) .gt. almostzero .and. prof(m) .gt. almostzero) then
!!               slope(m) = abs ( (log10 ( prof(m-1) ) - log10 ( prof(m) )) / dlograd )
               slope(m) = (log10 ( prof(m) ) - log10 ( prof(m-1) )) / dlograd
             else
               slope(m) = 1.0d-31
             endif
             if (abs ( slope(m) ) .lt. almostzero) slope(m) = 1.0d-33
           enddo
           slope(1) = slope(2)
           slope(npts) = slope(npts-1)

! Smooth the slopes by median filtering.

           nwin = 3

           do m=1,npts
             nm = 0
             do j=max(m-nwin,1),m
               nm = nm + 1
               slo(nm) = slope(j)
             enddo
             if (nm .ge. 1) then
               slom(m) = selectk ( (nm + 1) / 2, nm, slo )
             else
               slom(m) = slope(m)
             endif
           enddo
           do m=1,npts
             slope(m) = slom(m)
           enddo
         endif

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         datamin = funmin
         datamax = funmax

         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.pro'//cvalue(1:icva-1)//'.fits'
           ion = lastc ( outname )
         endif
         profname = outname(1:ion-5)//'.dat'
         ipn = lastc ( profname )
         if (cverbose .eq. '-verb2') then
           write (iotty,'(a)') '   Saving profile in '''//profname(1:ipn)//''''
           write (iotty,'(a,1x,2(1pe14.7))') '   Minmax values in the image:', funmin, funmax
         endif
         open ( 15, file=profname(1:ipn), status='unknown' )
         write (15,'(a)') '# MODFITS = Modify FITS Image & Header = '//compda
         write (15,'(a)') '# Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
         write (15,'(a)') '# Using '//clibname(lb:7)//' library version'//cfitsversion//' by William D Pence.'
         write (15,'(a)') '#'
         if (lazimuthaverage) then
         write (15,'(a)') '# AZIMUTHALLY-AVERAGED RADIAL PROFILE EXTRACTED FROM IMAGE: '//filename(1:fnlen1)
         write (15,'(a)') '# ABOUT CENTRAL POSITION (X1,Y1): '//arg2(1:ia2)//' (PIXELS)'
         write (15,'(a,1pe14.7)') '# TOTAL VALUE: ', totint
         write (15,'(a,1pe14.7)') '# TOTAL STDEV: ', totrms
         else
         write (15,'(a)') '# RADIAL PROFILE EXTRACTED FROM IMAGE: '//filename(1:fnlen1)
         write (15,'(a)') '# ALONG STRAIGHT LINE (X1,Y1,X2,Y2): '//arg2(1:ia2)//' (PIXELS)'
         endif
         write (15,'(a)') '#'
         write (15,'(a)') '# TABULATED QUANTITIES:'
         write (15,'(a)') '#'
         write (15,'(a)') '#   N ........ Radial coordinate in pixels'
         write (15,'(a)') '#   ARCSEC ... Radial coordinate in acseconds'
         write (15,'(a)') '#   ARCSECM .. Radial coordinate in acseconds (zero in the middle)'
         if (lazimuthaverage) then
         write (15,'(a)') '#   VALUE .... Azimuthally-averaged image value at radial coordinate'
         write (15,'(a)') '#   STDEV .... Azimuthal standard deviations of image values'
         else
         write (15,'(a)') '#   VALUE .... Image pixel value at radial coordinate'
         endif
         write (15,'(a)') '#   SLOPE .... Logarithmic slope of the radial profile'
         write (15,'(a)') '#'
         if (lazimuthaverage) then
           write (15,'(a)') '#   N      ARCSEC        ARCSECM         VALUE          STDEV          SLOPE'
         else
           write (15,'(a)') '#   N      ARCSEC        ARCSECM         VALUE          SLOPE'
         endif
         write (15,'(a)') '#'
         do k=1,npts-1
           arg(k) = min ( max ( arg(k), -1.0d37 ), 1.0d37 )
           argz(k) = min ( max ( argz(k), -1.0d37 ), 1.0d37 )
           prof(k) = min ( max ( prof(k), -1.0d37 ), 1.0d37 )
           slope(k) = min ( max ( slope(k), -1.0d37 ), 1.0d37 )
           if (lazimuthaverage) then
             profrms(k) = min ( max ( profrms(k), -1.0d37 ), 1.0d37 )
             write (15,'(i5,1x,5(1x,1pe14.7))') k, arg(k), argz(k), prof(k), profrms(k), -slope(k)
           else
             write (15,'(i5,1x,4(1x,1pe14.7))') k, arg(k), argz(k), prof(k), -slope(k)
           endif
         enddo

         deallocate ( arg, argz, prof, profrms, slope, slo, slom )
         close ( 15 )

         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'spread' .or. maction(1:iact) .eq. 'spr') then

         i1 = nint ( value1 )
         j1 = nint ( value2 )
         i2 = nint ( value3 )
         j2 = nint ( value4 )

         if (i1 .eq. i2) then
           if (cverbose .eq. '-verb2')
     &         write (iotty,'(a)') '   Spreading Y-line of pixels along X-axis: '//arg2(1:ia2)

           if (j1 .lt. 1 .or. j2 .gt. ny) then
             write (iotty,'(/a)') '   MODFITS: ERROR: Line coordinates outside of image dimensions.'
             stop 99
           endif
     
           do j=j1,j2
             funb(1,j) = funa(i1,j)
           enddo
           do i=2,nx
             do j=j1,j2
               funb(i,j) = funb(1,j)
             enddo
           enddo

           if (lnoname) then
             outname = filename(isp1:fnlen1-5)//'.yx.fits'
           endif
         endif

         if (j1 .eq. j2) then
           if (cverbose .eq. '-verb2')
     &         write (iotty,'(a)') '   Spreading X-line of pixels along Y-axis: '//arg2(1:ia2)

           if (i1 .lt. 1 .or. i2 .gt. nx) then
             write (iotty,'(/a)') '   MODFITS: ERROR: Line coordinates outside of image dimensions.'
             stop 99
           endif

           do i=i1,i2
             funb(i,1) = funa(i,j1)
           enddo
           do j=2,ny
             do i=i1,i2
               funb(i,j) = funb(i,1)
             enddo
           enddo

           if (lnoname) then
             outname = filename(isp1:fnlen1-5)//'.xy.fits'
           endif
         endif

         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'fixpix' .or. maction(1:iact) .eq. 'fix') then

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funa(i,j) .gt. funmax) funmax = funa(i,j)
             if (funa(i,j) .lt. funmin) funmin = funa(i,j)
             funb(i,j) = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'(''     Old minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         beam = value3
         if (beam .gt. almostzero) then
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Replacing area of pixels using a Gaussian: '//arg2(1:ia2)
         else
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Replacing area of pixels within a circle: '//arg2(1:ia2)
         endif
         i0 = nint ( value1 )
         j0 = nint ( value2 )
         if (i0 .lt. 1 .or. i0 .gt. nx .or. j0 .lt. 1 .or. j0 .gt. ny) then
           write (iotty,'(/a,2i5)') '   MODFITS: ERROR: Pixel outside of image: ', i0, j0
           stop 99
         endif
         beam = value3
         peak = value4

         hwhm = abs ( beam ) / 2.0d0
         sigm2g = hwhm**2 / log ( 4.0d0 )
         if (beam .gt. almostzero) then
           radmax2 = (1.3d0 * hwhm)**2
         else
           radmax2 = beam**2
         endif
         nw = nint ( radmax2 / (dx * dy) )
         
         do j=max(j0-nw,1),min(j0+nw,ny)
           dely = (dble ( j - j0 ) * dy)**2
           do i=max(i0-nw,1),min(i0+nw,nx)
             delx = (dble ( i - i0 ) * dx)**2
             rad2 = delx + dely
             if (peak .gt. almostzero) then
               if (rad2 .le. radmax2) then
                 argex = min ( rad2 / (2.0d0 * sigm2g), 50.0d0 )
                 funb(i,j) = max ( funb(i,j), peak * exp ( -argex ) )
               endif
             else
               if (rad2 .le. beam**2) then
                 funb(i,j) = abs ( peak )
               endif
             endif
           enddo
         enddo

         funmin =  1.0d+30
         funmax = -1.0d+30
         do j=1,ny
           do i=1,nx
             if (funb(i,j) .gt. funmax) funmax = funb(i,j)
             if (funb(i,j) .lt. funmin) funmin = funb(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') write (iotty,'( ''     New maximum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmax
         if (cverbose .eq. '-verb2') write (iotty,'( ''     New minimum: '',$)')
         if (cverbose .eq. '-verb2') write (iotty,*) funmin
         datamin = funmin
         datamax = funmax
         if (lnoname) then
           outname = filename(isp1:fnlen1-5)//'.fix.fits'
         endif
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'expand') then
        
         inrem = int ( log10 ( dble ( max ( nint ( value ), 1 ) ) ) ) + 1
         write (cnrem,'(i7)') nint ( value )
         if (cverbose .eq. '-verb2') then
           write (iotty,'(a)') '   Expanding image from its edges by '//cnrem(7-inrem+1:7)//' pixels'
           write (cx,'(i7)') nx
           write (cy,'(i7)') ny
           ncx = int ( log10 ( dble ( nx ) ) ) + 1
           ncy = int ( log10 ( dble ( ny ) ) ) + 1
           write (iotty,'(a)') '   Old image dimensions: '//cx(7-ncx+1:7)//' '//cy(7-ncy+1:7)//' pixels'
         endif

         nbw = nint ( value )
         if (nbw .gt. 0) then
           nx1 = nx + 2 * nbw
           ny1 = ny + 2 * nbw
           
           call expandit ( 'outflows', nx, ny, nx1, ny1, nbw, inmx, inmx, inmx, inmx, funb, funa )
       
           nx = nx1
           ny = ny1
           crpix1 = crpix1 + dble ( nbw )
           crpix2 = crpix2 + dble ( nbw )
         endif
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') then
           write (cx,'(i7)') nx
           write (cy,'(i7)') ny
           ncx = int ( log10 ( dble ( nx ) ) ) + 1
           ncy = int ( log10 ( dble ( ny ) ) ) + 1
           write (iotty,'(a)') '   New image dimensions: '//cx(7-ncx+1:7)//' '//cy(7-ncy+1:7)//' pixels'
           write (iotty,'(a,1x,2(1pe14.7))') '   Minmax values in the image:', funmin, funmax
         endif

         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.exp'//cvalue(1:icva-1)//'.fits'
         endif
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (maction(1:iact) .eq. 'shrink') then
        
         inrem = int ( log10 ( dble ( max ( nint(value), 1 ) ) ) ) + 1
         write (cnrem,'(i7)') nint(value)
         if (cverbose .eq. '-verb2') then
           write (iotty,'(a)') '   Shrinking image by '//cnrem(7-inrem+1:7)//' pixels'
           write (cx,'(i7)') nx
           write (cy,'(i7)') ny
           ncx = int ( log10 ( dble ( nx ) ) ) + 1
           ncy = int ( log10 ( dble ( ny ) ) ) + 1
           write (iotty,'(a)') '   Old image dimensions: '//cx(7-ncx+1:7)//' '//cy(7-ncy+1:7)//' pixels'
         endif

         nbw = nint ( value )
         if (nbw .gt. 0) then
           nx1 = nx - 2 * nbw
           ny1 = ny - 2 * nbw
           
           call shrinkit ( nx, ny, nx1, ny1, nbw, inmx, inmx, funb, funa )
       
           nx = nx1
           ny = ny1
           crpix1 = crpix1 - dble ( nbw )
           crpix2 = crpix2 - dble ( nbw )
         endif
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j)
           enddo
         enddo
         if (cverbose .eq. '-verb2') then
           write (cx,'(i7)') nx
           write (cy,'(i7)') ny
           ncx = int ( log10 ( dble ( nx ) ) ) + 1
           ncy = int ( log10 ( dble ( ny ) ) ) + 1
           write (iotty,'(a)') '   New image dimensions: '//cx(7-ncx+1:7)//' '//cy(7-ncy+1:7)//' pixels'
           write (iotty,'(a,1x,2(1pe14.7))') '   Minmax values in the image:', funmin, funmax
         endif

         if (lnoname) then
           ic1 = index ( cvalue(1:icva-1), '.' )
           if (ic1 .gt. 0) cvalue(ic1:ic1) = 'p'
           outname = filename(isp1:fnlen1-5)//'.shr'//cvalue(1:icva-1)//'.fits'
         endif
         changed = .true.
       endif

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       if (ltiming) then
         if (maction(1:5) .eq. 'clip' .or. maction(1:3) .eq. 'cli') then
           cput_modfits = cput_modfits + timer ( 'cpu', cpu_modfits )
           walt_modfits = walt_modfits + timer ( 'wal', wal_modfits )
           write (60,'(a)') ' modfits '//arg1(1:ia1)//' '//arg2(1:ia2)//' '//arg3(1:ia3)//' '//arg4(1:ia4)
     &                            //' '//arg5(1:ia5)//' '  //arg6(1:ia6)//' '//arg7(1:ia7)

! Total cpu usage, expressed in hours, minutes, and seconds.

           cputot = cput_modfits
           cpuhours = int ( cputot / 3600.0d0 )
           cpumins  = int ( cputot / 60.0d0 ) -  60 * cpuhours
           cpusecs  = int ( cputot + 0.5d0 ) - 3600 * cpuhours - 60 * cpumins

! Final wall clock time, expressed in hours, minutes, and seconds.

           wctot = walt_modfits
           wchours = int ( wctot / 3600.0d0 )
           wcmins  = int ( wctot / 60.0d0 ) -  60 * wchours
           wcsecs  = int ( wctot + 0.5d0 ) - 3600 * wchours - 60 * wcmins
          
           write (60,'(a)') ' ------------------------------------------------------'
           write (60,'(a,2(1pe9.2,a))') ' Execution CPU & WALL times:', cputot, ' sec ', wctot, ' sec'
           write (60,'(a,i3.2,a,i2.2,a,i2.2,a,i3.2,a,i2.2,a,i2.2,a,f7.2,a)') ' ---------------------------'
     &              , cpuhours,':',cpumins,':',cpusecs, ' ----', wchours,':',wcmins,':',wcsecs, ' ---'
          
           call when ( lunix, ctime, cdate, ndate, 4 )
          
           write (60,'(a)') ' MODFITS finished: '//cdate//' '//ctime
           close ( 60 )
         endif
       endif

       if (changed) then
         ion = lastc ( outname )
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Writing output image '''//outname(1:ion)//''''

         call wfits ( cfitsversion, inmx, nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, funb, dx, dy, object
     &              , crval1, crval2, outname(1:ion), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history, maction(1:iact) )
                                                                    
         deallocate ( funa, funb, level )
         
         if (cverbose .eq. '-verb2') write (iotty,'(/a)') ' Done.'
       else
         if (cverbose .eq. '-verb2') write (iotty,'(/a)') ' No changes, nothing to do.'
       endif

!!       stop !<-- commented out because it would lead to run-time messages about denormalized values, when using gfortran.
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine sizemeasure 
       
     &            ( iotty, nextr, nxmin, nxmax, nymin, nymax, inmx, nx, ny, image, segmn, dx, dy, momxco, momyco, afwhm, bfwhm
     &            , atheta, afoot, bfoot, equivrad, elongation, nonzero, almostzero )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       integer       n, nx, ny, l, m, nextr, irc, iotty, inmx
       integer       nxmin(*), nxmax(*), nymin(*), nymax(*)
 
       real*8        dx, dy, pi, sqpixel, sourcearea, dblm, dbll, convtodegrees, sqreightlogtwo, pixel, intsxx, intsyy, intsxy
     &             , xcorsc, ycorsc, momp2, momm2, sqrmo, asig, bsig, intsd, ints0m, ints0, almostzero, intx, inty, int0, int0m
     &             , intd, asize, bsize
     
       real*8        image(inmx,inmx), segmn(nx,ny), momxco(*), momyco(*), afwhm(*), bfwhm(*), atheta(*), equivrad(*), elongation(*)
     &             , nonzero(*), afoot(*), bfoot(*)

       real*8 , allocatable :: momxx(:), momyy(:), momxy(:)

       parameter   ( pi = 3.14159265358979d0, convtodegrees = 180.0d0 / pi )
!__________________________________________________________________________________________________________________________________
!                  
       sqreightlogtwo = sqrt ( 8.0d0 * log ( 2.0d0 ) )
       sqpixel = dx * dy
       pixel = sqrt ( sqpixel )
       asig = 0.0d0
       bsig = 0.0d0

       allocate ( momxx(nextr), momyy(nextr), momxy(nextr), stat=irc )

       if (irc .ne. 0) then
         write (iotty,'(/a)') '   MEASURE: ERROR: Trouble allocating memory (76).'
         stop 76
       endif

       do n=1,nextr
         intx = 0.0d0                        
         inty = 0.0d0                        
         int0 = 0.0d0                        
         if (nonzero(n) .gt. 9.0d0) then
           do m=nymin(n),nymax(n)              
             dblm = dble ( m )                 
             do l=nxmin(n),nxmax(n)            
               dbll = dble ( l )               
               if (n .eq. nint ( segmn(l,m) )) then
                 if (image(l,m) .gt. almostzero) then
                   intd = image(l,m) * sqpixel
                   int0 = int0 + intd
                   intx = intx + dbll * intd
                   inty = inty + dblm * intd
                 endif
               endif
             enddo
           enddo
           int0m = int0 + almostzero
           momxco(n) = intx / int0m
           momyco(n) = inty / int0m
         else
           momxco(n) = (nxmin(n) + nxmax(n)) / 2.0d0
           momyco(n) = (nymin(n) + nymax(n)) / 2.0d0
         endif
       enddo
       
       do n=1,nextr
         ints0 = 0.0d0
         intsxx = 0.0d0
         intsyy = 0.0d0
         intsxy = 0.0d0           
         if (nonzero(n) .gt. 9.0d0) then
           do m=nymin(n),nymax(n)
             dblm = dble ( m )
             do l=nxmin(n),nxmax(n)
               dbll = dble ( l )       
               if (n .eq. nint ( segmn(l,m) )) then
                 intsd = segmn(l,m) * sqpixel
                 ints0 = ints0 + intsd
                 xcorsc = dbll - momxco(n)
                 ycorsc = dblm - momyco(n)
                 intsxx = intsxx + xcorsc * xcorsc * intsd
                 intsyy = intsyy + ycorsc * ycorsc * intsd
                 intsxy = intsxy + xcorsc * ycorsc * intsd
               endif
             enddo
           enddo
         endif
         ints0m = ints0 + almostzero
         momxx(n) = intsxx / ints0m
         momyy(n) = intsyy / ints0m
         momxy(n) = intsxy / ints0m
       enddo

       do n=1,nextr
         if (nonzero(n) .gt. 9.0d0) then
           momp2 = (momxx(n) + momyy(n)) / 2.0d0
           momm2 = (momxx(n) - momyy(n)) / 2.0d0
           
           if (abs ( momxy(n) ) .lt. almostzero) then
             if (momxx(n) .gt. momyy(n)) then
               momp2 = momxx(n)
               momm2 = momyy(n)
               atheta(n) = 90.0d0
             else
               momp2 = momyy(n)
               momm2 = momxx(n)
               atheta(n) = 0.0d0
             endif
             sqrmo = sqrt ( momm2**2 + momxy(n)**2 )
             asig  = sqrt ( momp2 + sqrmo ) * pixel
             if (momp2 .gt. sqrmo) then 
               bsig = sqrt ( momp2 - sqrmo ) * pixel
             else
               asig = 0.0d0
               bsig = 0.0d0
             endif
           else       
             sqrmo = sqrt ( momm2**2 + momxy(n)**2 )
             if (momp2 + sqrmo .ge. 0.0d0) then 
               asig  = sqrt ( momp2 + sqrmo ) * pixel
             endif
             if (momp2 - sqrmo .ge. 0.0d0) then 
               bsig  = sqrt ( momp2 - sqrmo ) * pixel
             endif
             if (momp2 + sqrmo .lt. 0.0d0 .or. momp2 - sqrmo .lt. 0.0d0) then 
               asig = 0.0d0
               bsig = 0.0d0
             endif
             atheta(n) = 180.0d0 - atan2 ( sqrmo + momm2, momxy(n) ) * convtodegrees
           endif                                        
           
           asize = sqreightlogtwo * max ( asig, pixel )
           bsize = sqreightlogtwo * max ( bsig, pixel )
           
           afwhm(n) = asize
           bfwhm(n) = bsize
           elongation(n) = asize / bsize
  
           sourcearea = nonzero(n) * sqpixel
           equivrad(n) = max ( sqrt ( sourcearea / pi ), pixel )
  
           afoot(n) = 2.0d0 * afwhm(n)
           bfoot(n) = 2.0d0 * bfwhm(n)
         else
           asize = sqreightlogtwo * pixel
           bsize = asize
           elongation(n) = 1.0d0
           afwhm(n) = asize
           bfwhm(n) = bsize
           equivrad(n) = asize
           afoot(n) = 2.0d0 * afwhm(n)
           bfoot(n) = afwhm(n)
         endif
       enddo
       
       deallocate ( momxx, momyy, momxy )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine rfits 
       
     &            ( inmx, nx, ny, bunit, ctype1, ctype2, rp1, rp2, xr, yr, fun, dx, dy, object, ra, dec, fname, iotty, iolog
     &            , creator, beam, funmin, funmax, blank, rot1, rot2, cd11, cd12, cd21, cd22, equinox, bzero, bscale, wave
     &            , datamn, datamx, history, cverbose, rcore )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       simple, extend, anynull, myfile
       character*(*) ctype1, ctype2, creator, fname, bunit, object, history, cverbose
       character*3   cnan
       character*25  cdatetime
       character*80  comment
       integer       inmx, i, j, nx, ny, stat, blksize, group, naxis, pcount, gcount, bitpix, unit, iotty, iolog, maxdim
     &             , rw, nfound, blank, stat05, stat06, stat07, stat08, stat09, stat10, stat11, stat12, stat13, stat14, stat15
     &             , stat16, stat17, stat18, stat19, stat20, stat21, stat22, stat23, stat24, stat25, stat26, stat27, stat28, stat29
     &             , stat30, stat31, stat32, stat0
       parameter   ( maxdim = 2 )
       integer       naxes(maxdim)
       real*8        dx, dy, ra, dec, beam, funmin, funmax, xr, yr, rp1, rp2, rot1, rot2, epoch, equinox, bzero, bscale, wave
     &             , dxdeg, dydeg, rcore, cd11, cd12, cd21, cd22, almostzero, as2deg, datamx, datamn, sign1, arg1, arg2, argm
     &             , sinrota, cosrota, pi, angle
       real*8        fun(inmx,*)
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
         creator = 'MODFITS'
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

       call ftgkyd ( unit, 'RCORE' , rcore , comment, stat )

       if (stat .ne. 0) then
         rcore = 0.0d0
         stat = 0
       endif

       stat0 = stat05 + stat06 + stat07 + stat09 + stat10 + stat11 + stat12 + stat13 + stat14 + stat15 + stat16 + stat17 + stat18
     &       + stat19 + stat20 + stat21 + stat22 + stat23 + stat24 + stat25 + stat26 + stat27 + stat28 + stat30 + stat31

       if (cverbose .eq. '-verb2') then
!!         if (stat0  .gt. 0) write (iotty,'()')
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

       call ftg2dd ( unit, group, dble ( blank ), inmx, naxes(1), naxes(2), fun, anynull, stat )

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
       
     &            ( cfitsversion, inmx, nx, ny, bunit, ctype1, ctype2, rp1, rp2, xr, yr, fun, dx, dy, object, ra, dec, fname, cdate
     &            , ctime, creator, beam, blank, rot1, rot2, cd11, cd12, cd21, cd22, equinox, bzero, bscale, wave, datamn, datamx
     &            , history, maction )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       simple, extend
       character*(*) fname, ctype1, ctype2, object, creator, ctime, cdate, bunit, cfitsversion, history, maction
       character*19  cdatetime
       integer       inmx, i, j, nx, ny, status, blocksze, group, naxis, pcount, gcount, bitpix, unit, naxes(2), blank, ica
     &             , itilde
       real*8        xr, yr, dx, dy, rp1, rp2, ra, dec, beam, bzero, bscale, as2deg, rot1, rot2, equinox, datamn
     &             , datamx, dxdeg, dydeg, wave, cd11, cd12, cd21, cd22, almostzero, pi, angle
       real*8        fun(inmx,*)
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
cc!      bzero  = 0.5d0 * (datamin + datamax)
cc!      bscale = dabs ( datamax - bzero ) / 32767.0d0

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

       call ftpkys ( unit, 'CREATOR', creator,'Alexander Menshchikov, DAp IRFU CEA Saclay', status )
       call ftpcom ( unit,' ',status)
       call ftpkys ( unit, 'DATE'  , cdatetime, 'creation date and time', status )

! Write all parameters.

       dxdeg = - dx * as2deg
       dydeg = dy * as2deg

! CDELT* are defined correctly: determine CD1* and CD2*.

       angle = rot2 * pi / 180.0d0
       cd11 = dxdeg * cos ( angle )
       cd12 = abs ( dydeg ) * sign ( 1.0d0, dxdeg ) * sin ( angle )
       cd21 = - abs ( dxdeg ) * sign ( 1.0d0, dydeg ) * sin ( angle )
       cd22 = dydeg * cos ( angle )
       if (abs ( cd12 ) .lt. almostzero) cd12 = 0.0d0
       if (abs ( cd21 ) .lt. almostzero) cd21 = 0.0d0

       ica = 1
       itilde = index ( maction, '~' )
       if (itilde .gt. 0) then
         bunit = maction(itilde+1:)
       endif

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

       call ftp2dd ( unit, group, inmx, nx, ny, fun, status )

! Close the file.

       call ftclos ( unit, status )

! Free the unit number.

       call ftfiou ( unit, status )

! Check for any error, and if so print out error messages

       if (status .gt. 0) call printerr ( status )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||