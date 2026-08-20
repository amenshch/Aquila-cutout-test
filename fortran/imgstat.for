
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       program imgstat
!__________________________________________________________________________________________________________________________________
!
! Compute various statistical quantities for a FITS image.
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       logical       lfnam1, lfnam2, lunix, linsave, lsaveinmask, lstd, lmedimage, lmeaimage, lhilog, llolog
     &             , lmedian, lpdf, lpdflin, lpdflog, lpdfpos, ltiming, lmedapprox, lringwin, lskeimage, lkurimage, lmodimage, lplus
     &             , lmaxmin, ldone, inellipse
                                                                                                                        
       integer       infu, iotty, iolog, nfu

       parameter   ( infu = 2000, iotty = 6, iolog = 0, nfu = 10 )
                                                                                                                        
       character*3   dot, cwavenum
       character*4   cpmx
       character*6   cfitsversion, cverbose                                                                           
       character*7   clibname, cx, cy, cdistpc
       character*8   ctime, cfpptime
       character*9   ccpu, cwal, cccpu, ccwal
       character*10  cdate
       character*21  compda
       character*80  object, creator, ctype1, ctype2, coption, bunit, history
       character*500 filename1, filename2, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, args, imagename, maskname, outname
     &             , outfits

       integer       firstb, lastc, i, j, k, l, ndate, ia1, ia2, ia3, ia4, ia5, ia6, ia7, ia8, ia9, blank, nx1, ny1, nx2, ny2, iof
     &             , nbins, nxy, fnlen1, fnlen2, ion, nto, m, n, npx, irc, ico, nedge, im1, jm1, ip1, jp1, lb, nm, kk, npp1, lenimg
     &             , lenmsk, npoints, nbinslo, nbinshi, ncx, ncy, imax2, jmax2, npmx, n1, n2, m1, m2, wchours, wcmins !,ir, nrmax
     &             , wcsecs, cpuhours, cpumins, cpusecs, ndc, ndw, nxco, nyco, nn, mm, iml, ipl, jml, jpl, icd, iar, npxstat

       real          fitsvers
       real*8        timer, selectk, medapprox, modeapprox, wave, ra, dec, beam, sigma1, funmin1, funmax1, funmin2, funmax2, fmean2
     &             , fmean, sigmas, variance, rnpx, totsum2, totsumedge, maskedgemean, crpix1, crpix2, crval1, crval2, sigma, funhi
     &             , sigma2, dx1, dy1, dx2, dy2, funlo, crota1, crota2, equinox, bzero, bscale, datamin, datamax, mom2, mom3
     &             , mom4, skewness2, kurtosis2, deltaf, fmedian2, cd11, cd12, cd21, cd22, pdfbin, pdfmin, sq_arcsecs_per_sterad
     &             , convftot, pdflo, rc2, xcentr, xau, funmin2nozero, almostzero, funmin, funmax, pc, amu, muH2, Msun, bin, distpc
     &             , convmass, fwhm, funhalf, ycentr, isumpct, numtotal, numpart, rw2, rx2, ry2, cpu, wal, cpu_sum, wal_sum, cputot
     &             , wctot, xpoints, stdmean, stdstdm, afwhm, bfwhm, atheta, equivrad, elongation, mxco, myco, fhmax, dbli, dblj
     &             , rxy2, afwhmx, bfwhmx, fmode2, pixvalmax, pixvalmin, ramx, maskedgemedian !!, foomin1, foomin2, foomax1, foomax2
     &             , xfwhmx, yfwhmx, xcentrx, ycentrx, fwhmxy, funhalfx, funhalff, xfwfoot, yfwfoot, fwhmxyf, xcentrf, ycentrf

       integer, allocatable :: cimask(:,:)
       real*8, allocatable :: fun1(:,:), fun2(:,:), funstd(:,:), funm(:,:), funpdf(:,:), funone(:), funint(:), plane(:), funx(:)
     &                      , funbin(:,:)  !!, rad(:), npts(:), funmean(:), fun1mean(:,:)

       parameter   ( sq_arcsecs_per_sterad = 3282.80635d0 * 3600.0d0**2, pc = 3.085678d18, amu = 1.6605402D-24, muH2 = 2.8d0
     &             , xau = 1.4959787e+13, Msun = 1.9891d33, almostzero = 1.0d-20, ltiming = .true., dot = '•' )

       external      firstb, lastc, when, osystem, ftvers, getfitshead, rfits, wfits, halfmaxsizes, momentsizes, inellipse
     &             , selectk, modeapprox, medapprox, timer, showprogress, xyhalfmaxsizes
!__________________________________________________________________________________________________________________________________
!
! Determine operating system.

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
       call getarg ( 8, arg8 )
       call getarg ( 9, arg9 )

       ia1 = lastc ( arg1 )
       ia2 = lastc ( arg2 )
       ia3 = lastc ( arg3 )
       ia4 = lastc ( arg4 )
       ia5 = lastc ( arg5 )
       ia6 = lastc ( arg6 )
       ia7 = lastc ( arg7 )
       ia8 = lastc ( arg8 )
       ia9 = lastc ( arg9 )

       cverbose = '-verb2'
       cwavenum = ''
       lpdflin = .false.
       lpdfpos = .false.
       npoints = 0
       nbinslo = 0

       if (ia9 .gt. 0 .and. arg9(1:5) .eq. '-verb') then
         cverbose = arg9(1:6)
         if (ia9 .eq. 9) cwavenum = arg9(7:9)
       endif
       if (ia8 .gt. 0 .and. arg8(1:5) .eq. '-verb') then
         cverbose = arg8(1:6)
         if (ia8 .eq. 9) cwavenum = arg8(7:9)
       endif
       if (ia7 .gt. 0 .and. arg7(1:5) .eq. '-verb') then
         cverbose = arg7(1:6)
         if (ia7 .eq. 9) cwavenum = arg7(7:9)
       endif
       if (ia6 .gt. 0 .and. arg6(1:5) .eq. '-verb') then
         cverbose = arg6(1:6)
         if (ia6 .eq. 9) cwavenum = arg6(7:9)
       endif
       if (ia5 .gt. 0 .and. arg5(1:5) .eq. '-verb') then
         cverbose = arg5(1:6)
         if (ia5 .eq. 9) cwavenum = arg5(7:9)
       endif
       if (ia4 .gt. 0 .and. arg4(1:5) .eq. '-verb') then
         cverbose = arg4(1:6)
         if (ia4 .eq. 9) cwavenum = arg4(7:9)
       endif
       if (ia3 .gt. 0 .and. arg3(1:5) .eq. '-verb') then
         cverbose = arg3(1:6)
         if (ia3 .eq. 9) cwavenum = arg3(7:9)
       endif

       if (cverbose .eq. '-verb2') then
         write (iotty,'( )')
         if (arg1(1:ia1) .eq. ':') write (iotty,'(a$)') '  '
         write (iotty,'(a$)') ' IMGSTAT '//dot//' Compute Image Statistics '//dot//' '//compda
       endif
       if (arg1(1:ia1) .eq. ':') stop
       if (cverbose .eq. '-verb2') write (iotty,'()')
       if (cverbose .eq. '-verb2') write (iotty,'( a)') ' Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
       if (cverbose .eq. '-verb2') write (iotty,'( a)') ' Using '//clibname(lb:7)//' library version'//cfitsversion
     &                                                //' by William D Pence.'
       if (ia1 .eq. 0) then
         write (iotty,'( )')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USAGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')
         write (iotty,'(a)') ' imgstat [<action>] <imagein> [<maskimage>] [<distpc>] [-o <imageout>] [-verb{0|1|2}]' 
         write (iotty,'( )')   
         write (iotty,'(a)') ' This utility computes various statistical quantities in a FITS image;'
         write (iotty,'(a)') ' <imagein> and <imageout> are the names of the input and output images. If the'
         write (iotty,'(a)') ' optional mask is supplied in <maskimage>, then the quantities can be computed'
         write (iotty,'(a)') ' inside or outside the masked area. If the mask image name is not specified,'
         write (iotty,'(a)') ' all pixels of the input image are used.'
         write (iotty,'( )')   
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ <action> ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')   
         write (iotty,'(a)') ' save    [-nomed] ........... compute statistics inside a mask and save'
         write (iotty,'(a)') ' savein  [-nomed] ........... compute statistics inside a mask and save'
         write (iotty,'(a)') ' stdev <npts> [+]<N> ........ make N*stdev image using a window with R=<npts>'
         write (iotty,'(a)') ' median|medbin [+|-]<npts> .. median-filter image using a window with R=<npts>'
         write (iotty,'(a)') ' mode [+|-]<npts> ........... mode-filter image using a window with R=<npts>'
         write (iotty,'(a)') ' mean [+|-]<npts> ........... mean-filter image using a window with R=<npts>'
         write (iotty,'(a)') ' skewness <npts> ............ image of skewness using a window with R=<npts>'
         write (iotty,'(a)') ' kurtosis <npts> ............ image of kurtosis using a window with R=<npts>'
         write (iotty,'(a)') ' pdf [+|-]<pdfbin> .......... compute pixel distribution function (histogram)'
         write (iotty,'(a)') ' maxmin <max> <min> ......... compute statistics within the limits and save'
         write (iotty,'( )')   
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ NOTES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')   
         write (iotty,'(a)') ' Sliding filtering window is circular if <npts> > 0 and ring if <npts> < 0;'
         write (iotty,'(a)') ' bin size <pdfbin> for creating PDFs is logarithmic if > 0 or linear if < 0;'
         write (iotty,'(a)') ' use +<pdfbin> to create PDFs with only positive values. To compute mass from'
         write (iotty,'(a)') ' an image of surface densities, distance (pc) must be supplied in the optional'
         write (iotty,'(a)') ' <distpc> parameter.'
         write (iotty,'(a)')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         stop 99
       endif

       args = arg1(1:ia1)//' '//arg2(1:ia2)//' '//arg3(1:ia3)//' '//arg4(1:ia4)//' '//arg5(1:ia5)//' '//arg6(1:ia6)//' '
     &      //arg7(1:ia7)//' '//arg8(1:ia8)//' '//arg9(1:ia9)
       iar = lastc ( args )

       if (cverbose .eq. '-verb1')       !!cverbose .eq. '-verb0' .or. 
     &    write (iotty,'(a)') ' IMGSTAT: '//args(1:iar)

       fnlen1 = 0
       fnlen2 = 0
       lstd = .false.
       lmodimage = .false.
       lmedimage = .false.
       lmeaimage = .false.
       lskeimage = .false.
       lkurimage = .false.
       lpdf = .false.
       linsave = .false.
       lsaveinmask = .false.
       lmedian = .true.
       lmedapprox = .false.
       lringwin = .false.
       lplus = .false.
       lmaxmin = .false.

       inquire ( file=arg1(1:ia1), exist=lfnam1 )                      
       if (.not.lfnam1) then
         if (arg1(max(ia1-4,1):ia1) .ne. '.fits') then
           filename1 = arg1(1:ia1)//'.fits'
           fnlen1 = lastc ( filename1 )
           inquire ( file=filename1(1:fnlen1), exist=lfnam1 )
         endif
       endif

       if (arg1(1:ia1) .ne. 'save' .and. arg1(1:ia1) .ne. 'savein' .and. arg1(1:ia1) .ne. 'stdev' .and. 
     &     arg1(1:ia1) .ne. 'median' .and. arg1(1:ia1) .ne. 'medbin' .and. arg1(1:ia1) .ne. 'mode' .and.
     &     arg1(1:ia1) .ne. 'mean' .and. arg1(1:ia1) .ne. 'skewness' .and. arg1(1:ia1) .ne. 'kurtosis' .and.
     &     arg1(1:ia1) .ne. 'maxmin' .and. .not.lfnam1) then
         write (iotty,'(/a)') '   IMGSTAT: ERROR: Wrong action or file name: '''//arg1(1:ia1)//''''
         stop 99
       endif

       if (arg1(1:ia1) .eq. 'save') then
         linsave = .true.
       endif
       if (arg1(1:ia1) .eq. 'savein') then
         lsaveinmask = .true.
         linsave = .true.
       endif
       if (arg2(1:ia2) .eq. '-nomed') then
         lmedian = .false.     
       endif
       if (arg1(1:ia1) .eq. 'stdev') then
         lstd = .true.
       endif
       if (arg1(1:ia1) .eq. 'median') then
         lmedimage = .true.
       endif
       if (arg1(1:ia1) .eq. 'medbin') then
         lmedimage = .true.
         lmedapprox = .true.
       endif
       if (arg1(1:ia1) .eq. 'mode') then
         lmodimage = .true.
       endif
       if (arg1(1:ia1) .eq. 'mean') then
         lmeaimage = .true.
       endif
       if (arg1(1:ia1) .eq. 'skewness') then
         lskeimage = .true.
       endif
       if (arg1(1:ia1) .eq. 'kurtosis') then
         lkurimage = .true.
       endif
       if (arg1(1:ia1) .eq. 'pdf') then
         lpdf = .true.
       endif
       if (arg1(1:ia1) .eq. 'maxmin') then
         lmaxmin = .true.
         linsave = .true.
       endif
       coption = arg1(1:ia1)
       ico = lastc ( coption )
       if (lstd) then       
         ia3 = lastc ( arg3 )
         if (ia3 .eq. 0) then
           write (iotty,'(/a)') '   IMGSTAT: ERROR: Too few command line parameters.'
           stop 99
         endif
         coption = coption(1:ico)//' '//arg2(1:ia2)//' '//arg3(1:ia3)
         ico = lastc ( coption )
         read (arg2,*) xpoints
         npoints = nint ( xpoints )
         read (arg3,*,err=1122) sigmas
         if (arg3(1:1) .eq. '+') lplus = .true.
         goto 1133
 1122    write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble reading SIGMAS, verify it was specified after XPOINTS.'
         stop 99
 1133    arg1 = arg4
         arg2 = arg5
         arg3 = arg6
         arg4 = arg7
         arg5 = arg8
         arg6 = arg9
       elseif (lmodimage .or. lmedimage .or. lmeaimage .or. lskeimage .or. lkurimage) then       
         ia3 = lastc ( arg3 )
         if (ia3 .eq. 0) then
           write (iotty,'(/a)') '   IMGSTAT: ERROR: Too few command line parameters.'
           stop 99
         endif
         coption = coption(1:ico)//' '//arg2(1:ia2)
         ico = lastc ( coption )
         read (arg2,*) xpoints
         if (xpoints .lt. 0.0d0) then
           xpoints = abs ( xpoints )
           lringwin = .true.
         endif
         npoints = nint ( xpoints )
         arg1 = arg3
         arg2 = arg4
         arg3 = arg5
         arg4 = arg6
         arg5 = arg7
         arg6 = arg8
         arg7 = arg9
       elseif (lpdf) then       
         ia3 = lastc ( arg3 )
         if (ia3 .eq. 0) then
           write (iotty,'(/a)') '   IMGSTAT: ERROR: Too few command line parameters.'
           stop 99
         endif
         coption = coption(1:ico)//' '//arg2(1:ia2)
         ico = lastc ( coption )
         read (arg2,*) pdfbin
         lpdflin = .false.
         lpdflog = .false.
         lpdfpos = .false.
         if (arg2(1:1) .eq. '-') then
           lpdflin = .true.
         else
           lpdflog = .true.
           if (arg2(1:1) .eq. '+') lpdfpos = .true.
         endif
         pdfbin = abs ( pdfbin )
         arg1 = arg3
         arg2 = arg4
         arg3 = arg5
         arg4 = arg6
         arg5 = arg7
         arg6 = arg8
         arg7 = arg9
       elseif (lmaxmin) then
         read (arg2,*) pixvalmax
         read (arg3,*) pixvalmin
         arg1 = arg4
         arg2 = arg5
         arg3 = arg6
         arg4 = arg7
         arg5 = arg8
         arg6 = arg9
       else
         if (linsave .or. lpdf) then
           if (lmedian) then
             arg1 = arg2
             arg2 = arg3
             arg3 = arg4
             arg4 = arg5
             arg5 = arg6
             arg6 = arg7
             arg7 = arg8
             arg8 = arg9
           else
             arg1 = arg3
             arg2 = arg4
             arg3 = arg5
             arg4 = arg6
             arg5 = arg7
             arg6 = arg8
             arg7 = arg9
           endif
         endif
       endif
       ia1 = lastc ( arg1 )
       ia2 = lastc ( arg2 )

       filename1 = arg1(1:ia1)
       fnlen1 = lastc ( filename1 )

! Check if the input file exists and open it.

       inquire ( file=filename1(1:fnlen1), exist=lfnam1 )                      
       if (.not.lfnam1) then
         if (filename1(max(fnlen1-4,1):fnlen1) .ne. '.fits') then
           filename1 = filename1(1:fnlen1)//'.fits'
           fnlen1 = lastc ( filename1 )
           inquire ( file=filename1(1:fnlen1), exist=lfnam1 )
         endif
       endif
       imagename = filename1(1:fnlen1)
       lenimg = fnlen1

       lfnam2 = .false.
       if (ia2 .gt. 0 .and. arg2(1:5) .ne. '-verb') then
         filename2 = arg2(1:ia2)
         fnlen2 = lastc ( filename2 )
         inquire ( file=filename2(1:fnlen2), exist=lfnam2 )
         if (.not.lfnam2) then
           if (filename2(max(fnlen2-4,1):fnlen2) .ne. '.fits') then
             filename2 = filename2(1:fnlen2)//'.fits'
             fnlen2 = lastc ( filename2 )
             inquire ( file=filename2(1:fnlen2), exist=lfnam2 )
           endif
         endif
       endif

       ia3 = lastc ( arg3 )
       ia4 = lastc ( arg4 )
       ia5 = lastc ( arg5 )
       ia6 = lastc ( arg6 )

       if (lfnam2) then
         maskname = filename2(1:fnlen2)
         lenmsk = fnlen2
         if (arg3(1:ia3) .ne. '') then
           read (arg3(1:ia3),*,err=1110) distpc
           cdistpc = arg3(1:7)
           goto 1112
1110       continue
         endif
       else
         maskname = ' '
         lenmsk = 6
         if (arg2(1:ia2) .ne. '') then
           read (arg2(1:ia2),*,err=1111) distpc
           cdistpc = arg2(1:7)
           goto 1112
1111       continue
         endif
       endif
       distpc = 140.0d0
       write (cdistpc,'(i7)') int ( distpc )
       cdistpc = adjustl ( cdistpc )
1112   continue
       icd = lastc ( cdistpc )
       
       outname = ''
       ion = 0
       if (arg2(1:2) .eq. '-o') then
         outname = arg3
       endif
       if (arg3(1:2) .eq. '-o') then 
         outname = arg4
       endif
       if (arg4(1:2) .eq. '-o') then 
         outname = arg5
       endif
       if (arg5(1:2) .eq. '-o') then 
         outname = arg6
       endif

       if (outname .ne. '') then
         ion = lastc ( outname )
         if (outname(max(ion-4,1):ion) .ne. '.fits') then
           outname = outname(1:ion)//'.fits'
           ion = lastc ( outname )
         endif
       endif

       if (.not.lfnam1) then
         write (iotty,'(/a)') '   IMGSTAT: ERROR: FITS file '''//filename1(1:fnlen1)//''' not found.'
         stop 1
       endif

! Determine numbers of pixels in the FITS image.

       call getfitshead ( imagename(1:lenimg), nx1, ny1, dx1, dy1, bunit )
       
       write (cx,'(i7)') nx1
       write (cy,'(i7)') ny1
       ncx = int ( log10 ( dble ( nx1 ) ) ) + 1
       ncy = int ( log10 ( dble ( ny1 ) ) ) + 1

       if (cverbose .eq. '-verb2') write (iotty,'(/a)') '   Reading ('//cx(7-ncx+1:7)//' x '//cy(7-ncy+1:7)//') '''
     &                                                //imagename(1:lenimg)//''' '

       allocate ( fun1(nx1,ny1), stat=irc )
   
       if (irc .ne. 0) then
         write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (10).'
         stop 10
       endif

       call rfits ( nx1, ny1, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, fun1, dx1, dy1, object, ra, dec
     &            , imagename(1:lenimg), iotty, 0, creator, beam, funmin1, funmax1, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &            , equinox, bzero, bscale, wave, datamin, datamax, history, cverbose )

       if (cverbose .eq. '-verb2') write (iotty,'(a,1x,2(1pe14.7))'  ) '   Minmax values:', funmin1, funmax1
       if (cverbose .eq. '-verb2') write (iotty,'(a,1x,2(1pe14.7),a)') '     Pixel sizes:', dx1, dy1, ' arcsec'

! Second file (mask image) if exists.       

       if (lfnam2) then

! Determine numbers of pixels in the FITS image.

         call getfitshead ( maskname(1:lenmsk), nx2, ny2, dx2, dy2, bunit )

         write (cx,'(i7)') nx2
         write (cy,'(i7)') ny2
         ncx = int ( log10 ( dble ( nx2 ) ) ) + 1
         ncy = int ( log10 ( dble ( ny2 ) ) ) + 1
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Reading ('//cx(7-ncx+1:7)//' x '//cy(7-ncy+1:7)//') '''
     &                                                 //maskname(1:lenmsk)//''' '

         allocate ( fun2(nx2,ny2), stat=irc )
       
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (20).'
           stop 20
         endif

         call rfits ( nx2, ny2, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, fun2, dx2, dy2, object, ra, dec
     &              , maskname(1:lenmsk), iotty, 0, creator, beam, funmin2, funmax2, blank, crota1, crota2, cd11, cd12, cd21
     &              , cd22, equinox, bzero, bscale, wave, datamin, datamax, history, cverbose )
       
         if (cverbose .eq. '-verb2') write (iotty,'(a,1x,2(1pe14.7))'  ) '   Minmax values:', funmin2, funmax2
         if (cverbose .eq. '-verb2') write (iotty,'(a,1x,2(1pe14.7),a)') '     Pixel sizes:', dx2, dy2, ' arcsec'

         if (nx1 .ne. nx2 .or. ny1 .ne. ny2 .or. abs ( dx1 - dx2 ) .gt. 1.0d-3 * dx1 .or. abs ( dy1 - dy2 ) .gt. 1.0d-3 * dy1) then
           write (iotty,'(/a)') '   IMGSTAT: ERROR: Images have different parameters (look in headers):'
           write (iotty,*) nx1, nx2, ny1, ny2, abs ( dx1 - dx2 ), abs ( dy1 - dy2 )
           stop 99
         endif
       else

         allocate ( fun2(nx1,ny1), stat=irc )
       
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (30).'
           stop 30
         endif

         do j=1,ny1
           do i=1,nx1
             fun2(i,j) = 1.0d0
           enddo
         enddo
       endif

       if (lmaxmin) then
         do j=1,ny1
           do i=1,nx1
             if (fun1(i,j) .gt. pixvalmax .or. fun1(i,j) .lt. pixvalmin) then
               fun2(i,j) = 0.0d0
             endif
           enddo
         enddo
       endif

       convftot = 1.0e6 * dx1 * dy1 / sq_arcsecs_per_sterad  !<-- intensity in MJy/sr to flux in Jy
       convmass = dx1 * dy1 * (distpc * xau)**2 * amu * muH2 / Msun
       
! Compute statistics only if production of stdev or median images was not requested.

       if (.not.lstd .and. .not.lmodimage .and. .not.lmedimage .and. .not.lmeaimage .and. .not.lskeimage .and. .not.lkurimage) then
         npxstat = 0
         totsum2 = 0.0d0
         funmin2 = 1.0d+30
         funmax2 =-1.0d+30
         funmin2nozero = 1.0d+30
         fmean2 = 0.0d0
         fmedian2 = 0.0d0
         fmode2 = 0.0d0
         imax2 = 0
         jmax2 = 0
         do j=1,ny1
           do i=1,nx1
             if (fun2(i,j) .gt. almostzero) then
               npxstat = npxstat + 1
               totsum2 = totsum2 + fun1(i,j)
               funmin2 = min ( fun1(i,j), funmin2 )
               if (fun1(i,j) .ge. funmax2) then
                 imax2 = i
                 jmax2 = j
                 funmax2 = fun1(i,j)
               endif
               if (fun1(i,j) .gt. almostzero) funmin2nozero = min ( fun1(i,j), funmin2nozero )
             endif
           enddo
         enddo
         if (npxstat .eq. 0) then
           funmin2 = 0.0d0
           funmax2 = 0.0d0
           funmin2nozero = 0.0d0
         else
           fmean2 = totsum2 / dble ( npxstat )
         endif

         sigma2 = 0.0d0
         nxy = 0
         do j=1,ny1
           do i=1,nx1
             if (fun2(i,j) .gt. almostzero) then
               nxy = nxy + 1
               sigma2 = sigma2 + (fun1(i,j) - fmean2)**2
             endif
           enddo
         enddo
         if (nxy .gt. 1) sigma2 = sqrt ( sigma2 / dble ( nxy - 1 ) )

         nedge = 0
         totsumedge = 0.0d0
         maskedgemean = 0.0d0
         maskedgemedian = 0.0d0
!!         if (lsaveinmask) then
         if (lfnam2) then
           allocate ( funx(nx1*ny1), stat=irc )
           if (irc .ne. 0) then
             write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (130).'
             stop 130
           endif
           do j=1,ny1
             do i=1,nx1
               if (fun2(i,j) .gt. almostzero) then
                 im1 = max ( i - 1, 1 )
                 jm1 = max ( j - 1, 1 )
                 ip1 = min ( i + 1, nx1 )
                 jp1 = min ( j + 1, ny1 )
                 if (fun2(im1,j  ) .lt. almostzero .or. fun2(ip1,j  ) .lt. almostzero .or.
     &               fun2(i  ,jm1) .lt. almostzero .or. fun2(i  ,jp1) .lt. almostzero .or.
     &               fun2(im1,jm1) .lt. almostzero .or. fun2(im1,jp1) .lt. almostzero .or.
     &               fun2(ip1,jm1) .lt. almostzero .or. fun2(ip1,jp1) .lt. almostzero) then
                   nedge = nedge + 1
                   totsumedge = totsumedge + fun1(i,j)
                   funx(nedge) = fun1(i,j)
                 endif
               endif
             enddo
           enddo
           if (nedge .gt. 0) maskedgemean = totsumedge / dble ( nedge )
           if (nedge .ne. 0) then
             maskedgemedian = selectk ( (nedge + 1) / 2, nedge, funx )
           endif
           deallocate ( funx )
         endif
         
! Find the peak sizes in two independent ways.

         call halfmaxsizes ( nx1, ny1, fun2, fun1, dx1, dy1, almostzero, afwhmx, bfwhmx, fwhm, xcentr, ycentr )
         
         funhalf = fun1(nint(xcentr),nint(ycentr)) / 2.0d0
         
         call xyhalfmaxsizes ( nx1, ny1, fun2, fun1, dx1, dy1, almostzero, xfwhmx, yfwhmx, fwhmxy, xcentrx, ycentrx )
         
         funhalfx = fun1(nint(xcentrx),nint(ycentrx)) / 2.0d0

         call xyfootsizes ( nx1, ny1, fun2, fun1, dx1, dy1, almostzero, xfwfoot, yfwfoot, fwhmxyf, xcentrf, ycentrf )
         
         funhalff = fun1(nint(xcentrf),nint(ycentrf)) / 2.0d0

         call momentsizes ( 1, nx1, 1, ny1, nx1, ny1, fun2, fun1, dx1, dy1, mxco, myco, afwhm, bfwhm, atheta, afwhmx, bfwhmx
     &                    , equivrad, elongation, almostzero )
     
         mxco = min ( max ( mxco, 1.0d0 ), dble ( nx1 ) )
         myco = min ( max ( myco, 1.0d0 ), dble ( ny1 ) )
         nxco = nint ( mxco )
         nyco = nint ( myco )

         if (fun1(nxco,nyco) .gt. almostzero) fhmax = -1.0d99
         if (fun1(nxco,nyco) .lt.-almostzero) fhmax = 1.0d99
         do j=1,ny1
           dblj = dble ( j )
           do i=1,nx1
             dbli = dble ( i )
             if (fun2(i,j) .gt. almostzero) then
               if (.not.inellipse ( dbli, dblj, dx1, dy1, mxco, myco, afwhm, bfwhm, atheta )) then
                 if (fun1(nxco,nyco) .gt. almostzero) fhmax = max ( fhmax, fun1(i,j) )
                 if (fun1(nxco,nyco) .lt.-almostzero) fhmax = min ( fhmax, fun1(i,j) )
               endif
             endif
           enddo
         enddo

! Skip this fragment for now...

!!         if (fhmax .eq. -1.0d100) then
!!           nrmax = min ( nx1 - 1, ny1 - 1 ) / 2 + 1
!!          
!!           allocate ( rad(nrmax+1), npts(nrmax+1), funmean(nrmax+1), fun1mean(nx1,ny1), stat=irc )
!!          
!!           if (irc .ne. 0) then
!!             write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (78).'
!!             stop 78
!!           endif
!!          
!!           do ir=1,nrmax+1
!!             npts(ir) = 0.0d0
!!             funmean(ir) = 0.0d0
!!             rad(ir) = (dble ( ir - 1 ) + 0.5d0) * sqrt ( dx1 * dy1 )
!!           enddo
!!          
!!           do ir=1,nrmax
!!             foomin1 = 2.0d0 * rad(ir)
!!             foomin2 = 2.0d0 * rad(ir+1)
!!             foomax1 = foomin1
!!             foomax2 = foomin2
!!               
!!             do m=1,ny1
!!               dblj = dble ( m )
!!               do l=1,nx1
!!                 dbli = dble ( l )
!!               
!!                 if (inellipse ( dbli, dblj, dx1, dy1, mxco, myco, foomax2, foomin2, atheta )) then
!!                   if (.not.inellipse ( dbli, dblj, dx1, dy1, mxco, myco, foomax1, foomin1, atheta )) then
!!                     npts(ir) = npts(ir) + 1.0d0
!!                     funmean(ir) = funmean(ir) + fun1(l,m)
!!                   endif
!!                 endif
!!               enddo
!!             enddo
!!           enddo
!!          
!!           do ir=1,nrmax
!!             if (npts(ir) .gt. 0.0d0) funmean(ir) = funmean(ir) / npts(ir)
!!           enddo
!!          
!!           do ir=1,nrmax
!!             foomin1 = 2.0d0 * rad(ir)
!!             foomin2 = 2.0d0 * rad(ir+1)
!!             foomax1 = foomin1
!!             foomax2 = foomin2
!!           
!!             do m=1,ny1
!!               dblj = dble ( m )
!!               do l=1,nx1
!!                 dbli = dble ( l )
!!                 
!!                 if (inellipse ( dbli, dblj, dx1, dy1, mxco, myco, foomax2, foomin2, atheta )) then
!!                   if (.not.inellipse ( dbli, dblj, dx1, dy1, mxco, myco, foomax1, foomin1, atheta )) then
!!                     fun1mean(l,m) = fun1mean(l,m) + funmean(ir)
!!                   endif
!!                 endif
!!               enddo
!!             enddo
!!           enddo
!!          
!!           outfits = imagename(1:lenimg-5)//'.avg.fits'
!!           iof = lastc ( outfits )
!!          
!!           if (cverbose .eq. '-verb2') 
!!     &       write (iotty,'(a,f4.1,a)') '   Writing averaged image '''//outfits(1:iof)//''''
!!          
!!             call wfits ( cfitsversion, nx1, ny1, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, fun1mean, dx1, dy1
!!     &                  , object, crval1, crval2, outfits(1:iof), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12
!!     &                  , cd21, cd22, equinox, bzero, bscale, wave, datamin, datamax, history )
!!          
!!           deallocate ( rad, npts, funmean, fun1mean )
!!         endif
       
! Calculation of the median value; do it only if requested explicitly, to save time.

         if (lmedian) then
           allocate ( funx(nx1*ny1), stat=irc )
           if (irc .ne. 0) then
             write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (131).'
             stop 131
           endif
           npx = 0
           do j=1,ny1
             do i=1,nx1
               if (fun2(i,j) .gt. almostzero) then
                 npx = npx + 1
                 funx(npx) = fun1(i,j)
               endif
             enddo
           enddo
           if (npx .ne. 0) then
             fmedian2 = selectk ( (npx + 1) / 2, npx, funx )
             fmode2 = modeapprox ( npx, funx, funmin2, funmax2 )
!!             fmedian2 = medapprox ( npx, funx, funmin2, funmax2 )
           endif
           deallocate ( funx )

           if (lpdf) then
             lhilog = .false.
             llolog = .false.
             pdflo = funmin2
             if (funmin2 .lt. -1.5d0 * pdfbin) then
               pdflo = -0.5d0 * pdfbin
             elseif (funmin2 .lt. 0.5d0 * pdfbin) then
               pdflo = funmin2nozero
             endif
             
             if (lpdflin) then
               if (pdflo .lt. 0.0d0) then
                 nbinslo = max ( nint ( (abs ( funmin2 ) - abs ( pdflo ) ) / pdfbin ) + 1, 1 )
                 nbinshi = max ( nint ( (abs ( funmax2 ) - abs ( pdflo ) ) / pdfbin ) + 1, 1 )
                 nbins = nbinslo + nbinshi
               else
                 nbins = max ( nint ( (funmax2 - funmin2) / pdfbin ) + 1, 2 )
               endif
             else
               if (pdflo .lt. 0.0d0) then
                 nbinslo = max ( nint ( (log10 ( abs ( funmin2 ) ) - log10 ( abs ( pdflo ) )) / pdfbin ) + 1, 1 )
                 nbinshi = max ( nint ( (log10 ( funmax2 ) - log10 ( abs ( pdflo ) )) / pdfbin ) + 1, 1 )
                 nbins = nbinslo + nbinshi
               else
                 nbins = max ( nint ( (log10 ( funmax2 ) - log10 ( pdflo ) ) / pdfbin ) + 1, 2 )
               endif
             endif

             allocate ( funone(nbins+1), funint(nbins+1), plane(nbins+1), stat=irc )
       
             if (irc .ne. 0) then
               write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (10).'
               stop 11
             endif

             do k=1,nbins
               funint(k) = 0.0d0
               funone(k) = 0.0d0
             enddo

             if (lpdflin) then
               do k=nbinslo,1,-1
                 plane(k) = -(abs ( pdflo ) + dble ( nbinslo - k ) * pdfbin)
               enddo
               do k=nbinslo+1,nbins
                 plane(k) = abs ( pdflo ) + dble ( k - nbinslo - 1 ) * pdfbin
               enddo
             else
               if (pdflo .lt. 0.0d0) then
                 do k=nbinslo,1,-1
                   plane(k) = -10.0d0**( log10 ( abs ( pdflo ) ) + dble ( nbinslo - k ) * pdfbin )
                 enddo
                 do k=nbinslo+1,nbins
                   plane(k) = 10.0d0**( log10 ( abs ( pdflo ) ) + dble ( k - nbinslo - 1 ) * pdfbin )
                 enddo
               else
                 do k=1,nbins
                   plane(k) = 10.0d0**( log10 ( pdflo ) + dble ( k - 1 ) * pdfbin )
                 enddo
               endif
             endif

             allocate ( funpdf(nx1,ny1), stat=irc )
             if (irc .ne. 0) then
               write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (55).'
               stop 55
             endif

             if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing pixel distribution (bin = '//coption(5:ico)//')'

             open (111,file=imagename(1:lenimg)//'.pd.dat',status='unknown')
             write (111,'(a)') '# IMGSTAT = Compute Image Statistics = '//compda
             write (111,'(a)') '# Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
             write (111,'(a)') '# Using '//clibname(lb:7)//' library version'//cfitsversion//' by William D Pence.'
             write (111,'(a)') '#'
             write (111,'(a)') '# PIXEL DISTRIBUTION FROM THE IMAGE: '//imagename(1:lenimg)
             write (111,'(a)') '# FOR A BIN SIZE: '//coption(5:ico)
             write (111,'(a)') '#'
             write (111,'(a)') '# TABULATED QUANTITIES:'
             write (111,'(a)') '#'
             write (111,'(a)') '#   N ......... Bin number'
             write (111,'(a)') '#   BINLOW .... Lower edge of the bin'
             write (111,'(a)') '#   BINCOO .... Bin-centered coordinate'
             write (111,'(a)') '#   NPIXELS ... Number of pixels per bin'
             write (111,'(a)') '#   NPIXBIN ... Product NPIXELS x BINCOO'
             write (111,'(a)') '#   FMEAN ..... Mean value within the bin'
             write (111,'(a)') '#   SIGMA ..... Standard deviation within the bin'
             write (111,'(a)') '#'
             write (111,'(a)') '#      N     BINLOW      BINCOO     NPIXELS     NPIXBIN     FMEAN       SIGMA'
             write (111,'(a)') '#'

             npx = 0
             do k=1,nbins-1
               funint(k) = (plane(k) + plane(k+1)) / 2.0d0
               do j=1,ny1
                 do i=1,nx1
                   if (fun2(i,j) .gt. almostzero) then
                     if (fun1(i,j) .gt. plane(k)) then
                       if (fun1(i,j) .le. plane(k+1)) then
                         npx = npx + 1
                         funone(k) = funone(k) + 1.0d0
                       endif
                     endif
                   endif
                 enddo
               enddo
             enddo

             if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Saving PDF data in '''//imagename(1:lenimg)//'.pd.dat'''

             funmax = 0.0d0
             funmin = 1.0d99
             do k=1,nbins-1
               if (lpdflin .or. funone(k) .gt. almostzero) then
                 if (funone(k) .lt. funmin) funmin = funone(k)
                 if (funone(k) .gt. funmax) funmax = funone(k)
               endif
             enddo
             pdfmin = funint(1)
             do k=1,nbins-1
               if (lpdflin .or. .not.lpdfpos .or. funint(k) .gt. almostzero) then
                 if (funone(k) .gt. almostzero) pdfmin = funint(k)
                 exit
               endif
             enddo
             
             do j=1,ny1
               do i=1,nx1
                 funpdf(i,j) = 0.0d0
               enddo
             enddo
             kk = 0
             do k=1,nbins-1
               funone(k) = max ( funone(k), funmin )
               nxy = 0
               fmean = 0.0d0
               do j=1,ny1
                 do i=1,nx1
                   if (fun2(i,j) .gt. almostzero) then
                     if (fun1(i,j) .gt. plane(k)) then
                       if (fun1(i,j) .le. plane(k+1)) then
                         funpdf(i,j) = funone(k)
                         nxy = nxy + 1
                         fmean = fmean + fun1(i,j)
                       endif
                     endif
                   endif
                 enddo
               enddo
               if (nxy .gt. 0) fmean = fmean / dble ( nxy )
               nxy = 0
               sigma = 0.0d0
               do j=1,ny1
                 do i=1,nx1
                   if (fun2(i,j) .gt. almostzero) then
                     if (fun1(i,j) .gt. plane(k)) then
                       if (fun1(i,j) .le. plane(k+1)) then
                         nxy = nxy + 1
                         sigma = sigma + (fun1(i,j) - fmean)**2
                       endif
                     endif
                   endif
                 enddo
               enddo
               if (nxy .gt. 1) sigma = sqrt ( sigma / dble ( nxy - 1 ) )
               if (funint(k) .ge. pdfmin .and. funone(k) .gt. 3.0d0) then
                 kk = kk + 1
                 write (111,'(i8,6(1pe12.4))') kk, plane(k), funint(k), funone(k), abs ( funint(k) ) * funone(k), fmean, sigma
               endif
             enddo

             close (111)
             
             outfits = imagename(1:lenimg-5)//'.pd.fits'
             iof = lastc ( outfits )

             if (cverbose .eq. '-verb2') 
     &         write (iotty,'(a,f4.1,a)') '   Writing PDF image '''//outfits(1:iof)//''''

             call wfits ( cfitsversion, nx1, ny1, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, funpdf, dx1, dy1, object
     &                  , crval1, crval2, outfits(1:iof), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12, cd21
     &                  , cd22, equinox, bzero, bscale, wave, datamin, datamax, history )

             deallocate ( funpdf )
             deallocate ( funone, funint, plane )
           endif
         endif
       
! Higher-order statistical moments.       
       
         mom2 = 0.0d0
         mom3 = 0.0d0
         mom4 = 0.0d0
         nxy = 0
         
         do j=1,ny1
           do i=1,nx1
             if (fun2(i,j) .gt. almostzero) then
               nxy = nxy + 1
               deltaf = fun1(i,j) - fmean2
               mom2 = mom2 + deltaf**2
               mom3 = mom3 + deltaf**3
               mom4 = mom4 + deltaf**4
             endif
           enddo
         enddo
         
         if (nxy .gt. 0) then
           rnpx = dble ( nxy )
           mom2 = mom2 / rnpx
           mom3 = mom3 / rnpx
           mom4 = mom4 / rnpx
         endif
         
         skewness2 = 0.0d0
         kurtosis2 = 0.0d0
         if (mom2 .gt. 0.0d0) then
           skewness2 = mom3 / sqrt ( mom2**3 )
           kurtosis2 = mom4 / mom2**2 - 3.0d0
         endif
         
         if (cverbose .eq. '-verb2') write (iotty,'(3(a,1pe15.7))') '   Maximum:', funmax2  , '      Mean:', fmean2 
         if (lmedian) then
         if (cverbose .eq. '-verb2') write (iotty,'(2(a,1pe15.7))') '   Minimum:', funmin2, '    Median:', fmedian2
         if (cverbose .eq. '-verb2') write (iotty,'(2(a,1pe15.7))') '      Mode:', fmode2
         else
         if (cverbose .eq. '-verb2') write (iotty,'(a,1pe15.7,a )') '   Minimum:', funmin2, '    Median:  not requested'
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(2(a,1pe15.7))') '    TotSum:', totsum2  , '     StDev:', sigma2                         
         if (cverbose .eq. '-verb2') write (iotty,'(2(a,1pe15.7))') '  Skewness:', skewness2, '  Kurtosis:', kurtosis2
         
         if (lfnam2) then
           if (cverbose .eq. '-verb2') write (iotty,'(2(a,1pe15.7))') '  EdgeMean:', maskedgemean, '  EdgeMedi:', maskedgemedian
           
           if (lsaveinmask .and. nxy .lt. 9) then
             if (cverbose .eq. '-verb2' ) then
               write (iotty,'(/a,i3)') '   IMGSTAT: WARNING: Number of pixels in the mask is too small:', nxy
             endif
           endif
         endif

         if (cverbose .eq. '-verb2') then
           write (iotty,'(a,1pe15.7,a)') '   Minimum:', funmin2nozero, ' over positive non-zero pixels'
           write (iotty,'(2(a,1pe15.7,a),1pe14.7,a)') '   TotFlux:', totsum2 * convftot, '  ', 'Mass:', totsum2 * convmass
     &                                              , ' Msun (distance '//cdistpc(1:icd)//' pc, convmass:', convmass, ')'
           write (cx,'(f7.2)') xcentr
           write (cy,'(f7.2)') ycentr
           ncx = int ( log10 ( xcentr ) ) + 4
           ncy = int ( log10 ( ycentr ) ) + 4
           write (iotty,'(a,2(1pe14.7),a,e14.7,a)') '   Maximum: '//cx(7-ncx+1:7)//' '//cy(7-ncy+1:7)//'  max/2 absize:', afwhmx
     &                                            , bfwhmx, ' arcsec', fwhm, ' mean'
           write (cx,'(f7.2)') mxco
           write (cy,'(f7.2)') myco
           ncx = int ( log10 ( mxco ) ) + 4
           ncy = int ( log10 ( myco ) ) + 4
           write (iotty,'(a,2(1pe14.7),a,e14.7,a)') '  Centroid: '//cx(7-ncx+1:7)//' '//cy(7-ncy+1:7)//'  from moments:', afwhm
     &                                            , bfwhm, ' arcsec', atheta, ' deg'
           write (iotty,'(2(a,1pe14.7),a)') '  Half-max:', funhalf, '  from moments:', fhmax, ' (size levels)'
           write (cx,'(f7.2)') xcentrx
           write (cy,'(f7.2)') ycentrx
           ncx = int ( log10 ( xcentrx ) ) + 4
           ncy = int ( log10 ( ycentrx ) ) + 4
           write (iotty,'(a,2(1pe14.7),a,2(e14.7),a)') '   Maximum: '//cx(7-ncx+1:7)//' '//cy(7-ncy+1:7)//'  max/2 xysize:', xfwhmx
     &                                            , yfwhmx, ' arcsec', xfwhmx * xAU / pc * distpc, yfwhmx * xAU / pc * distpc ,' pc'
           write (cx,'(f7.2)') xcentrf
           write (cy,'(f7.2)') ycentrf
           ncx = int ( log10 ( xcentrf ) ) + 4
           ncy = int ( log10 ( ycentrf ) ) + 4
           write (iotty,'(a,2(1pe14.7),a,2(e14.7),a)') '   max/100 xysize:'
     &                                , xfwfoot, yfwfoot, ' arcsec', xfwfoot * xAU / pc * distpc, yfwfoot * xAU / pc * distpc ,' pc'
           write (iotty,'(a,1pe14.7)') '  Number of pixels used in the calculations: ', real ( npxstat )
         endif
          
         if (linsave) then
           if (.not.lfnam2) then
             if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Saving image statistics in ''.+imgstat1'''
             open  (11,file='.+imgstat1'//cwavenum,status='unknown')     
             write (11,'(9(f56.19,1x))') funmin2, funmax2, fmean2, fmedian2, totsum2, sigma2, skewness2, kurtosis2, funmin2nozero
             write (11,'(9(f56.19,1x))') afwhmx, bfwhmx, fwhm
             write (11,'(9(f56.19,1x))') afwhm, bfwhm, atheta
             write (11,'(9(f56.19,1x))') funhalf, fhmax
             write (11,'(9(f56.19,1x))') totsum2 * convftot, totsum2 * convmass
             write (11,'(9(f56.19,1x))') xfwhmx, yfwhmx, xfwhmx * xAU / pc * distpc, yfwhmx * xAU / pc * distpc
             write (11,'(9(f56.19,1x))') xfwfoot, yfwfoot, xfwfoot * xAU / pc * distpc, yfwfoot * xAU / pc * distpc
             write (11,'(3i15,1x,2(f10.5))') npxstat, nx1, ny1, dx1, dy1
             close (11)
           else
             if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Saving image statistics in ''.+imgstat2'''
             open  (22,file='.+imgstat2'//cwavenum,status='unknown')     
             write (22,'(11(f56.19,1x))') funmin2, funmax2, fmean2, fmedian2, totsum2, sigma2, skewness2, kurtosis2, funmin2nozero
     &                                  , maskedgemean, maskedgemedian
             write (22,'(9(f56.19,1x))') afwhmx, bfwhmx, fwhm
             write (22,'(9(f56.19,1x))') afwhm, bfwhm, atheta
             write (22,'(9(f56.19,1x))') funhalf, fhmax
             write (22,'(9(f56.19,1x))') totsum2 * convftot, totsum2 * convmass
             write (22,'(9(f56.19,1x))') xfwhmx, yfwhmx, xfwhmx * xAU / pc * distpc, yfwhmx * xAU / pc * distpc
             write (22,'(9(f56.19,1x))') xfwfoot, yfwfoot, xfwfoot * xAU / pc * distpc, yfwfoot * xAU / pc * distpc
             write (22,'(3i15,1x,2(f10.5))') npxstat, nx1, ny1, dx1, dy1
             close (22)
           endif

           if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') funmin2
           if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') funmax2
           if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') fmean2
           if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') fmedian2
           if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') totsum2
           if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') sigma2
           if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') skewness2
           if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') kurtosis2
           if (lfnam2) then
             if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') maskedgemean
             if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') maskedgemedian
           endif
           if (cverbose .eq. '-verb2') write (iotty, '(f56.19)') funmin2nozero
         endif

         if (lmaxmin .and. outname .ne. '') then
           if (cverbose .eq. '-verb2') 
     &       write (iotty,'(/a,f4.1,a)') '   Writing implicitly used mask '''//outname(1:ion)//''''
     
           call wfits ( cfitsversion, nx1, ny1, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, fun2, dx1, dy1, object
     &                , crval1, crval2, outname(1:ion), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &                , equinox, bzero, bscale, wave, datamin, datamax, history )
         endif
       endif
                        
! Computing and writing stdev image.

       if (lstd) then

         if (.not.lfnam2) then 
           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Mask image for STDEV image calculation not supplied.'
         endif

         allocate ( funstd(nx1,ny1), stat=irc )
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (50).'
           stop 50
         endif

         if (.not.allocated ( fun2 )) then
           allocate ( fun2(nx1,ny1), stat=irc )
           if (irc .ne. 0) then
             write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (70).'
             stop 70
           endif
           do j=1,ny1
             do i=1,nx1
               fun2(i,j) = 1.0d0    
             enddo
           enddo
         endif

         write (cpmx,'(i4.4)') 2 * npoints + 1
         npmx = int ( log10 ( dble ( 2 * npoints + 1 ) ) ) + 1
         ramx = dble ( 2 * npoints + 1 ) * sqrt ( dx1 * dy1 )
!!         write (camx,'(f4.1)') dble ( 2 * npoints + 1 ) * sqrt ( dx1 * dy1 )
!!         namx = int ( log10 ( dble ( 2 * npoints + 1 ) * sqrt ( dx1 * dy1 ) ) ) + 1

         if (cverbose .eq. '-verb2') then
           write (iotty,'(a,1pe9.2,a)') '   Computing stdev in a circular window of D = '//cpmx(4-npmx+1:4)//' pixels = ', ramx
     &                                , ' arcsec'
         endif

         if (xpoints .le. 2.0d0) then
           rw2 = (xpoints + 0.5d0)**2
         else
           rw2 = (xpoints)**2
         endif
         npp1 = npoints + 1

         allocate ( cimask(-npoints:npoints,-npoints:npoints), stat=irc )
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (160).'
           stop 160
         endif

! The circular mask.
         
!!         if (cverbose .eq. '-verb2') write (iotty,'()')
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
!!           if (cverbose .eq. '-verb2') write (iotty,'(1000i2)') (cimask(k,l),k=-npoints,min(-npoints+38,npoints))
         enddo
         
         if (ltiming .and. cverbose .eq. '-verb2') then
           cpu = 0.0d0
           wal = 0.0d0
           cpu_sum = 0.0d0
           wal_sum = 0.0d0
!!           if (iotty .gt. 0) write (iotty,'(/a)') ' STARTED: '//cdate//' '//ctime
!!           if (iolog .gt. 0) write (iolog,'(/a)') ' STARTED: '//cdate//' '//ctime
         endif
         if (ltiming) then
           cpu = timer ( 'cpu', 0.0d0 )
           wal = timer ( 'wal', 0.0d0 )
         endif

         isumpct = 0.0d0
         numtotal = dble ( ny1 )
         numpart = 0.0d0
       
         if (cverbose .eq. '-verb2') then
           if (iotty .gt. 0) then
             write (iotty,'(11x,a)') '__________________________________________________'
             write (iotty,'(a)', advance='no') ' Progress: '
             if (iotty .gt. 0) endfile   ( iotty, err=666 )
 666         continue 
             if (iotty .gt. 0) backspace ( iotty )
           endif       
           if (iolog .gt. 0) then
             write (iolog,'(11x,a)') '__________________________________________________'
             write (iolog,'(a)', advance='no') ' Progress: '
           endif
         endif

         do j=1,ny1
           n1 = max ( j - npoints, 1 )
           n2 = min ( j + npoints, ny1 )
           if (cverbose .eq. '-verb2') then
             numpart = dble ( j )
             call showprogress ( iotty, iolog, lunix, numpart, numtotal, isumpct )
           endif

           do i=1,nx1
             m1 = max ( i - npoints, 1 )
             m2 = min ( i + npoints, nx1 )
             funstd(i,j) = 0.0d0
             
! First compute the mean value.

             nto = 0
             fmean = 0.0d0
             do m=m1,m2
               do n=n1,n2
                 if (fun2(m,n) .gt. almostzero) then
                   if (cimask(m-i,n-j) .ge. 1) then
                     nto = nto + 1
                     fmean = fmean + fun1(m,n)
                   endif
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
                 if (fun2(m,n) .gt. almostzero) then
                   if (cimask(m-i,n-j) .ge. 1) then
                     variance = variance + (fun1(m,n) - fmean)**2
                   endif
                 endif
               enddo
             enddo
             if (nto .ge. 2) then
               sigma1 = sqrt ( variance / dble ( nto - 1 ) )
             else
               sigma1 = 0.0d0
             endif
             if (fun2(i,j) .gt. almostzero) then
               funstd(i,j) = sigmas * sigma1
             endif
           enddo
         enddo

         outfits = outname(1:ion)
         if (ion .eq. 0) then
           outfits = imagename(1:lenimg)
           iof = lastc ( outfits )
           if (outfits(max(iof-4,1):iof) .ne. '.fits') then
             outfits = outfits(1:iof)//'.fits'
           endif                            
           iof = lastc ( outfits )
           if (outfits(max(iof-3,1):iof) .eq. '.std') then
             outfits = outfits(1:iof)//'.fits'
           elseif (outfits(max(iof-8,1):iof) .eq. '.std'//cpmx//'.fits') then
             outfits = outfits(1:iof)
           else
             outfits = outfits(1:max(iof-5,1))//'.std'//cpmx//'.fits'
           endif
         endif
         iof = lastc ( outfits )
         if (outfits(max(iof-4,1):iof) .ne. '.fits') then
           outfits = outfits(1:iof)//'.fits'
         else
           outfits = outfits(1:iof-5)//'.fits'
         endif                            
         iof = lastc ( outfits )

         if (cverbose .eq. '-verb2') 
     &     write (iotty,'(/a,f4.1,a)') '   Writing', sigmas, '*sigma stdev image '''//outfits(1:iof)//''''

         call wfits ( cfitsversion, nx1, ny1, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, funstd, dx1, dy1, object
     &              , crval1, crval2, outfits(1:iof), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history )

         if (lplus) then
           funhi = 0.0d0
           do j=1,ny1
             do i=1,nx1
               funhi = max ( funhi, fun1(i,j) )
             enddo
           enddo
           funlo = funhi / 1.0d4
           bin = 7.0d-2
           nbins = max ( nint ( (log10 ( funhi ) - log10 ( funlo ) ) / bin ) + 1, 2 )

           allocate ( plane(nbins+1), funbin(nx1,ny1), stat=irc )
       
           if (irc .ne. 0) then
             write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (12).'
             stop 12
           endif
           do k=1,nbins+1
             plane(k) = 10.0d0**( log10 ( funlo ) + dble ( k - 1 ) * bin )
           enddo
           write (iotty,'(/a)') '   Computing stdev distribution within intensity bins'
  
           open (111,file=outfits(1:iof)//'.bin.dat',status='unknown')
           write (111,'(a)') '# IMGSTAT = Compute Image Statistics = '//compda
           write (111,'(a)') '# Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
           write (111,'(a)') '# Using '//clibname(lb:7)//' library version'//cfitsversion//' by William D Pence.'
           write (111,'(a)') '#'
           write (111,'(a)') '# STDEV DISTRIBUTION FROM IMAGE: '//imagename(1:lenimg)
           write (111,'(a)') '#'
           write (111,'(a,1pe11.4)') '# FUNHI:', funhi
           write (111,'(a,1pe11.4)') '# FUNLO:', funlo
           write (111,'(a,1pe11.4)') '# LGBIN:', bin
           write (111,'(a,i5)') '# NBINS:', nbins
           write (111,'(a)') '#'
           write (111,'(a)') '# TABULATED QUANTITIES:'
           write (111,'(a)') '#'
           write (111,'(a)') '#   N ......... Bin number'
           write (111,'(a)') '#   BINLOW .... Lower edge of the bin'
           write (111,'(a)') '#   BINCOO .... Bin-centered coordinate'
           write (111,'(a)') '#   NPIXELS ... Number of pixels per bin'
           write (111,'(a)') '#   STDMEAN ... Mean stdev value within the bin'
           write (111,'(a)') '#   STDSTDM ... STDEV of STDMEAN within the bin'
           write (111,'(a)') '#'
           write (111,'(a)') '#      N    BINLOW      BINCOO      NPIXELS     STDMEAN    STDSTDM'
           write (111,'(a)') '#'
           kk = 0
           do k=1,nbins
             nto = 0
             stdmean = 0.0d0
             do j=1,ny1
               do i=1,nx1
                 if (fun2(i,j) .gt. almostzero) then
                   if (fun1(i,j) .gt. plane(k)) then
                     if (fun1(i,j) .le. plane(k+1)) then
                       nto = nto + 1
                       stdmean = stdmean + funstd(i,j)
                     endif
                   endif
                 endif
               enddo
             enddo
             if (nto .gt. 1) stdmean = stdmean / dble ( nto )
             variance = 0.0d0
             do j=1,ny1
               do i=1,nx1
                 if (fun2(i,j) .gt. almostzero) then
                   if (fun1(i,j) .gt. plane(k)) then
                     if (fun1(i,j) .le. plane(k+1)) then
                       nto = nto + 1
                       variance = variance + (funstd(i,j) - stdmean)**2
                     endif
                   endif
                 endif
               enddo
             enddo
             if (nto .gt. 1) then
               stdstdm = sqrt ( variance / dble ( nto - 1 ) )
             else
               stdstdm = 0.0d0
             endif
             if (nto .gt. 1) then
               kk = kk + 1
               write (111,'(i8,6(1pe12.4))') kk, plane(k), (plane(k) + plane(k+1)) / 2.0d0, dble ( nto ), stdmean, stdstdm
             endif
           enddo
           close (111)

           do j=1,ny1
             do i=1,nx1
               funbin(i,j) = 0.0d0
             enddo
           enddo
           do k=1,nbins
             do j=1,ny1
               do i=1,nx1
                 if (fun2(i,j) .gt. almostzero) then
                   if (fun1(i,j) .gt. plane(k)) then
                     if (fun1(i,j) .le. plane(k+1)) then
                       funbin(i,j) = plane(k)
                     endif
                   endif
                 endif
               enddo
             enddo
           enddo
           outfits = outfits(1:iof-5)//'.bin.fits'
           iof = lastc ( outfits )
           
           if (cverbose .eq. '-verb2') 
     &       write (iotty,'(/a,f4.1,a)') '   Writing binned intensity image '''//outfits(1:iof)//''''

           call wfits ( cfitsversion, nx1, ny1, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, funbin, dx1, dy1, object
     &                , crval1, crval2, outfits(1:iof), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &                , equinox, bzero, bscale, wave, datamin, datamax, history )

           deallocate ( plane, funbin )
         endif
         deallocate ( funstd, cimask )
       endif
       
! Computing and writing median- or mean- filtered image.

       if (lmodimage .or. lmedimage .or. lmeaimage .or. lskeimage .or. lkurimage) then

         if (.not.lfnam2) then 
           if (lmodimage .and. cverbose .eq. '-verb2') write (iotty,'(a)') '   Mask image for mode filtering not supplied.'
           if (lmedimage .and. cverbose .eq. '-verb2') write (iotty,'(a)') '   Mask image for median filtering not supplied.'
           if (lmeaimage .and. cverbose .eq. '-verb2') write (iotty,'(a)') '   Mask image for mean filtering not supplied.'
           if (lskeimage .and. cverbose .eq. '-verb2') write (iotty,'(a)') '   Mask image for skewness image not supplied.'
           if (lkurimage .and. cverbose .eq. '-verb2') write (iotty,'(a)') '   Mask image for kurtosis image not supplied.'
         endif

         allocate ( funm(nx1,ny1), funx(nx1*ny1), stat=irc )
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (150).'
           stop 150
         endif
         do j=1,ny1
           do i=1,nx1
             funm(i,j) = fun1(i,j)  
           enddo
         enddo
         if (.not.allocated ( fun2 )) then
           allocate ( fun2(nx1,ny1), stat=irc )
           if (irc .ne. 0) then
             write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (170).'
             stop 170
           endif
           do j=1,ny1
             do i=1,nx1
               fun2(i,j) = 1.0d0    
             enddo
           enddo
         endif
         do i=1,nx1*ny1
           funx(i) = 0.0d0    
         enddo

         write (cpmx,'(i4.4)') 2 * npoints + 1
         npmx = int ( log10 ( dble ( 2 * npoints + 1 ) ) ) + 1
         ramx = dble ( 2 * npoints + 1 ) * sqrt ( dx1 * dy1 )
!!         write (camx,'(f4.1)') dble ( 2 * npoints + 1 ) * sqrt ( dx1 * dy1 )
!!         namx = int ( log10 ( dble ( 2 * npoints + 1 ) * sqrt ( dx1 * dy1 ) ) ) + 1

         if (cverbose .eq. '-verb2') then
           if (lmodimage) 
     &         write (iotty,'(a,1pe9.2,a)') '   Mode filtering in a circular window of D = '//cpmx(4-npmx+1:4)//' pixels =', ramx
     &                                    , ' arcsec'
           if (lmedimage .and. lmedapprox) 
     &         write (iotty,'(a,1pe9.2,a)') '   Median filtering in a circular window of D = '//cpmx(4-npmx+1:4)//' pixels =', ramx
     &                                    , ' arcsec'
               write (iotty,'(a         )') '   Using an approximate binning method'
           if (lmedimage .and. .not.lmedapprox) 
     &         write (iotty,'(a,1pe9.2,a)') '   Median filtering in a circular window of D = '//cpmx(4-npmx+1:4)//' pixels =', ramx
     &                                    , ' arcsec'
               write (iotty,'(a         )') '   Using an accurate sorting method'
           if (lmeaimage) 
     &         write (iotty,'(a,1pe9.2,a)') '   Mean filtering in a circular window of D = '//cpmx(4-npmx+1:4)//' pixels =', ramx
     &                                    , ' arcsec'
           if (lskeimage) 
     &         write (iotty,'(a,1pe9.2,a)') '   Skewness in a circular window of D = '//cpmx(4-npmx+1:4)//' pixels =', ramx
     &                                    , ' arcsec'
           if (lkurimage) 
     &         write (iotty,'(a,1pe9.2,a)') '   Kurtosis in a circular window of D = '//cpmx(4-npmx+1:4)//' pixels =', ramx
     &                                    , ' arcsec'
         endif

         if (xpoints .le. 2.0d0) then
           rw2 = (xpoints + 0.5d0)**2
           rc2 = 0.0d0
         else
           rw2 = (xpoints)**2
           rc2 = (xpoints - 1.0d0)**2
         endif
         if (.not.lringwin) rc2 = 0.0d0
         npp1 = npoints + 1

         allocate ( cimask(-npoints:npoints,-npoints:npoints), stat=irc )
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   IMGSTAT: ERROR: Trouble allocating memory (160).'
           stop 160
         endif

! The circular mask could also serve as weights for weighted median.
         
         if (cverbose .eq. '-verb2') write (iotty,'()')
         do l=-npoints,npoints
           ry2 = (dble ( l ))**2
           do k=-npoints,npoints
             rx2 = (dble ( k ))**2
             rxy2 = rx2 + ry2
             if (rxy2 .le. rw2 .and. rxy2 .ge. rc2) then
!!               cimask(k,l) = npoints + 1 - nint ( sqrt ( rx2 + ry2 ) )  !<-- arbitrary (linear circular) weights
               cimask(k,l) = 1 
             else
               cimask(k,l) = 0
             endif
           enddo
           if (cverbose .eq. '-verb2') write (iotty,'(1000i2)') (cimask(k,l),k=-npoints,min(-npoints+38,npoints))
         enddo

         if (ltiming .and. cverbose .eq. '-verb2') then
           cpu = 0.0d0
           wal = 0.0d0
           cpu_sum = 0.0d0
           wal_sum = 0.0d0
!!           if (iotty .gt. 0) write (iotty,'(/a)') ' STARTED: '//cdate//' '//ctime
!!           if (iolog .gt. 0) write (iolog,'(/a)') ' STARTED: '//cdate//' '//ctime
         endif
         if (ltiming) then
           cpu = timer ( 'cpu', 0.0d0 )
           wal = timer ( 'wal', 0.0d0 )
         endif

         isumpct = 0.0d0
         numtotal = dble ( ny1 )
         numpart = 0.0d0
       
         if (cverbose .eq. '-verb2') then
           if (iotty .gt. 0) then
             write (iotty,'(11x,a)') '__________________________________________________'
             write (iotty,'(a)', advance='no') ' Progress: '
             if (iotty .gt. 0) endfile   ( iotty, err=888 )
 888         continue 
             if (iotty .gt. 0) backspace ( iotty )
           endif       
           if (iolog .gt. 0) then
             write (iolog,'(11x,a)') '__________________________________________________'
             write (iolog,'(a)', advance='no') ' Progress: '
           endif
         endif

! Mode filtering using an approximate binning method.

         if (lmodimage) then
           do j=1,ny1
             n1 = max ( j - npoints, 1 )
             n2 = min ( j + npoints, ny1 )
             if (cverbose .eq. '-verb2') then
               numpart = dble ( j )
               call showprogress ( iotty, iolog, lunix, numpart, numtotal, isumpct )
             endif
             do i=1,nx1
               if (fun2(i,j) .gt. almostzero) then
                 m1 = max ( i - npoints, 1 )
                 m2 = min ( i + npoints, nx1 )
                 nm = 0
                 funmin = 1.0d99
                 funmax = -1.0d99
                 do m=m1,m2
                   do n=n1,n2
                     if (fun2(m,n) .gt. almostzero .and. cimask(m-i,n-j) .ge. 1) then
                       nm = nm + 1
                       funx(nm) = fun1(m,n)
                       funmin = min ( funmin, funx(nm) )
                       funmax = max ( funmax, funx(nm) )
                     endif
                   enddo
                 enddo
                 if (nm .ne. 0) then
                   funm(i,j) = modeapprox ( nm, funx, funmin, funmax )
                 endif
               endif
             enddo
           enddo
         endif

! Median filtering using an approximate binning method.

         if (lmedimage .and. lmedapprox) then
           do j=1,ny1
             n1 = max ( j - npoints, 1 )
             n2 = min ( j + npoints, ny1 )
             if (cverbose .eq. '-verb2') then
               numpart = dble ( j )
               call showprogress ( iotty, iolog, lunix, numpart, numtotal, isumpct )
             endif
             do i=1,nx1
               if (fun2(i,j) .gt. almostzero) then
                 m1 = max ( i - npoints, 1 )
                 m2 = min ( i + npoints, nx1 )
                 nm = 0
                 funmin = 1.0d99
                 funmax = -1.0d99
                 do m=m1,m2
                   do n=n1,n2
                     if (fun2(m,n) .gt. almostzero .and. cimask(m-i,n-j) .ge. 1) then
                       nm = nm + 1
                       funx(nm) = fun1(m,n)
                       funmin = min ( funmin, funx(nm) )
                       funmax = max ( funmax, funx(nm) )
                     endif
                   enddo
                 enddo
                 if (nm .ne. 0) then
                   funm(i,j) = medapprox ( nm, funx, funmin, funmax )
                 endif
               endif
             enddo
           enddo
         endif
         
! Median filtering using an accurate sorting method.

         if (lmedimage .and. .not.lmedapprox) then

           do j=1,ny1
!!             n1 = max ( j - npoints, 1 )
!!             n2 = min ( j + npoints, ny1 )
             if (cverbose .eq. '-verb2') then
               numpart = dble ( j )
               call showprogress ( iotty, iolog, lunix, numpart, numtotal, isumpct )
             endif

             do i=1,nx1
               if (fun2(i,j) .gt. almostzero) then

!!                 m1 = max ( i - npoints, 1 )
!!                 m2 = min ( i + npoints, nx1 )
!!                 nm = 0

!!                 do m=m1,m2
!!                   do n=n1,n2
!!                     if (fun2(m,n) .gt. almostzero .and. cimask(m-i,n-j) .ge. 1) then
!!                       nm = nm + 1
!!                       funx(nm) = fun1(m,n)
!!                     endif
!!                   enddo
!!                 enddo

! Next spiral scanning of the sliding window is to automatically shrink the window radius when it comes to the mask edges.
! This avoids jumps of the filtered areas adjacent to the edges of the masked area, it makes the result smooth at egdes.

                 nm = 1
                 funx(nm) = fun1(i,j)

                 do l=1,npoints
                   
                   ldone = .false.

                   do nn=min(j+l,ny1),max(j-l+1,1),-1
                     ipl = min ( i + l, nx1 )
                     if (fun2(ipl,nn) .gt. almostzero .and. cimask(ipl-i,nn-j) .ge. 1) then
                       nm = nm + 1
                       funx(nm) = fun1(ipl,nn)
                     else
                       ldone = .true.
                     endif
                   enddo
                   do mm=min(i+l,nx1),max(i-l+1,1),-1
                     jml = max ( j - l, 1 )
                     if (fun2(mm,jml) .gt. almostzero .and. cimask(mm-i,jml-j) .ge. 1) then
                       nm = nm + 1
                       funx(nm) = fun1(mm,jml)
                     else
                       ldone = .true.
                     endif
                   enddo
                   do nn=max(j-l,1),min(j+l-1,ny1),1
                     iml = max ( i - l, 1 )
                     if (fun2(iml,nn) .gt. almostzero .and. cimask(iml-i,nn-j) .ge. 1) then
                       nm = nm + 1
                       funx(nm) = fun1(iml,nn)
                     else
                       ldone = .true.
                     endif
                   enddo
                   do mm=max(i-l,1),min(i+l-1,nx1),1
                     jpl = min ( j + l, ny1 )
                     if (fun2(mm,jpl) .gt. almostzero .and. cimask(mm-i,jpl-j) .ge. 1) then
                       nm = nm + 1
                       funx(nm) = fun1(mm,jpl)
                     else
                       ldone = .true.
                     endif
                   enddo
                   if (ldone) exit
                 enddo

                 if (nm .ne. 0) then
                   funm(i,j) = selectk ( (nm + 1) / 2, nm, funx )    !<-- SELECTK is slower by a factor of 1.5 than MEDAPPROX
                 endif
               endif
             enddo
           enddo
         endif

! Mean filtering.

         if (lmeaimage) then
           do j=1,ny1
             n1 = max ( j - npoints, 1 )
             n2 = min ( j + npoints, ny1 )
             if (cverbose .eq. '-verb2') then
               numpart = dble ( j )
               call showprogress ( iotty, iolog, lunix, numpart, numtotal, isumpct )
             endif
             do i=1,nx1
               if (fun2(i,j) .gt. almostzero) then
                 m1 = max ( i - npoints, 1 )
                 m2 = min ( i + npoints, nx1 )
                 nm = 0
                 funx(1) = 0.0d0
                 do m=m1,m2
                   do n=n1,n2
                     if (fun2(m,n) .gt. almostzero .and. cimask(m-i,n-j) .ge. 1) then
                       nm = nm + 1
                       funx(1) = funx(1) + fun1(m,n)
                     endif
                   enddo
                 enddo
                 if (nm .ne. 0) then
                   funm(i,j) = funx(1) / dble ( nm )
                 endif
               endif
             enddo
           enddo
         endif
         
! Skewness and kurtosis images.

         if (lskeimage .or. lkurimage) then
           do j=1,ny1
             n1 = max ( j - npoints, 1 )
             n2 = min ( j + npoints, ny1 )
             if (cverbose .eq. '-verb2') then
               numpart = dble ( j )
               call showprogress ( iotty, iolog, lunix, numpart, numtotal, isumpct )
             endif
             do i=1,nx1
               funm(i,j) = 0.0d0
               if (fun2(i,j) .gt. almostzero) then
                 m1 = max ( i - npoints, 1 )
                 m2 = min ( i + npoints, nx1 )
                 nm = 0
                 fmean = 0.0d0
                 do m=m1,m2
                   do n=n1,n2
                     if (fun2(m,n) .gt. almostzero .and. cimask(m-i,n-j) .ge. 1) then
                       nm = nm + 1
                       fmean = fmean + fun1(m,n)
                     endif
                   enddo
                 enddo
                 if (nm .ne. 0) then
                   fmean = fmean / dble ( nm )
                 endif
                 mom2 = 0.0d0
                 mom3 = 0.0d0
                 mom4 = 0.0d0
                 nm = 0
                 do m=m1,m2
                   do n=n1,n2
                     if (fun2(m,n) .gt. almostzero .and. cimask(m-i,n-j) .ge. 1) then
                       nm = nm + 1
                       deltaf = fun1(m,n) - fmean
                       mom2 = mom2 + deltaf**2
                       mom3 = mom3 + deltaf**3
                       mom4 = mom4 + deltaf**4
                     endif
                   enddo
                 enddo
                 if (nm .gt. 0) then
                   rnpx = dble ( nm )
                   mom2 = mom2 / rnpx
                   mom3 = mom3 / rnpx
                   mom4 = mom4 / rnpx
                 endif
                 if (mom2 .gt. 0.0d0) then
                   if (lskeimage) funm(i,j) = mom3 / sqrt ( mom2**3 )
                   if (lkurimage) funm(i,j) = mom4 / mom2**2 - 3.0d0
                 endif
               endif
             enddo
           enddo
         endif
         
         outfits = outname(1:ion)
         if (ion .eq. 0) then
           outfits = imagename(1:lenimg)
           iof = lastc ( outfits )
           if (outfits(max(iof-4,1):iof) .ne. '.fits') then
             outfits = outfits(1:iof)//'.fits'
           endif                            
           iof = lastc ( outfits )
           if (lmodimage) then
             if (outfits(max(iof-3,1):iof) .eq. '.mod') then
               outfits = outfits(1:iof)//'.fits'
             elseif (outfits(max(iof-8,1):iof) .eq. '.mod'//cpmx//'.fits') then
               outfits = outfits(1:iof)
             else
               outfits = outfits(1:max(iof-5,1))//'.mod'//cpmx//'.fits'
             endif
           endif
           if (lmedimage) then
             if (outfits(max(iof-3,1):iof) .eq. '.med') then
               outfits = outfits(1:iof)//'.fits'
             elseif (outfits(max(iof-8,1):iof) .eq. '.med'//cpmx//'.fits') then
               outfits = outfits(1:iof)
             else
               outfits = outfits(1:max(iof-5,1))//'.med'//cpmx//'.fits'
             endif
           endif
           if (lmeaimage) then
             if (outfits(max(iof-3,1):iof) .eq. '.mea') then
               outfits = outfits(1:iof)//'.fits'
             elseif (outfits(max(iof-8,1):iof) .eq. '.mea'//cpmx//'.fits') then
               outfits = outfits(1:iof)
             else
               outfits = outfits(1:max(iof-5,1))//'.mea'//cpmx//'.fits'
             endif
           endif
           if (lskeimage) then
             if (outfits(max(iof-3,1):iof) .eq. '.ske') then
               outfits = outfits(1:iof)//'.fits'
             elseif (outfits(max(iof-8,1):iof) .eq. '.ske'//cpmx//'.fits') then
               outfits = outfits(1:iof)
             else
               outfits = outfits(1:max(iof-5,1))//'.ske'//cpmx//'.fits'
             endif
           endif
           if (lkurimage) then
             if (outfits(max(iof-3,1):iof) .eq. '.kur') then
               outfits = outfits(1:iof)//'.fits'
             elseif (outfits(max(iof-8,1):iof) .eq. '.kur'//cpmx//'.fits') then
               outfits = outfits(1:iof)
             else
               outfits = outfits(1:max(iof-5,1))//'.kur'//cpmx//'.fits'
             endif
           endif
         endif
         iof = lastc ( outfits )
         if (outfits(max(iof-4,1):iof) .ne. '.fits') then
           outfits = outfits(1:iof)//'.fits'
         else
           outfits = outfits(1:iof-5)//'.fits'
         endif                            
         iof = lastc ( outfits )

         if (lmodimage .and. cverbose .eq. '-verb2') write (iotty,'(/a)') '   Writing mode-filtered '''//outfits(1:iof)//''''
         if (lmedimage .and. cverbose .eq. '-verb2') write (iotty,'(/a)') '   Writing median-filtered '''//outfits(1:iof)//''''
         if (lmeaimage .and. cverbose .eq. '-verb2') write (iotty,'(/a)') '   Writing mean-filtered '''//outfits(1:iof)//''''
         if (lskeimage .and. cverbose .eq. '-verb2') write (iotty,'(/a)') '   Writing skewness image '''//outfits(1:iof)//''''
         if (lkurimage .and. cverbose .eq. '-verb2') write (iotty,'(/a)') '   Writing kurtosis image '''//outfits(1:iof)//''''

         call wfits ( cfitsversion, nx1, ny1, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, funm, dx1, dy1, object
     &              , crval1, crval2, outfits(1:iof), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history )

         deallocate ( funm, funx, cimask )
       endif

       deallocate ( fun1 )
  
       if ((lmodimage .or. lmedimage .or. lmeaimage .or. lskeimage .or. lkurimage .or. lstd) .and. ltiming .and.
     &     cverbose .eq. '-verb2') then
         cpu_sum = cpu_sum + timer ( 'cpu', cpu )
         wal_sum = wal_sum + timer ( 'wal', wal )

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
    
!!         if (iotty .gt. 0) write (iotty,'(/a)') ' ENDED: '//cdate//' '//ctime//' CPU:'//ccpu//' s '//cccpu(ndc:)
!!     &                                      //' WALL:'//cwal//' s '//ccwal(ndw:)
!!         if (iolog .gt. 0) write (iolog,'(/a)') ' ENDED: '//cdate//' '//ctime//' CPU:'//ccpu//' s '//cccpu(ndw:)
!!     &                                      //' WALL:'//cwal//' s '//ccwal(ndw:)
!!         if (iolog .gt. 0) close ( iolog )
       endif

       if (cverbose .eq. '-verb2') write (iotty,'(/a)') ' Done.'
       
!!       stop !<-- commented out because it would lead to run-time messages about denormalized values, when using gfortran.
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine halfmaxsizes
       
     &            ( nx, ny, mask, image, dx, dy, almostzero, afwhm, bfwhm, fwhm, xcentr, ycentr )
!__________________________________________________________________________________________________________________________________
!
! Direct measurements of the minimum and maximum sizes (at peak half-maximum) for an arbitrary intensity distribution.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       
       logical       lcond1, lcond2, lcond3, lcond4, lcond5, lcond6, lcond7, lcond8

       integer       nn, i, j, ii, jj, imax, jmax, nx, ny, nedge, im1, jm1, ip1, jp1
 
       real*8        dx, dy, dblii2, dbljj2, fwhm, afwhm, bfwhm, imhalf, funmax, almostzero, radius2, drad2, radhalf, hwhm, sigm2
     &             , dlogim, sqpixel, xhalf, yhalf, xcentr, ycentr, imageijabs, image(nx,ny), mask(nx,ny)
!__________________________________________________________________________________________________________________________________
!                  
       funmax = -1.0d99
       imax = 0
       jmax = 0
       do j=1,ny
         do i=1,nx
           if (mask(i,j) .gt. almostzero) then
             if (abs ( image(i,j) ) .ge. funmax) then
               imax = i
               jmax = j
               funmax = abs ( image(i,j) )
             endif
           endif
         enddo
       enddo
       
       hwhm = 0.0d0
       afwhm = 0.0d0
       bfwhm = 0.0d0
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
             
             if (imageijabs .ge. imhalf .and. mask(i,j) .gt. 0.9d0) then
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
           afwhm = fwhm
           bfwhm = fwhm
         endif
       endif
       
       if (xcentr .lt. almostzero) xcentr = dble ( nx ) / 2.0d0       
       if (ycentr .lt. almostzero) ycentr = dble ( ny ) / 2.0d0       

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine xyhalfmaxsizes
       
     &            ( nx, ny, mask, image, dx, dy, almostzero, xfwhm, yfwhm, fwhm, xcentr, ycentr )
!__________________________________________________________________________________________________________________________________
!
! Direct measurements of the minimum and maximum sizes (at peak half-maximum) for an arbitrary intensity distribution.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       
       logical       lcond1, lcond2, lcond3, lcond4

       integer       nn, i, j, ii, jj, imax, jmax, nx, ny, nedge, im1, jm1, ip1, jp1
 
       real*8        dx, dy, dblii2, dbljj2, fwhm, xfwhm, yfwhm, imhalf, funmax, almostzero, radius2, drad2, radhalf, hwhm, sigm2
     &             , dlogim, sqpixel, xhalf, yhalf, xcentr, ycentr, imageijabs, image(nx,ny), mask(nx,ny)
!__________________________________________________________________________________________________________________________________
!                  
       funmax = -1.0d99
       imax = 0
       jmax = 0
       do j=1,ny
         do i=1,nx
           if (mask(i,j) .gt. almostzero) then
             if (abs ( image(i,j) ) .ge. funmax) then
               imax = i
               jmax = j
               funmax = abs ( image(i,j) )
             endif
           endif
         enddo
       enddo
       
       if (imax .eq. nx/2 .and. jmax .ne. ny/2) jmax = ny/2
       if (imax .eq. nx/2+1 .and. jmax .ne. ny/2+1) jmax = ny/2+1
       if (jmax .eq. ny/2 .and. imax .ne. nx/2) imax = nx/2
       if (jmax .eq. ny/2+1 .and. imax .ne. nx/2+1) imax = nx/2+1
       
       xfwhm = 0.0d0
       yfwhm = 0.0d0
       xcentr = 0.0d0
       ycentr = 0.0d0
       imhalf = abs ( funmax ) / 2.0d0
       sqpixel = dx * dy
       nedge = 0
       hwhm = 0.0d0
       
       if (imax .gt. 0 .and. jmax .gt. 0) then
         
         j = jmax
         jj = j - jmax
         do i=1,nx
           ii = i - imax
           imageijabs = abs ( image(i,j) )
           
           if (imageijabs .ge. imhalf .and. mask(i,j) .gt. 0.9d0) then
             im1 = max ( i - 1, 1 )
             ip1 = min ( i + 1, nx )
             
             lcond1 = abs ( image(im1,j  ) ) .lt. imhalf .and. abs ( image(im1,j  ) ) .ge. 0.0d0
             lcond2 = abs ( image(ip1,j  ) ) .lt. imhalf .and. abs ( image(ip1,j  ) ) .ge. 0.0d0
             
             if (lcond1 .or. lcond2) then
               nedge = nedge + 1
               dblii2 = dble ( ii )**2
               dbljj2 = dble ( jj )**2
               radius2 = (dblii2 + dbljj2) * sqpixel
               dlogim = log ( imageijabs ) - log ( abs ( imhalf ) )
               radhalf = 0.0d0
               xhalf = 0.0d0
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
               endif
               if (lcond2) then
                 nn = nn + 1
                 drad2 = (dble ( ip1 - imax )**2 + dbljj2) * sqpixel - radius2
                 sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(ip1,j) ) ))
                 radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
                 xhalf = xhalf + dble ( ip1 )
               endif
               
               if (nn .gt. 0) then
                 radhalf = radhalf / dble ( nn )
                 xhalf = xhalf / dble ( nn )
               endif
               hwhm = hwhm + radhalf
               xcentr = xcentr + xhalf
             endif
           endif
         enddo
         
         if (nedge .gt. 0) then
           hwhm = hwhm / dble ( nedge )
           xcentr = xcentr / dble ( nedge )
           fwhm = max ( 2.0d0 * hwhm, sqrt ( sqpixel ) )
           xfwhm = fwhm
         endif

         nedge = 0
         hwhm = 0.0d0
         i = imax
         ii = i - imax
         do j=1,ny
           jj = j - jmax
           imageijabs = abs ( image(i,j) )
           
           if (imageijabs .ge. imhalf .and. mask(i,j) .gt. 0.9d0) then
             jm1 = max ( j - 1, 1 )
             jp1 = min ( j + 1, ny )
             
             lcond3 = abs ( image(i,jm1) ) .lt. imhalf .and. abs ( image(i,jm1) ) .ge. 0.0d0
             lcond4 = abs ( image(i,jp1) ) .lt. imhalf .and. abs ( image(i,jp1) ) .ge. 0.0d0
             
             if (lcond3 .or. lcond4) then
               nedge = nedge + 1
               sqpixel = dx * dy
               dblii2 = dble ( ii )**2
               dbljj2 = dble ( jj )**2
               radius2 = (dblii2 + dbljj2) * sqpixel
               dlogim = log ( imageijabs ) - log ( abs ( imhalf ) )
               radhalf = 0.0d0
               yhalf = 0.0d0
               nn = 0
                     
! Gaussian interpolation to obtain accurate values of full width at half-maximum for Gaussian-like shapes.
! From two equations: ln(G1) - ln(G2) = -x1^2 + x2^2 / (2 sigma^2) and ln(G1) - ln(Gh) = -x1^2 + xh^2 / (2 sigma^2):
! (2 sigma^2) = (x2^2 - x1^2) / (ln(G1) - ln(G2)) and xh = sqrt ( x1^2 + (2 sigma^2) * (ln(G1) - ln(Gh)) )
                 
               if (lcond3) then
                 nn = nn + 1
                 drad2 = (dblii2 + dble ( jm1 - jmax )**2) * sqpixel - radius2
                 sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(i,jm1) ) ))
                 radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
                 yhalf = yhalf + dble ( jm1 )
               endif
               if (lcond4) then
                 nn = nn + 1
                 drad2 = (dblii2 + dble ( jp1 - jmax )**2) * sqpixel - radius2
                 sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(i,jp1) ) ))
                 radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
                 yhalf = yhalf + dble ( jp1 )
               endif
               
               if (nn .gt. 0) then
                 radhalf = radhalf / dble ( nn )
                 yhalf = yhalf / dble ( nn )
               endif
               hwhm = hwhm + radhalf
               ycentr = ycentr + yhalf
             endif
           endif
         enddo
       endif

       if (nedge .gt. 0) then
         hwhm = hwhm / dble ( nedge )
         ycentr = ycentr / dble ( nedge )
         fwhm = max ( 2.0d0 * hwhm, sqrt ( sqpixel ) )
         yfwhm = fwhm
       endif
       
       fwhm = (xfwhm + yfwhm) / 2.0d0
       
       if (xcentr .lt. almostzero) xcentr = dble ( nx ) / 2.0d0       
       if (ycentr .lt. almostzero) ycentr = dble ( ny ) / 2.0d0       

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine xyfootsizes
       
     &            ( nx, ny, mask, image, dx, dy, almostzero, xfwhm, yfwhm, fwhm, xcentr, ycentr )
!__________________________________________________________________________________________________________________________________
!
! Direct measurements of the minimum and maximum sizes (at peak half-maximum) for an arbitrary intensity distribution.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       
       logical       lcond1, lcond2, lcond3, lcond4

       integer       nn, i, j, ii, jj, imax, jmax, nx, ny, nedge, im1, jm1, ip1, jp1
 
       real*8        dx, dy, dblii2, dbljj2, fwhm, xfwhm, yfwhm, imhalf, funmax, almostzero, radius2, drad2, radhalf, hwhm, sigm2
     &             , dlogim, sqpixel, xhalf, yhalf, xcentr, ycentr, imageijabs, image(nx,ny), mask(nx,ny)
!__________________________________________________________________________________________________________________________________
!                  
       funmax = -1.0d99
       imax = 0
       jmax = 0
       do j=1,ny
         do i=1,nx
           if (mask(i,j) .gt. almostzero) then
             if (abs ( image(i,j) ) .ge. funmax) then
               imax = i
               jmax = j
               funmax = abs ( image(i,j) )
             endif
           endif
         enddo
       enddo
       
       if (imax .eq. nx/2 .and. jmax .ne. ny/2) jmax = ny/2
       if (imax .eq. nx/2+1 .and. jmax .ne. ny/2+1) jmax = ny/2+1
       if (jmax .eq. ny/2 .and. imax .ne. nx/2) imax = nx/2
       if (jmax .eq. ny/2+1 .and. imax .ne. nx/2+1) imax = nx/2+1
       
       xfwhm = 0.0d0
       yfwhm = 0.0d0
       xcentr = 0.0d0
       ycentr = 0.0d0
!!       imax = nint ( xcentr )
!!       jmax = nint ( ycentr )
       imhalf = abs ( funmax ) / 100.0d0
       sqpixel = dx * dy
       nedge = 0
       hwhm = 0.0d0
       
       if (imax .gt. 0 .and. jmax .gt. 0) then
         
         j = jmax
         jj = j - jmax
         do i=1,nx
           ii = i - imax
           imageijabs = abs ( image(i,j) )
           
           if (imageijabs .ge. imhalf .and. mask(i,j) .gt. 0.9d0) then
             im1 = max ( i - 1, 1 )
             ip1 = min ( i + 1, nx )
             
             lcond1 = abs ( image(im1,j  ) ) .lt. imhalf .and. abs ( image(im1,j  ) ) .ge. 0.0d0
             lcond2 = abs ( image(ip1,j  ) ) .lt. imhalf .and. abs ( image(ip1,j  ) ) .ge. 0.0d0
             
             if (lcond1 .or. lcond2) then
               nedge = nedge + 1
               dblii2 = dble ( ii )**2
               dbljj2 = dble ( jj )**2
               radius2 = (dblii2 + dbljj2) * sqpixel
               dlogim = log ( imageijabs ) - log ( abs ( imhalf ) )
               radhalf = 0.0d0
               xhalf = 0.0d0
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
               endif
               if (lcond2) then
                 nn = nn + 1
                 drad2 = (dble ( ip1 - imax )**2 + dbljj2) * sqpixel - radius2
                 sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(ip1,j) ) ))
                 radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
                 xhalf = xhalf + dble ( ip1 )
               endif
               
               if (nn .gt. 0) then
                 radhalf = radhalf / dble ( nn )
                 xhalf = xhalf / dble ( nn )
               endif
               hwhm = hwhm + radhalf
               xcentr = xcentr + xhalf
             endif
           endif
         enddo
         
         if (nedge .gt. 0) then
           hwhm = hwhm / dble ( nedge )
           xcentr = xcentr / dble ( nedge )
           fwhm = max ( 2.0d0 * hwhm, sqrt ( sqpixel ) )
           xfwhm = fwhm
         endif

         nedge = 0
         hwhm = 0.0d0
         i = imax
         ii = i - imax
         do j=1,ny
           jj = j - jmax
           imageijabs = abs ( image(i,j) )
           
           if (imageijabs .ge. imhalf .and. mask(i,j) .gt. 0.9d0) then
             jm1 = max ( j - 1, 1 )
             jp1 = min ( j + 1, ny )
             
             lcond3 = abs ( image(i,jm1) ) .lt. imhalf .and. abs ( image(i,jm1) ) .ge. 0.0d0
             lcond4 = abs ( image(i,jp1) ) .lt. imhalf .and. abs ( image(i,jp1) ) .ge. 0.0d0
             
             if (lcond3 .or. lcond4) then
               nedge = nedge + 1
               sqpixel = dx * dy
               dblii2 = dble ( ii )**2
               dbljj2 = dble ( jj )**2
               radius2 = (dblii2 + dbljj2) * sqpixel
               dlogim = log ( imageijabs ) - log ( abs ( imhalf ) )
               radhalf = 0.0d0
               yhalf = 0.0d0
               nn = 0
                     
! Gaussian interpolation to obtain accurate values of full width at half-maximum for Gaussian-like shapes.
! From two equations: ln(G1) - ln(G2) = -x1^2 + x2^2 / (2 sigma^2) and ln(G1) - ln(Gh) = -x1^2 + xh^2 / (2 sigma^2):
! (2 sigma^2) = (x2^2 - x1^2) / (ln(G1) - ln(G2)) and xh = sqrt ( x1^2 + (2 sigma^2) * (ln(G1) - ln(Gh)) )
                 
               if (lcond3) then
                 nn = nn + 1
                 drad2 = (dblii2 + dble ( jm1 - jmax )**2) * sqpixel - radius2
                 sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(i,jm1) ) ))
                 radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
                 yhalf = yhalf + dble ( jm1 )
               endif
               if (lcond4) then
                 nn = nn + 1
                 drad2 = (dblii2 + dble ( jp1 - jmax )**2) * sqpixel - radius2
                 sigm2 = drad2 / (log ( imageijabs ) - log ( abs ( image(i,jp1) ) ))
                 radhalf = radhalf + sqrt ( radius2 + sigm2 * dlogim )
                 yhalf = yhalf + dble ( jp1 )
               endif
               
               if (nn .gt. 0) then
                 radhalf = radhalf / dble ( nn )
                 yhalf = yhalf / dble ( nn )
               endif
               hwhm = hwhm + radhalf
               ycentr = ycentr + yhalf
             endif
           endif
         enddo
       endif

       if (nedge .gt. 0) then
         hwhm = hwhm / dble ( nedge )
         ycentr = ycentr / dble ( nedge )
         fwhm = max ( 2.0d0 * hwhm, sqrt ( sqpixel ) )
         yfwhm = fwhm
       endif
       
       fwhm = (xfwhm + yfwhm) / 2.0d0
       
       if (xcentr .lt. almostzero) xcentr = dble ( nx ) / 2.0d0       
       if (ycentr .lt. almostzero) ycentr = dble ( ny ) / 2.0d0       

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine momentsizes
       
     &            ( nxmin, nxmax, nymin, nymax, nx, ny, mask, image, dx, dy, momxco, momyco, afwhm, bfwhm, atheta, afwhmx
     &            , bfwhmx, equivrad, elongation, almostzero )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       integer       nx, ny, l, m, nxmin, nxmax, nymin, nymax
 
       real*8        dx, dy, pi, sqpixel, sourcearea, dblm, dbll, convtodegrees, sqreightlogtwo, intsxx, intsyy, intsxy
     &             , xcorsc, ycorsc, momp2, momm2, sqrmo, asig, bsig, intsd, ints0m, ints0, almostzero, intx, inty, int0, int0m
     &             , intd, nonzero, afwhmx, bfwhmx
     
       real*8        image(nx,ny), mask(nx,ny), momxco, momyco, afwhm, bfwhm, atheta, equivrad, elongation, momxx, momyy, momxy

       parameter   ( pi = 3.14159265358979d0, convtodegrees = 180.0d0 / pi, sqreightlogtwo = sqrt ( 8.0d0 * log ( 2.0d0 ) ) )
!__________________________________________________________________________________________________________________________________
!                  
       sqpixel = dx * dy
       intx = 0.0d0                        
       inty = 0.0d0                        
       int0 = 0.0d0
       nonzero = 0.0d0
       asig = 0.0d0
       bsig = 0.0d0
       
       do m=nymin,nymax
         dblm = dble ( m )                 
         do l=nxmin,nxmax            
           dbll = dble ( l )               
           
           if (mask(l,m) .gt. 0.9d0) then
             intd = image(l,m) * sqpixel
             int0 = int0 + intd
             intx = intx + dbll * intd
             inty = inty + dblm * intd
             if (image(l,m) .gt. almostzero) then
               nonzero = nonzero + 1.0d0
             endif
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
           
           if (mask(l,m) .gt. 0.9d0) then
             intsd = image(l,m) * sqpixel
             ints0 = ints0 + intsd
             xcorsc = dbll - momxco
             ycorsc = dblm - momyco
             intsxx = intsxx + xcorsc * xcorsc * intsd
             intsyy = intsyy + ycorsc * ycorsc * intsd
             intsxy = intsxy + xcorsc * ycorsc * intsd
           endif
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

       bfwhmx = bfwhmx / sqrt ( elongation )
       afwhmx = bfwhmx * elongation
       
       sourcearea = nonzero * sqpixel
       equivrad = sqrt ( sourcearea / pi )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine rfits 
       
     &            ( nx, ny, bunit, ctype1, ctype2, rp1, rp2, xr, yr, fun, dx, dy, object, ra, dec, fname, iotty, iolog, creator
     &            , beam, funmin, funmax, blank, rot1, rot2, cd11, cd12, cd21, cd22, equinox, bzero, bscale, wave, datamn, datamx
     &            , history, cverbose )
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       simple, extend, anynull, myfile
       character*(*) ctype1, ctype2, creator, fname, bunit, object, history, cverbose
       character*3   cnan
       character*25  cdatetime
       character*80  comment
       integer       i, j, nx, ny, stat, blksize, group, naxis, pcount, gcount, bitpix, unit, iotty, iolog, maxdim
     &             , rw, nfound, blank, stat05, stat06, stat07, stat08, stat09, stat10, stat11, stat12, stat13, stat14, stat15
     &             , stat16, stat17, stat18, stat19, stat20, stat21, stat22, stat23, stat24, stat25, stat26, stat27, stat28, stat29
     &             , stat30, stat31, stat32, stat0
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
         creator = 'IMGSTAT'
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

       if (cverbose .eq. '-verb2') then
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
             if (fun(i,j) .lt. datamn) datamn = fun(i,j)
             if (fun(i,j) .gt. datamx) datamx = fun(i,j)
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

       call ftpkys ( unit, 'CREATOR', creator,'Alexander Menshchikov, DAp IRFU CEA Saclay', status )
       call ftpcom ( unit,' ',status)
       call ftpkys ( unit, 'DATE'  , cdatetime, 'creation date and time', status )

! Write all parameters.

       dxdeg = - dx * as2deg
       dydeg = dy * as2deg

       call ftpkyd ( unit, 'BZERO'  , bzero  ,13 , 'zero point in scaling equation', status )
       call ftpkyd ( unit, 'BSCALE' , bscale ,13 , 'linear factor in scaling equation' , status )
       call ftpkyd ( unit, 'DATAMAX', datamx ,13 , 'maximum data value', status )
       call ftpkyd ( unit, 'DATAMIN', datamn ,13 , 'minimum data value', status )
!!       if (blank .gt. 0)
!!     & call ftpkyj ( unit, 'BLANK'  , blank      , 'value used for undefined array elements', status )
       call ftpkys ( unit, 'BUNIT'  , bunit      , 'physical units of the array values', status )
       call ftpkys ( unit, 'CTYPE1' , ctype1     , 'name of the coordinate axis', status )
       call ftpkys ( unit, 'CTYPE2' , ctype2     , 'name of the coordinate axis', status )
       call ftpkyd ( unit, 'CRPIX1' , rp1    ,13 , 'coordinate system reference pixel', status )
       call ftpkyd ( unit, 'CRPIX2' , rp2    ,13 , 'coordinate system reference pixel', status )
       call ftpkyd ( unit, 'CROTA1' , rot1   ,13 , 'coordinate system rotation angle', status )
       call ftpkyd ( unit, 'CROTA2' , rot2   ,13 , 'coordinate system rotation angle', status )
       call ftpkyd ( unit, 'CD1_1'  , cd11   ,13 , 'linear projection matrix', status )
       call ftpkyd ( unit, 'CD1_2'  , cd12   ,13 , 'linear projection matrix', status )
       call ftpkyd ( unit, 'CD2_1'  , cd21   ,13 , 'linear projection matrix', status )
       call ftpkyd ( unit, 'CD2_2'  , cd22   ,13 , 'linear projection matrix', status )
       call ftpkyd ( unit, 'CRVAL1' , xr     ,13 , 'coordinate value at reference pixel', status )
       call ftpkyd ( unit, 'CRVAL2' , yr     ,13 , 'coordinate value at reference pixel', status )
       call ftpkyd ( unit, 'CDELT1' , dxdeg  ,13 , 'coordinate increment along axis', status )
       call ftpkyd ( unit, 'CDELT2' , dydeg  ,13 , 'coordinate increment along axis', status )
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