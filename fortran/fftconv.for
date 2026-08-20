
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       program fftconv
!__________________________________________________________________________________________________________________________________
!
! Convolution using a Fast Fourier Transform algorithm RLFT3 from Press, Teukolsky, Vetterling, & Flannery (1992):
! "Numerical recipes in FORTRAN. The art of scientific computing", Cambridge University Press, 2nd edition.
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!__________________________________________________________________________________________________________________________________
!                 
       implicit      none

       logical       lfname, lunix, lnormalize, lnoscaling, lsavekern, lnoimage, lpsfname
                                  
       character*3   dot
       character*4   cpsfbeam
       character*6   cfitsversion, cverbose
       character*7   clibname, cx, cy
       character*8   ctime, cfpptime, ckerntype, cdowhat
       character*10  cdate
       character*21  compda
       character*80  object, creator, ctype1, ctype2, comm(10), bunit, history
       character*1000 filename, filenamec, outname, psfname, arg1, arg2, arg3, arg4, arg5, arg6, arg7

       integer       lastc, lasts, i, j, fnl, fnlc, irc, ndate, nxo, nyo, ixo, iyo, ia1, ia2, ia3, ia4, ia5, ia6, ia7, isp1, iotty
     &             , iolog, ion, lb, blank, nxbig, nybig, nx, ny, i0, j0, nn, jx, psfl, nxp, nyp, ifit, ncx, ncy, idot1, idot2
     &             , idot3, icom, ii, jj, inam

       real          fitsvers
       real*8        wave, crpix1, crpix2, crval1, crval2, ra, dec, correct, dx, dy, eps, pi, beam, fract, pixfactor, funmin
     &             , funmax, crota1, crota2, equinox, bzero, bscale, datamin, datamax, almostzero, beamrfits, expo, cd11
     &             , cd12, cd21, cd22, delx, dely, rad2, argex, sigm2m, hwhm, psfbeam, dxk, dyk

       real*8    , allocatable :: fun(:,:), kern(:,:), kernw(:,:), psf(:,:)

       parameter   ( eps = 0.01d0, pi = 3.14159265358979d0, almostzero = 1.0d-20, iotty = 6, iolog = 0, dot = '•' )

       external      rfits, wfits, when, osystem, ftvers, lastc, lasts, getfitshead, fftconvol
!__________________________________________________________________________________________________________________________________
!
! Determine operating system.

       call osystem ( lunix )
       creator = 'FFT Convolution'
       call when ( lunix, ctime, cdate, ndate, 4 )
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

! Get command line parameters.

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
       ion = 0
       psfl = 0

       if (cverbose .eq. '-verb2') then
         write (iotty,'( )')
         if (arg1(1:ia1) .eq. ':') write (iotty,'(a$)') '  '
         write (iotty,'(a$)') ' FFTCONV '//dot//' FFTransform & Convolve FITS Images '//dot//' '//compda
       endif
       if (arg1(1:ia1) .eq. ':') stop
       if (cverbose .eq. '-verb2') write (iotty,'()')
       if (cverbose .eq. '-verb2') write (iotty,'( a)') ' Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
       if (cverbose .eq. '-verb2') write (iotty,'( a)') ' Using RLFT3 from Numerical Recipes for F77 by William H Press et al.'
       if (cverbose .eq. '-verb2') write (iotty,'( a)') ' Using '//clibname(lb:7)//' library version'//cfitsversion
     &                                                //' by William D Pence.'
       if (ia2 .eq. 0) then
         write (iotty,'( )')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USAGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')
         write (iotty,'(a)') ' fftconv <fwhm> {<expo>[,<frac>]|<psfn>} [<imagein>] [<OPTION>] [-o <imageout>] [-verb{0|1|2}]'
         write (iotty,'( )')   
         write (iotty,'(a)') ' This utility convolves a FITS image <imagein> creating <imageout> as output.'  
         write (iotty,'(a)') ' If <imagein> is not supplied, then only a kernel image is created. An integer'
         write (iotty,'(a)') ' number instead of <imagein> is the length (in pixels) of a vertical filament'
         write (iotty,'(a)') ' with the kernel width and profile to be saved.'
         write (iotty,'( )')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')   
         write (iotty,'(a)') ' <fwhm> .................. full width at half-maximum of the convolution kernel'
         write (iotty,'(a)') ' <expo>[,<frac>]|<psfn> .. exponent for power-law kernels or kernel image name'
         write (iotty,'(a)') ' <imagein> ............... name of an image to convolve or (if name not given)'
         write (iotty,'(a)') '                           save an image of the convolution kernel'
         write (iotty,'( )')   
         write (iotty,'(a)') ' --------------------------------- <OPTION> ----------------------------------'
         write (iotty,'( )')   
         write (iotty,'(a)') ' normalize ...... scale convolved image <imageout> to original peak value'
         write (iotty,'(a)') ' noscaling ...... do not rescale <imageout> for total flux conservation'
         write (iotty,'( )')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ NOTES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')
         write (iotty,'(a)') ' <fwhm> =-0, <expo> =-0: compute and save an image of Fourier amplitudes'
         write (iotty,'(a)') ' <fwhm> > 0, <expo> = 0: kernel is a Gaussian with a width <fwhm>'
         write (iotty,'(a)') ' <fwhm> > 0, <expo> < 0: kernel is a power law with <expo>, divided by <fwhm>'
         write (iotty,'(a)') ' <fwhm> > 0, <expo> > 0: kernel is a Gauss-Moffat with <fwhm>, <expo>[,<frac>]'
         write (iotty,'( )')
         write (iotty,'(a)') ' <fwhm> < 0: kernel is a cylinder with a width |<fwhm>|'
         write (iotty,'(a)') ' <fwhm> = 0: kernel is an exponential exp(-r^2)'
         write (iotty,'( )')
         write (iotty,'(a)') ' <frac> is the (optional) proportion of Gaussian and Moffat:'
         write (iotty,'(a)') ' kernel = (1-frac) * gauss + frac * moffat'
         write (iotty,'(a)') ' if <expo> > 0 and <frac> is not given, then frac = 1 and kernel = moffat'
         write (iotty,'( )')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         stop 99
       endif
       
       if (arg1(1:ia1) .eq. '-0' .and. arg2(1:ia2) .eq. '-0') then
         cdowhat = 'fftransf'
       else
         cdowhat = 'convolve'
       endif
       psfname = 'none'
       lpsfname = .false.
       read (arg1(1:ia1),*,err=6) beam
       goto 7
   6   continue
         write (iotty,'(/a)') '   FFTCONV: ERROR: Trouble reading BEAM from this string: '''//arg1(1:ia1)//'''.'
         stop 99
   7   continue
       if (arg2(1:1) .ne. '/') then  !<-- to avoid a weird special case: reading a slash would give 0...
         icom = index ( arg2(1:ia2), ',' )
         if (icom .eq. 0) then
           read (arg2(1:ia2),*,err=8) expo
           fract = 1.0d0
           arg2 = arg2(1:ia2)//',1'
           icom = ia2 + 1
           ia2 = ia2 + 2
         else
           read (arg2(1:icom-1),*,err=8) expo
           read (arg2(icom+1:ia2),*,err=8) fract
           if (fract .lt. almostzero) then
             ia2 = icom-1
             icom = 0
           endif
         endif
         goto 10
       endif
   8   continue
       icom = 0
       expo = 0.0d0
       fract = 0.0d0
       read (arg2(1:ia2),'(a)',err=9) psfname
       psfl = lastc ( psfname )
       inquire ( file=psfname(1:psfl), exist=lpsfname )
       if (.not.lpsfname) then
         psfname = psfname(1:psfl)//'.fits'
         psfl = lastc ( psfname )
         inquire ( file=psfname(1:psfl), exist=lpsfname )
       endif
       goto 10
   9   continue
         write (iotty,'(/a)') '   FFTCONV: ERROR: Trouble reading EXPO or PSFNAME from this string: '''//arg2(1:ia2)//'''.'
         stop 99
  10   continue
       if (psfname .eq. 'none') then
         if (beam .gt. almostzero) then
           if ((expo .gt. -almostzero .and. expo .lt. almostzero) .or. fract .lt. almostzero) then
             ckerntype = 'gaussian'
           elseif (expo .lt. -almostzero) then
             ckerntype = 'powerlaw'
           else
             ckerntype = 'moffatfn'
           endif
         else
           if (beam .lt. -almostzero) then
             ckerntype = 'cylinder'
           else
             ckerntype = 'exponent'
           endif
         endif
         pixfactor = 1.0d0
       else
         ckerntype = 'psfimage'
         inam = index ( arg2(1:ia2), '.r' )
         ifit = index ( arg2(inam+1:ia2), '.' )
         cpsfbeam = arg2(inam+2:inam+ifit-1)
         inam = index ( cpsfbeam, 'p' )
         cpsfbeam(inam:inam) = '.'
         read (cpsfbeam,*) psfbeam
         pixfactor = beam / psfbeam
         psfbeam = beam
       endif
       beam = abs ( beam )

       filename = arg3(1:ia3)
       fnl = lastc ( filename )
       lfname = .false.
       lnoimage = .true.
       nn = 0

! If the 3rd parameter is a number (not a filename), then it is the length (in pixels) of a vertical linear filament.

       if (ia3 .gt. 0) then
         if (arg3(1:1) .ne. '/') then  !<-- to avoid a weird special case: reading a slash would give 0...
           read (arg3(1:ia3),*,err=998) nn
           lfname = .false.
           lnoimage = .true.
           goto 999
         endif
 998     continue
         nn = -1

! Check if the input file exists and open it.

         inquire ( file=filename(1:fnl), exist=lfname )
         if (.not.lfname .and. filename(max(fnl-4,1):fnl) .ne. '.fits') then
           filename = filename(1:fnl)//'.fits'
           fnl = lastc ( filename )
           inquire ( file=filename(1:fnl), exist=lfname )
         endif
         if (lfname) lnoimage = .false.
       endif
 999   continue
       
!!       if (lnoimage .and. abs ( beam ) .lt. almostzero) then
!!         write (iotty,'(/a)') '   FFTCONV: ERROR: With NOIMAGE and PSFNAME, convolution BEAM must not be zero.'
!!         stop 99
!!       endif

       lnormalize = .false.
       lnoscaling = .false.
       if (arg4(1:ia4) .eq. 'normalize') lnormalize = .true.
       if (arg4(1:ia4) .eq. 'noscaling') lnoscaling = .true.
       lsavekern = .false.
       nxp = 1
       nyp = 1

       if (cverbose .eq. '-verb1') then
         write (iotty,'(a)') ' FFTCONV: '//arg1(1:ia1)//' '//arg2(1:ia2)//' '//arg3(1:ia3)//' '//arg4(1:ia4)//' '//arg5(1:ia5)
     &                                   //' '//arg6(1:ia6)//' '//arg7(1:ia7)
       endif
       if (.not.lpsfname .and. ckerntype .eq. 'psfimage') then
         write (iotty,'(/a)') '   FFTCONV: ERROR: Image '''//psfname(1:psfl)//''' not found.'
         stop 0
       endif
       if (.not.lfname .and. nn .eq. -1) then
         write (iotty,'(/a)') '   FFTCONV: ERROR: Image '''//filename(1:fnl)//''' not found.'
         stop 1
       endif
       
       if (.not.lnoimage) then
         isp1 = lasts ( filename ) + 1
         if (cdowhat .eq. 'convolve') then
           if (ckerntype .eq. 'psfimage') then
             ifit = index ( psfname(1:psfl), '.fits' ) - 1
             inam = lasts ( arg2(1:ia2) )
!!             filenamec = 'kern.'//ckerntype//'.k'//arg1(1:ia1)//'as.'//arg2(inam+1:ifit)//'.fits'
             filenamec = filename(isp1:fnl-5)//'.c'//arg1(1:ia1)//'as.'//arg2(inam+1:ifit)//'.fits'
           else
             idot1 = index ( arg1(1:ia1), '.' )
             idot2 = index ( arg2(1:ia2), '.' )
             idot3 = index ( arg2(idot2+1:ia2), '.' ) + idot2
             if (idot1 .gt. 0) arg1(idot1:idot1) = 'p'
             if (idot2 .gt. 0) arg2(idot2:idot2) = 'p'
             if (idot3 .gt. 0) arg2(idot3:idot3) = 'p'
             if (ckerntype .eq. 'moffatfn') then
               if (icom .eq. 0) then
                 filenamec = filename(isp1:fnl-5)//'.c'//arg1(1:ia1)//'x'//arg2(1:ia2)//'as.fits'
               else
                 filenamec = filename(isp1:fnl-5)//'.c'//arg1(1:ia1)//'x'//arg2(1:icom-1)//'as~'//arg2(icom+1:ia2)//'.fits'
               endif
             else
               filenamec = filename(isp1:fnl-5)//'.c'//arg1(1:ia1)//'as.fits'
             endif
           endif
         else
           filenamec = filename(isp1:fnl-5)//'.fourier.amp.fits'
         endif
         fnlc = lastc ( filenamec )
         outname = filenamec(1:fnlc)
         if (arg4(1:2) .eq. '-o') outname = arg5
         if (arg5(1:2) .eq. '-o') outname = arg6
         ion = lastc ( outname )
         if (outname(ion-4:ion) .ne. '.fits') then
           outname = outname(1:ion)//'.fits'
           ion = lastc ( outname )
         endif
       endif

       if (ckerntype .eq. 'psfimage') then
         call getfitshead ( psfname(1:psfl), nxp, nyp, dx, dy, bunit )
       endif
       if (.not.lnoimage) then
         call getfitshead ( filename(1:fnl), nxo, nyo, dx, dy, bunit )
       else
         lsavekern = .true.
         dx = 0.1d0
         dy = 0.1d0
         dxk = dx
         dyk = dy
         if (ckerntype .eq. 'psfimage') then
           nxo = nxp
           nyo = nyp
         else if (ckerntype .eq. 'moffatfn') then
           nxo = 2 * nint ( 5.0d0 * beam / dx + exp ( 7.0d0 / abs ( expo )**(0.363636d0) ) ) + 1
           nyo = 2 * nint ( 5.0d0 * beam / dy + exp ( 7.0d0 / abs ( expo )**(0.363636d0) ) ) + 1
         else if (ckerntype .eq. 'powerlaw') then
           nxo = 2 * nint ( 5.0d0 * beam / dx + exp ( 7.0d0 / abs ( expo )**(0.363636d0) ) ) + 1
           nyo = 2 * nint ( 5.0d0 * beam / dy + exp ( 7.0d0 / abs ( expo )**(0.363636d0) ) ) + 1
         else if (ckerntype .ne. 'exponent') then
           nxo = 2 * nint ( 5.0d0 * beam / dx ) + 1
           nyo = 2 * nint ( 5.0d0 * beam / dy ) + 1
         else
           nxo = 11
           nyo = 11
         endif
         if (nn .eq. 0) then
           nxbig = nxo
           nybig = nyo
         else
           nxbig = 2801
           nybig = 2801
         endif
         blank = 0
         if (cverbose .eq. '-verb2') write (iotty,'()')
       endif

! The next error message informs that fftconv cannot handle huge images, because the fortran loop integer value would
! exceed the limiting value of ~2e9, which would lead to code crashes.
       
       if (nxo .gt. 32760 .or. nyo .gt. 32760) then
         write (iotty,'(/a)') '   FFTCONV: ERROR: Image '''//filename(1:fnl)//''' has too large dimensions (> 32760).'
         stop 3
       endif

       do i=1,10
         comm(i) = 'no comment'
       enddo

       allocate ( fun(nxo,nyo), kern(nxo,nyo), psf(nxp,nyp), stat=irc )

       if (irc .ne. 0) then
         write (iotty,'(/a,2i5)') '   FFTCONV: ERROR: Trouble allocating memory (10) for arrays of sizes: ', nxo, nyo
         stop 10
       endif

       if (lpsfname) then
         write (cx,'(i7)') nxp
         write (cy,'(i7)') nyp
         ncx = int ( log10 ( dble ( nxp ) ) ) + 1
         ncy = int ( log10 ( dble ( nyp ) ) ) + 1

         if (cverbose .eq. '-verb2' .and. .not.lnoimage) write (iotty,'()')
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Reading ('//cx(7-ncx+1:7)//' x '//cy(7-ncy+1:7)//') '''
     &                                                 //psfname(1:psfl)//''''

         call rfits ( nxp, nyp, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, psf, dxk, dyk, object, ra, dec
     &              , psfname(1:psfl), iotty, 0, creator, beamrfits, funmin, funmax, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history, cverbose )
           
         if (cverbose .eq. '-verb2') write (iotty,'(a,1x,2(1pe14.7))') '   Minmax values in the image:', funmin, funmax
       endif

       if (.not.lnoimage) then
         write (cx,'(i7)') nxo
         write (cy,'(i7)') nyo
         ncx = int ( log10 ( dble ( nxo ) ) ) + 1
         ncy = int ( log10 ( dble ( nyo ) ) ) + 1

         if (cverbose .eq. '-verb2' .and. .not.lpsfname) write (iotty,'()')
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Reading ('//cx(7-ncx+1:7)//' x '//cy(7-ncy+1:7)//') '''
     &                                                 //filename(1:fnl)//''''

         call rfits ( nxo, nyo, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, fun, dx, dy, object, ra, dec
     &              , filename(1:fnl), iotty, 0, creator, beamrfits, funmin, funmax, blank, crota1, crota2, cd11, cd12, cd21
     &              , cd22, equinox, bzero, bscale, wave, datamin, datamax, history, cverbose )
           
         if (cverbose .eq. '-verb2') write (iotty,'(a,1x,2(1pe14.7))') '   Minmax values in the image:', funmin, funmax
       endif
         
       if (dx .le. 1.0d-10 .or. dy .le. 1.0d-10) then
         write (iotty,'(/a,2(1pe15.7))') '   FFTCONV: ERROR: Wrong pixel size(s):', dx, dy
         stop 99
       endif
       
! To compute and save just an image of Fourier amplitudes, we need to redefine (increase) the image sizes.

       if (cdowhat .eq. 'fftransf') then
         ixo = 0
         iyo = 0
         do i=1,14
           nx = mod ( 2**i, nxo )
           if (nx .lt. ixo) exit
           ixo = nx
         enddo
         do j=1,14
           ny = mod ( 2**j, nyo )
           if (ny .lt. iyo) exit
           iyo = ny
         enddo
         nx = 2**i
         ny = 2**j

! Save the image temporarily, to redefine its dimensions.

         ixo = nxo
         iyo = nyo
         do j=1,iyo
           do i=1,ixo
             kern(i,j) = fun(i,j)
           enddo
         enddo
         
         deallocate ( fun )

         nxo = max ( nx, ny, 1 )
         nyo = nxo
              
         allocate ( fun(nxo,nyo), stat=irc )
  
         if (irc .ne. 0) then
           write (iotty,'(/a,2i5)') '   FFTCONV: ERROR: Trouble allocating memory (11) for an array of sizes: ', nxo, nyo
           stop 11
         endif

         do j=1,nyo
           do i=1,nxo
             fun(i,j) = kern(min(i,ixo),min(j,iyo))
           enddo
         enddo
       endif

       if (lnoimage) then
         crpix1 = (nxo - 1) / 2 + 1
         crpix2 = (nyo - 1) / 2 + 1
         crval1 = 0.0d0
         crval2 = 0.0d0
         wave = 100.0d0
         object = ckerntype//' kernel'
         history = 'no history'
         ctype1 = 'RA---TAN'
         ctype2 = 'DEC--TAN'
         bzero = 0.0d0
         bscale = 1.0d0
         bunit = 'MJy/sr'
         equinox = 2.0d3
         crota1 = 0.0d0
         crota2 = 0.0d0
         cd11 = -dx
         cd12 = 0.0d0
         cd21 = 0.0d0
         cd22 = dy
       endif
       
       call fftconvol ( iotty, nxo, nyo, dx, dy, fun, beam, expo, nxp, nyp, psf, fract, psfbeam, cdowhat, ckerntype, lnoimage
     &                , lnormalize, lnoscaling, cverbose )

       if (lsavekern) then
         if (ckerntype .eq. 'psfimage') then
           nxo = nxp
           nyo = nyp
         endif
         i0 = nxbig / 2 - nxo / 2
         j0 = nybig / 2 - nyo / 2

         allocate ( kernw(nxbig,nybig), stat=irc )

         if (irc .ne. 0) then
           write (iotty,'(/a,2i5)') '   FFTCONV: ERROR: Trouble allocating memory (20) for an array of sizes: ', nxbig, nybig
           stop 20
         endif

         do j=1,nybig
           do i=1,nxbig
             kernw(i,j) = 0.0d0
           enddo
         enddo
         
         do jx=1,2*nn+1
           do j=1,nybig
             jj = j - j0 - (nn + 1) + jx
             do i=1,nxbig
               ii = i - i0
               if (ii .ge. 1 .and. ii .le. nxo .and. jj .ge. 1 .and. jj .le. nyo) then
                 kernw(i,j) = max ( kernw(i,j), fun(ii,jj) )
               endif
             enddo
           enddo
         enddo

         if (ckerntype .eq. 'psfimage') then
           ifit = index ( arg2(1:ia2), '.fits' ) - 1
           if (ifit .eq. -1) ifit = ia2
           inam = lasts ( arg2(1:ia2) )
           filenamec = 'kern.'//ckerntype//'.k'//arg1(1:ia1)//'as.'//arg2(inam+1:ifit)//'.fits'
         else
           idot1 = index ( arg1(1:ia1), '.' )
           idot2 = index ( arg2(1:ia2), '.' )
           idot3 = index ( arg2(icom+1:ia2), '.' ) + icom
           if (idot1 .gt. 0) arg1(idot1:idot1) = 'p'
           if (idot2 .gt. 0) arg2(idot2:idot2) = 'p'
           if (idot3 .gt. 0) arg2(idot3:idot3) = 'p'
           if (icom .eq. 0) then
             filenamec = 'kern.gaus'//ckerntype(1:4)//'.k'//arg1(1:ia1)//'x'//arg2(1:ia2)//'as.fits'
           else
             filenamec = 'kern.gaus'//ckerntype(1:4)//'.k'//arg1(1:ia1)//'x'//arg2(1:icom-1)//'as~'//arg2(icom+1:ia2)//'.fits'
           endif
         endif
         fnlc = lastc ( filenamec )

         if (cverbose .eq. '-verb2')
     &      write (iotty,'(a)') '   Writing kernel '''//filenamec(1:fnlc)//''''
     
         dxk = dxk * pixfactor
         dyk = dyk * pixfactor
     
         call wfits ( cfitsversion, nxbig, nybig, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, kernw, dxk, dyk, object
     &              , crval1, crval2, filenamec(1:fnlc), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history )

         deallocate ( kernw )

         if (ckerntype .eq. 'moffatfn') then
           i0 = nxbig / 2 + 1

           allocate ( kernw(nxbig,nybig), stat=irc )

           if (irc .ne. 0) then
             write (iotty,'(/a,2i5)') '   FFTCONV: ERROR: Trouble allocating memory (20) for an array of sizes: ', nxbig, nybig
             stop 20
           endif

           hwhm = beam / 2.0d0
           correct = (2.0d0**(1.0d0 / dble ( expo )) - 1.0d0)
           sigm2m = hwhm**2

           do j=1,nybig
             do i=1,nxbig
               kernw(i,j) = 0.0d0
             enddo
           enddo
           do jx=1,2*nn+1
             j0 = nybig / 2 + 1 - (nn + 1) + jx
             do j=1,nybig
               dely = (dble ( j - j0 ) * dy)**2
               do i=1,nxbig
                 delx = (dble ( i - i0 ) * dx)**2
                 rad2  = delx + dely
                 argex = rad2 / sigm2m
                 kernw(i,j) = 1.0d0 / (1.0d0 + correct * argex)**expo
               enddo
             enddo
           enddo
           idot1 = index ( arg1(1:ia1), '.' )
           idot2 = index ( arg2(1:ia2), '.' )
           idot3 = index ( arg2(icom+1:ia2), '.' ) + icom
           if (idot1 .gt. 0) arg1(idot1:idot1) = 'p'
           if (idot2 .gt. 0) arg2(idot2:idot2) = 'p'
           if (idot3 .gt. 0) arg2(idot3:idot3) = 'p'
           if (icom .eq. 0) then
             filenamec = 'kern.'//ckerntype//'.k'//arg1(1:ia1)//'x'//arg2(1:ia2)//'as.fits'
           else
             filenamec = 'kern.'//ckerntype//'.k'//arg1(1:ia1)//'x'//arg2(1:icom-1)//'as~'//arg2(icom+1:ia2)//'.fits'
           endif
           fnlc = lastc ( filenamec )

           if (cverbose .eq. '-verb2')
     &        write (iotty,'(a)') '   Writing kernel '''//filenamec(1:fnlc)//''''
     
           call wfits ( cfitsversion, nxbig, nybig, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, kernw, dxk, dyk, object
     &                , crval1, crval2, filenamec(1:fnlc), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12, cd21
     &                , cd22, equinox, bzero, bscale, wave, datamin, datamax, history )

           deallocate ( kernw )
         endif
       endif

       if (cdowhat .eq. 'fftransf') then
         do j=1,nyo
           do i=1,nxo
             fun(i,j) = fun(i,j) / sqrt ( dble ( nxo * nyo ))
           enddo
         enddo
         dx = 1.0d0 / ( dx * dble ( nxo ) )
         dy = 1.0d0 / ( dy * dble ( nyo ) )
       endif

       if (.not.lnoimage) then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Writing output image '''//outname(1:ion)//''''

         call wfits ( cfitsversion, nxo, nyo, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, fun, dx, dy, object, crval1
     &              , crval2, outname(1:ion), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &              , equinox, bzero, bscale, wave, datamin, datamax, history )
       endif

       deallocate ( fun, kern, psf )

       if (cverbose .eq. '-verb2') write (iotty,'(/a)') ' Done.'
       
!!       stop !<-- commented out because it would lead to run-time messages about denormalized values, when using gfortran.
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine rfits 
       
     &            ( nx, ny, bunit, ctype1, ctype2, rp1, rp2, xr, yr, fun, dx, dy, object, ra, dec, fname, iotty, iolog
     &            , creator, beam, funmin, funmax, blank, rot1, rot2, cd11, cd12, cd21, cd22, equinox, bzero, bscale, wave
     &            , datamn, datamx, history, cverbose )
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

       call ftghpr ( unit, maxdim, simple, bitpix, naxis, naxes, pcount, gcount, extend, stat )

! Read all parameters.

       call ftgkys ( unit, 'CREATOR', creator, comment, stat )
       
       if (stat .eq. 0 .and. (comment(1:21) .eq. 'Alexander Menshchikov' .or. comment(1:22) .eq. 'Alexander Men''shchikov')) then
         myfile = .true.
       else
         myfile = .false.
         stat = 0
         creator = 'FFTCONV'
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

       call ftg2dd ( unit, group, dble ( blank ), nx, naxes(1), naxes(2), fun, anynull, stat )

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
       cnan='NAN'
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
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       simple, extend
       character*(*) fname, ctype1, ctype2, object, creator, ctime, cdate, bunit, cfitsversion, history
       character*19  cdatetime
       integer       i, j, nx, ny, status, blocksze, group, naxis, pcount, gcount, bitpix, unit, naxes(2), blank
       real*8        xr, yr, dx, dy, rp1, rp2, ra, dec, beam, bzero, bscale, as2deg, rot1, rot2, equinox, datamn
     &             , datamx, dxdeg, dydeg, wave, cd11, cd12, cd21, cd22, almostzero, pi, angle
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

       call ftpkys ( unit, 'CREATOR', creator,'Alexander Men''shchikov, DAp IRFU CEA Saclay', status )
       call ftpcom ( unit,' ',status)
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