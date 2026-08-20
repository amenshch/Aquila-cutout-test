
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       program operate
!__________________________________________________________________________________________________________________________________
!
! Operate on two FITS images.
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay  
!__________________________________________________________________________________________________________________________________
!                 
       implicit      none

       integer       iotty
       parameter   ( iotty = 6 )

       logical       lfnama, lfnamb, lunix
                                          
       character*3   cnan, dot
       character*6   cfitsversion, cverbose         
       character*7   clibname, cnexp, cx, cy, cnite
       character*8   ctime, cfpptime
       character*10  cdate
       character*21  compda
       character*80  object, bunit, bunit1, bunit2, ctype1, ctype2, ctype11, ctype12, ctype21, ctype22, object1, object2, creator
     &             , history, history1, history2
       character*500 filenama, filenamb, outname, arg1, arg2, arg3, arg4, arg5, arg6

       integer       i, j, ion, fnlen1, fnlen2, ndate, nx, ny, nx1, ny1, nx2, ny2, firstb, lastc, lasts, ia1, im1, jm1, ip1
     &             , jp1, ia2, ia3, ia4, ia5, ia6, isp1, isp2, blank, blank1, blank2, irc, lb, nexpand, iter, inexp, itermax, ncx
     &             , ncy, ind, nt, l, inite

       real          fitsvers
       real*8        dx, dy, dx1, dy1, dx2, dy2, ra, dec, beam, bbody, almostzero, beam1, beam2, tem1, tem2, dtem
     &             , equinox, equinox1, equinox2, ra1, ra2, dec1, dec2, bzero, bscale, bzero1, bzero2, bscale1, bscale2
     &             , funmin, funmax, crpix11, crpix12, crval11, crval12, datamin1, datamax1, crpix1, crpix2, crval1, crval2, bbmax
     &             , crota11, crota12, crota1, crota2, datamin2, datamax2, datamin, datamax, crpix21, crpix22, crval21, crval22
     &             , crota21, crota22, cd11, cd12, cd21, cd22, wave, wave1, wave2, cd111, cd112, cd121, cd122, cd211, cd212, cd221
     &             , cd222, muH2, amu, speedolight, dust2gas, frequency, frequency0, beta, opacity, opacity0, factor, tdustx, bbmin
     &             , blackbody, bbody1, bbody2

       real*8, allocatable :: funa(:,:), funb(:,:), func(:,:), tdust(:), bbodyx(:)

       parameter   ( almostzero = 1.0d-30, dot = '•', nt = 30000 )
       parameter   ( muH2 = 2.8d0, speedolight = 2.99792458d10, amu = 1.6605402D-24, frequency0 = 1.0d12, opacity0 = 10.0d0 !! 9.306
     &             , dust2gas = 1.0d-2, beta = 2.0d0 )

       external      firstb, lastc, rfits, wfits, when, planck, osystem, ftvers, lasts, getfitshead
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
       cnan = 'NAN'

! Get command line parameters (file names).

       call getarg ( 1, arg1 )
       call getarg ( 2, arg2 )
       call getarg ( 3, arg3 )
       call getarg ( 4, arg4 )
       call getarg ( 5, arg5 )
       call getarg ( 6, arg6 )

       ia1 = lastc ( arg1 )
       ia2 = lastc ( arg2 )
       ia3 = lastc ( arg3 )
       ia4 = lastc ( arg4 )
       ia5 = lastc ( arg5 )
       ia6 = lastc ( arg6 )

       cverbose = '-verb2'
       if (ia6 .gt. 0 .and. arg6(1:ia6-1) .eq. '-verb') cverbose = arg6(1:ia6)
       if (ia5 .gt. 0 .and. arg5(1:ia5-1) .eq. '-verb') cverbose = arg5(1:ia5)
       if (ia4 .gt. 0 .and. arg4(1:ia4-1) .eq. '-verb') cverbose = arg4(1:ia4)
       if (ia3 .gt. 0 .and. arg3(1:ia3-1) .eq. '-verb') cverbose = arg3(1:ia3)

       if (cverbose .eq. '-verb2') then
         write (iotty,'( )')
         if (arg1(1:ia1) .eq. ':') write (iotty,'(a$)') '  '
         write (iotty,'(a$)') ' OPERATE '//dot//' Operate on Two FITS Images '//dot//' '//compda
       endif
       if (arg1(1:ia1) .eq. ':') stop
       if (cverbose .eq. '-verb2') write (iotty,'()')
       if (cverbose .eq. '-verb2') write (iotty,'( a)') ' Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
       if (cverbose .eq. '-verb2') write (iotty,'( a)') ' Using '//clibname(lb:7)//' library version'//cfitsversion
     &                                                //' by William D Pence.'
       if (ia3 .eq. 0 .or. 
     &     (arg2(1:ia2) .ne. '+' .and. arg2(1:ia2) .ne. '-' .and. arg2(1:ia2) .ne. 'x' .and. arg2(1:ia2) .ne. '~'.and.
     &      arg2(1:ia2) .ne. '%' .and. arg2(1:ia2) .ne. 's' .and. arg2(1:ia2) .ne. 'o' .and. arg2(1:ia2) .ne. 'l' .and. 
     &      arg2(1:ia2) .ne. 'g' .and. arg2(1:ia2) .ne. 'r' .and. arg2(1:ia2) .ne. 'e' .and. arg2(1:ia2) .ne. 'p' .and.
     &      arg2(1:ia2) .ne. 'h' .and. arg2(1:1) .ne. 'd' .and. arg2(1:1) .ne. 't' .and. arg2(1:1) .ne. 'q' .and.
     &      arg2(1:1) .ne. 'i')) then
         write (iotty,'( )')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USAGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')
         write (iotty,'(a)') ' operate <image1> <operator> <image2> [-o <image>] [-verb{0|1|2}]'
         write (iotty,'( )')   
         write (iotty,'(a)') ' This simple utility allows various modifications of two input FITS images'
         write (iotty,'(a)') ' <image1> and <image2>, producing an output with an optional name <image>.'
         write (iotty,'( )')   
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OPERATOR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')   
         write (iotty,'(a)') ' + .......... summation: image = image1 + image2'
         write (iotty,'(a)') ' - .......... subtraction: image = image1 - image2'
         write (iotty,'(a)') ' x .......... multiplication: image = image1 * image2'
         write (iotty,'(a)') ' % .......... division: image = image1 / image2'
         write (iotty,'(a)') ' ~ or \~ .... relative difference: image = (image1 - image2) / image2'
         write (iotty,'(a)') ' s .......... squared difference: image = (image1 - image2)^2'
         write (iotty,'(a)') ' o .......... if { image1 < image2} then image = 0, else image = image1'
         write (iotty,'(a)') ' l .......... minimization: image = min { image1, image2 }'
         write (iotty,'(a)') ' g .......... maximization: image = max { image1, image2 }'
         write (iotty,'(a)') ' r .......... replacement of non-zero pixels: image = image1 ~> image2'
         write (iotty,'(a)') ' e .......... extend masks of image1 by adding connected pixels of image2'
         write (iotty,'(a)') ' p .......... expand masks of areas image1 > image2 until image1 = image2'
         write (iotty,'(a)') ' h .......... replacement of FITS header: image = image2 (header1)'
         write (iotty,'(a)') ' d[<wave>] .. compute surfdensity = factor * intens1 / Bnu(temps2)'
         write (iotty,'(a)') ' t[<wave>] .. compute temperature = Invert { factor * intens1 / surfdens2 }'
         write (iotty,'(a)') ' i[<wave>] .. compute intensities = surfdens1 * Bnu(temps2) / factor'
         write (iotty,'(a)') ' q[<Td>] .... equalization of Tdust: image = image1 * Bnu(Td) / Bnu(temps2)'
         write (iotty,'(a)')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FITS KEYWORDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')   
         write (iotty,'(a)') ' Supported keywords: HEADER CREATOR OBJECT CTYPE1 CTYPE2 CRPIX1 CRPIX2'
         write (iotty,'(a)') ' CRVAL1 CRVAL2 CROTA1 CROTA2 CDELT1 CDELT2 CD1_1 CD1_2 CD2_1 CD2_2 RA DEC'
         write (iotty,'(a)') ' NAXIS1 NAXIS2 BZERO BSCALE BLANK DATAMIN DATAMAX EQUINOX WAVE HISTORY'
         write (iotty,'( )')   
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         if (ia2 .gt. 0 .and. ia3 .gt. 0) then
           write (iotty,'(/a)') '   OPERATE: ERROR: Unknown operator: '''//arg2(1:ia2)//''''
         endif
         stop 99
       endif
       
       if (cverbose .eq. '-verb1')   !! cverbose .eq. '-verb0' .or. 
     &    write (iotty,'(a)') ' OPERATE: '//arg1(1:ia1)//' '//arg2(1:ia2)//' '//arg3(1:ia3)//' '//arg4(1:ia4)//' '//arg5(1:ia5)
     &                                    //' '//arg6(1:ia6)

       filenama = arg1(1:ia1)
       fnlen1 = lastc ( filenama )
       
       filenamb = arg3(1:ia3)
       fnlen2 = lastc ( filenamb )

! Check if the input file exists and open it.

       inquire ( file=filenama(1:fnlen1), exist=lfnama )
       if (.not.lfnama .and. filenama(fnlen1-4:fnlen1) .ne. '.fits') then
         filenama = filenama(1:fnlen1)//'.fits'
         fnlen1 = lastc ( filenama )
         inquire ( file=filenama(1:fnlen1), exist=lfnama )
       endif

! Check if the input file exists and open it.

       inquire ( file=filenamb(1:fnlen2), exist=lfnamb )
       if (.not.lfnamb .and. filenamb(fnlen2-4:fnlen2) .ne. '.fits') then
         filenamb = filenamb(1:fnlen2)//'.fits'
         fnlen2 = lastc ( filenamb )
         inquire ( file=filenamb(1:fnlen2), exist=lfnamb )
       endif
                             
       isp1 = lasts ( filenama ) + 1
       isp2 = lasts ( filenamb ) + 1

       if (arg2(1:ia2) .eq. 'q') then
         tdustx = 15.0d0
         arg2 = 'q15'
         ia2 = lastc ( arg2 )
       endif

       outname = filenama(isp1:fnlen1-5)//'.'//arg2(1:ia2)//'.'//filenamb(isp2:fnlen2)
       if (arg2(1:2) .eq. '-o') outname = arg3
       if (arg3(1:2) .eq. '-o') outname = arg4
       if (arg4(1:2) .eq. '-o') outname = arg5
       ion = lastc ( outname )
       if (outname(ion-4:ion) .ne. '.fits') then
         outname = outname(1:ion)//'.fits'
         ion = lastc ( outname )
       endif

       if (arg2(1:1) .eq. 'd' .or. arg2(1:1) .eq. 't' .or. arg2(1:1) .eq. 'i') then
         if (cverbose .eq. '-verb2') write (iotty,'(/a,f7.3,a)') ' OPERATE: '//filenama(isp1:fnlen1-5)//' '//arg2(1:ia2)
     &                                                                       //' '//filenamb(isp2:fnlen2-5)
     &                                                                       //' '//'-> using kappa0 ='
     &                                                                       , opacity0, ' cm^2/g at 300 um'
       endif

       if (.not.lfnama) then
         write (iotty,'(/a)') '   OPERATE: ERROR: File '''//filenama(1:fnlen1-5)//''' not found.'
         stop 1
       endif
       if (.not.lfnamb) then
         write (iotty,'(/a)') '   OPERATE: ERROR: File '''//filenamb(1:fnlen2-5)//''' not found.'
         stop 2
       endif

! Determine numbers of pixels in FITS images. 
! Signal GETFITSHEAD to ignore bad values of CDELT1, CDELT2 in the second image, if transferring header.

       if (arg2(1:ia2) .eq. 'h') then
         dx2 = -1.0d100
         dy2 = -1.0d100
       endif

       call getfitshead ( filenama(1:fnlen1), nx1, ny1, dx1, dy1, bunit )
       call getfitshead ( filenamb(1:fnlen2), nx2, ny2, dx2, dy2, bunit )

       if ((nx1 .ne. nx2 .or. ny1 .ne. ny2) .and. arg2(1:ia2) .ne. 'h') then
         write (iotty, '(/a)') '   OPERATE: ERROR: Images have different numbers of pixels:'
         write (iotty, '(a )') '                   '//filenama(1:fnlen1)
         write (iotty, '(a )') '                   '//filenamb(1:fnlen2)
         stop 3
       endif

       allocate ( funa(nx1,ny1), tdust(nt), bbodyx(nt), stat=irc )
       if (irc .ne. 0) then
         write (iotty,'(/a)') '   OPERATE: ERROR: Trouble allocating memory (10).'
         stop 10
       endif
       allocate ( funb(nx2,ny2), stat=irc )
       if (irc .ne. 0) then
         write (iotty,'(/a)') '   OPERATE: ERROR: Trouble allocating memory (20).'
         stop 20
       endif
       
       write (cx,'(i7)') nx1
       write (cy,'(i7)') ny1
       ncx = int ( log10 ( dble ( nx1 ) ) ) + 1
       ncy = int ( log10 ( dble ( ny1 ) ) ) + 1

       if (cverbose .eq. '-verb2') write (iotty,'(/a)') '   Reading ('//cx(7-ncx+1:7)//' x '//cy(7-ncy+1:7)//') '''
     &                                                //filenama(1:fnlen1)//''''

       call rfits ( nx1, ny1, bunit1, ctype11, ctype12, crpix11, crpix12, crval11, crval12, funa, dx1, dy1, object1, ra1, dec1
     &            , filenama(1:fnlen1), iotty, 0, creator, beam1, funmin, funmax, blank1, crota11, crota12, cd111, cd112, cd121
     &            , cd122, equinox1, bzero1, bscale1, wave1, datamin1, datamax1, history1, cverbose )

       if (cverbose .eq. '-verb2') write (iotty,'(a,1x,2(1pe14.7))') '   Minmax values in the image:', funmin, funmax

       write (cx,'(i7)') nx2
       write (cy,'(i7)') ny2
       ncx = int ( log10 ( dble ( nx2 ) ) ) + 1
       ncy = int ( log10 ( dble ( ny2 ) ) ) + 1

       if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Reading ('//cx(7-ncx+1:7)//' x '//cy(7-ncy+1:7)//') '''
     &                                               //filenamb(1:fnlen2)//''''

       call rfits ( nx2, ny2, bunit2, ctype21, ctype22, crpix21, crpix22, crval21, crval22, funb, dx2, dy2, object2, ra2, dec2
     &            , filenamb(1:fnlen2), iotty, 0, creator, beam2, funmin, funmax, blank2, crota21, crota22, cd211, cd212, cd221
     &            , cd222, equinox2, bzero2, bscale2, wave2, datamin2, datamax2, history2, cverbose )

       if (cverbose .eq. '-verb2') write (iotty,'(a,1x,2(1pe14.7))') '   Minmax values in the image:', funmin, funmax

       if (beam1 .eq. beam2) then
         beam = beam1
       else
         beam = 0.0d0
       endif
       if (arg2(1:ia2) .ne. 'h') then
         nx = nx1
         ny = ny1
       else
         nx = nx2
         ny = ny2
       endif
       bunit = bunit1
       ctype1 = ctype11
       ctype2 = ctype12
       dx = dx1
       dy = dy1
       crpix1 = crpix11
       crpix2 = crpix12
       crval1 = crval11
       crval2 = crval12
       blank = blank1
       crota1 = crota11
       crota2 = crota12
       cd11 = cd111
       cd12 = cd112
       cd21 = cd121
       cd22 = cd122
       equinox = equinox1
       object = object1
       ra = ra1
       dec = dec1
       bzero = bzero1
       bscale = bscale1
       history = history1
       wave = wave1
       if (wave1 .lt. almostzero) then
         wave = wave2
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. '+') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Adding images: image1 + image2'
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
               read (cnan,*) funb(i,j)
             else
               funb(i,j) = funa(i,j) + funb(i,j)
             endif
           enddo
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. '-') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Subtracting images: image1 - image2' 
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
               read (cnan,*) funb(i,j)
             else
               funb(i,j) = funa(i,j) - funb(i,j)
             endif
           enddo
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. 'x') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Multiplying images: image1 * image2'
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
               read (cnan,*) funb(i,j)
             else
               funb(i,j) = funa(i,j) * funb(i,j)
             endif
           enddo
         enddo
         ind = index ( bunit2, '(H2/cm^2)/(MJy/sr)' )
         if (ind .gt. 0) bunit = 'H2/cm^2'
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. '%') then 
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Dividing images: image1 / image2'
         do j=1,ny
           do i=1,nx
             if (.not.isnan ( funa(i,j) ) .and. .not.isnan ( funb(i,j) ) .and. abs ( funb(i,j) ) .ge. almostzero) then
               funb(i,j) = funa(i,j) / funb(i,j)
             else if (isnan ( funa(i,j) ) .or. abs ( funb(i,j) ) .eq. 0.0d0) then
               read (cnan,*) funb(i,j)
             endif
           enddo
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. '~') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Relative differencing images: (image1 - image2) / |image2|'
         do j=1,ny
           do i=1,nx
             if (abs ( funb(i,j) ) .ge. almostzero .and. abs ( funa(i,j) / (funb(i,j) + 1.0d-30) ) .le. 1.0d10) then
               funb(i,j) = (funa(i,j) - funb(i,j)) / abs ( funb(i,j) )
             else
               if (abs ( funb(i,j) ) .lt. almostzero) read (cnan,*) funb(i,j)
               if (isnan ( funa(i,j) )) then 
                 funb(i,j) = 0.0d0
               endif
             endif
           enddo
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. 's') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Squared differencing: (image1 - image2)^2'
         do j=1,ny
           do i=1,nx
             if (abs ( funb(i,j) ) .ge. almostzero) then
               funb(i,j) = (funa(i,j) - funb(i,j))**2
             else
               funb(i,j) = 0.0d0
             endif
           enddo
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. 'o') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Zeroing image1 where image1 < image2'
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
               read (cnan,*) funb(i,j)
             else
               if(funa(i,j) .gt. funb(i,j)) then
                 funb(i,j) = funa(i,j)
               else  
                 funb(i,j) = 0.0d0
               endif
             endif
           enddo
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. 'l') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Minimizing images: min ( image1, image2 )'
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
               read (cnan,*) funb(i,j)
             else
               funb(i,j) = min ( funa(i,j), funb(i,j) )
             endif
           enddo
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. 'g') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Maximizing images: max ( image1, image2 )'
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
               read (cnan,*) funb(i,j)
             else
               funb(i,j) = max ( funa(i,j), funb(i,j) )
             endif
           enddo
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. 'r') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Implanting non-zero pixels: image1 -> image2'
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
               read (cnan,*) funb(i,j)
             else
               if(funa(i,j) .gt. almostzero) then
                 funb(i,j) = funa(i,j)
               endif
             endif
           enddo
         enddo
       endif 
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. 'e') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Extending masks of image1 adding non-zero pixels of image2'
         allocate ( func(nx1,ny1), stat=irc )
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   OPERATE: ERROR: Trouble allocating memory (30).'
           stop 30
         endif
         do j=1,ny
           do i=1,nx
             func(i,j) = funa(i,j)
           enddo
         enddo
         itermax = 500
         do iter=1,itermax
           nexpand = 0
           do j=1,ny
             jm1 = max (j - 1, 1)
             jp1 = min (j + 1, ny)
             do i=1,nx
               im1 = max (i - 1, 1)
               ip1 = min (i + 1, nx)
               if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
                 read (cnan,*) funb(i,j)
               else
                 if(funb(i,j) .gt. almostzero) then
                   if(func(i,j) .lt. almostzero) then
                     if (func(im1,j) .gt. almostzero) then
                       funa(i,j) = funa(im1,j)
                       nexpand = nexpand + 1
                     else if (func(ip1,j) .gt. almostzero) then
                       funa(i,j) = funa(ip1,j)
                       nexpand = nexpand + 1
                     else if (func(i,jm1) .gt. almostzero) then
                       funa(i,j) = funa(i,jm1)
                       nexpand = nexpand + 1
                     else if (func(i,jp1) .gt. almostzero) then
                       funa(i,j) = funa(i,jp1)
                       nexpand = nexpand + 1
                     endif
                   endif
                 endif
               endif
             enddo
           enddo
           do j=1,ny
             do i=1,nx
               func(i,j) = funa(i,j)
             enddo
           enddo
           if (cverbose .eq. '-verb2') then
             inexp = int ( log10 ( dble ( max ( nexpand, 1 ) ) ) ) + 1
             write (cnexp,'(i7)') nexpand
             inite = int ( log10 ( dble ( max ( iter, 1 ) ) ) ) + 1
             write (cnite,'(i7)') iter
             write (iotty,'(a)') '   Extended '//cnexp(7-inexp+1:7)//' pixels at iteration '//cnite(7-inite+1:7)
           endif
           if (iter .eq. itermax) then
             write (iotty,'(/a)') '   OPERATE: ERROR: No convergence in extending mask'
             stop 99
           endif
           if (nexpand .eq. 0) exit
         enddo
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j)
           enddo
         enddo
         deallocate ( func )
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. 'p') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Expanding areas of image1 > image2 until image1 = image2'
         allocate ( func(nx1,ny1), stat=irc )
         if (irc .ne. 0) then
           write (iotty,'(/a)') '   OPERATE: ERROR: Trouble allocating memory (30).'
           stop 30
         endif
         do j=1,ny
           do i=1,nx
             func(i,j) = funa(i,j)
           enddo
         enddo
         itermax = 500
         do iter=1,itermax
           nexpand = 0
           do j=1,ny
             jm1 = max (j - 1, 1)
             jp1 = min (j + 1, ny)
             do i=1,nx
               im1 = max (i - 1, 1)
               ip1 = min (i + 1, nx)
               if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
                 read (cnan,*) funb(i,j)
               else
                 if(abs ( func(i,j) ) .lt. almostzero) then
                   if (abs ( func(im1,j) ) .gt. almostzero .and. funa(im1,j) .gt. funb(im1,j)) then
                     funa(i,j) = funa(im1,j)
                     nexpand = nexpand + 1
                   else if (abs ( func(ip1,j) ) .gt. almostzero .and. funa(ip1,j) .gt. funb(ip1,j)) then
                     funa(i,j) = funa(ip1,j)
                     nexpand = nexpand + 1
                   else if (abs ( func(i,jm1) ) .gt. almostzero .and. funa(i,jm1) .gt. funb(i,jm1)) then
                     funa(i,j) = funa(i,jm1)
                     nexpand = nexpand + 1
                   else if (abs ( func(i,jp1) ) .gt. almostzero .and. funa(i,jp1) .gt. funb(i,jp1)) then
                     funa(i,j) = funa(i,jp1)
                     nexpand = nexpand + 1
                   else if (abs ( func(im1,jp1) ) .gt. almostzero .and. funa(im1,jp1) .gt. funb(im1,jp1)) then
                     funa(i,j) = funa(im1,jp1)
                     nexpand = nexpand + 1
                   else if (abs ( func(ip1,jp1) ) .gt. almostzero .and. funa(ip1,jp1) .gt. funb(ip1,jp1)) then
                     funa(i,j) = funa(ip1,jp1)
                     nexpand = nexpand + 1
                   else if (abs ( func(ip1,jm1) ) .gt. almostzero .and. funa(ip1,jm1) .gt. funb(ip1,jm1)) then
                     funa(i,j) = funa(ip1,jm1)
                     nexpand = nexpand + 1
                   else if (abs ( func(im1,jm1) ) .gt. almostzero .and. funa(im1,jm1) .gt. funb(im1,jm1)) then
                     funa(i,j) = funa(im1,jm1)
                     nexpand = nexpand + 1
                   endif
                 endif
               endif
             enddo
           enddo
           do j=1,ny
             do i=1,nx
               func(i,j) = funa(i,j)
             enddo
           enddo
           if (cverbose .eq. '-verb2') then
             inexp = int ( log10 ( dble ( max ( nexpand, 1 ) ) ) ) + 1
             write (cnexp,'(i7)') nexpand
             inite = int ( log10 ( dble ( max ( iter, 1 ) ) ) ) + 1
             write (cnite,'(i7)') iter
             write (iotty,'(a)') '   Expanded '//cnexp(7-inexp+1:7)//' pixels at iteration '//cnite(7-inite+1:7)
           endif
           if (iter .eq. itermax) then
             write (iotty,'(/a)') '   OPERATE: ERROR: No convergence in expanding mask'
             stop 99
           endif
           if (nexpand .eq. 0) exit
         enddo
         do j=1,ny
           do i=1,nx
             funb(i,j) = funa(i,j)
           enddo
         enddo
         deallocate ( func )
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:1) .eq. 'd') then
         if (cverbose .eq. '-verb2') ! NOTE: image1 is intensities I_lambda (MJy/sr) and image2 is temperatures (K)
     &      write (iotty,'(a)') '   Computing surface densities: surfdens = intens1 / (Bnu(tempers2) * kappa * eta * mu * mH)'
         if (ia2 .gt. 1) then
           read (arg2(2:ia2),*) wave
         else
           wave = wave1
           if (wave .lt. almostzero) then
             write (iotty,'(/a)') ' OPERATE: ERROR: Keyword WAVE wrong or missing in '''//filenama(1:fnlen1)//'''.'
             stop 4
           endif
           if (cverbose .eq. '-verb2') then
             write (iotty,'(a,i4,a)') ' OPERATE: WAVELENGTH =', nint ( wave ), ' µm ASSUMED FOR THE IMAGE ~> IS IT CORRECT?'
             write (iotty,'(a)'     ) ' OPERATE: IF NOT, CORRECT "WAVE" IN '''//filenama(1:fnlen1)//''''
           endif
         endif

         if (cverbose .eq. '-verb2') then
           write (iotty,'(a,f7.3,a)') '   NOTE: Using reference opacity: kappa0 =', opacity0, ' cm^2/g at 300 um'
         endif

! Inu = Bnu(Td) * opacity * dust2gas * mu * mH * NH2  ~>  NH2 = Inu / (Bnu(Td) * opacity * dus2gas * mu * mH)
! Inu / Bnu(Td) = image1 / Bnu(image2) = opacity * dus2gas * mu * mH * NH2

         frequency = speedolight * 1.0d4 / wave
         opacity = opacity0 * (frequency / frequency0)**beta
         factor = 1.0d-17 / (opacity * dust2gas * muH2 * amu)
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) ) .or. funb(i,j) .lt. almostzero) then
!!               read (cnan,*) funb(i,j)
               funb(i,j) = 0.0d0
             else
               funb(i,j) = max ( min ( funb(i,j), 100.0d0 ), 3.0d0 )
               call planck ( funb(i,j), 1.0d4 / wave, bbody )
               if (isnan ( bbody ) .or. abs ( bbody ) .le. almostzero) then
                 read (cnan,*) funb(i,j)
               else
                 funb(i,j) = factor * funa(i,j) / bbody
               endif
             endif
           enddo
         enddo
         bunit = 'H2/cm^2'
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:1) .eq. 't') then
         if (cverbose .eq. '-verb2') ! NOTE: image1 is intensities I_lambda (MJy/sr) and image2 is surface densities NH2 (1/cm^2)
     &      write (iotty,'(a)') '   Computing temperatures: Invert { intens1 / (surfdens2 * kappa * eta * mu * mH) }'
         if (ia2 .gt. 1) then
           read (arg2(2:ia2),*) wave
         else
           wave = wave1
           if (wave .lt. almostzero) then
             write (iotty,'(/a)') ' OPERATE: ERROR: Keyword WAVE wrong or missing in '''//filenama(1:fnlen1)//'''.'
             stop 4
           endif
           if (cverbose .eq. '-verb2') then
             write (iotty,'(a,i4,a)') ' OPERATE: WAVELENGTH =', nint ( wave ), ' µm ASSUMED FOR THE IMAGE ~> IS IT CORRECT?'
             write (iotty,'(a)'     ) ' OPERATE: IF NOT, CORRECT "WAVE" IN '''//filenama(1:fnlen1)//''''
           endif
         endif

         if (cverbose .eq. '-verb2') then
           write (iotty,'(a,f7.3,a)') '   NOTE: Using reference opacity: kappa0 =', opacity0, ' cm^2/g at 300 um'
         endif

! Inu = Bnu(Td) * opacity * dust2gas * mu * mH * NH2  ~>  NH2 = Inu / (Bnu(Td) * opacity * dus2gas * mu * mH)
! Inu / Bnu(Td) = image1 / Bnu(image2) = opacity * dus2gas * mu * mH * NH2

         frequency = speedolight * 1.0d4 / wave
         opacity = opacity0 * (frequency / frequency0)**beta
         factor = 1.0d-17 / (opacity * dust2gas * muH2 * amu)
         tem1 = log ( 3.0d0 )
         tem2 = log ( 2000.0d0 )
         dtem = (tem2 - tem1) / dble ( nt - 1 )
         call planck ( exp ( tem1 ), 1.0d4 / wave, bbody1 )
         call planck ( exp ( tem2 ), 1.0d4 / wave, bbody2 )

         bbmax = -99.0
         bbmin = 99.0
         do l=1,nt
           tdust(l) = tem1 + dtem * dble ( l - 1 )
           call planck ( exp ( tdust(l) ), 1.0d4 / wave, bbodyx(l) )
           bbodyx(l) = log ( bbodyx(l))
           bbmin = min ( bbmin, bbodyx(l) )
           bbmax = max ( bbmax, bbodyx(l) )
         enddo
                 
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
!!               read (cnan,*) funb(i,j)
               funb(i,j) = 0.0d0
             else
               if (funa(i,j) .le. almostzero .or. funb(i,j) .le. almostzero) then
!!                 read (cnan,*) funb(i,j)
                 funb(i,j) = 0.0d0
               else
                 blackbody = factor * funa(i,j) / max ( funb(i,j), almostzero )
                 blackbody = log ( min ( max ( blackbody, bbody1 ), bbody2 ) )

                 do l=2,nt
                   if (blackbody .ge. bbodyx(l-1) .and. blackbody .le. bbodyx(l)) then
                     tdustx = tdust(l-1) + (blackbody - bbodyx(l-1)) * (tdust(l) - tdust(l-1)) / (bbodyx(l) - bbodyx(l-1))
                     goto 333
                   endif
                 enddo
                 
                 write (iotty,'(/a,2i4,10(1pe10.2))') '   OPERATE: ERROR: Intensity interval not found:', 
     &                  i,j, funa(i,j), funb(i,j), exp ( blackbody ), factor * funa(i,j) / max ( funb(i,j), almostzero )
                 stop 33
 333             continue
                 tdustx = exp ( tdustx )
                 if (isnan ( tdustx ) .or. tdustx .le. almostzero) then
                   read (cnan,*) funb(i,j)
                 else
                   funb(i,j) = tdustx
                 endif
               endif
             endif
           enddo
         enddo
         bunit = 'Kelvin'
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:1) .eq. 'i') then
         if (cverbose .eq. '-verb2') ! NOTE: image2 is temperatures (K) and image1 is surfdensities (1/cm^2)
     &      write (iotty,'(a)') '   Computing intensities = surfdens1 * Bnu(tempers2) * kappa * eta * mu * mH'
         if (ia2 .gt. 1) then
           read (arg2(2:ia2),*) wave
         else
           wave = wave1
           if (wave .lt. almostzero) then
             write (iotty,'(/a)') ' OPERATE: ERROR: Keyword WAVE wrong or missing in '''//filenama(1:fnlen1)//'''.'
             stop 4
           endif
           if (cverbose .eq. '-verb2') then
             write (iotty,'(a,i4,a)') ' OPERATE: WAVELENGTH =', nint ( wave ), ' µm ASSUMED FOR THE IMAGE ~> IS IT CORRECT?'
             write (iotty,'(a)'     ) ' OPERATE: IF NOT, CORRECT "WAVE" IN '''//filenama(1:fnlen1)//''''
           endif
         endif

         if (cverbose .eq. '-verb2') then
           write (iotty,'(a,f7.3,a)') '   NOTE: Using reference opacity: kappa0 =', opacity0, ' cm^2/g at 300 um'
         endif

! Inu = Bnu(Td) * opacity * dust2gas * mu * mH * NH2  ~>  NH2 = Inu / (Bnu(Td) * opacity * dus2gas * mu * mH)
! Inu / Bnu(Td) = image1 / Bnu(image2) = opacity * dus2gas * mu * mH * NH2

         frequency = speedolight * 1.0d4 / wave
         opacity = opacity0 * (frequency / frequency0)**beta
         factor = 1.0d-17 / (opacity * dust2gas * muH2 * amu)
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) ) .or. funb(i,j) .lt. almostzero) then
!!               read (cnan,*) funb(i,j)
               funb(i,j) = 0.0d0
             else
               funb(i,j) = max ( min ( funb(i,j), 100.0d0 ), 3.0d0 )
               call planck ( funb(i,j), 1.0d4 / wave, bbody )
               if (isnan ( bbody ) .or. abs ( bbody ) .le. almostzero) then
                 read (cnan,*) funb(i,j)
               else
                 funb(i,j) = funa(i,j) * bbody / factor 
               endif
             endif
           enddo
         enddo
         bunit = 'MJy/sr'
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:1) .eq. 'q') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Equalizing temperatures: image1 * Bnu(Td) / Bnu(image2)'
         read (arg2(2:ia2),*,err=555) tdustx
         goto 444
 555     continue         
           write (iotty,'(/a)') '   OPERATE: ERROR: Trouble reading Tdust from operand: '//arg2(1:ia2)
           stop 5
 444     continue
         if (wave1 .lt. almostzero) then
           write (iotty,'(/a)') '   OPERATE: ERROR: Keyword WAVE wrong or missing in '''//filenama(1:fnlen1)//'''.'
           stop 4
         endif
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   The image will correspond to a constant temperature of '
     &                                                 //arg2(2:ia2)//'K'
         call planck ( tdustx, 1.0d4 / wave1, bbody1 )
         do j=1,ny
           do i=1,nx
             if (isnan ( funa(i,j) ) .or. isnan ( funb(i,j) )) then
               read (cnan,*) funb(i,j)
             else
               call planck ( funb(i,j), 1.0d4 / wave1, bbody2 )
               if (isnan ( bbody2 ) .or. abs ( bbody2 ) .le. almostzero) then
                 read (cnan,*) funb(i,j)
               else
                 funb(i,j) = funa(i,j) * bbody1 / bbody2
               endif
             endif
           enddo
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
       if (arg2(1:ia2) .eq. 'h') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Replacing header: image1 -> image2'
         
! Preserve some header keywords of the second file (being corrected), they may have some meaningful (important) information.
         
         if (wave2 .gt. almostzero) wave = wave2
         if (object2 .ne. ' ') object = object2
         if (equinox2 .gt. almostzero) equinox = equinox2
         if (history2 .ne. ' ') history = history2
       endif
!__________________________________________________________________________________________________________________________________
!
! Get the current date and time.

       call when ( lunix, ctime, cdate, ndate, 4 )

       if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Writing output image '''//outname(1:ion)//''''

       call wfits ( cfitsversion, nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, funb, dx, dy, object, crval1
     &            , crval2, outname(1:ion), cdate, ctime, creator, beam, blank, crota1, crota2, cd11, cd12, cd21, cd22
     &            , equinox, bzero, bscale, wave, datamin, datamax, history )

       deallocate ( funa, funb, tdust, bbodyx )
       
       if (cverbose .eq. '-verb2') write (iotty,'(/a)') ' Done.'

!!       stop !<-- commented out because it would lead to run-time messages about denormalized values, when using gfortran.
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
       integer       i, j, nx, ny, stat, blksize, group, naxis, pcount, gcount, bitpix, unit, iotty, iolog, maxdim, rw
     &             , nfound, blank, stat05, stat06, stat07, stat08, stat09, stat10, stat11, stat12, stat13, stat14, stat15, stat16
     &             , stat17, stat18, stat19, stat20, stat21, stat22, stat23, stat24, stat25, stat26, stat27, stat28, stat29, stat30
     &             , stat31, stat32, stat0
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
         creator = 'OPERATE'
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
       wave = 0.0d0
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

       call ftpkys ( unit, 'CREATOR', creator,'Alexander Menshchikov, DAp IRFU CEA Saclay', status )
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