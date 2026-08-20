!__________________________________________________________________________________________________________________________________
!
! Various tools that can be linked to other programs.
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!__________________________________________________________________________________________________________________________________
!                 
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine osystem 
       
     &            ( lsystem )
!__________________________________________________________________________________________________________________________________
!
! This subroutine determines the operating system it runs on - whether it's UNIX-like or OS/2 (DOS)-like.
!__________________________________________________________________________________________________________________________________
!
       implicit none
       logical  lsystem
!__________________________________________________________________________________________________________________________________
!
       lsystem = .false.
       open ( 666, file='/dev/null', status='unknown', err=66 )
       lsystem = .true.
   66  close ( 666 )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine round1d ( vect, imax, isigfig )
!__________________________________________________________________________________________________________________________________
!
!  written by: David Clarke; modified by Alexander Men'shchikov
!
!  PURPOSE: This subroutine sets the number of significant figures of the elements in the input vector to the specified number.
!
!  INPUT VARIABLES:
!    vect        vector to be "rounded off"
!    imax        dimension of "vect"
!    isigfig     number of significant figures on output
!
!  OUTPUT VARIABLES:
!    vect        rounded off input vector
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       i, imax, isigfig, i1, i2
       real*8        q1, q2, q3, tiny
       real*8        vect(*)
       parameter     ( tiny = 1.0d-99 )
!__________________________________________________________________________________________________________________________________
!
       if (isigfig .le. 0) return
       q1 = 0.0d0
       do i=1,imax
         q1 = dmax1 ( q1, dabs(vect(i)) )
       enddo
       if (q1 .le. tiny) return
       q1 = dlog10 ( q1 )
       q2 = sign ( 0.5d0, q1 )
       i1 = isigfig - idint ( q1 + q2 + 0.5d0 )
       i2 = min0 ( -1, isigfig - 13 )
       q1 = 10.0d0**i1
       q2 = 10.0d0**i2
       do i=1,imax
         q3 = sign  ( q2, vect(i) )
         vect(i) = dnint ( vect(i) * q1 + q3 ) / q1
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine expandit 
       
     &            ( ctype, nx, ny, nxb, nyb, nbw, nxmx, nymx, nxmxb, nymxb, image, imagewide )
!__________________________________________________________________________________________________________________________________
!
! Expand image into borders of width 'nbw'.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) ctype
       integer       i, j, nx, ny, nxb, nyb, nbw, nxmxb, nymxb, nxmx, nymx, k, ni, nj
       real*8        image(nxmx,nymx), imagewide(nxmxb,nymxb), value
!__________________________________________________________________________________________________________________________________
!
! Initialise a larger image for no border effects in the area we are interested in.

       do j=1,nyb
         do i=1,nxb
           imagewide(i,j) = 0.0d0
         enddo
       enddo

! Shift the image to the center of the new image.

       do j=1,ny
         do i=1,nx
           imagewide(i+nbw,j+nbw) = image(i,j)
         enddo
       enddo

       if (ctype .eq. 'periodic') then
              
! Assume periodicity to extend the image in x- and y- directions.

         do i=1,nbw
           do j=1+nbw,ny+nbw
             imagewide(i,j) = imagewide(i+nx,j)
           enddo
         enddo
         do i=1+nbw+nx,nxb
           do j=1+nbw,ny+nbw
             imagewide(i,j) = imagewide(i-nx,j)
           enddo
         enddo
         do j=1,nbw
           do i=1+nbw,nx+nbw
             imagewide(i,j) = imagewide(i,j+ny)
           enddo
         enddo
         do j=1+nbw+ny,nyb
           do i=1+nbw,nx+nbw
             imagewide(i,j) = imagewide(i,j-ny)
           enddo
         enddo

! In a similar fashion, take care of the corners.

         do j=1,nbw
           do i=1,nbw
             imagewide(i,j) = imagewide(i+nx,j+ny)
           enddo
         enddo
         do j=1+nbw+ny,nyb
           do i=1+nbw+nx,nxb
             imagewide(i,j) = imagewide(i-nx,j-ny)
           enddo
         enddo
         do i=1,nbw
           do j=1+nbw+ny,nyb
             imagewide(i,j) = imagewide(i,j-ny)
           enddo
         enddo
         do j=1,nbw
           do i=1+nbw+nx,nxb
             imagewide(i,j) = imagewide(i-nx,j)
           enddo
         enddo
       else

! Assume "outflowing" boundary to extend the image in the x- and y- directions.

         do i=1,nbw
           do j=1+nbw,ny+nbw
             imagewide(i,j) = imagewide(nbw+1,j)
           enddo
         enddo
         do i=1+nbw+nx,nxb
           do j=1+nbw,ny+nbw
             imagewide(i,j) = imagewide(nbw+nx,j)
           enddo
         enddo
         do j=1,nbw
           do i=1+nbw,nx+nbw
             imagewide(i,j) = imagewide(i,nbw+1)
           enddo
         enddo
         do j=1+nbw+ny,nyb
           do i=1+nbw,nx+nbw
             imagewide(i,j) = imagewide(i,nbw+ny)
           enddo
         enddo

! In a similar fashion, take care of the corners. 
                                 
         ni = 10
         nj = 10
         value = 0.0d0
         k = 0
         do j=1+nbw,1+nbw+nj
           do i=1+nbw,1+nbw+ni
             k = k + 1
             value = value + imagewide(i,j)
           enddo
         enddo
         value = value / dble ( k )
         do j=1,nbw
           do i=1,nbw
             imagewide(i,j) = value
           enddo
         enddo

         value = 0.0d0
         k = 0
         do j=1+nbw+ny-nj,1+nbw+ny
           do i=1+nbw+nx-ni,1+nbw+nx
             k = k + 1
             value = value + imagewide(i,j)
           enddo
         enddo
         value = value / dble ( k )
         do j=1+nbw+ny,nyb
           do i=1+nbw+nx,nxb
             imagewide(i,j) = value
           enddo
         enddo

         value = 0.0d0
         k = 0
         do i=1+nbw,1+nbw+ni
           do j=1+nbw+ny-nj,1+nbw+ny
             k = k + 1
             value = value + imagewide(i,j)
           enddo
         enddo
         value = value / dble ( k )
         do i=1,nbw
           do j=1+nbw+ny,nyb
             imagewide(i,j) = value
           enddo
         enddo

         value = 0.0d0
         k = 0
         do j=1+nbw,1+nbw+nj
           do i=1+nbw+nx-ni,1+nbw+nx
             k = k + 1
             value = value + imagewide(i,j)
           enddo
         enddo
         value = value / dble ( k )
         do j=1,nbw
           do i=1+nbw+nx,nxb
             imagewide(i,j) = value
           enddo
         enddo                                                                                                                            
       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine shrinkit 
       
     &            ( nx, ny, nxs, nys, mm, nxmx, nymx, image, imagesmall )
!__________________________________________________________________________________________________________________________________
!
! Shrink image into itself by shifting quadrants and averaging.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       i, j, nx, ny, nxs, nys, mm, nxh, nyh, nxsh, nysh, nxmx, nymx
       real*8        image(nxmx,nymx), imagesmall(nxmx,nymx)
!__________________________________________________________________________________________________________________________________
!
       nxh = nx / 2
       nyh = ny / 2
       nxsh = nxs / 2
       nysh = nys / 2
       
       do j=1,nys
         do i=1,nxs
           imagesmall(i,j) = 0.0d0
         enddo
       enddo

! Shift the image to the center of the new image.

       do j=1,nyh
         do i=1,nxh
           imagesmall(i,j) = imagesmall(i,j) + image(i,j)
         enddo
       enddo
       do j=1,nyh
         do i=nxh+1,nx
           if (i-2*mm .ge. 1) then
             imagesmall(i-2*mm,j) = imagesmall(i-2*mm,j) + image(i,j)
           endif
         enddo
       enddo
       do j=nyh+1,ny
         do i=1,nxh
           if (j-2*mm .ge. 1) then
             imagesmall(i,j-2*mm) = imagesmall(i,j-2*mm) + image(i,j)
           endif
         enddo
       enddo
       do j=nyh+1,ny
         do i=nxh+1,nx
           if (i-2*mm .ge. 1 .and. j-2*mm .ge. 1) then
             imagesmall(i-2*mm,j-2*mm) = imagesmall(i-2*mm,j-2*mm) + image(i,j)
           endif
         enddo
       enddo

       do j=1,nys
         do i=max(nxsh-mm+1,1),min(nxsh+mm,nxs)
           imagesmall(i,j) = imagesmall(i,j) / 2.0d0
         enddo
       enddo
       do j=max(nysh-mm+1,1),min(nysh+mm,nys)
         do i=1,nxs
           imagesmall(i,j) = imagesmall(i,j) / 2.0d0
         enddo
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function inellipse

     &          ( x, y, dx, dy, x0, y0, amajor, bminor, thepa )
!__________________________________________________________________________________________________________________________________
!
! This subroutine determines whether a given point (x,y) is inside an ellipse centered on (x0,y0).
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       inellipse
       real*8        x, y, x0, y0, amajor, bminor, thepa, radx, rady, angle, dex, dey, sina, cosa, delx, dely, theta
     &             , pi, pio2, pi2, dx, dy, radx2, rady2, rad2, radell2, sinth2
       parameter   ( pi = 3.14159265358979d0, pio2 = pi / 2.0d0, pi2 = pi * 2.0d0 )
!__________________________________________________________________________________________________________________________________
!
       inellipse = .false.
       
       radx = amajor / 2.0d0  !<-- major radius of the current ellipse
       rady = bminor / 2.0d0  !<-- minor radius of the current ellipse
       angle = (90.0d0 - thepa) * (pi / 180.0d0)  !<-- angle of the ellipse's major axis in radians
         
       dex = (x - x0) * dx
       dey = (y - y0) * dy
       sina = sin ( angle )
       cosa = cos ( angle )
       delx = dex * cosa - dey * sina
       dely = dex * sina + dey * cosa
       
       if (dely .gt. 0.0d0) then
         if (delx .eq. 0.0d0) then
           theta = pio2
         else
           if (delx .gt. 0.0d0) then
             theta = atan ( dely / delx )
           else
             theta = pi - atan ( dely / abs ( delx ) )
           endif
         endif
       else
         if (delx .eq. 0.0d0) then
           theta = pi + pio2
         else
           if (delx .gt. 0.0d0) then
             theta = pi2 - atan ( abs ( dely / delx ) )
           else
             theta = pi + atan ( abs ( dely / delx ) )
           endif
         endif
       endif
       
       radx2 = radx**2
       rady2 = rady**2
       sinth2 = sin ( theta )**2
       radell2 = radx2 * rady2 / (rady2 + (radx2 - rady2) * sinth2 + 1.0d-30)
       rad2  = delx**2 + dely**2
       
       if (rad2 .le. radell2) inellipse = .true.

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine findmin
       
     &            ( xn, yn, xk, yk, nx, ny, image, imin, jmin, imagemin )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       nx, ny, imin, jmin, l, npt, ic, jc
       real*8        xn, yn, xk, yk, dx, dy, dbll, imagemin
       real*8        image(nx,ny)
!__________________________________________________________________________________________________________________________________
!
!!!       imin = nint ( xn )
!!!       jmin = nint ( yn )
!!!       imagemin = image(imin,jmin)

       npt = 3 * nint ( max ( abs ( xk - xn ), abs ( yk - yn ), 1.0d0 ) )
       dx = (xk - xn) / dble ( npt - 1 )
       dy = (yk - yn) / dble ( npt - 1 )
       imagemin = 1.0d+90
       imin = 0
       jmin = 0
       do l=1,npt
         dbll = dble ( l - 1 )
         ic = nint ( xn + dx * dbll )
         jc = nint ( yn + dy * dbll )
         if (image(ic,jc) .lt. imagemin) then
           imin = ic
           jmin = jc
           imagemin = image(imin,jmin)
         endif
       enddo
       if (imin .eq. 0 .or. jmin .eq. 0) then
         write (*,'(/a,4(1pe10.2),i5)') '   FINDMIN: ERROR: Wrong input?: ', xn, yn, xk, yk, npt
         stop 99     
       endif
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine skipcomm 
       
     &            ( iunit, filename, nheader )
!__________________________________________________________________________________________________________________________________
!
!                                     MCM3D, Alexander Men'shchikov, 12/07/2007
!
! Find data in the character line being read from file (by skipping comment lines) and count the number of header comment lines.
!
! INPUT PARAMETERS:
!
!    iunit    - file unit number
!    filename - file name
!
! OUTPUT PARAMETERS:
!
!    nheader  - number of header comment lines in the file
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       iunit, i, lastc, nheader, nhead, lfn, lenc
       character*(*) filename
       character*1   ccomm, first
       character*2   clasto
       character*80  clast
       external      lastc
!__________________________________________________________________________________________________________________________________
!
       nhead = 0
       first = '?'
       clasto = ''
       do
          nhead = nhead + 1
          read (iunit,'(a)', end=98, err=99) clast
          ccomm = clast(1:1)

!!          if (first .eq. '?' .and. ccomm .ne. ';' .and. ccomm .ne. ':' .and. ccomm .ne. '*' .and. ccomm .ne. '#' .and. 
!!     &        ccomm .ne. '!') exit

          if ((first .eq. ';' .and. ccomm .ne. ';') .or. (first .eq. ':' .and. ccomm .ne. ':') .or. 
     &        (first .eq. '*' .and. ccomm .ne. '*') .or. (first .eq. '#' .and. ccomm .ne. '#') .or. 
     &        (first .eq. '!' .and. ccomm .ne. '!' .and. ccomm .ne. '#')) exit

!!          if (first .ne. '?' .and. first .ne. ccomm) then
!!            write (*,'(/a)') '   SKIPCOMM: ERROR: Comment symbol differs within header: '''//first//''' /= '''//ccomm//'''.'
!!            stop 99
!!          endif
          if (first .eq. '?') first = ccomm
          clasto = clast(1:2)
       enddo

 98    backspace ( iunit )
       nheader = nhead - 1
 
       lenc = lastc ( clasto )
       if (lenc .gt. 1) then
         do
           nhead = nheader
           nheader = nheader - 1
           backspace ( iunit )
           read (iunit,'(a)', end=98, err=99) clast
           backspace ( iunit )
           lenc = lastc ( clast )
           if (lenc .eq. 1) exit
         enddo
       endif

       return

 99    lfn = lastc ( filename )
       write (*,'(/a,i3)') '   SKIPCOMM: ERROR: Trouble reading '''// filename(1:lfn) //''' at line:', i

       stop 99
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine linescolumns 
       
     &            ( iounit, filename, nheader, nlines )
!__________________________________________________________________________________________________________________________________
!
!                                     MCM3D, Alexander Men'shchikov, 20/08/2007
!
! Find actual number of data columns and data lines in the file, based on just the first data line. File is assumed to be
! positioned at the first data line, after all possible header comments. At the end of this routine, the file is positioned back
! at the first data line.
!
! INPUT PARAMETERS:
!
!    iounit   - file unit number
!    filename - file name
!    nheader  - number of header comment lines in the file
!
! OUTPUT PARAMETERS:
!
!    nlines   - number of data lines in the file
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) filename
!!       character*20000 inline
       integer       lastc, iounit, nheader, nlines, maxlines, i, nl, lfn, irc     !! ncolumns
       parameter   ( maxlines = 100000 )
       character*1, allocatable :: cc(:)
       external      lastc
!__________________________________________________________________________________________________________________________________
!
       allocate ( cc(maxlines), stat=irc )
       
       if (irc .ne. 0) then
         write (iounit,'(/a)') '   LINESCOLUMNS: ERROR: Trouble allocating memory (10).'
         stop 10
       endif

       do i=1,maxlines
         cc(i) = ' '
       enddo
                         
! Find number of lines.

       nl = 0
       do
         nl = nl + 1
         read (iounit,'(a)',end=91,err=98) cc(nl)
       enddo
  91   nlines = nl - 1

! Don't count as data lines those which are footer comments (at the end of the file).
                      
       nl = 0
       do i=nlines,1,-1
         if (cc(i) .ne. ';' .and. cc(i) .ne. ':' .and. cc(i) .ne. '#' .and. cc(i) .ne. '!' .and. cc(i) .ne. '*') exit
         nl = nl + 1
       enddo
       nlines = nlines - nl
       rewind ( iounit )
       do i=1,nheader
         read (iounit,'(a)')  !! inline
       enddo

       deallocate ( cc )       

       return

  98   lfn = lastc ( filename )
       write (*,'(/a)') '   LINESCOLUMNS: ERROR: Trouble reading ''' // filename(1:lfn) // '''.'

       deallocate ( cc )       

       stop
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine getnumpix
       
     &            ( fname, nx, ny )
!__________________________________________________________________________________________________________________________________
!
! Kept here for compatibility with GETSOURCES v1. Determine numbers of pixels in a FITS image.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*80  record
       character*(*) fname
       integer       unit, status, readwrite, nkeys, blocksize, nspace, lfn, lastc, i, nx, ny
       external      ftgiou, ftopen, ftghsp, printerr, ftgrec, lastc
!__________________________________________________________________________________________________________________________________
!
       status = 0
       readwrite = 0
       lfn = lastc ( fname )
       call ftgiou ( unit, status )
       call ftopen ( unit, fname(1:lfn), readwrite, blocksize, status )
       call ftghsp ( unit, nkeys, nspace, status )
       if (status .gt. 0) then
         call printerr ( status )
         stop 99
       endif
       nx = 0
       ny = 0
       do i=1,nkeys
         call ftgrec ( unit, i, record, status )
         if (index ( record, 'NAXIS1' ) .gt. 0) read (record(10:31),*) nx
         if (index ( record, 'NAXIS2' ) .gt. 0) read (record(10:31),*) ny
       enddo
       if (nx .eq. 0 .or. ny .eq. 0) then
         write (*,'(/a)') '   GETNUMPIX: ERROR: Trouble getting image dimensions in ''' // fname(1:lfn) // '''.'
         stop 99
       endif

       return
       end
 
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine getfitshead
       
     &            ( fname, nx, ny, dx, dy, bunit )
!__________________________________________________________________________________________________________________________________
!
! Obtain some header info from a FITS image.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       lgotdx, lgotdy, lnocheckdx, lnocheckdy
       character*80  record, bunit
       character*(*) fname
       integer       unit, status, readwrite, nkeys, blocksize, nspace, lfn, lastc, i, nx, ny
       real*8        dx, dy
       external      lastc, ftgiou, ftopen, ftghsp, printerr, ftgrec
!__________________________________________________________________________________________________________________________________
!
       status = 0
       readwrite = 0
       lfn = lastc ( fname )
       call ftgiou ( unit, status )
       call ftopen ( unit, fname(1:lfn), readwrite, blocksize, status )
       call ftghsp ( unit, nkeys, nspace, status )
       if (status .gt. 0) then
         call printerr ( status )
         stop 99
       endif
       nx = 0
       ny = 0
       lnocheckdx = dx .ge. 1.0d99
       lnocheckdy = dy .ge. 1.0d99
       dx = 0.0d0
       dy = 0.0d0
       lgotdx = .false.
       lgotdy = .false.
       do i=1,nkeys
         call ftgrec ( unit, i, record, status )
         if (index ( record, 'NAXIS1' ) .eq. 1) read (record(10:31),*) nx
         if (index ( record, 'NAXIS2' ) .eq. 1) read (record(10:31),*) ny
         if (index ( record, 'CDELT1' ) .eq. 1) then
           read (record(10:31),*) dx
           lgotdx = .true.
         endif
         if (index ( record, 'CDELT2' ) .eq. 1) then
           read (record(10:31),*) dy
           lgotdy = .true.
         endif
         if (index ( record, 'BUNIT' ) .eq. 1) read (record(10:31),*) bunit
       enddo
       if (lgotdx) dx = abs ( dx ) * 3600.0d0
       if (lgotdy) dy = abs ( dy ) * 3600.0d0
       if (nx .eq. 0 .or. ny .eq. 0) then
         write (*,'(/a)') '   GETFITSHEAD: ERROR: Trouble getting numbers of pixels in ''' // fname(1:lfn) // '''.'
         stop 99
       endif
       if ((.not.lnocheckdx .and. lgotdx .and. dx .lt. 1.0d-99) .or. (.not.lnocheckdy .and. lgotdy .and. dy .lt. 1.0d-99)) then
         write (*,'(/a)') '   GETFITSHEAD: ERROR: Trouble getting sizes of pixels in ''' // fname(1:lfn) // '''.'
         stop 99
       endif

       return
       end
 
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine ds9region

     &            ( nunit, dx, dy, nhead, nhcols, lhl, excat, headline, nextr, cnwave, obsbeam, ncolumnx, ncolumny, ncolumnxas
     &            , ncolumnyas, ncolumng, ncolumngm, ncolumnfm, ncolumna, ncolumnb, ncolumnt, ncolumnfp, ncolumnft, ncolumnr
     &            , ncolumnfa, ncolumnfb, sizefact, ncolumnfpe, ncolumnfte, thickness, sigmin, goodmin, snr1min, snr2min, elongmax
     &            , lfootprints, regcol, phy, fk5 )
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       logical       lfootprints

       character*2   sharp
       character*5   cattype
       character*6   cj
       character*7   cnextr
       character*11  cafwhm, cbfwhm, cfoota, cfootb, ctheta
       character*13  cxpixm, cypixm
       character*300 headreg(48)
       character*15000 cline
       character*(*) cnwave, regcol, excat, phy, fk5, headline(*)
     
       integer       j, k, ia, ib, il, ii, lastc, lhl, nhead, nextr, nxe, nye, nunit, ncolumnx, ncolumny, ncolumna, ncolumnb, ncx
     &             , ncolumnt, ncolumnfp, ncolumnft, ncolumnr, ncolumnfpe, ncolumnfte, no, m, nxco, nyco, nhcols, ncolumng
     &             , ncolumnfm, flagmon, is, ix, iy, ncolumnfa, ncolumnfb, ncolumngm, ic, id, ihc, it, ncolumnxas, ncolumnyas, ncj
     &             , icl
     
       real*8        dx, dy, goodmin, x0o, y0o, x0oas, y0oas, majdiam, mindiam, angle, radx, rady, obsbeam, thickness, sigmono
     &             , sigmin, fxpbest, fxperro, fxtbest, fxterro, snr1min, snr2min, sizefact, goodness, goodmono, fmajdiam, fmindiam
     &             , elongmax
!__________________________________________________________________________________________________________________________________
!
       if (.not.lfootprints) then
         if (phy .eq. 'phy') open ( 15, file=excat//'.ell.phy.reg', status='unknown' )
         if (fk5 .eq. 'fk5') open ( 17, file=excat//'.ell.fk5.reg', status='unknown' )
       else
         if (phy .eq. 'phy') open ( 16, file=excat//'.foo.phy.reg', status='unknown' )
         if (fk5 .eq. 'fk5') open ( 18, file=excat//'.foo.fk5.reg', status='unknown' )
       endif
  
       headreg( 1) = '# SYMBOLIC COLOR NAMES RECOGNIZED BY DS9 (https://www.tcl.tk/man/tcl/TkCmd/colors.html'
     &             //' Tcl8.6.11/Tk8.6.11 Documentation > Tk Commands > colors)'
       headreg( 2) = '# '
       headreg( 3) = '# AliceBlue AntiqueWhite AntiqueWhite1 AntiqueWhite2 AntiqueWhite3 AntiqueWhite4 agua aquamarine'
     &             //' aquamarine1 aquamarine2 aquamarine3 aquamarine4 azure azure1'
       headreg( 4) = '# azure3 azure4 beige bisque bisque1 bisque2 bisque3 bisque4 black BlanchedAlmond blue blue1 blue2'
     &             //' blue3 blue4 BlueViolet brown brown1 brown2 brown3 brown4 azure2 '
       headreg( 5) = '# burlywood burlywood1 burlywood2 burlywood3 burlywood4 CadetBlue CadetBlue1 CadetBlue2 CadetBlue3'
     &             //' CadetBlue4 chartreuse chartreuse1 chartreuse2 chartreuse3'
       headreg( 6) = '# chartreuse4 chocolate chocolate1 chocolate2 chocolate3 chocolate4 coral coral1 coral2 coral3'
     &             //' coral4 CornflowerBlue cornsilk cornsilk1 cornsilk2 cornsilk3'
       headreg( 7) = '# cornsilk4 crymson cyan cyan1 cyan2 cyan3 cyan4 dark blue dark cyan dark goldenrod DarkBlue'
     &             //' DarkCyan DarkGoldenrod DarkGoldenrod1 DarkGoldenrod2 DarkGoldenrod3'
       headreg( 8) = '# DarkGoldenrod4 DarkGray DarkGreen DarkKhaki DarkMagenta DarkOliveGreen DarkOliveGreen1'
     &             //' DarkOliveGreen2 DarkOliveGreen3 DarkOliveGreen4 DarkOrange DarkOrange1'
       headreg( 9) = '# DarkOrange2 DarkOrange3 DarkOrange4 DarkOrchid DarkOrchid1 DarkOrchid2 DarkOrchid3 DarkOrchid4'
     &             //' DarkRed DarkSalmon DarkSeaGreen DarkSeaGreen1 DarkSeaGreen2'
       headreg(10) = '# DarkSeaGreen3 DarkSeaGreen4 DarkSlateBlue DarkSlateGray DarkSlateGray1 DarkSlateGray2'
     &             //' DarkSlateGray3 DarkSlateGray4 DarkSlateGrey DarkTurquoise DarkViolet'
       headreg(11) = '# DeepPink DeepPink1 DeepPink2 DeepPink3 DeepPink4 DeepSkyBlue DeepSkyBlue1 DeepSkyBlue2'
     &             //' DeepSkyBlue3 DeepSkyBlue4 DimGray DodgerBlue DodgerBlue1 DodgerBlue2'
       headreg(12) = '# DodgerBlue3 DodgerBlue4 firebrick firebrick1 firebrick2 firebrick3 firebrick4 FloralWhite'
     &             //' ForestGreen fuchsia gainsboro GhostWhite gold gold1 gold2 gold3 gold4'
       headreg(13) = '# goldenrod goldenrod1 goldenrod2 goldenrod3 goldenrod4 gray gray0 gray1 gray2 gray3 gray4 gray5'
     &             //' gray6 gray7 gray8 gray9 gray10 gray11 gray12 gray13 gray14 gray15'
       headreg(14) = '# gray16 gray17 gray18 gray19 gray20 gray21 gray22 gray23 gray24 gray25 gray26 gray27 gray28'
     &             //' gray29 gray30 gray31 gray32 gray33 gray34 gray35 gray36 gray37 gray38'
       headreg(15) = '# gray39 gray40 gray41 gray42 gray43 gray44 gray45 gray46 gray47 gray48 gray49 gray50 gray51'
     &             //' gray52 gray53 gray54 gray55 gray56 gray57 gray58 gray59 gray60 gray61'
       headreg(16) = '# gray62 gray63 gray64 gray65 gray66 gray67 gray68 gray69 gray70 gray71 gray72 gray73 gray74'
     &             //' gray75 gray76 gray77 gray78 gray79 gray80 gray81 gray82 gray83 gray84'
       headreg(17) = '# gray85 gray86 gray87 gray88 gray89 gray90 gray91 gray92 gray93 gray94 gray95 gray96 gray97'
     &             //' gray98 gray99 gray100 green green1 green2 green3 green4 GreenYellow'
       headreg(18) = '# honeydew honeydew1 honeydew2 honeydew3 honeydew4 HotPink HotPink1 HotPink2 HotPink3 HotPink4'
     &             //' IndianRed IndianRed1 IndianRed2 IndianRed3 IndianRed4 indigo ivory'
       headreg(19) = '# ivory1 ivory2 ivory3 ivory4 khaki khaki1 khaki2 khaki3 khaki4 lavender LavenderBlush'
     &             //' LavenderBlush1 LavenderBlush2 LavenderBlush3 LavenderBlush4 LawnGreen'
       headreg(20) = '# LemonChiffon LemonChiffon1 LemonChiffon2 LemonChiffon3 LemonChiffon4 LightBlue LightBlue1'
     &             //' LightBlue2 LightBlue3 LightBlue4 LightCoral LightCyan LightCyan1'
       headreg(21) = '# LightCyan2 LightCyan3 LightCyan4 LightGoldenrod LightGoldenrod1 LightGoldenrod2 LightGoldenrod3'
     &             //' LightGoldenrod4 LightGoldenrodYellow LightGray LightGreen'
       headreg(22) = '# LightGrey LightPink LightPink1 LightPink2 LightPink3 LightPink4 LightSalmon LightSalmon1'
     &             //' LightSalmon2 LightSalmon3 LightSalmon4 LightSeaGreen LightSkyBlue'
       headreg(23) = '# LightSkyBlue1 LightSkyBlue2 LightSkyBlue3 LightSkyBlue4 LightSlateBlue LightSlateGray'
     &             //' LightSlateGrey LightSteelBlue LightSteelBlue1 LightSteelBlue2'
       headreg(24) = '# LightSteelBlue3 LightSteelBlue4 LightYellow LightYellow1 LightYellow2 LightYellow3 LightYellow4'
     &             //' lime LimeGreen linen magenta magenta1 magenta2 magenta3 magenta4'
       headreg(25) = '# maroon maroon1 maroon2 maroon3 maroon4 red MediumAquamarine MediumBlue MediumOrchid'
     &             //' MediumOrchid1 MediumOrchid2 MediumOrchid3 MediumOrchid4 MediumPurple'
       headreg(26) = '# MediumPurple1 MediumPurple2 MediumPurple3 MediumPurple4 MediumSeaGreen MediumSlateBlue'
     &             //' MediumSpringGreen MediumTurquoise MediumVioletRed MidnightBlue MintCream'
       headreg(27) = '# MistyRose MistyRose1 MistyRose2 MistyRose3 MistyRose4 moccasin NavajoWhite NavajoWhite1'
     &             //' NavajoWhite2 NavajoWhite3 NavajoWhite4 navy NavyBlueOldLace olive'
       headreg(28) = '# OliveDrab OliveDrab1 OliveDrab2 OliveDrab3 OliveDrab4 orange orange1 orange2 orange3 orange4'
     &             //' OrangeRed OrangeRed1 OrangeRed2 OrangeRed3 OrangeRed4 orchid orchid1'
       headreg(29) = '# orchid2 orchid3 orchid4 red PaleGoldenrod PaleGreen PaleGreen1 PaleGreen2 PaleGreen3 PaleGreen4'
     &             //' PaleTurquoise PaleTurquoise1 PaleTurquoise2 PaleTurquoise3'
       headreg(30) = '# PaleTurquoise4 PaleVioletRed PaleVioletRed1 PaleVioletRed2 PaleVioletRed3 PaleVioletRed4'
     &             //' PapayaWhip PeachPuff PeachPuff1 PeachPuff2 PeachPuff3 PeachPuff4 peru'
       headreg(31) = '# pink pink1 pink2 pink3 pink4 plum plum1 plum2 plum3 plum4 PowderBlue purple purple1 purple2'
     &             //' purple3 purple4 red red1 red2 red3 red4 RosyBrown RosyBrown1'
       headreg(32) = '# RosyBrown2 RosyBrown3 RosyBrown4 RoyalBlue RoyalBlue1 RoyalBlue2 RoyalBlue3 RoyalBlue4'
     &             //' SaddleBrown salmon salmon1 salmon2 salmon3 salmon4 SandyBrown SeaGreen'
       headreg(33) = '# SeaGreen1 SeaGreen2 SeaGreen3 SeaGreen4 seashell seashell1 seashell2 seashell3 seashell4 sienna'
     &             //' sienna1 sienna2 sienna3 sienna4 silver SkyBlue SkyBlue1 SkyBlue2'
       headreg(34) = '# SkyBlue3 SkyBlue4 SlateBlue SlateBlue1 SlateBlue2 SlateBlue3 SlateBlue4 SlateGray SlateGray1'
     &             //' SlateGray2 SlateGray3 SlateGray4 SlateGrey snow snow1 snow2 snow3'
       headreg(35) = '# snow4 SpringGreen SpringGreen1 SpringGreen2 SpringGreen3 SpringGreen4 steel blue SteelBlue'
     &             //' SteelBlue1 SteelBlue2 SteelBlue3 SteelBlue4 tan tan1 tan2 tan3 tan4'
       headreg(36) = '# teal thistle thistle1 thistle2 thistle3 thistle4 tomato tomato1 tomato2 tomato3 tomato4'
     &             //' turquoise turquoise1 turquoise2 turquoise3 turquoise4 violet VioletRed'
       headreg(37) = '# VioletRed1 VioletRed2 VioletRed3 VioletRed4 wheat wheat1 wheat2 wheat3 wheat4 white WhiteSmoke'
     &             //' yellow yellow1 yellow2 yellow3 yellow4 YellowGreen'
       headreg(38) = '#'
       headreg(39) = '# Region file format: DS9 version 4.1'
       headreg(40) = '#'
       headreg(41) = '# GETSF COLORS ARE DEFINED BY THEIR HEXADECIMAL VALUES:'
       headreg(42) = '#'
       headreg(43) = '# red #f0194b, green #3ce04b, yellow #fff019, blue #7383ff, orange #ff9241, purple #a14eff,'
     &             //' cyan #66ffff, magenta #ff52ff, lime #bcf60c, pink #ffb5b5,'
       headreg(44) = '# teal #00a0a0, lavender #e6beff, brown #ba7334, beige #f2f2ca, maroon #aa0000, mint #aaffc3,'
     &             //' olive #a5b500, apricot #ffd8b1, gray #909090, white #ffffff, black #000000'
       headreg(45) = '#'
       headreg(47) = '#'
       headreg(48) = 'global color=#000000 width=5 dashlist=8 3 font="helvetica 10 normal roman" select=1'
     &             //' highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
          
       do m=1,nhead
         ihc = lastc ( headline(m) )
         ia = index ( headline(m), '! TABULATED QUANTITIES' )
         if (ia .gt. 0) then
           do ii=1,37
             il = lastc ( headreg(ii) )
             if (.not.lfootprints) then
               if (phy .eq. 'phy') write (15,'(a)') headreg(ii)(1:il)
               if (fk5 .eq. 'fk5') write (17,'(a)') headreg(ii)(1:il)
             else
               if (phy .eq. 'phy') write (16,'(a)') headreg(ii)(1:il)
               if (fk5 .eq. 'fk5') write (18,'(a)') headreg(ii)(1:il)
             endif
           enddo
           if (headline(1)(2:3) .eq. '__') then
             cline = headline(1)
           else if (headline(2)(2:3) .eq. '__') then
             cline = headline(2)
           else if (headline(3)(2:3) .eq. '__') then
             cline = headline(3)
           else
             cline = headline(4)
           endif
           id = lastc ( cline )
           write (cnextr,'(i7)') nextr
           ncx = int ( log10 ( dble ( nextr ) ) ) + 1
           if (.not.lfootprints) then
             if (phy .eq. 'phy') write (15,'(a)') '#'//cline(2:id)
             if (fk5 .eq. 'fk5') write (17,'(a)') '#'//cline(2:id)
             if (phy .eq. 'phy') write (15,'(a)') '#'
             if (fk5 .eq. 'fk5') write (17,'(a)') '#'
             if (phy .eq. 'phy') write (15,'(a)') '# HALF-MAX ELLIPSES, TOTAL OF '//cnextr(7-ncx+1:7)//' SOURCES FROM '''//excat
     &                      //''' (COMMENTED OUT ARE THOSE DEEMED NOT ACCEPTABLE)'
             if (fk5 .eq. 'fk5') write (17,'(a)') '# HALF-MAX ELLIPSES, TOTAL OF '//cnextr(7-ncx+1:7)//' SOURCES FROM '''//excat
     &                      //''' (COMMENTED OUT ARE THOSE DEEMED NOT ACCEPTABLE)'
             if (phy .eq. 'phy') write (15,'(a)') '#'
             if (fk5 .eq. 'fk5') write (17,'(a)') '#'
             if (phy .eq. 'phy') write (15,'(a,f6.2,6(a,f4.1))') '# ACCEPTABILITY CRITERIA:  NWAVE = '//cnwave//'  OBSBEAM ='
     &                     , obsbeam, '  GOODMIN >=', goodmin, '  SIGMIN >=', sigmin, '  SNR1MIN >', snr1min, '  SNR2MIN >', snr2min
     &                     , '  ELONGMAX <', elongmax
             if (fk5 .eq. 'fk5') write (17,'(a,f6.2,6(a,f4.1))') '# ACCEPTABILITY CRITERIA:  NWAVE = '//cnwave//'  OBSBEAM ='
     &                     , obsbeam, '  GOODMIN >=', goodmin, '  SIGMIN >=', sigmin, '  SNR1MIN >', snr1min, '  SNR2MIN >', snr2min
     &                     , '  ELONGMAX <', elongmax
           else
             if (phy .eq. 'phy') write (16,'(a)') '#'//cline(2:id)
             if (fk5 .eq. 'fk5') write (18,'(a)') '#'//cline(2:id)
             if (phy .eq. 'phy') write (16,'(a)') '#'
             if (fk5 .eq. 'fk5') write (18,'(a)') '#'
             if (phy .eq. 'phy') write (16,'(a)') '# HALF-MAX ELLIPSES, TOTAL OF '//cnextr(7-ncx+1:7)//' SOURCES FROM '''//excat
     &                      //''' (COMMENTED OUT ARE THOSE DEEMED NOT ACCEPTABLE)'
             if (fk5 .eq. 'fk5') write (18,'(a)') '# HALF-MAX ELLIPSES, TOTAL OF '//cnextr(7-ncx+1:7)//' SOURCES FROM '''//excat
     &                      //''' (COMMENTED OUT ARE THOSE DEEMED NOT ACCEPTABLE)'
             if (phy .eq. 'phy') write (16,'(a)') '#'
             if (fk5 .eq. 'fk5') write (18,'(a)') '#'
             if (phy .eq. 'phy') write (16,'(a,f6.2,6(a,f4.1))') '# ACCEPTABILITY CRITERIA:  NWAVE = '//cnwave//'  OBSBEAM ='
     &                     , obsbeam, '  GOODMIN >=', goodmin, '  SIGMIN >=', sigmin, '  SNR1MIN >', snr1min, '  SNR2MIN >', snr2min
     &                     , '  ELONGMAX <', elongmax
             if (fk5 .eq. 'fk5') write (18,'(a,f6.2,6(a,f4.1))') '# ACCEPTABILITY CRITERIA:  NWAVE = '//cnwave//'  OBSBEAM ='
     &                     , obsbeam, '  GOODMIN >=', goodmin, '  SIGMIN >=', sigmin, '  SNR1MIN >', snr1min, '  SNR2MIN >', snr2min
     &                     , '  ELONGMAX <', elongmax
           endif
           
           exit
         endif
         if (.not.lfootprints) then
           if (phy .eq. 'phy') write (15,'(a)') '#'//headline(m)(2:ihc)
           if (fk5 .eq. 'fk5') write (17,'(a)') '#'//headline(m)(2:ihc)
         else
           if (phy .eq. 'phy') write (16,'(a)') '#'//headline(m)(2:ihc)
           if (fk5 .eq. 'fk5') write (18,'(a)') '#'//headline(m)(2:ihc)
         endif
       enddo
       
       do ii=38,45
         il = lastc ( headreg(ii) )
         if (.not.lfootprints) then
           if (phy .eq. 'phy') write (15,'(a)') headreg(ii)(1:il)
           if (fk5 .eq. 'fk5') write (17,'(a)') headreg(ii)(1:il)
         else
           if (phy .eq. 'phy') write (16,'(a)') headreg(ii)(1:il)
           if (fk5 .eq. 'fk5') write (18,'(a)') headreg(ii)(1:il)
         endif
       enddo
       if (.not.lfootprints) then
         if (phy .eq. 'phy') write (15,'(a)') 'global color='//regcol//' width=3 dashlist=8 3 font="helvetica 10 normal roman"'
     &                 //' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
         if (phy .eq. 'phy') write (15,'(a)') 'physical'
         if (fk5 .eq. 'fk5') write (17,'(a)') 'global color='//regcol//' width=3 dashlist=8 3 font="helvetica 10 normal roman"'
     &                 //' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
         if (fk5 .eq. 'fk5') write (17,'(a)') 'fk5'
       else
         if (phy .eq. 'phy') write (16,'(a)') 'global color='//regcol//' width=3 dashlist=8 3 font="helvetica 10 normal roman"'
     &                 //' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
         if (phy .eq. 'phy') write (16,'(a)') 'physical'
         if (fk5 .eq. 'fk5') write (18,'(a)') 'global color='//regcol//' width=3 dashlist=8 3 font="helvetica 10 normal roman"'
     &                 //' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
         if (fk5 .eq. 'fk5') write (18,'(a)') 'fk5'
       endif
       
       rewind ( nunit )
       do k=1,nhead
         read (nunit,'(a)') headline(k)
         lhl = lastc ( headline(k) )
         if (lhl .gt. 5) nhcols = k
       enddo
       lhl = lastc ( headline(nhcols) )
       
       do j=1,nextr
       
         read (nunit,'(a)') cline
         icl = lastc ( cline )
         
         if (cline(1:1) .eq. '!' .or. cline(1:1) .eq. '#') goto 88
         
         read (cline,*) no
         
         read (cline(ncolumnx:),*) x0o
         read (cline(ncolumny:),*) y0o
         read (cline(ncolumnxas:),*) x0oas
         read (cline(ncolumnyas:),*) y0oas
         
         goodness = 1.0d99
         if (ncolumng .gt. 0) then
           read (cline(ncolumng:),*) goodness
         endif
         goodmono = 1.0d99
         if (ncolumngm .gt. 0) then
           read (cline(ncolumngm:),*) goodmono
         endif
         flagmon = 0
         if (ncolumnfm .gt. 0) then
           read (cline(ncolumnfm:),*) flagmon
         endif
         
         read (cline(ncolumna:), *) majdiam
         read (cline(ncolumnb:), *) mindiam
         read (cline(ncolumnt:), *) angle
         read (cline(ncolumnfp:),*) fxpbest
         read (cline(ncolumnft:),*) fxtbest
         read (cline(ncolumnr:), *) sigmono
       
         if (ncolumnfa .gt. 0) then
           read (cline(ncolumnfa:), *) fmajdiam
         else
           fmajdiam = 2.3d0 * majdiam
         endif
         if (ncolumnfb .gt. 0) then
           read (cline(ncolumnfb:), *) fmindiam
         else
           fmindiam = 2.3d0 * mindiam
         endif
       
         majdiam = majdiam * sizefact
         mindiam = mindiam * sizefact
         fmajdiam = fmajdiam * sizefact
         fmindiam = fmindiam * sizefact
         
         if (cattype .eq. 'gauss') angle = - angle - 90.0d0
                
         if (cattype .ne. 'truth') then
           read (cline(ncolumnfpe:),*) fxperro
           read (cline(ncolumnfte:),*) fxterro
         else
           fxperro = 0.0d0
           fxterro = 0.0d0
         endif
         
         if (thickness .gt. 0.0d0) then
           nxe = int ( max ( majdiam, mindiam ) / dx * 2.0d0 ) + 10 !<-- number of pixels in the ellipse window in the X-dir
           nye = int ( max ( majdiam, mindiam ) / dy * 2.0d0 ) + 10 !<-- number of pixels in the ellipse window in the Y-dir
         else
           nxe = int ( obsbeam / dx * 2.0d0 ) + 2
           nye = int ( obsbeam / dy * 2.0d0 ) + 2
         endif
         nxco = nint ( x0o )  !<-- ellipse center's X-coordinate expressed in pixels
         nyco = nint ( y0o )  !<-- ellipse center's Y-coordinate expressed in pixels
       
         radx = majdiam / 2.0d0  !<-- major radius of the current ellipse
         rady = mindiam / 2.0d0  !<-- minor radius of the current ellipse
         
         write (cj,'(i6)') j
         ncj = int ( log10 ( dble ( j ) ) ) + 1
       
         if (j .gt. 0) then
         
           if (abs ( sigmono ) .ge. sigmin .and. 
     &         abs ( goodmono ) .ge. goodmin .and. 
     &         fxpbest .gt. snr1min * fxperro .and.
     &         fxtbest .gt. snr2min * fxterro .and.
     &         fmajdiam .lt. elongmax * fmindiam .and.
     &         fmajdiam .gt. 1.15d0 * majdiam) then
             sharp = ''
             is = 0
           else
             sharp = '# '
             is = 2
           endif

           ix = int ( log10 ( max ( x0o, 1.0d0 ) ) ) + 1
           iy = int ( log10 ( max ( y0o, 1.0d0 ) ) ) + 1
           write (cxpixm,'(f11.4)') x0o
           write (cypixm,'(f11.4)') y0o
           cxpixm = cxpixm(6-ix+1:11)
           cypixm = cypixm(6-iy+1:11)
           ix = lastc ( cxpixm )
           iy = lastc ( cypixm )
           
           it = int ( log10 ( max ( angle, 1.0d0 ) ) ) + 1
           write (ctheta,'(f11.4)') angle
           ctheta = ctheta(6-it+1:11)
           it = lastc ( ctheta )
           
           ia = int ( log10 ( max ( majdiam / 2.0d0, 1.0d0 ) ) ) + 1
           ib = int ( log10 ( max ( mindiam / 2.0d0, 1.0d0 ) ) ) + 1
           write (cafwhm,'(f11.4)') majdiam / 2.0d0
           write (cbfwhm,'(f11.4)') mindiam / 2.0d0
           cafwhm = cafwhm(6-ia+1:11)
           cbfwhm = cbfwhm(6-ib+1:11)
           ia = lastc ( cafwhm )
           ib = lastc ( cbfwhm )
           
           ic = int ( log10 ( max ( fmajdiam / 2.0d0, 1.0d0 ) ) ) + 1
           id = int ( log10 ( max ( fmindiam / 2.0d0, 1.0d0 ) ) ) + 1
           write (cfoota,'(f11.4)') fmajdiam / 2.0d0
           write (cfootb,'(f11.4)') fmindiam / 2.0d0
           cfoota = cfoota(6-ic+1:11)
           cfootb = cfootb(6-id+1:11)
           ic = lastc ( cfoota )
           id = lastc ( cfootb )
           
           if (.not.lfootprints) then
             if (phy .eq. 'phy') write (15,'(a)') sharp(1:is)//'ellipse('//cxpixm(1:ix)//','//cypixm(1:iy)//','//cbfwhm(1:ib)
     &                      //'",'//cafwhm(1:ia)//'",'//ctheta(1:it)//') # '//cj(6-ncj+1:6)
           else
             if (phy .eq. 'phy') write (16,'(a)') sharp(1:is)//'ellipse('//cxpixm(1:ix)//','//cypixm(1:iy)//','//cfootb(1:id)
     &                      //'",'//cfoota(1:ic)//'",'//ctheta(1:it)//') # '//cj(6-ncj+1:6)
           endif
         
           ix = int ( log10 ( max ( abs ( x0oas ), 1.0d0 ) ) ) + 1
           iy = int ( log10 ( max ( abs ( y0oas ), 1.0d0 ) ) ) + 1
           write (cxpixm,'(f13.7)') x0oas
           write (cypixm,'(f13.7)') y0oas
           cxpixm = cxpixm(4-ix+1:13)
           cypixm = cypixm(4-iy+1:13)
           ix = lastc ( cxpixm )
           iy = lastc ( cypixm )
           if (cxpixm(1:1) .eq. ' ') cxpixm = cxpixm(2:ix)
           if (cypixm(1:1) .eq. ' ') cypixm = cypixm(2:iy)
           ix = lastc ( cxpixm )
           iy = lastc ( cypixm )
         
           if (.not.lfootprints) then
             if (fk5 .eq. 'fk5') write (17,'(a)') sharp(1:is)//'ellipse('//cxpixm(1:ix)//','//cypixm(1:iy)//','//cbfwhm(1:ib)
     &                      //'",'//cafwhm(1:ia)//'",'//ctheta(1:it)//') # '//cj(6-ncj+1:6)
           else
             if (fk5 .eq. 'fk5') write (18,'(a)') sharp(1:is)//'ellipse('//cxpixm(1:ix)//','//cypixm(1:iy)//','//cfootb(1:id)
     &                      //'",'//cfoota(1:ic)//'",'//ctheta(1:it)//') # '//cj(6-ncj+1:6)
           endif
         endif
 88      continue     
       enddo
       
       do ii=47,48
         il = lastc ( headreg(ii) )
         if (.not.lfootprints) then
           if (phy .eq. 'phy') write (15,'(a)') headreg(ii)(1:il)
           if (fk5 .eq. 'fk5') write (17,'(a)') headreg(ii)(1:il)
         else
           if (phy .eq. 'phy') write (16,'(a)') headreg(ii)(1:il)
           if (fk5 .eq. 'fk5') write (18,'(a)') headreg(ii)(1:il)
         endif
       enddo
       
       if (.not.lfootprints) then
         if (phy .eq. 'phy') write (15,'(a)') 'physical'
         if (fk5 .eq. 'fk5') write (17,'(a)') 'fk5'
       else
         if (phy .eq. 'phy') write (16,'(a)') 'physical'
         if (fk5 .eq. 'fk5') write (18,'(a)') 'fk5'
       endif
       
       rewind ( nunit )
       do k=1,nhead
         read (nunit,'(a)') headline(k)
         lhl = lastc ( headline(k) )
         if (lhl .gt. 5) nhcols = k
       enddo
       lhl = lastc ( headline(nhcols) )
       
       do j=1,nextr
       
         read (nunit,'(a)') cline
         
         if (cline(1:1) .eq. '!' .or. cline(1:1) .eq. '#') goto 99
         
         read (cline,*) no
         read (cline(ncolumnx:),*) x0o
         read (cline(ncolumny:),*) y0o
         read (cline(ncolumnxas:),*) x0oas
         read (cline(ncolumnyas:),*) y0oas
       
         goodness = 1.0d99
         if (ncolumng .gt. 0) then
           read (cline(ncolumng:),*) goodness
         endif
         goodmono = 1.0d99
         if (ncolumngm .gt. 0) then
           read (cline(ncolumngm:),*) goodmono
         endif
         flagmon = 0
         if (ncolumnfm .gt. 0) then
           read (cline(ncolumnfm:),*) flagmon
         endif
         
         read (cline(ncolumna:), *) majdiam
         read (cline(ncolumnb:), *) mindiam
         read (cline(ncolumnt:), *) angle
         read (cline(ncolumnfp:),*) fxpbest
         read (cline(ncolumnft:),*) fxtbest
         read (cline(ncolumnr:), *) sigmono
       
         if (ncolumnfa .gt. 0) then
           read (cline(ncolumnfa:), *) fmajdiam
         else
           fmajdiam = 2.3d0 * majdiam
         endif
         if (ncolumnfb .gt. 0) then
           read (cline(ncolumnfb:), *) fmindiam
         else
           fmindiam = 2.3d0 * mindiam
         endif
       
         majdiam = majdiam * sizefact
         mindiam = mindiam * sizefact
         fmajdiam = fmajdiam * sizefact
         fmindiam = fmindiam * sizefact
         
         if (cattype .eq. 'gauss') angle = - angle - 90.0d0
                
         if (cattype .ne. 'truth') then
           read (cline(ncolumnfpe:),*) fxperro
           read (cline(ncolumnfte:),*) fxterro
         else
           fxperro = 0.0d0
           fxterro = 0.0d0
         endif
         
         if (thickness .gt. 0.0d0) then
           nxe = int ( max ( majdiam, mindiam ) / dx * 2.0d0 ) + 10 !<-- number of pixels in the ellipse window in the X-dir
           nye = int ( max ( majdiam, mindiam ) / dy * 2.0d0 ) + 10 !<-- number of pixels in the ellipse window in the Y-dir
         else
           nxe = int ( obsbeam / dx * 2.0d0 ) + 2
           nye = int ( obsbeam / dy * 2.0d0 ) + 2
         endif
         nxco = nint ( x0o )  !<-- ellipse center's X-coordinate expressed in pixels
         nyco = nint ( y0o )  !<-- ellipse center's Y-coordinate expressed in pixels
       
         radx = majdiam / 2.0d0  !<-- major radius of the current ellipse
         rady = mindiam / 2.0d0  !<-- minor radius of the current ellipse
       
         write (cj,'(i6)') j
         ncj = int ( log10 ( dble ( j ) ) ) + 1

         if (j .gt. 0) then
       
           if (abs ( sigmono ) .ge. sigmin .and. 
     &         abs ( goodmono ) .ge. goodmin .and. 
     &         fxpbest .gt. snr1min * fxperro .and.
     &         fxtbest .gt. snr2min * fxterro .and.
     &         fmajdiam .lt. elongmax * fmindiam .and.
     &         fmajdiam .gt. 1.15d0 * majdiam) then
             sharp = ''
             is = 0
           else
             sharp = '# '
             is = 2
           endif
       
           ix = int ( log10 ( max ( x0o, 1.0d0 ) ) ) + 1
           iy = int ( log10 ( max ( y0o, 1.0d0 ) ) ) + 1
           write (cxpixm,'(f11.4)') x0o
           write (cypixm,'(f11.4)') y0o
           cxpixm = cxpixm(6-ix+1:11)
           cypixm = cypixm(6-iy+1:11)
           ix = lastc ( cxpixm )
           iy = lastc ( cypixm )
           
           it = int ( log10 ( max ( angle, 1.0d0 ) ) ) + 1
           write (ctheta,'(f11.4)') angle
           ctheta = ctheta(6-it+1:11)
           it = lastc ( ctheta )
           
           ia = int ( log10 ( max ( majdiam / 2.0d0, 1.0d0 ) ) ) + 1
           ib = int ( log10 ( max ( mindiam / 2.0d0, 1.0d0 ) ) ) + 1
           write (cafwhm,'(f11.4)') majdiam / 2.0d0
           write (cbfwhm,'(f11.4)') mindiam / 2.0d0
           cafwhm = cafwhm(6-ia+1:11)
           cbfwhm = cbfwhm(6-ib+1:11)
           ia = lastc ( cafwhm )
           ib = lastc ( cbfwhm )
           
           ic = int ( log10 ( max ( fmajdiam / 2.0d0, 1.0d0 ) ) ) + 1
           id = int ( log10 ( max ( fmindiam / 2.0d0, 1.0d0 ) ) ) + 1
           write (cfoota,'(f11.4)') fmajdiam / 2.0d0
           write (cfootb,'(f11.4)') fmindiam / 2.0d0
           cfoota = cfoota(6-ic+1:11)
           cfootb = cfootb(6-id+1:11)
           ic = lastc ( cfoota )
           id = lastc ( cfootb )
           
           if (.not.lfootprints) then
             if (phy .eq. 'phy') write (15,'(a)') sharp(1:is)//'ellipse('//cxpixm(1:ix)//','//cypixm(1:iy)//','//cbfwhm(1:ib)
     &                      //'",'//cafwhm(1:ia)//'",'//ctheta(1:it)//') # '//cj(6-ncj+1:6)
           else
             if (phy .eq. 'phy') write (16,'(a)') sharp(1:is)//'ellipse('//cxpixm(1:ix)//','//cypixm(1:iy)//','//cfootb(1:id)
     &                      //'",'//cfoota(1:ic)//'",'//ctheta(1:it)//') # '//cj(6-ncj+1:6)
           endif
       
           ix = int ( log10 ( max ( abs ( x0oas ), 1.0d0 ) ) ) + 1
           iy = int ( log10 ( max ( abs ( y0oas ), 1.0d0 ) ) ) + 1
           write (cxpixm,'(f13.7)') x0oas
           write (cypixm,'(f13.7)') y0oas
           cxpixm = cxpixm(4-ix+1:13)
           cypixm = cypixm(4-iy+1:13)
           ix = lastc ( cxpixm )
           iy = lastc ( cypixm )
           if (cxpixm(1:1) .eq. ' ') cxpixm = cxpixm(2:ix)
           if (cypixm(1:1) .eq. ' ') cypixm = cypixm(2:iy)
           ix = lastc ( cxpixm )
           iy = lastc ( cypixm )
       
           if (.not.lfootprints) then
             if (fk5 .eq. 'fk5') write (17,'(a)') sharp(1:is)//'ellipse('//cxpixm(1:ix)//','//cypixm(1:iy)//','//cbfwhm(1:ib)
     &                      //'",'//cafwhm(1:ia)//'",'//ctheta(1:it)//') # '//cj(6-ncj+1:6)
           else
             if (fk5 .eq. 'fk5') write (18,'(a)') sharp(1:is)//'ellipse('//cxpixm(1:ix)//','//cypixm(1:iy)//','//cfootb(1:id)
     &                      //'",'//cfoota(1:ic)//'",'//ctheta(1:it)//') # '//cj(6-ncj+1:6)
           endif
         endif
99       continue
       enddo
      
       if (.not.lfootprints) then
         close ( 15 )
         close ( 17 )
       else
         close ( 16 )
         close ( 18 )
       endif
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine showprogress ( iotty, iolog, lunix, icount, ntotal, isumpct )
!__________________________________________________________________________________________________________________________________
!
! Display progress bar in percentage of current work; each point means that 2 % of work has been done (50 points in total)
!__________________________________________________________________________________________________________________________________
!
       logical       lshow, lunix
       character*3   dot
       integer       i, ishow, iotty, iolog
       real*8        icount, ntotal, isumpct, ipcent
       parameter   ( dot = '•' )
       external      isleep
!__________________________________________________________________________________________________________________________________
!
       ipcent = icount * 100.0d0 / ntotal
       ishow = int ( ipcent / 2.0d0 - isumpct )
       ishow = max ( 0, min ( ishow, int ( 50.0d0 - isumpct ) ) )
       lshow = mod ( int ( ipcent ), 2 ) .eq. 0
     
       if (lshow) then
         do i=1,ishow
           if (iotty .gt. 0) write (iotty,'(a)', advance='no') dot
           if (iolog .gt. 0) write (iolog,'(a)', advance='no') dot
           if (lunix) then
             if (iotty .gt. 0) endfile ( iotty, err=888 )
 888         continue 
             if (iotty .gt. 0) backspace ( iotty )
           endif
           isumpct = isumpct + 1.0d0
         enddo
       endif
     
       if (abs ( icount - ntotal ) .lt. ntotal * 1.0d-14) then
          call isleep ( 1 )
          if (iotty .gt. 0) write (iotty,'(a)', advance='no') ' done'
          if (iolog .gt. 0) write (iolog,'(a)', advance='no') ' done'
          if (lunix) then
            if (iotty .gt. 0) endfile ( iotty, err=999 )
 999        continue 
            if (iotty .gt. 0) backspace ( iotty )
          endif
          call isleep ( 1 )
          if (iotty .gt. 0) write (iotty,'()')
          if (iolog .gt. 0) write (iolog,'()')
       endif
     
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine fdelete 
       
     &            ( fname, status )
!__________________________________________________________________________________________________________________________________
!
! A simple little routine to delete a FITS file
!
! May contain fragments from "The FITSIO Cookbook" v1.1 (April 1995) by William D Pence (HEASARC).
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) fname
       integer       status, unit, blocksize
       external      ftgiou, ftopen, ftdelt, ftcmsg, ftfiou
!__________________________________________________________________________________________________________________________________
!
! Simply return if status is greater than zero

       if (status .gt. 0) return

! Get an unused Logical Unit Number to use to open the FITS file

       call ftgiou ( unit, status )

! Try to open the file, to see if it exists

       call ftopen ( unit, fname, 1, blocksize, status )

       if (status .eq. 0) then

! File was opened; so now delete it

         call ftdelt ( unit, status )

       elseif (status .eq. 103) then

! File doesn't exist, so just reset status to zero and clear errors

         status = 0
         call ftcmsg
       else

! There was some other error opening the file; delete the file anyway

         status = 0
         call ftcmsg
         call ftdelt ( unit, status )
       endif

! Free the unit number for later reuse

       call ftfiou ( unit, status )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine printerr 
       
     &            ( status )
!__________________________________________________________________________________________________________________________________
!
! Print out the FITSIO error messages to the user
!
! May contain fragments from "The FITSIO Cookbook" v1.1 (April 1995) by William D Pence (HEASARC).
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*30  errtext
       character*80  errmess
       integer       status
       external      ftgerr, ftgmsg
!__________________________________________________________________________________________________________________________________
!
! Check if status is OK (no error); if so, simply return

       if (status .le. 0) return

! Get the text string which describes the error

       call ftgerr ( status, errtext )
       write (*,*) '   FITSIO Error Status =',status,': ',errtext

! Read and print out all the error messages on the FITSIO stack

       call ftgmsg ( errmess )
       do while (errmess .ne. ' ')
         write (*,*) errmess
         call ftgmsg ( errmess )
       enddo

       return
       end
                                                                                                                                   
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine when 
       
     &            ( lunix, ctime, cdate, ndate, iform )
!__________________________________________________________________________________________________________________________________
!
!  OUTPUT VARIABLES:
!
!    ctime       time returned as hh:mm:ss
!    cdate       date returned as specified by "iform"
!    ndate       number of characters set in "cdate"
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       lunix
       character*(*) ctime, cdate
       character*8   ttime, tdate
       character*12  timda
       character*24  udate
       integer*2     hrs, mins, secs, hsecs, year, month, day
       integer       iform, ndate, i, lastc
       external      lastc
!_______
!
! >>>>> UNCOMMENT THE NEXT LINE TO COMPILE ON UNIX:
!_______
!!!not4gfortran       !external fdate
      !external fdate
       external gettim, getdat
!__________________________________________________________________________________________________________________________________
!
! Date and time: "ttime" is in the format hh:mm:ss, "tdate" is in the format mm/dd/yy, "udate" is in a system-dependent format

       udate = ' '
       hrs   = 0
       mins  = 0
       secs  = 0
       hsecs = 0
       year  = 0
       month = 0
       day   = 0

       if (lunix) then
!_______
!
! >>>>> UNCOMMENT THE NEXT LINE TO COMPILE ON UNIX:
!_______
          call fdate ( udate )

          ttime      = udate(12:19)
          tdate(4:5) = udate( 9:10)
          tdate(7:8) = udate(23:24)
          tdate(3:3) = '/'
          tdate(6:6) = '/'
          if (udate(5:7) .eq. 'Jan') tdate(1:2) = '01'
          if (udate(5:7) .eq. 'Feb') tdate(1:2) = '02'
          if (udate(5:7) .eq. 'Mar') tdate(1:2) = '03'
          if (udate(5:7) .eq. 'Apr') tdate(1:2) = '04'
          if (udate(5:7) .eq. 'May') tdate(1:2) = '05'
          if (udate(5:7) .eq. 'Jun') tdate(1:2) = '06'
          if (udate(5:7) .eq. 'Jul') tdate(1:2) = '07'
          if (udate(5:7) .eq. 'Aug') tdate(1:2) = '08'
          if (udate(5:7) .eq. 'Sep') tdate(1:2) = '09'
          if (udate(5:7) .eq. 'Oct') tdate(1:2) = '10'
          if (udate(5:7) .eq. 'Nov') tdate(1:2) = '11'
          if (udate(5:7) .eq. 'Dec') tdate(1:2) = '12'

       else
!_______
!
! >>>>> COMMENT OUT THE NEXT 2 LINES TO COMPILE ON UNIX:
!_______
!!          call gettim ( hrs, mins, secs, hsecs )
!!          call getdat ( year, month, day )

          write (timda, *) hrs
          if (timda(11:11) .eq. ' ') timda(11:11)='0'
          ttime(1:2) = timda(11:12)
          ttime(3:3) = ':'
          write (timda, *) mins
          if (timda(11:11) .eq. ' ') timda(11:11)='0'
          ttime(4:5) = timda(11:12)
          ttime(6:6) = ':'
          write (timda, *) secs
          if (timda(11:11) .eq. ' ') timda(11:11)='0'
          ttime(7:8) = timda(11:12)
          write (timda, *) month
          if (timda(11:11) .eq. ' ') timda(11:11)='0'
          tdate(1:2) = timda(11:12)
          tdate(3:3) = '/'
          write (timda, *) day
          if (timda(11:11) .eq. ' ') timda(11:11)='0'
          tdate(4:5) = timda(11:12)
          tdate(6:6) = '/'
          write (timda, *) year
          tdate(7:8) = timda(11:12)
       endif

! Convert to desired format.

       if (iform .eq. 1) then
         ctime(1:8) = ttime
         cdate(1:8) = tdate
         cdate(1:2) = tdate(4:5)
         cdate(4:5) = tdate(1:2)
         ndate      = 8
       endif

       if (iform .eq. 2) then
         ctime(1:8) = ttime
         cdate(1:8) = tdate
         ndate      = 8
       endif

       if (iform .eq. 3) then
         ctime(1:8) = ttime
         if (tdate(1:2) .eq. '01') cdate(1:8)  = 'January '
         if (tdate(1:2) .eq. '02') cdate(1:9)  = 'February '
         if (tdate(1:2) .eq. '03') cdate(1:6)  = 'March '
         if (tdate(1:2) .eq. '04') cdate(1:6)  = 'April '
         if (tdate(1:2) .eq. '05') cdate(1:4)  = 'May '
         if (tdate(1:2) .eq. '06') cdate(1:5)  = 'June '
         if (tdate(1:2) .eq. '07') cdate(1:5)  = 'July '
         if (tdate(1:2) .eq. '08') cdate(1:7)  = 'August '
         if (tdate(1:2) .eq. '09') cdate(1:10) = 'September '
         if (tdate(1:2) .eq. '10') cdate(1:8)  = 'October '
         if (tdate(1:2) .eq. '11') cdate(1:9)  = 'November '
         if (tdate(1:2) .eq. '12') cdate(1:9)  = 'December '
         i = lastc ( cdate )
         cdate(i+2:i+3) = tdate(4:5)
         cdate(i+4:i+7) = ', 20'
         cdate(i+8:i+9) = tdate(7:8)
         ndate          = i + 9
       endif

       if (iform .eq. 4) then
         ctime(1:8)  = ttime
         cdate(1:8)  = tdate
         cdate(1:2)  = tdate(4:5)
         cdate(4:5)  = tdate(1:2)
         cdate(9:10) = tdate(7:8)
         cdate(7:8)  = '20'
         ndate       = 10
       endif

       ctime(3:3) = ':'
       ctime(6:6) = ':'
       if (cdate(1:1) .eq. ' ') cdate(1:1) = '0'
       if (cdate(4:4) .eq. ' ') cdate(4:4) = '0'

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function timer 
       
     &          ( ttype, otime )
!__________________________________________________________________________________________________________________________________
!
!  PURPOSE: Timer for measuring execution time and CPU time. This code has been extracted from the main code of ZEUS-3D into this 
!  function in order to simplify timing of various parts of the code. 
!
!  INPUT VARIABLES:
!    ttype       'cpu' - to measure CPU time, 'wal' - to measure wall clock time
!    otime       previously taken measurement that will be subtracted from the new system timer value, giving time interval
!                (otime can be zero).
!
!  OUTPUT VARIABLES:
!    timer       time interval or initial timer value (for otime = 0.0d0).
!__________________________________________________________________________________________________________________________________
!
       character*3   ttype
       character*8   crdate
       character*10  crtime
       real          cptime
       real*8        timer, otime, ryear, rmon, rday, rhour, rmin, rsec, rtime
!__________________________________________________________________________________________________________________________________
!
       if (ttype .eq. 'cpu') then
         call cpu_time ( cptime )
         rtime = dble ( cptime )
       else
         call date_and_time ( crdate, crtime )
         read (crdate(3:4),*) ryear
         read (crdate(5:6),*) rmon
         read (crdate(7:8),*) rday
         read (crtime(1:2),*) rhour
         read (crtime(3:4),*) rmin
         read (crtime(5:10),*) rsec
         rtime = 31536000.0d0 * ryear + 2629822.96584d0 * rmon + 86400.0d0 * rday + 3600.0d0 * rhour + 60.0d0 * rmin + rsec
       endif
       timer = rtime - otime

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function lastc 
       
     &          ( string )
!__________________________________________________________________________________________________________________________________
!
! Returns position of a last non-blank character in a string.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) string
       integer       lastc, ifinal, i, j
!__________________________________________________________________________________________________________________________________
!
       ifinal = len ( string )
       do i=1,ifinal
         j = ifinal + 1 - i
         if (string(j:j) .ne. ' ') goto 10
       enddo
       j = 0
10     continue
       lastc = j

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function firstb 
       
     &          ( string )
!__________________________________________________________________________________________________________________________________
!
! Returns position of the first blank (" ") in a string.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) string
       integer       firstb, ifinal, i
!__________________________________________________________________________________________________________________________________
!
       ifinal = len ( string )
       do i=1,ifinal
         if (string(i:i) .eq. ' ') goto 10
       end do
       i = 0
10     continue
       firstb = i

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function lastd 
       
     &          ( string )
!__________________________________________________________________________________________________________________________________
!
! Returns position of last dot (".") character in a string.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) string
       integer       lastd, ifinal, i, j
!__________________________________________________________________________________________________________________________________
!
       ifinal = len ( string )
       do i=1,ifinal
         j = ifinal + 1 - i
         if (string(j:j) .eq. '.') goto 10
       end do
       j = 0
10     continue
       lastd = j

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function lasts ( string )
!__________________________________________________________________________________________________________________________________
!
! Returns position of last slash ("/") character in a string.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) string
       integer       lasts, ifinal, i, j
!__________________________________________________________________________________________________________________________________
!
       ifinal = len ( string )
       do i=1,ifinal
         j = ifinal + 1 - i
         if (string(j:j) .eq. '/') goto 10
       end do
       j = 0
10     continue
       lasts = j

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function lastbs 
       
     &          ( string )
!__________________________________________________________________________________________________________________________________
!
! Returns position of last backslash ("\") character in a string.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) string
       integer       lastbs, ifinal, i, j
!__________________________________________________________________________________________________________________________________
!
       ifinal = len ( string )
       do i=1,ifinal
         j = ifinal + 1 - i
         if (string(j:j) .eq. '\') goto 10
       end do
       j = 0
10     continue
       lastbs = j

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function lastp
       
     &          ( string )
!__________________________________________________________________________________________________________________________________
!
! Returns position of last plus ("+") character in a string.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) string
       integer       lastp, ifinal, i, j
!__________________________________________________________________________________________________________________________________
!
       ifinal = len ( string )
       do i=1,ifinal
         j = ifinal + 1 - i
         if (string(j:j) .eq. '+') goto 10
       end do
       j = 0
10     continue
       lastp = j

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function last
       
     &          ( sub, string )
!__________________________________________________________________________________________________________________________________
!
! Returns position of the last instance of a substring in a string.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) sub, string
       integer       last, isub, ifinal, i, j, je
!__________________________________________________________________________________________________________________________________
!
       isub = len ( sub ) - 1
       ifinal = len ( string ) - isub
       do i=1,ifinal
         j = ifinal + 1 - i
         je = j + isub
         if (string(j:je) .eq. sub) goto 10
       enddo
       j = 0
10     continue
       last = j

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine sortarr 
       
     &            ( n, x, y )
!__________________________________________________________________________________________________________________________________
!
! Sorts array x(j) in ascending order and rearranges array y(j) according to x(j).
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       i, j, n
       real*8        x(*), y(*), a, b
!__________________________________________________________________________________________________________________________________
!
       do j=2,n
         a = y(j)
         b = x(j)
         do i=j-1,1,-1
           if (x(i) .le. b) then
             goto 10
           else
             y(i+1) = y(i)
             x(i+1) = x(i)
           end if
         enddo
         i = 0
 10      continue
         y(i+1) = a
         x(i+1) = b
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function medapprox
      
     &          ( npts, xdata, leftend, rightend )
!__________________________________________________________________________________________________________________________________
!
! Adapted by A. Men'shchikov from a code for the binapprox algorithm by Ryan J. Tibshirani (2008). 
! Webpage: http://www.stat.cmu.edu/~ryantibs/median/ Paper: http://www.stat.cmu.edu/~ryantibs/papers/median.pdf
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       nbins
       parameter   ( nbins = 1000 )
       integer       i, j, nb, kmed, npts, leftcount, ntotcount, ncounts(0:nbins)
       real*8        medapprox, scalefac, leftend, rightend, xdata(npts)           !! , xmean, sigma, dnpts
!__________________________________________________________________________________________________________________________________
!
!!       dnpts = dble ( npts )
!!       xmean = sum ( xdata ) / dnpts
!!       sigma = sqrt ( dot_product ( xdata - xmean, xdata - xmean ) / dnpts )
!!       leftend  = xmean - sigma
!!       rightend = xmean + sigma

       if (leftend .eq. rightend) then
         medapprox = leftend
         return
       endif

       scalefac = dble ( nbins ) / (rightend - leftend)
       leftcount = 0
       do i=0,nbins
         ncounts(i) = 0
       enddo

       do i=1,npts
         if (xdata(i) .lt. leftend) then
           leftcount = leftcount + 1
         else if (xdata(i) .lt. rightend) then
           nb = int ( (xdata(i) - leftend ) * scalefac )
           ncounts(nb) = ncounts(nb) + 1
         endif
       enddo

       if (mod ( npts, 2 ) .eq. 1) then
       
! Odd number 'npts' of data set elements, find the bin that contains the median.
       
         kmed = (npts + 1) / 2
         ntotcount = leftcount
         do i=0,nbins
           ntotcount = ntotcount + ncounts(i)
           if (ntotcount .ge. kmed) then
             medapprox = dble ( 2 * i + 1 ) / (2.0d0 * scalefac) + leftend
             exit
           endif
         enddo
       else
       
! Even number 'npts' of data set elements, find the bins that contain the medians.
       
         kmed = npts / 2
         ntotcount = leftcount
         do i=0,nbins
           ntotcount = ntotcount + ncounts(i)
           if (ntotcount .ge. kmed) then
             j = i
   1         continue
             if (ntotcount .eq. kmed) then
               j = j + 1
               ntotcount = ntotcount + ncounts(j)
               goto 1
             endif
             medapprox = dble ( i + j + 1 ) / (2.0d0 * scalefac) + leftend
             exit
           endif
         enddo
       endif

       medapprox = 0.0d0 !<-- this statement should never be executed!

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function modeapprox
      
     &          ( npts, xdata, leftend, rightend )
!__________________________________________________________________________________________________________________________________
!
! Adapted by A. Men'shchikov from a code for the binapprox algorithm by Ryan J. Tibshirani (2008). 
! Webpage: http://www.stat.cmu.edu/~ryantibs/median/ Paper: http://www.stat.cmu.edu/~ryantibs/papers/median.pdf
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       nbins
       parameter   ( nbins = 1000 )
       integer       i, nb, npts, leftcount, ncountsmax, ncounts(0:nbins)
       real*8        modeapprox, scalefac, leftend, rightend, xdata(npts)
!__________________________________________________________________________________________________________________________________
!
       modeapprox = leftend
       if (leftend .eq. rightend) then
         return
       endif

       scalefac = dble ( nbins ) / (rightend - leftend)
       leftcount = 0
       do i=0,nbins
         ncounts(i) = 0
       enddo

       do i=1,npts
         if (xdata(i) .lt. leftend) then
           leftcount = leftcount + 1
         else if (xdata(i) .lt. rightend) then
           nb = int ( (xdata(i) - leftend ) * scalefac )
           ncounts(nb) = ncounts(nb) + 1
         endif
       enddo

! Find the mode value.

       ncountsmax = 0
       do i=0,nbins
         if (ncounts(i) .gt. ncountsmax) then
           ncountsmax = ncounts(i)
           modeapprox = dble ( 2 * i + 1 ) / (2.0d0 * scalefac) + leftend
         endif
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine traceskels

     &            ( nx, ny, skeletons, filename, iotty, nextr, npixfil, np, xpix, ypix, xpixok, ypixok, pixd, pixdist, cverbose
     &            , almostzero, lnobranches )
!__________________________________________________________________________________________________________________________________
!
! Find coordinates of skeletons.
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       lendpoint, lused, lnobranches, ldebug, looping
       character*(*) filename, cverbose
       integer       k, i, j, im1, ip1, jm1, jp1, ntry, ico1, jco1, nx, ny, npixconnmin, npixconn4, npixconn8, npix4, npix8, numfil
     &             , l, n0, im, jm, ii, jj, m, n, lm1, nextr, npixconn, iotty, ia, ja, nzero, ico, jco, numpix4, numpix8
       integer       npixfil(*), np(*)
       real*8        xcoo, ycoo, almostzero, delta, deltamin, dist
       real*8        skeletons(nx,ny), pixdist(nx,ny), xpix(*), ypix(*), xpixok(nextr,*), ypixok(nextr,*), pixd(*)
!__________________________________________________________________________________________________________________________________
!
       lnobranches = .true.

       do k=1,nextr
         ntry = 0
 80      continue
         xcoo = 0.0d0
         ycoo = 0.0d0
         npixfil(k) = 0
         npixconnmin = 8
         ico1 = 0
         jco1 = 0

         do j=1,ny
           jm1 = max ( j - 1, 1 )
           jp1 = min ( j + 1, ny )
           do i=1,nx
             im1 = max ( i - 1, 1 )
             ip1 = min ( i + 1, nx )
             numfil = nint ( skeletons(i,j) )

             if (k .eq. numfil) then
               if (ico1 .eq. 0) ico1 = i
               if (jco1 .eq. 0) jco1 = j
               npixfil(k) = npixfil(k) + 1
               npixconn4 = 0
               npixconn8 = 0
               npix4 = 0
               npix8 = 0
               if (skeletons(im1,j) .gt. almostzero) then
                 npixconn4 = npixconn4 + 1
                 npix4 = 1
               endif
               if (skeletons(ip1,j) .gt. almostzero) then
                 npixconn4 = npixconn4 + 1
                 npix4 = 2
               endif
               if (skeletons(i,jm1) .gt. almostzero) then
                 npixconn4 = npixconn4 + 1
                 npix4 = 3
               endif
               if (skeletons(i,jp1) .gt. almostzero) then
                 npixconn4 = npixconn4 + 1
                 npix4 = 4
               endif
               if (skeletons(im1,jm1) .gt. almostzero) then
                 npixconn8 = npixconn8 + 1
                 npix8 = 1
               endif
               if (skeletons(ip1,jp1) .gt. almostzero) then
                 npixconn8 = npixconn8 + 1
                 npix8 = 2
               endif
               if (skeletons(im1,jp1) .gt. almostzero) then
                 npixconn8 = npixconn8 + 1
                 npix8 = 3
               endif
               if (skeletons(ip1,jm1) .gt. almostzero) then
                 npixconn8 = npixconn8 + 1
                 npix8 = 4
               endif
               npixconn = npixconn4 + npixconn8
               lendpoint = .false.
               if ((npixconn4 .le. 1 .and. npixconn8 .le. 1) .and. npixconn .le. npixconnmin) then
                 if (npixconn4 .eq. 1 .and. npixconn8 .eq. 1) then
                   if (npix4 .eq. 1 .and. (npix8 .eq. 1 .or. npix8 .eq. 3)) lendpoint = .true.
                   if (npix4 .eq. 2 .and. (npix8 .eq. 2 .or. npix8 .eq. 4)) lendpoint = .true.
                   if (npix4 .eq. 3 .and. (npix8 .eq. 1 .or. npix8 .eq. 4)) lendpoint = .true.
                   if (npix4 .eq. 4 .and. (npix8 .eq. 3 .or. npix8 .eq. 2)) lendpoint = .true.
                 else
                   lendpoint = .true.
                 endif
               endif
               if (lendpoint) then
                 xcoo = dble ( i )
                 ycoo = dble ( j )
               endif
               npixconnmin = min ( npixconnmin, npixconn )
             endif
           enddo
         enddo
         
! Fix a problem and try again, if the first pixel of the current filament has not been found (using the above logic).      
         
         if (nint ( xcoo ) .eq. 0 .or. nint ( ycoo ) .eq. 0) then
           if (ico1 .gt. 0 .and. jco1 .gt. 0) then
             if (skeletons(ico1,jco1) .gt. almostzero) then
               skeletons(ico1,jco1) = 0.0d0
               np(k) = np(k) - 1
             endif
             if (skeletons(ico1-1,jco1) .gt. almostzero) then
               skeletons(ico1-1,jco1) = 0.0d0
               np(k) = np(k) - 1
             endif
             if (skeletons(ico1+1,jco1) .gt. almostzero) then
               skeletons(ico1+1,jco1) = 0.0d0
               np(k) = np(k) - 1
             endif
             if (skeletons(ico1,jco1-1) .gt. almostzero) then
               skeletons(ico1,jco1-1) = 0.0d0
               np(k) = np(k) - 1
             endif
             if (skeletons(ico1,jco1+1) .gt. almostzero) then
               skeletons(ico1,jco1+1) = 0.0d0
               np(k) = np(k) - 1
             endif
             if (skeletons(ico1-1,jco1-1) .gt. almostzero) then
               skeletons(ico1-1,jco1-1) = 0.0d0
               np(k) = np(k) - 1
             endif
             if (skeletons(ico1+1,jco1+1) .gt. almostzero) then
               skeletons(ico1+1,jco1+1) = 0.0d0
               np(k) = np(k) - 1
             endif
             if (skeletons(ico1-1,jco1+1) .gt. almostzero) then
               skeletons(ico1-1,jco1+1) = 0.0d0
               np(k) = np(k) - 1
             endif
             if (skeletons(ico1+1,jco1-1) .gt. almostzero) then
               skeletons(ico1+1,jco1-1) = 0.0d0
               np(k) = np(k) - 1
             endif
           endif
           ntry = ntry + 1
           if (ntry .le. 2) then  !<-- Do no more than 3 attempts to find the first pixel.
             goto 80
           else
             if (cverbose .eq. '-verb2' .and. iotty .gt. 0) then
               write (iotty,'(a,i5)') '   TRACESKELS: Trouble finding first pixel of skeleton:', k
               write (iotty,'(a)'   ) '   TRACESKELS: in '//filename
             endif
             lnobranches = .false.
!!!             stop 99
           endif
         endif
         
! The first pixel has been defined, now search for the last pixel of the current filament.
         
         ico = nint ( xcoo )
         jco = nint ( ycoo )
         l = 0
         n = 0
         lm1 = 0
         delta = 0.0d0
 111     continue     
         l = l + 1
         if (l .eq. 1) then
           n0 = 0  !<-- At the first pixel, a one-pixel window (zero pixel radius)
         else
           n0 = 1  !<-- At all other pixels, a 9-pixel window (one pixel radius)
         endif
         deltamin = 1.0d+30
         im = 0
         jm = 0
         do jj=max(jco-n0,1),min(jco+n0,ny)
           do ii=max(ico-n0,1),min(ico+n0,nx)
             if (k .eq. nint ( skeletons(ii,jj) )) then

! The pixel (ii,jj) in this small window belongs to the filament k. 
! Check if it has already been considered (used) in this search for the last pixel.
               
               lused = .false.
               if (l .gt. 1) then
                 do m=1,l-1
                   if (ii .eq. nint ( xpix(m) ) .and. jj .eq. nint ( ypix(m) )) then
                     lused = .true.
                     exit
                   endif
                 enddo
               endif
               
! If the pixel (ii,jj) in this small window has not been considered (used), find its distance from the previously accepted pixel.
! Find the minimum distance between all pixels of the window and the previously accepted pixel.
               
               if (.not.lused) then
                 if (l .gt. 1) then
                   delta = sqrt ( (dble ( ii ) - xpix(lm1))**2 + (dble ( jj ) - ypix(lm1))**2 )
                 endif
                 if (delta .lt. deltamin) then
                   deltamin = delta
                   im = ii
                   jm = jj
                 endif
               endif
             endif
           enddo
         enddo

         if (im .gt. 0 .and. jm .gt. 0) then

! The minimum distance within the small 9-pixel window has been found. 
! Consider the pixel (im,jm) as the next current pixel of the skeleton.
! Store the distance of the current pixel from the first pixel of the skeleton.
! Create an image of the distances of each pixel from the first one.

           xpix(l) = dble ( im )
           ypix(l) = dble ( jm )
           if (l .eq. 1) then
             pixd(l) = 1.0d0
           else
             pixd(l) = pixd(lm1) + deltamin
           endif
           pixdist(im,jm) = pixd(l)
           ico = im
           jco = jm
           lm1 = l
           n = l
           if (l .lt. np(k)) goto 111  !<-- turn to the next pixel of a skeleton.
         else
         
! The minimum distance has not been determined as all pixels in the 9-pixel window have been used.
! There are further pixels belonging to the skeleton, because of its branching.
! Skip back over the previously accepted pixels until the branch is found and continue.

           l = l - 1
           n = n - 1
           if (n .le. 1) then
             if (cverbose .eq. '-verb2' .and. iotty .gt. 0) then
               write (iotty,'(a,4i5)') '   TRACESKELS: Trouble finding last pixel of skeleton:', k, l, ico, jco
               write (iotty,'(a)'    ) '   TRACESKELS: in '//filename
             endif
             lnobranches = .false.
!!             stop 99
           else
             ico = nint ( xpix(n) )
             jco = nint ( ypix(n) )
             lm1 = n
             goto 111
           endif
         endif

! The pixel (im,jm) is accepted as the next current pixel of the skeleton k.
         
         do l=1,np(k)
           xpixok(k,l) = xpix(l)
           ypixok(k,l) = ypix(l)
         enddo

! Identify and separate skeleton branches.
         
         ldebug = .false.

         do l=2,np(k)
           dist = sqrt ( (xpixok(k,l) - xpixok(k,l-1))**2 + (ypixok(k,l) - ypixok(k,l-1))**2 )
           if (dist .gt. 1.45d0) then
             nzero = np(k) - (l - 1)
             ico = nint ( xpixok(k,l) )
             jco = nint ( ypixok(k,l) )
             looping = .false.
             
! Complicated case of a looping skeleton: branch is just one pixel and it has two 4-connected pixels of the same skeleton.

             if (l .eq. np(k)) then
               im1 = max ( ico - 1, 1 )
               ip1 = min ( ico + 1, nx )
               jm1 = max ( jco - 1, 1 )
               jp1 = min ( jco + 1, ny )
               numpix4 = 0
               numpix8 = 0
               if (skeletons(im1,jco) .gt. almostzero) numpix4 = numpix4 + 1
               if (skeletons(ip1,jco) .gt. almostzero) numpix4 = numpix4 + 1
               if (skeletons(ico,jm1) .gt. almostzero) numpix4 = numpix4 + 1
               if (skeletons(ico,jp1) .gt. almostzero) numpix4 = numpix4 + 1
               if (skeletons(im1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
               if (skeletons(ip1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
               if (skeletons(im1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
               if (skeletons(ip1,jm1) .gt. almostzero) numpix8 = numpix8 + 1
               if (numpix4 .ge. 2 .and. (numpix8 .eq. 2 .or. numpix8 .eq. 3)) looping = .true.
             endif

             if (cverbose .eq. '-verb2' .and. iotty .gt. 0) then
               ia = nint ( xpixok(k,np(k)) )
               ja = nint ( ypixok(k,np(k)) )
               write (iotty,'(a,5(i4,a))') '   TRACESKELS: Separating branch of skeleton ', k
     &                                   , ' at (', ico, ',', jco, ')-(', ia, ',', ja, ')'
               if (ldebug) then
                 do n=1,l-1
                   ia = nint ( xpixok(k,n) )
                   ja = nint ( ypixok(k,n) )
                   write (iotty,'(a,3(i4,a))') '   TRACESKELS: Skeleton pixel: ', n, ' (', ia, ',', ja, ')'
                 enddo
                 do n=l,np(k)
                   ia = nint ( xpixok(k,n) )
                   ja = nint ( ypixok(k,n) )
                   write (iotty,'(a,3(i4,a))') '   TRACESKELS: Branch pixel: ', n, ' (', ia, ',', ja, ')'
                 enddo
                 write (iotty,'(a$)') '   TRACESKELS: Removing branch pixels:'
               endif
             endif

             if (looping) then
               skeletons(ico,jco) = 0.0d0
               skeletons(im1,jco) = 0.0d0
               skeletons(ip1,jco) = 0.0d0
               skeletons(ico,jm1) = 0.0d0
               skeletons(ico,jp1) = 0.0d0
               skeletons(im1,jm1) = 0.0d0
               skeletons(ip1,jp1) = 0.0d0
               skeletons(im1,jp1) = 0.0d0
               skeletons(ip1,jm1) = 0.0d0
               if (cverbose .eq. '-verb2' .and. iotty .gt. 0) then
                 write (iotty,'(a,3(i4,a))') '   TRACESKELS: Looping skeleton identified: zeroing all 9 pixels at ('
     &                                     , ico, ',', jco, ')'
               endif
             else
               do n=l,min(l+2,np(k))
                 ia = nint ( xpixok(k,n) )
                 ja = nint ( ypixok(k,n) )
                 skeletons(ia,ja) = 0.0d0
                 if (cverbose .eq. '-verb2' .and. iotty .gt. 0) then
                   if (ldebug) then
                     write (iotty,'(a,(i4,a)$)') ' (', ia, ',', ja, ')'
                   endif
                 endif
               enddo
             endif
             if (cverbose .eq. '-verb2' .and. iotty .gt. 0) then
               if (ldebug) then
                 write (iotty,'()')
               endif
             endif
             lnobranches = .false.
           endif
         enddo
       enddo 

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine convert84
 
     &            ( inmx, nx, ny, fun, worker, filename, cverbose, iotty, almostzero, lskeletons )
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!
! To clean clusters, I use TintFill which works for 4-connected clusters only (in my implementation).
! This subroutine converts all 8-connected pixels to 4-connected before calling TintFill.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       lskeletons
       character*(*) filename, cverbose
       integer       i, j, im1, ip1, jm1, jp1, nx, ny, nump4a, nump8a, nump4b, nump8b, numpix4, numpix8, iotty, iz, jz, izm1, izp1
     &             , jzm1, jzp1, inmx
       real*8        almostzero
       real*8        fun(inmx,inmx), worker(nx,ny)
!__________________________________________________________________________________________________________________________________
!
       if (cverbose .eq. '-verb2') then
         write (iotty,'(a)') '   Converting clusters of 8-connected pixels to 4-connected'
       endif
       
       do j=1,ny
         jm1 = max ( j - 1, 1 )
         jp1 = min ( j + 1, ny )
         do i=1,nx
           im1 = max ( i - 1, 1 )
           ip1 = min ( i + 1, nx )

           if (fun(i,j) .gt. almostzero) then

! Up-right diagonal direction.
       
             if (fun(im1,jm1) .gt. almostzero .and. fun(i,jm1) + fun(im1,j) .lt. almostzero ) then
               iz = i
               jz = jm1
               izm1 = max ( iz - 1, 1 )
               izp1 = min ( iz + 1, nx )
               jzm1 = max ( jz - 1, 1 )
               jzp1 = min ( jz + 1, ny )
               nump4a = 0
               nump8a = 0
               if (fun(izm1,jz  ) .gt. almostzero) nump4a = nump4a + 1
               if (fun(izp1,jz  ) .gt. almostzero) nump4a = nump4a + 1
               if (fun(iz  ,jzm1) .gt. almostzero) nump4a = nump4a + 1
               if (fun(iz  ,jzp1) .gt. almostzero) nump4a = nump4a + 1
               if (fun(izm1,jzm1) .gt. almostzero) nump8a = nump8a + 1
               if (fun(izp1,jzp1) .gt. almostzero) nump8a = nump8a + 1
               if (fun(izm1,jzp1) .gt. almostzero) nump8a = nump8a + 1
               if (fun(izp1,jzm1) .gt. almostzero) nump8a = nump8a + 1
               iz = im1
               jz = j
               izm1 = max ( iz - 1, 1 )
               izp1 = min ( iz + 1, nx )
               jzm1 = max ( jz - 1, 1 )
               jzp1 = min ( jz + 1, ny )
               nump4b = 0
               nump8b = 0
               if (fun(izm1,jz  ) .gt. almostzero) nump4b = nump4b + 1
               if (fun(izp1,jz  ) .gt. almostzero) nump4b = nump4b + 1
               if (fun(iz  ,jzm1) .gt. almostzero) nump4b = nump4b + 1
               if (fun(iz  ,jzp1) .gt. almostzero) nump4b = nump4b + 1
               if (fun(izm1,jzm1) .gt. almostzero) nump8b = nump8b + 1
               if (fun(izp1,jzp1) .gt. almostzero) nump8b = nump8b + 1
               if (fun(izm1,jzp1) .gt. almostzero) nump8b = nump8b + 1
               if (fun(izp1,jzm1) .gt. almostzero) nump8b = nump8b + 1
             
               if (nump4a .le. nump4b .and. nump8a .le. nump8b) then
                 worker(i,jm1) = worker(i,j)
               else if (nump4a .gt. nump4b .or. nump8a .gt. nump8b) then
                 worker(im1,j) = worker(i,j)
               else
                 write (iotty,'(/a,6i5)') '   CONVERT84: ERROR: Trouble converting 8- to 4-connected (up-right): '
     &                                  , i, j, nump4a, nump4b, nump8a, nump8b
                 write (iotty,'(a)'     ) '   CONVERT84: Image: '//filename
                 stop 99
               endif
             endif

! Up-left diagonal direction.
           
             if (fun(ip1,jm1) .gt. almostzero .and. fun(i,jm1) + fun(ip1,j) .lt. almostzero ) then
               iz = i
               jz = jm1
               nump4a = 0
               nump8a = 0
               izm1 = max ( iz - 1, 1 )
               izp1 = min ( iz + 1, nx )
               jzm1 = max ( jz - 1, 1 )
               jzp1 = min ( jz + 1, ny )
               if (fun(izm1,jz  ) .gt. almostzero) nump4a = nump4a + 1
               if (fun(izp1,jz  ) .gt. almostzero) nump4a = nump4a + 1
               if (fun(iz  ,jzm1) .gt. almostzero) nump4a = nump4a + 1
               if (fun(iz  ,jzp1) .gt. almostzero) nump4a = nump4a + 1
               if (fun(izm1,jzm1) .gt. almostzero) nump8a = nump8a + 1
               if (fun(izp1,jzp1) .gt. almostzero) nump8a = nump8a + 1
               if (fun(izm1,jzp1) .gt. almostzero) nump8a = nump8a + 1
               if (fun(izp1,jzm1) .gt. almostzero) nump8a = nump8a + 1
               iz = ip1
               jz = j
               izm1 = max ( iz - 1, 1 )
               izp1 = min ( iz + 1, nx )
               jzm1 = max ( jz - 1, 1 )
               jzp1 = min ( jz + 1, ny )
               nump4b = 0
               nump8b = 0
               if (fun(izm1,jz  ) .gt. almostzero) nump4b = nump4b + 1
               if (fun(izp1,jz  ) .gt. almostzero) nump4b = nump4b + 1
               if (fun(iz  ,jzm1) .gt. almostzero) nump4b = nump4b + 1
               if (fun(iz  ,jzp1) .gt. almostzero) nump4b = nump4b + 1
               if (fun(izm1,jzm1) .gt. almostzero) nump8b = nump8b + 1
               if (fun(izp1,jzp1) .gt. almostzero) nump8b = nump8b + 1
               if (fun(izm1,jzp1) .gt. almostzero) nump8b = nump8b + 1
               if (fun(izp1,jzm1) .gt. almostzero) nump8b = nump8b + 1
             
               if (nump4a .le. nump4b .and. nump8a .le. nump8b) then
                 worker(i,jm1) = worker(i,j)
               else if (nump4a .gt. nump4b .or. nump8a .gt. nump8b) then
                 worker(ip1,j) = worker(i,j)
               else
                 write (iotty,'(/a,6i5)') '   CONVERT84: ERROR: Trouble converting 8- to 4-connected (up-left): '
     &                                  , i, j, nump4a, nump4b, nump8a, nump8b
                 write (iotty,'(a)'     ) '   CONVERT84: Image: '//filename
                 stop 99
               endif
             endif
               
! Additionally remove some pixels of groups of four (only for skeletons).
             
             if (lskeletons) then
               numpix4 = 0
               numpix8 = 0
               if (fun(ip1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
               if (fun(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
               if (fun(ip1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
               if (numpix4 .eq. 2 .and. numpix8 .eq. 1) then
                 iz = i
                 jz = jp1
                 izm1 = max ( iz - 1, 1 )
                 izp1 = min ( iz + 1, nx )
                 jzm1 = max ( jz - 1, 1 )
                 jzp1 = min ( jz + 1, ny )
                 numpix4 = 0
                 numpix8 = 0
                 if (fun(izm1,jz  ) .gt. almostzero) numpix4 = numpix4 + 1
                 if (fun(izp1,jz  ) .gt. almostzero) numpix4 = numpix4 + 1
                 if (fun(iz  ,jzm1) .gt. almostzero) numpix4 = numpix4 + 1
                 if (fun(iz  ,jzp1) .gt. almostzero) numpix4 = numpix4 + 1
                 if (fun(izm1,jzm1) .gt. almostzero) numpix8 = numpix8 + 1
                 if (fun(izp1,jzp1) .gt. almostzero) numpix8 = numpix8 + 1
                 if (fun(izm1,jzp1) .gt. almostzero) numpix8 = numpix8 + 1
                 if (fun(izp1,jzm1) .gt. almostzero) numpix8 = numpix8 + 1
               
                 if (numpix4 .eq. 2 .and. numpix8 .eq. 1) then
                   worker(i,jp1) = 0.0d0
                   if (cverbose .eq. '-verb2') then
                     write (iotty,'(a,2i5,4i2)') '   Additionally removed one pixel of four (a):', i, jp1, numpix4, numpix8
                   endif
                 else
                   iz = i
                   jz = j
                   izm1 = max ( iz - 1, 1 )
                   izp1 = min ( iz + 1, nx )
                   jzm1 = max ( jz - 1, 1 )
                   jzp1 = min ( jz + 1, ny )
                   numpix4 = 0
                   numpix8 = 0
                   if (fun(izm1,jz  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (fun(izp1,jz  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (fun(iz  ,jzm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (fun(iz  ,jzp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (fun(izm1,jzm1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (fun(izp1,jzp1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (fun(izm1,jzp1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (fun(izp1,jzm1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (numpix4 .eq. 2 .and. numpix8 .le. 2) then
                     worker(i,j) = 0.0d0
                     if (cverbose .eq. '-verb2') then
                       write (iotty,'(a,2i5,4i2)') '   Additionally removed one pixel of four (b):', i, j, numpix4, numpix8
                     endif
                   endif
                 endif
               endif
               
               numpix4 = 0
               numpix8 = 0
               if (fun(im1,j  ) .gt. almostzero) numpix4 = numpix4 + 1
               if (fun(i  ,jp1) .gt. almostzero) numpix4 = numpix4 + 1
               if (fun(im1,jp1) .gt. almostzero) numpix8 = numpix8 + 1
               if (numpix4 .eq. 2 .and. numpix8 .eq. 1) then
                 iz = i
                 jz = jp1
                 izm1 = max ( iz - 1, 1 )
                 izp1 = min ( iz + 1, nx )
                 jzm1 = max ( jz - 1, 1 )
                 jzp1 = min ( jz + 1, ny )
                 numpix4 = 0
                 numpix8 = 0
                 if (fun(izm1,jz  ) .gt. almostzero) numpix4 = numpix4 + 1
                 if (fun(izp1,jz  ) .gt. almostzero) numpix4 = numpix4 + 1
                 if (fun(iz  ,jzm1) .gt. almostzero) numpix4 = numpix4 + 1
                 if (fun(iz  ,jzp1) .gt. almostzero) numpix4 = numpix4 + 1
                 if (fun(izm1,jzm1) .gt. almostzero) numpix8 = numpix8 + 1
                 if (fun(izp1,jzp1) .gt. almostzero) numpix8 = numpix8 + 1
                 if (fun(izm1,jzp1) .gt. almostzero) numpix8 = numpix8 + 1
                 if (fun(izp1,jzm1) .gt. almostzero) numpix8 = numpix8 + 1
               
                 if (numpix4 .eq. 2 .and. numpix8 .eq. 1) then
                   worker(i,jp1) = 0.0d0
                   if (cverbose .eq. '-verb2') then
                     write (iotty,'(a,2i5,4i2)') '   Additionally removed one pixel of four (c):', i, jp1, numpix4, numpix8
                   endif
                 else
                   iz = i
                   jz = j
                   izm1 = max ( iz - 1, 1 )
                   izp1 = min ( iz + 1, nx )
                   jzm1 = max ( jz - 1, 1 )
                   jzp1 = min ( jz + 1, ny )
                   numpix4 = 0
                   numpix8 = 0
                   if (fun(izm1,jz  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (fun(izp1,jz  ) .gt. almostzero) numpix4 = numpix4 + 1
                   if (fun(iz  ,jzm1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (fun(iz  ,jzp1) .gt. almostzero) numpix4 = numpix4 + 1
                   if (fun(izm1,jzm1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (fun(izp1,jzp1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (fun(izm1,jzp1) .gt. almostzero) numpix8 = numpix8 + 1
                   if (fun(izp1,jzm1) .gt. almostzero) numpix8 = numpix8 + 1
               
                   if (numpix4 .eq. 2 .and. numpix8 .le. 2) then
                     worker(i,j) = 0.0d0
                     if (cverbose .eq. '-verb2') then
                       write (iotty,'(a,2i5,4i2)') '   Additionally removed one pixel of four (d):', i, j, numpix4, numpix8
                     endif
                   endif
                 endif
               endif
             endif
           endif
         enddo
       enddo
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      
       function selectk

     &          ( k, n, arr )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN. Recoded in a better language by A. Men'shchikov. Used to compute median value.
!
! Returns the k-th smallest value in the array arr(1:n). The input array will be rearranged to have this value in location arr(k),
! with all smaller elements moved to arr(1:k-1) (in arbitrary order) and all larger elements in arr(k+1..n) (in arbitrary order).
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       k, n, i, ir, j, l, mid, l1, irl
       real*8        selectk, arr(n), a, temp
!__________________________________________________________________________________________________________________________________
!
       l = 1
       ir = n
 1     continue
       irl = ir - l
       if (irl .le. 1) then           ! Active partition contains 1 or 2 elements.
         if (irl .eq. 1) then         ! Active partition contains 2 elements.
           if (arr(ir) .lt. arr(l)) then
             temp = arr(l)
             arr(l) = arr(ir)
             arr(ir) = temp
           endif
         endif
         selectk = arr(k)
         goto 9
       else
         l1 = l + 1
         mid = (l + ir) / 2           ! Choose median of left, center, and right elements as partitioning element a.
         temp = arr(mid)              ! Also rearrange so that arr(l) <= arr(l+1), arr(ir) >= arr(l+1).
         arr(mid) = arr(l1)
         arr(l1) = temp
         if (arr(l) .gt. arr(ir)) then
           temp = arr(l)
           arr(l) = arr(ir)
           arr(ir) = temp
         endif
         if (arr(l1) .gt. arr(ir)) then
           temp = arr(l1)
           arr(l1) = arr(ir)
           arr(ir) = temp
         endif
         if (arr(l) .gt. arr(l1)) then
           temp = arr(l)
           arr(l) = arr(l1)
           arr(l1) = temp
         endif

! Initialize pointers for partitioning.

         a = arr(l1)                ! Partitioning element.
         i = l1
         j = ir
 3       continue
         i = i + 1
         do while (arr(i) .lt. a)   ! Scan up to find element > a.
           i = i + 1
         enddo
         j = j - 1
         do while (arr(j) .gt. a)   ! Scan down to find element < a.
           j = j - 1
         enddo
         if (j .ge. i) then         ! Pointers crossed. Exit with partitioning complete.
           temp = arr(i)  
           arr(i) = arr(j)          ! Exchange elements.
           arr(j) = temp
           goto 3
         endif

! Pointers crossed. Exit with partitioning complete.

         arr(l1) = arr(j)           ! Insert partitioning element.
         arr(j) = a
         if (j .ge. k) ir = j - 1   ! Keep active the partition that contains the k-th element.
         if (j .le. k) l = i
       endif
       goto 1

 9     return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine fftconvol
       
     &            ( iotty, nxo, nyo, dx, dy, fun, beam, expo, nxp, nyp, psf, fract, psfbeam, cdowhat, ckerntype, lnoimage
     &            , lnormalize, lnoscaling, cverbose )
!__________________________________________________________________________________________________________________________________
!
! FFT convolution of an image, encapsulated in a subroutine.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       lnoimage, lnormalize, lnoscaling
       character*2   cfix
       character*6   cverbose, cnxo, cnyo, cnx, cny, cnxout, cnyout
       character*8   ckerntype, cdowhat

       integer       i, j, irc, nx, ny, ij, nxo, nyo, iotty, nxout, nyout, mxo, myo, nx0, ny0, inxo, inyo, inx, iny, inxout, inyout
     &             , imin, imax, jmin, jmax, insp, jn, jm, jp, im, ip, nx2, ny2, ixo, iyo, nxmx, nxp, nyp, i0, j0, ii, jj

       real*8        xout, yout, xfnyq, yfnyq, dx, dy, dxf, dyf, totfx1, totfx2, dilute, hwhm, rad2, argex, beamx, fnorm, beam
     &             , funmin, funmax, almostzero, constant, expo, correct, sigm2g, sigm2m, delx2, dely2, fract, expargex, omexpargex
     &             , psfbeam, fun(nxo,nyo), psf(nxp,nyp)
       
       real*8    , allocatable :: funi(:), funikern(:), funx(:), fuu(:,:), funn(:,:), kern(:,:)
       complex*16, allocatable :: speq(:), speqkern(:), spec(:), speckern(:)

       parameter   ( almostzero = 1.0d-20 )
       
       external      rlft3, equivalent, equivmult
!__________________________________________________________________________________________________________________________________
!
! Find powers of 2 in X and Y coordinates to use as array dimensions for the FFT.

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
     
       nxmx = max ( nx, ny, 1 )
       insp = nxmx + 1
       
       allocate ( funn(insp,insp), kern(insp,insp), stat=irc )

       if (irc .ne. 0) then
         write (iotty,'(/a,2i5)') '   FFTCONVOL: ERROR: Trouble allocating memory (10) for 2 arrays of sizes: ', insp, insp
         stop 10
       endif
       
       do j=1,insp
         do i=1,insp
           funn(i,j) = 0.0d0
           kern(i,j) = 0.0d0
         enddo
       enddo

! Compute total flux under the image for rescaling after the convolution. 
! The offsets 'nx0' and 'ny0' below are to avoid border effects in convolution and thus to reduce the error in rescaling map.
         
       nx0 = nint ( 2.0d0 * beam / dx )
       ny0 = nint ( 2.0d0 * beam / dy )

       imin = max ( min ( 1 + nx0, nxo / 2 - 10 ), 1 )
       jmin = max ( min ( 1 + ny0, nyo / 2 - 10 ), 1 )
       imax = min ( max ( nxo - nx0, nxo / 2 + 11 ), nxo )
       jmax = min ( max ( nyo - ny0, nyo / 2 + 11 ), nyo )

!!! ABM: replaced on 2022-11-20
!!       imin = 2
!!       jmin = 2
!!       imax = nxo - 1
!!       jmax = nyo - 1

       fnorm = 1.0d0
       constant = 0.0d0
       funmax =-1.0d+30
       funmin = 1.0d+30
                                                             
       if (.not.lnoimage) then
         do j=jmin,jmax
           do i=imin,imax
             if (.not.isnan ( fun(i,j) )) then  
               funmax = max ( funmax, fun(i,j) )
               funmin = min ( funmin, fun(i,j) )
             endif
           enddo
         enddo
         if (abs ( funmax - funmin ) .lt. almostzero) cdowhat = 'constant'
       endif

       do j=1,nyo
         do i=1,nxo
           if (isnan ( fun(i,j) )) then  
             funn(i,j) = 0.0d0
           else
             funn(i,j) = fun(i,j)
           endif
         enddo
       enddo

       if (cdowhat .eq. 'constant') then
         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   NOTE: The input image is constant; no transform will be performed.'
       else

! Highest (Nyquist) frequency of the 2D sample.

         xfnyq = 1.0d0 / (2.0d0 * dx)
         yfnyq = 1.0d0 / (2.0d0 * dy)

! Find the lowest frequency.

         dxf = 2.0d0 * xfnyq / nx
         dyf = 2.0d0 * yfnyq / ny

! Define the upper limit of frequency point numbers for output.

         xout = 0.0d0
         yout = 0.0d0
         nxout = nint ( 2.0d0 * xout / dxf )
         nyout = nint ( 2.0d0 * yout / dyf )
         if (nxout .eq. 0 .or. nxout .gt. nx + 1) nxout = nx + 1
         if (nyout .eq. 0 .or. nyout .gt. ny + 1) nyout = ny + 1

! The number of visibility points should be odd for output.

         if (mod ( nxout, 2 ) .eq. 0) nxout = min ( nxout + 1, nx + 1 )
         if (mod ( nyout, 2 ) .eq. 0) nyout = min ( nyout + 1, ny + 1 )

         if (lnoimage) then
           nx = nxo
           ny = nyo
           nxout = nxo
           nyout = nyo
         endif

         inxo = int ( log10 ( dble ( nxo ) ) ) + 1
         write (cnxo,'(i5)') nxo
         inyo = int ( log10 ( dble ( nyo ) ) ) + 1
         write (cnyo,'(i5)') nyo
         inx = int ( log10 ( dble ( nx ) ) ) + 1
         write (cnx,'(i5)') nx
         iny = int ( log10 ( dble ( ny ) ) ) + 1
         write (cny,'(i5)') ny
         inxout = int ( log10 ( dble ( nxout ) ) ) + 1
         write (cnxout,'(i5)') nxout
         inyout = int ( log10 ( dble ( nyout ) ) ) + 1
         write (cnyout,'(i5)') nyout

         if (cverbose .eq. '-verb2') write (iotty,'(2x,a)') ' Original: '//cnxo(5-inxo+1:5)//' '//cnyo(5-inyo+1:5)
     &                                                    //' Transform: '//cnx(5-inx+1:5)//' '//cny(5-iny+1:5)
     &                                                    //' Output: '//cnxout(5-inxout+1:5)//' '//cnyout(5-inyout+1:5)
         totfx1 = 1.0d0

! Compute total flux under the input image.

         if (.not.lnoimage) then
           
! Subtract the minimum value of the input image to make the image entirely positive.
! This is important for accurate re-normalization of the total flux under the image.
! Fix on 200206: moved subtraction of constant down to always subtract minimum, except for noscaling.

           if (lnoscaling) then
             totfx1 = 1.0d0
           else
             constant = funmin
             do j=1,ny
               do i=1,nx
                 funn(i,j) = funn(i,j) - constant
               enddo
             enddo
             if (lnormalize) then
               totfx1 = 0.0d0
               do j=jmin,jmax
                 do i=imin,imax
                   totfx1 = max ( totfx1, funn(i,j) )
                 enddo
               enddo
             else
               totfx1 = 0.0d0
               do j=jmin,jmax
                 do i=imin,imax
                   totfx1 = totfx1 + funn(i,j)
                 enddo
               enddo
             endif
           endif
         endif
         
! Compute a convolution kernel.

         if (cdowhat .eq. 'convolve') then
     
           if (cverbose .eq. '-verb2') then
             write (cfix,'(i2)') int ( log10 ( max ( beam, 1.0d0 ) ) ) + 10
             if (ckerntype .eq. 'moffatfn')
     &           write (iotty,'(3x,a,f'//cfix//'.7,2(a,f5.2),a)') 'Preparing a', beam
     &                                                          , ' arcsec Gauss-Moffat kernel (fit=',fract, ' expo=',expo,')'
             if (ckerntype .eq. 'powerlaw')
     &           write (iotty,'(3x,a,2(f5.2,a))') 'Preparing a power-law kernel (expo=', expo, ', factor=', beam, ')'
             if (ckerntype .eq. 'gaussian') 
     &           write (iotty,'(3x,a,f'//cfix//'.7,a)') 'Preparing', beam, ' arcsec Gaussian kernel'
             if (ckerntype .eq. 'exponent')
     &           write (iotty,'(3x,a,f'//cfix//'.7,a)') 'Preparing', beam, ' arcsec exponential kernel'
             if (ckerntype .eq. 'cylinder') 
     &           write (iotty,'(3x,a,f'//cfix//'.7,a)') 'Preparing', beam, ' arcsec cylindrical kernel'
             if (ckerntype .eq. 'psfimage') then
               write (cfix,'(i2)') int ( log10 ( max ( psfbeam, 1.0d0 ) ) ) + 10
               write (iotty,'(3x,a,f'//cfix//'.7,a)') 'Using external', psfbeam, ' arcsec convolution kernel'
             endif
           endif

! sigma, parameter of a gaussian, is ~20% smaller (by factor sqrt((log 4)) = 1.177) than its HWHM - HALF width at half maximum.
! if one takes sigma defined in this way, one gets correct value 1/2 (half maximum) of the gaussian at x = HWHM

           hwhm = beam / 2.0d0
           sigm2g = hwhm**2 / log ( 4.0d0 )
           if (ckerntype .eq. 'moffatfn') then
             correct = (2.0d0**(1.0d0 / dble ( expo )) - 1.0d0)
             sigm2m = hwhm**2
           endif
           beamx = 0.0d0
           i0 = nxo / 2 - nxp / 2
           j0 = nyo / 2 - nyp / 2
           
           do j=1,nyo
             jj = j - (nyo / 2 + 1)
             dely2 = (dble ( jj ) * dy)**2
             do i=1,nxo
               ii = i - (nxo / 2 + 1)
               delx2 = (dble ( ii ) * dx)**2
               rad2  = delx2 + dely2
               if (ckerntype .eq. 'gaussian') then
                 argex = min ( rad2 / (2.0d0 * sigm2g), 50.0d0 )
                 kern(i,j) = exp ( -argex )
               endif
               if (ckerntype .eq. 'moffatfn') then
                 argex = rad2 / sigm2m
                 beamx = 1.0d0 / (1.0d0 + correct * argex)**expo
                 argex = min ( rad2 / (2.0d0 * sigm2g), 50.0d0 )
                 expargex = exp ( -argex )
                 omexpargex = 1.0d0 - expargex
                 kern(i,j) = max ( kern(i,j), fract * omexpargex * beamx + (1.0d0 - fract * omexpargex) * expargex )
               endif
               if (ckerntype .eq. 'powerlaw') then
                 beamx = (sqrt ( max ( rad2, almostzero ) / dx))**expo / beam
                 kern(i,j) = min ( beamx, 100.0d0 )
               endif
               if (ckerntype .eq. 'exponent') then
                 argex = min ( rad2, 100.0d0 )
                 kern(i,j) = exp ( -argex )
               endif
               if (ckerntype .eq. 'cylinder') then
                 if (sqrt ( rad2 ) .le. hwhm) then
                   kern(i,j) = 1.0d0
                 else
                   kern(i,j) = 0.0d0
                 endif
               endif
               if (ckerntype .eq. 'psfimage') then
                 kern(i,j) = psf(min(max(i-i0,1),nxp),min(max(j-j0,1),nyp))
               endif
               if (lnoimage) fun(i,j) = kern(i,j)
             enddo
           enddo
         endif
         
         if (lnoimage) then
           deallocate ( funn, kern )
           return
         endif

         if (cdowhat .eq. 'convolve' .or. cdowhat .eq. 'fftransf') then

! Initialize FFT arrays with zeros.
               
           allocate ( speq(insp), speqkern(insp), funi(insp*insp), funikern(insp*insp), stat=irc )
          
           if (irc .ne. 0) then
             write (iotty,'(/a,2i5)') '   FFTCONVOL: ERROR: Trouble allocating memory (20) for 4 arrays of sizes: ', insp, insp
             stop 20
           endif
          
           do j=1,insp*insp
             funi(j) = 0.0d0
             funikern(j) = 0.0d0
           enddo

           fnorm = dble ( nx * ny ) / 2.0d0
           do j=1,ny
             jn = (j - 1) * nx
             do i=1,nx
               ij = i + jn
               funi(ij) = funn(i,j)
               funikern(ij) = kern(i,j)
             enddo
           enddo
         endif
         
! Fast Fourier Transform to Foruier space.

         if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing 2D fast Fourier transform of the image'
       
         call rlft3 ( funi, speq, nx, ny, 1, 1 )

! Reconstruct the power spectrum. Magnitude of the visibility is given by abs ( F ) or sqrt ( real(F)**2 + imag(F)**2 ),
! and phase is given by atan ( imag(F) / real(F) ). Power spectrum is F**2.

         if (cdowhat .eq. 'fftransf') then
     
           allocate ( spec(insp/2*insp), speckern(insp/2*insp), stat=irc )
     
           if (irc .ne. 0) then
             write (iotty,'(/a,2i5)') '   FFTCONVOL: ERROR: Trouble allocating memory (21) for 2 arrays of sizes: ', insp/2, insp
             stop 21
           endif
           if (cdowhat .eq. 'fftransf') then
             do j=1,insp/2*insp
               spec(j) = ( 0.0d0, 0.0d0 )
               speckern(j) = ( 0.0d0, 0.0d0 )
             enddo
           endif
     
           call equivalent ( insp/2*insp, funi, spec )
     
           nx2 = nx / 2
           ny2 = ny / 2
           do j=1,ny2+1
             jn = (j - 1) * nx2
             jp = j + ny2
             do i=1,nx2
               ij = i + jn
               ip = i + nx2
               funn(ip,jp) = abs ( spec(ij) )
             enddo
             funn(nx+1,jp) = abs ( speq(j) )
           enddo
           do j=ny2+1,ny
             jn = (j - 1) * nx2
             jm = j - ny2
             do i=1,nx2
               ij = i + jn
               ip = i + nx2
               funn(ip,jm) = abs ( spec(ij) )
             enddo
             funn(nx+1,jm) = abs ( speq(j) )
           enddo

! Symmetry is used to fill in the other half of the plane in X.

           do j=1,ny2
             jm = ny2 - j + 1
             jp = ny2 + j + 1
             do i=1,nx2
               im = nx2 - i + 1
               ip = nx2 + i + 1
               funn(im,jm) = funn(ip,jp)
             enddo
           enddo
           do j=ny2+1,ny+1
             jm = ny - j + 2
             do i=1,nx2
               im = nx2 - i + 1
               ip = nx2 + i + 1
               funn(im,j) = funn(ip,jm)
             enddo
           enddo
           do j=1,nyo
             do i=1,nxo
               fun(i,j) = funn(i,j)
             enddo
           enddo

           deallocate ( spec, speckern )
         else

! Fast Fourier Transform to Fourier space.

           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing 2D fast Fourier transform of the beam'

           call rlft3 ( funikern, speqkern, nx, ny, 1, 1 )

! Now convolve (multiply in Fourier space).

           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Convolving the image by multiplying the transforms'
       
           allocate ( funx(insp*insp), stat=irc )
       
           if (irc .ne. 0) then
             write (iotty,'(/a,2i5)') '   FFTCONVOL: ERROR: Trouble allocating memory (25) for array sizes: ', insp, insp
             stop 25
           endif

           call equivmult ( insp, funi, funikern, funx )

           do j=1,insp
             speq(j) = speq(j) * speqkern(j)
           enddo
       
           deallocate ( funi, funikern )

! Fast Fourier Transform back to image space.

           if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Computing 2D fast Fourier transform back to images'
       
           call rlft3 ( funx, speq, nx, ny, 1, -1 )
       
           deallocate ( speq, speqkern )
         
           if (cverbose .eq. '-verb2') write (iotty,'(a,$)') '   Unscrambling transform'
         
           allocate ( fuu(insp,insp), stat=irc )
         
           if (irc .ne. 0) then
             write (iotty,'(/a,2i5)') '   FFTCONVOL: ERROR: Trouble allocating memory (30) for 2 arrays of sizes: ', insp, insp
             stop 30
           endif
       
           do j=1,ny
             jn = (j - 1) * nx
             do i=1,nx
               ij = i + jn
               fuu(i,j) = funx(ij) / fnorm
             enddo
           enddo
       
           deallocate ( funx )

! Further unscramble the output.
                    
           mxo = nxo / 2
           myo = nyo / 2
           do j=1,ny-myo
             jp = j + myo
             do i=1,nx-mxo
               ip = i + mxo
               funn(i,j) = fuu(ip,jp)
             enddo
           enddo
           do j=1,nyo-(ny-myo)
             jp = j + ny - myo
             do i=1,nxo-(nx-mxo)
               ip = i + nx - mxo
               funn(ip,jp) = fuu(i,j)
             enddo
           enddo
           do j=1,ny-myo
             jp = j + myo
             do i=1,nxo-(nx-mxo)
               ip = i + nx - mxo
               funn(ip,j) = fuu(i,jp)
             enddo
           enddo
           do j=1,nyo-(ny-myo)
             jp = j + ny - myo
             do i=1,nx-mxo
               ip = i + mxo
               funn(i,jp) = fuu(ip,j)
             enddo
           enddo
       
           deallocate ( fuu )

           do j=1,nyo
             do i=1,nxo
               fun(i,j) = funn(i,j)
             enddo
           enddo

! Compute total flux under the convolved image.

           if (lnoscaling) then
             totfx2 = 1.0d0
           else
             if (lnormalize) then
               totfx2 = 0.0d0
               do j=jmin,jmax
                 do i=imin,imax
                   totfx2 = max ( totfx2, fun(i,j) )
                 enddo
               enddo
             else
               totfx2 = 0.0d0
               do j=jmin,jmax
                 do i=imin,imax
                   totfx2 = totfx2 + fun(i,j)
                 enddo
               enddo
!!               write (*,*) totfx2
             endif
           endif
           dilute = totfx1 / totfx2
           
!!           write (*,*) dilute

           if (cverbose .eq. '-verb2') write (iotty,'(a,1pe10.3)') ' scaling factor:', dilute
           if (abs ( dilute ) .lt. almostzero) then
             write (iotty,'(/a,1pe11.3)') '   FFTCONVOL: ERROR: Scaling factor is suspicious:', dilute
             stop 99
           endif
           
! Rescale the convolved image in order to conserve the total flux.

           do j=1,nyo
             do i=1,nxo
               fun(i,j) = fun(i,j) * dilute + constant
             enddo
           enddo
         endif
       endif

       deallocate ( funn, kern )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine equivmult
       
     &            ( insp, spec, speckern, funx )
!__________________________________________________________________________________________________________________________________
!
! Replacement for the EQUIVALENCE statement, allowing ALLOCATED arrays.
!__________________________________________________________________________________________________________________________________
!
       implicit none
       integer       j, insp
       complex*16    spec(insp/2*insp), speckern(insp/2*insp), funx(insp*insp)
!__________________________________________________________________________________________________________________________________
!
       do j=1,insp/2*insp
         funx(j) = spec(j) * speckern(j)
       enddo
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine equivalent
       
     &            ( ndim, infun, outfun )
!__________________________________________________________________________________________________________________________________
!
! Replacement for the EQUIVALENCE statement, allowing ALLOCATED arrays.
!__________________________________________________________________________________________________________________________________
!
       implicit none
       integer       j, ndim
       complex*16    infun(ndim), outfun(ndim)
!__________________________________________________________________________________________________________________________________
!
       do j=1,ndim
         outfun(j) = infun(j)
       enddo
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine rlft3 
       
     &            ( data, speq, nn1, nn2, nn3, isign )
!__________________________________________________________________________________________________________________________________
!
! Fast Fourier Transform algorithm RLFT3 from Press, Teukolsky, Vetterling, & Flannery (1992):
! "Numerical recipes in FORTRAN. The art of scientific computing", Cambridge University Press, 2nd edition.
!
! Given a two- or three- dimensional real array 'data' whose dimensions are nn1, nn2, nn3 (where nn3 is 1 for the case of a
! two-dimensional array), this routine returns (for isign=1) the complex fast fourier transform as two complex arrays: On output,
! 'data' contains the zero and positive frequency values for the first frequency component, while 'speq' contains the Nyquist
! critical frequency values of the first frequency component. Second (and third) frequency components are store for zero, positive,
! and negative frequencies, in standard wrap-around order. For isign=-1, the inverse transform (times nn1*nn2*nn3/2 as a constant
! multiplicative factor) is performed, with output 'data' (viewed as a real array) deriving from input 'data' (viewed as complex
! and 'speq'. The dimensions nn1, nn2, nn3 must always be integer powers of 2.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       isign, nn1, nn2, nn3

! Note that 'data' is dimensioned as complex, its output format.

       complex*16    data(nn1/2,nn2,nn3), speq(nn2,nn3)
       integer       i1, i2, i3, j1, j2, j3, nn(3)
       real*8        theta, wi, wpi, wpr, wr, wtemp
       complex*16    c1, c2, h1, h2, w, wc
       external      fourn
!__________________________________________________________________________________________________________________________________
!
       c1 = dcmplx ( 0.5d0, 0.0d0 )
       c2 = dcmplx ( 0.0d0, -0.5d0 * isign )
       theta = 6.28318530717959d0 / dble ( isign * nn1 )
       wpr = -2.0d0 * sin ( 0.5d0 * theta )**2
       wpi = sin ( theta )
       nn(1) = nn1 / 2
       nn(2) = nn2
       nn(3) = nn3

! Case of forward transform.

       if (isign .eq. 1) then

! Here is where most of the computer time is spent.

         call fourn ( data, nn, 3, isign )

! Extend data periodically into 'speq'.

         do i3=1,nn3
           do i2=1,nn2
             speq(i2,i3) = data(1,i2,i3)
           enddo
         enddo
       endif

       do i3=1,nn3

! Zero frequency is its own reflection; otherwise locate corresponding negative frequency in wrap-around order.

         j3 = 1
         if (i3 .ne. 1) j3 = nn3 - i3 + 2

! Initialize trigonometric recurrence.

         wr = 1.0d0
         wi = 0.0d0

         do i1=1,nn1/4+1

           j1 = nn1 / 2 - i1 + 2

           do i2=1,nn2
             j2 = 1
             if (i2 .ne. 1) j2 = nn2 - i2 + 2

             if (i1 .eq. 1) then            
               wc = conjg ( speq(j2,j3) )
               h1 = c1 * (data(1,i2,i3) + wc)
               h2 = c2 * (data(1,i2,i3) - wc)
               data(1,i2,i3) = h1 + h2
               speq(j2,j3) = conjg ( h1 - h2 )
             else                            
               wc = conjg ( data(j1,j2,j3) )
               h1 = c1 * (data(i1,i2,i3) + wc)
               h2 = c2 * (data(i1,i2,i3) - wc)
               data(i1,i2,i3) = h1 + w * h2
               data(j1,j2,j3) = conjg ( h1 - w * h2 )
             endif
           enddo

! Do the recurrence.

           wtemp = wr
           wr = wr * wpr - wi * wpi + wr
           wi = wi * wpr + wtemp * wpi + wi
           w = dcmplx ( dble ( wr ), dble ( wi ) )
         enddo
       enddo
                         
! Case of reverse transform.

       if (isign .eq. -1) then

         call fourn ( data, nn, 3, isign )

       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine fourn 
       
     &            ( data, nn, ndim, isign )
!__________________________________________________________________________________________________________________________________
!
! Fast Fourier Transform algorithm RLFT3 from Press, Teukolsky, Vetterling, & Flannery (1992):
! "Numerical recipes in FORTRAN. The art of scientific computing", Cambridge University Press, 2nd edition.
!
! -132 -O3 -fp-model precise -fp-speculation off -warn all -watch source -traceback -assume noold_unit_star -heap-arrays \
! -parallel -module modules
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       isign, ndim, nn(ndim)
       real*8        data(*)
       integer       i1, i2, i2rev, i3, i3rev, ibit, idim, ifp1, ifp2, ip1, ip2, ip3, k1, k2, n, nprev, nrem, ntot
       real*8        tempi, tempr, theta, wi, wpi, wpr, wr, wtemp
!__________________________________________________________________________________________________________________________________
!
       ntot = 1
       do idim=1,ndim
         ntot = ntot * nn(idim)
       enddo
       nprev = 1

       do idim=1,ndim
         n     = nn(idim)
         nrem  = ntot / (n * nprev)
         ip1   = 2 * nprev
         ip2   = ip1 * n
         ip3   = ip2 * nrem
         i2rev = 1

         do i2=1,ip2,ip1
           if (i2 .lt. i2rev) then
             do i1=i2,i2+ip1-2,2
               do i3=i1,ip3,ip2
                 i3rev = i2rev + i3 - i2
                 tempr = data(i3  )
                 tempi = data(i3+1)
                 data(i3       ) = data(i3rev    )
                 data(i3 + 1   ) = data(i3rev + 1)
                 data(i3rev    ) = tempr
                 data(i3rev + 1) = tempi
               enddo
             enddo
           endif
           ibit = ip2 / 2
 1        continue
          if ((ibit .ge. ip1) .and. (i2rev .gt. ibit)) then
             i2rev = i2rev - ibit
             ibit  = ibit / 2
             goto 1
           endif
           i2rev = i2rev + ibit
         enddo
         
         ifp1 = ip1

 2       continue        
         if (ifp1 .lt. ip2) then
           ifp2  = 2 * ifp1
           theta = isign * 6.28318530717959d0 / (ifp2 / ip1)
           wpr   = -2.d0 * sin ( 0.5d0 * theta )**2
           wpi   = sin ( theta )
           wr    = 1.d0
           wi    = 0.d0
           do i3=1,ifp1,ip1
             do i1=i3,i3+ip1-2,2
               do i2=i1,ip3,ifp2
                 k1 = i2
                 k2 = k1 + ifp1
                 tempr = dble ( wr ) * data(k2  ) - dble ( wi ) * data(k2+1)
                 tempi = dble ( wr ) * data(k2+1) + dble ( wi ) * data(k2  )
                 data(k2  ) = data(k1  ) - tempr
                 data(k2+1) = data(k1+1) - tempi
                 data(k1  ) = data(k1  ) + tempr
                 data(k1+1) = data(k1+1) + tempi
               enddo
             enddo
             wtemp = wr
             wr    = wr * wpr - wi * wpi + wr
             wi    = wi * wpr + wtemp * wpi + wi
           enddo
           ifp1 = ifp2
           goto 2
         endif
         nprev = n * nprev
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine histogrm
       
     &            ( nf, nfun, ndim, nbinmx, xn, fun1, fun2, ibeg, iend, xtype )
!__________________________________________________________________________________________________________________________________
!
! Convert functions [nf..nfun] into histograms (for linear or logarithmic x-axis).
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) xtype
       integer       m, i, j, l, k, nf, nfun, ndim, nbinmx
       integer       ibeg(nfun), iend(nfun)
       real*8        fact, shft, xn(nbinmx,nfun), fun1(nbinmx,nfun,ndim), fun2(nbinmx,nfun,ndim)
!__________________________________________________________________________________________________________________________________
!
       shft = 0.0d0
       fact = 1.0d0
       do m=nf,nfun
         if (2 * iend(m) + 1 .gt. nbinmx) then
           return
         endif
         if (xtype .eq. 'lin') then
           shft = (xn(ibeg(m)+1,m) - xn(ibeg(m),m)) / 2.0d0
         else
           fact = sqrt ( xn(ibeg(m)+1,m) / xn(ibeg(m),m) )
         endif
         do i=iend(m),ibeg(m),-1
           l = 2 * i
           k = min ( i+1, iend(m) )
           xn (l  ,m) = xn (k,m)
           xn (l-1,m) = xn (i,m)
           do j=1,ndim
             fun1(l  ,m,j) = fun1(i,m,j)
             fun1(l-1,m,j) = fun1(i,m,j)
             fun2(l  ,m,j) = fun2(i,m,j)
             fun2(l-1,m,j) = fun2(i,m,j)
           enddo
         enddo
         iend(m) = 2 * iend(m)
         do j=1,ndim
           fun1(iend(m),m,j) = fun1(iend(m)-1,m,j)
           fun2(iend(m),m,j) = fun2(iend(m)-1,m,j)
         enddo
!!         iend(m) = iend(m) + 1
         if (xtype .eq. 'lin') then
           xn(iend(m),m) = xn(iend(m)-2,m) + 2.0d0 * shft
         else
           xn(iend(m),m) = xn(iend(m)-2,m) * fact**2
         endif
         do l=ibeg(m),iend(m)
           if (xtype .eq. 'lin') then
             xn(l,m) = xn(l,m) - shft
           else
             xn(l,m) = xn(l,m) / fact
           endif
         enddo
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine isleep
       
     &            ( iseconds )
!__________________________________________________________________________________________________________________________________
!
! For portability, I make my own sleep function.
!__________________________________________________________________________________________________________________________________
!
!!!not4gfortran       !use ifport, only: sleep
      !use ifport, only: sleep
       implicit      none
       integer       iseconds
      !external      sleep
!__________________________________________________________________________________________________________________________________
!
       call sleep ( iseconds )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine ps2pdf 
       
     &            ( psfile )
!__________________________________________________________________________________________________________________________________
!
! Converts PS file to PDF format by calling macos PS2PDF utility.
!__________________________________________________________________________________________________________________________________
!
!!!not4gfortran       !use ifport, only: system, ierrno
      !use ifport, only: system, ierrno
       implicit      none
       character*(*) psfile
       integer       returncode, ilc, lastc    !!, errnumber
       external      lastc
!__________________________________________________________________________________________________________________________________
!
       ilc = lastc ( psfile )
       returncode = system ( 'ps2pdf -sPAPERSIZE=a4 '//psfile//' '//psfile(1:ilc-3)//'.pdf' )

!!       returncode = system ( 'pstopdf '//psfile//' -o '//psfile(1:ilc-3)//'.pdf 2> /dev/null' )
       
       if (returncode .eq. 32512) then
         write (*,'(/a)') '  PS2PDF: System command (PS2PDF) not found.'
       endif
       
       if (returncode .eq. 0) then
         returncode = system ( '\rm '//psfile )
       else
         write (*,'(a,i5)') '  PS2PDF: Cannot convert the PS file to PDF - rc: ', returncode
!!         stop 999
       endif
         
!!!not4gfortran       if (returncode .eq. -1) then
!!!not4gfortran         errnumber = ierrno ( )
!!!not4gfortran         write (*,'(a,i3,a)') '   PStoPDF: Error converting '//psfile//' to PDF (error:', errnumber, ')'
!!!not4gfortran       else
!!!not4gfortran         returncode = system ( '\rm '//psfile )
!!!not4gfortran         if (returncode .eq. -1) then
!!!not4gfortran           errnumber = ierrno ( )
!!!not4gfortran           write (*,'(a,i3,a)') '   PStoPDF: Error deleting '//psfile//' (error:', errnumber, ')'
!!!not4gfortran         endif
!!!not4gfortran       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine trimpdf 
       
     &            ( pdffile )
!__________________________________________________________________________________________________________________________________
!
! Trims PDF file of wite margins by calling GhostScript's script (assuming Ghostscript installed and the script exists).
!__________________________________________________________________________________________________________________________________
!
!!!not4gfortran       !use ifport, only: system, ierrno
      !use ifport, only: system, ierrno
       implicit      none
       character*(*) pdffile
       integer       returncode, errnumber, ilc, lastc
       external      lastc
!__________________________________________________________________________________________________________________________________
!
       ilc = lastc ( pdffile )
       returncode = system ( 'trimpdf '//pdffile )
       if (returncode .eq. -1) then
         errnumber = ierrno ( )
         write (*,'(a,i3,a)') '   TRIMPDF: Error trimming PDF file '//pdffile//' (error:', errnumber, ')'
       endif
       if (returncode .ne. 0) stop 999
         
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine makedir 
       
     &            ( dirname )
!__________________________________________________________________________________________________________________________________
!__________________________________________________________________________________________________________________________________
!
!!!not4gfortran       !use ifport, only: system, ierrno
      !use ifport, only: system, ierrno
       implicit      none
       character*(*) dirname
       integer       returncode, errnumber, ilc, lastc
       external      lastc
!__________________________________________________________________________________________________________________________________
!
       ilc = lastc ( dirname )
       returncode = system ( '\mkdir -p '//dirname )
       if (returncode .eq. -1) then
         errnumber = ierrno ( )
         write (*,'(a,i3,a)') '   MAKEDIR: Error creating directory '//dirname//' (error:', errnumber, ')'
       endif
       if (returncode .ne. 0) stop 999
         
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine fitfun
       
     &            ( iotty, iolog, ma, npt, x, y, sig, a, ia, covar, alpha, chisq, iter, lgood, lconv, lsing, cmodel, modelfun
     &            , cwhat )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       lsingular, lgood, lconv, lsing
       character*1   cma, crc, crcx
       character*2   cnpt
       character*3   citer
       character*9   chead(8)
       character*(*) cmodel, cwhat
       integer       i, ma, maa, ia(ma), iconverging, iter, npt, lastc, iotty, iolog, itermax, nit, nh, nhm3, iverbose
       real*8        x(npt), y(npt), sig(npt), a(ma), covar(ma,ma), alpha(ma,ma), fit(1000), tau(1000), amin(10), amax(10), ao(10)
     &             , dfda(10), alamda, alamdao, chisq, ochisq, dchisq, xnu, freq0, wave0, kappa0, amuxmuH2, speedolight, dust2gas
     &             , da, damx, almostzero, pc, opacit0relerr, dustgasrelerr
       parameter   ( amuxmuH2 = 1.6605402D-24 * 2.8d0, speedolight = 2.99792458d10, pc = 3.085678d18, almostzero = 1.0d-30
     &             , damx = 0.2d0, itermax = 99 )
       common / copacity / freq0, wave0, kappa0, dust2gas, opacit0relerr, dustgasrelerr, iverbose
       external      lastc, modelfun, mrqmin
!__________________________________________________________________________________________________________________________________
!
       iter = 0
       chisq = 0.0d0
       dchisq = 0.0d0
       alamda = 0.0d0
       crc = ' '
       crcx = crc
       
! Compute a number of the degrees of freedom. For good results, CHI^2 should be approximately equal to XNU.

       maa = 0
       do i=1,ma
         maa = maa + ia(i)
       enddo
       xnu = npt - maa

! Set some reasonable limits on the parameter values.

       if (cmodel .eq. 'modbody') then
         chead(1) = ' TEMPDUST'
         amin(1) = 3.0d0
         amax(1) = 2000.0d0 
         if (a(5) .gt. 1.0001d0 * pc) then
           chead(2) = ' TOTLMASS'
           amin(2) = 1.0d25
           amax(2) = 1.0d40
         else
           chead(2) = '  SURFDEN'
           amin(2) = 1.0d15
           amax(2) = 1.0d30
         endif
         chead(3) = ' OPASLOPE'
         amin(3) = 0.1d0
         amax(3) = 1.0d1
         chead(4) = ' SLDANGLE'
         amin(4) = 1.0d-16
         amax(4) = 1.0d+01
         chead(5) = ' DISTANCE'
         amin(5) = 3.08572d20
         amax(5) = 6.17144d23
         chead(6) = '   CHISQ '
         chead(7) = '  DCHISQ '
         chead(8) = '  ALAMDA '
         nh = 8
         nhm3 = nh - 3
       endif
       if (cmodel .eq. 'thinbody') then
         chead(1) = ' TEMPDUST'
         amin(1) = 3.0d0
         amax(1) = 2000.0d0 
         if (a(4) .gt. 1.0001d0 * pc) then
           chead(2) = ' TOTLMASS'
           amin(2) = 1.0d25
           amax(2) = 1.0d40
         else
           chead(2) = '  SURFDEN'
           amin(2) = 1.0d15
           amax(2) = 1.0d30
         endif
         chead(3) = ' OPASLOPE'
         amin(3) = 0.1d0
         amax(3) = 1.0d1
         chead(4) = ' DISTANCE'
         amin(4) = 3.08572d20
         amax(4) = 6.17144d23
         chead(5) = '   CHISQ '
         chead(6) = '  DCHISQ '
         chead(7) = '  ALAMDA '
         nh = 7
         nhm3 = nh - 3
       endif

       if (cwhat .eq. 'finalfit' .and. iverbose .gt. 0) then
         if (iotty .gt. 0) write (iotty,'(/9a)') '   FITFUN:  IT ', (chead(i), i=1,nh)
         if (iolog .gt. 0) write (iolog,'(/9a)') '   FITFUN:  IT ', (chead(i), i=1,nh)
         if (iotty .gt. 0) write (iotty,'(/a,i4,a1,1p8e9.2)') '   FITFUN:', iter, crcx, (a(i), i=1,nhm3), chisq, dchisq, alamda
         if (iolog .gt. 0) write (iolog,'(/a,i4,a1,1p8e9.2)') '   FITFUN:', iter, crcx, (a(i), i=1,nhm3), chisq, dchisq, alamda
       endif

       do i=1,ma
         ao(i) = a(i)
       enddo
       alamda = -1.0d0
       lsingular = .false.

       call mrqmin ( x, y, sig, npt, a, ia, ma, covar, alpha, ma, chisq, modelfun, alamda, crc )

       if (crc .ne. ' ') then
         lsingular = .true.
         crcx = crc
       endif

       do i=1,ma
         da = (a(i) - ao(i)) / (ao(i) + 1.0d-30)
         if (abs ( da ) .gt. damx) a(i) = ao(i) * (1.0d0 + sign ( damx, da ))
       enddo

       do i=1,ma
         a(i) = min ( max ( a(i), amin(i) ), amax(i) )
         ao(i) = a(i)
       enddo

       iter = 1
       iconverging = 0

       if (cwhat .eq. 'finalfit' .and. iverbose .gt. 0) then
         if (iotty .gt. 0) write (iotty,'(a,i4,a1,1p8e9.2)') '   FITFUN:', iter, crcx, (a(i), i=1,nhm3), chisq, dchisq, alamda
         if (iolog .gt. 0) write (iolog,'(a,i4,a1,1p8e9.2)') '   FITFUN:', iter, crcx, (a(i), i=1,nhm3), chisq, dchisq, alamda
       endif
       
 10    continue
         iter = iter + 1
         ochisq = chisq
         alamdao = alamda

         call mrqmin ( x, y, sig, npt, a, ia, ma, covar, alpha, ma, chisq, modelfun, alamda, crc )
     
         if (crc .ne. ' ') then
           lsingular = .true.
           crcx = crc
         endif

         dchisq = chisq - ochisq
         if (ochisq .gt. almostzero) then
           dchisq = dchisq / ochisq
         endif

         if (alamda .lt. alamdao .or. (dchisq .gt. -1.0d-3 .and. dchisq .lt. -almostzero)) then
           iconverging = iconverging + 1
         else
           if (abs ( dchisq ) .gt. almostzero) iconverging = 0
         endif

         do i=1,ma
           a(i) = min ( max ( a(i), amin(i) ), amax(i) )
           ao(i) = a(i)
         enddo
         
         if (cwhat .eq. 'finalfit' .and. iverbose .eq. 2) then
           if (iotty .gt. 0) write (iotty,'(a,i4,a1,1p8e9.2)') '   FITFUN:', iter, crcx, (a(i), i=1,nhm3), chisq, dchisq, alamda
           if (iolog .gt. 0) write (iolog,'(a,i4,a1,1p8e9.2)') '   FITFUN:', iter, crcx, (a(i), i=1,nhm3), chisq, dchisq, alamda
         endif

         if (.not.lsingular .and. iter .lt. itermax .and. (alamda .lt. alamdao .or. alamda .lt. 1.0d3) .and. 
     &       (iconverging .lt. 4 .or. alamda .gt. 1.0d-5)) goto 10

         if (cwhat .eq. 'finalfit' .and. iverbose .ge. 1) then
           if (iotty .gt. 0) write (iotty,'(a,i4,a1,1p8e9.2)') '   FITFUN:', iter, crcx, (a(i), i=1,nhm3), chisq, dchisq, alamda
           if (iolog .gt. 0) write (iolog,'(a,i4,a1,1p8e9.2)') '   FITFUN:', iter, crcx, (a(i), i=1,nhm3), chisq, dchisq, alamda
         endif
       continue

       lsing = lsingular
       if (cwhat .ne. 'finalfit') then
         if (lsing .or. (alamda .ge. alamdao .and. alamda .ge. 1.0d3)) iter = itermax
       endif
       lgood = .not.lsingular .and. chisq / (xnu + 1.0d0) .le. 1.0d0
       lconv = .not.lsingular .and. iter .lt. itermax

       alamda = 0.0d0

       call mrqmin ( x, y, sig, npt, a, ia, ma, covar, alpha, ma, chisq, modelfun, alamda, crc )

       if (cwhat .eq. 'finalfit' .and. lsingular) then
         if (iotty .gt. 0) write (iotty,'(/a)') '   FITFUN: ('//crcx//') Fitting aborted in GAUSSJ with a singular matrix.'
         if (iolog .gt. 0) write (iolog,'(/a)') '   FITFUN: ('//crcx//') Fitting aborted in GAUSSJ with a singular matrix.'
       endif

       if (cwhat .eq. 'finalfit' .and. iverbose .eq. 2) then
         write (cma,'(i1)') nhm3
         if (iotty .gt. 0) write (iotty,'(/a,1p'//cma//'e9.2,a)') '   FITFUN: Abs±', (sqrt ( abs ( covar(i,i) ) ) 
     &                                                          , i=1,nhm3), ' <~ absolute uncertainties'
         if (iolog .gt. 0) write (iolog,'(/a,1p'//cma//'e9.2,a)') '   FITFUN: Abs±', (sqrt ( abs ( covar(i,i) ) ) 
     &                                                          , i=1,nhm3), ' <~ absolute uncertainties'
         if (iotty .gt. 0) write (iotty,'( a,1p'//cma//'e9.2,a)') '   FITFUN: Rel±', (sqrt ( abs ( covar(i,i) ) ) 
     &                                                          / (a(i) + 1.0d-30), i=1,nhm3), ' <~ relative uncertainties'
         if (iolog .gt. 0) write (iolog,'( a,1p'//cma//'e9.2,a)') '   FITFUN: Rel±', (sqrt ( abs ( covar(i,i) ) ) 
     &                                                          / (a(i) + 1.0d-30), i=1,nhm3), ' <~ relative uncertainties'
       endif

!  For good results, chi^2 should be approximately equal to xnu.

       if (cwhat .eq. 'finalfit' .and. iverbose .gt. 0) then
         nit = 1
         write (citer,'(i3)') iter
         if (iter .gt. 0) nit = int ( log10 ( dble ( iter ) ) ) + 1
         if (chisq / (xnu + 1.0d0) .le. 1.0d0) then
           if (iotty .gt. 0) write (iotty,'(/a,2(1pe9.2,a),a)') '   FITFUN: chi^2:', chisq, ' xnu:', xnu, ' <~ GOOD FIT: '
     &                                                        //citer(3-nit+1:3)// ' iterations'
           if (iolog .gt. 0) write (iolog,'(/a,2(1pe9.2,a),a)') '   FITFUN: chi^2:', chisq, ' xnu:', xnu, ' <~ GOOD FIT: '
     &                                                        //citer(3-nit+1:3)// ' iterations'
         else
           if (iotty .gt. 0) write (iotty,'(/a,2(1pe9.2,a),a)') '   FITFUN: chi^2:', chisq, ' xnu:', xnu, ' <~ BAD FIT: '
     &                                                        //citer(3-nit+1:3)// ' iterations'
           if (iolog .gt. 0) write (iolog,'(/a,2(1pe9.2,a),a)') '   FITFUN: chi^2:', chisq, ' xnu:', xnu, ' <~ BAD FIT: '
     &                                                        //citer(3-nit+1:3)// ' iterations'
         endif
       endif

! Output modbody values.

       do i=1,npt
         call modelfun ( x(i), a, fit(i), dfda, ma )
       enddo
       if (cmodel .eq. 'modbody') then
         do i=1,npt
           if (a(5) .gt. 1.0001d0 * pc) then 
             tau(i) = a(2) * dust2gas * kappa0 * (x(i) / freq0)**a(3) / a(4) / a(5)**2
           else
             tau(i) = a(2) * dust2gas * amuxmuH2 * kappa0 * (x(i) / freq0)**a(3)
           endif
         enddo
       endif

       write (cnpt,'(i2)') npt

       if (cwhat .eq. 'finalfit' .and. iverbose .eq. 2) then
         if (iotty .gt. 0) write (iotty,'(/a,'//cnpt//'(1pe8.1),a)') '   FITFUN:     Wavelengths:'
     &                   , (1.0d4 * speedolight / x(i), i=1,npt), ' µm'
         if (iolog .gt. 0) write (iolog,'(/a,'//cnpt//'(1pe8.1),a)') '   FITFUN:     Wavelengths:'
     &                   , (1.0d4 * speedolight / x(i), i=1,npt), ' µm'
         if (iotty .gt. 0) write (iotty,'( a,'//cnpt//'(1pe8.1),a)') '   FITFUN:     Frequencies:', (x(i), i=1,npt), ' Hz'
         if (iolog .gt. 0) write (iolog,'( a,'//cnpt//'(1pe8.1),a)') '   FITFUN:     Frequencies:', (x(i), i=1,npt), ' Hz'
         if (iotty .gt. 0) write (iotty,'( a,'//cnpt//'(1pe8.1),a)') '   FITFUN: Original fluxes:', (y(i), i=1,npt), ' Jy'
         if (iolog .gt. 0) write (iolog,'( a,'//cnpt//'(1pe8.1),a)') '   FITFUN: Original fluxes:', (y(i), i=1,npt), ' Jy'
         if (iotty .gt. 0) write (iotty,'( a,'//cnpt//'(1pe8.1),a)') '   FITFUN:  Modbody fluxes:', (fit(i), i=1,npt), ' Jy'
         if (iolog .gt. 0) write (iolog,'( a,'//cnpt//'(1pe8.1),a)') '   FITFUN:  Modbody fluxes:', (fit(i), i=1,npt), ' Jy'
         if (cmodel .eq. 'modbody') then
           if (iotty .gt. 0) write (iotty,'( a,'//cnpt//'(1pe8.1),a)') '   FITFUN:  Optical depths:', (tau(i), i=1,npt)
           if (iolog .gt. 0) write (iolog,'( a,'//cnpt//'(1pe8.1),a)') '   FITFUN:  Optical depths:', (tau(i), i=1,npt)
         endif
       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine modbody 
       
     &            ( freqx, a, y, dyda, na )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       na
       real*8        freqx, y, arg, k_B, hh, cc, bbody, taux, factor, freq0, wave0, kappa0, expargm1, exparg, dbbda1, dfada2, dfada3
     &             , dfada4, dfada5, freq, exptaux, hhokB, amuxmuH2, constant, dust2gas, a4bbody, pc, a(na), exptau, dyda(na)
     &             , iverbose
       parameter   ( k_B = 1.380658d-16, hh = 6.626d-27, cc = 2.99792458d+10, amuxmuH2 = 1.6605402D-24 * 2.8d0, hhokB = hh / k_B
     &             , constant = 1.0d23 * 2.0d0 * hh / cc**2, pc = 3.085678d18 )
       common / copacity / freq0, wave0, kappa0, dust2gas, iverbose
!__________________________________________________________________________________________________________________________________
!
! Construct a model function (modbody) used in fitting.
! a(1): temperature, a(2): mass or surface density, a(3): beta, a(4): solid angle, a(5): distance.
! Normally, beta and distance are kept fixed.

       arg = min ( hhokB * freqx / min ( max ( a(1), 3.0d0 ), 2.0d3 ), 1.0d2 )
       exparg = exp ( arg )
       expargm1 = exparg - 1.0d0
       bbody = constant * freqx**3 / expargm1
       dbbda1 = bbody / expargm1 * exparg * arg / a(1)
       factor = 1.0d0
       dfada2 = 0.0d0
       dfada3 = 0.0d0
       dfada4 = 0.0d0
       dfada5 = 0.0d0
       freq = freqx / freq0
       if (a(5) .gt. 1.0001d0 * pc) then 
         taux = a(2) * dust2gas * kappa0 * freq**a(3) / a(4) / a(5)**2
         if (taux .le. 1.0d2) then
           exptau = exp ( -taux )
           factor = 1.0d0 - exptau
           exptaux = exptau * taux
           dfada2 = exptaux / a(2)
           dfada3 = exptaux * log ( freq )
           dfada4 = (1.0d0 - exptaux) / a(4)
           dfada5 = -2.0d0 * exptaux / a(5)
         endif
       else
         taux = a(2) * dust2gas * amuxmuH2 * kappa0 * freq**a(3)
         if (taux .le. 1.0d2) then
           exptau = exp ( -taux )
           factor = 1.0d0 - exptau
           exptaux = exptau * taux
           dfada2 = exptaux / a(2)
           dfada3 = exptaux * log ( freq )
         endif
       endif

       a4bbody = a(4) * bbody
       y = a4bbody * factor

! Compute partial derivatives of the model function with respect to its parameters.

       dyda(1) = a(4) * factor * dbbda1
       dyda(2) = a4bbody * dfada2
       dyda(3) = a4bbody * dfada3
       dyda(4) = a4bbody * dfada4
       dyda(5) = a4bbody * dfada5

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine thinbody 

     &            ( freqx, a, y, dyda, na )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       na, iverbose
       real*8        freqx, y, arg, k_B, hh, cc, bbody, taux, freq0, wave0, kappa0, expargm1, exparg, dbbda1, dtada2, dtada3, dtada4
     &             , freq, hhokB, dust2gas, amuxmuH2, constant, pc, a(na), dyda(na)
       parameter   ( k_B = 1.380658d-16, hh = 6.626d-27, cc = 2.99792458d+10, amuxmuH2 = 1.6605402D-24 * 2.8d0, pc = 3.085678d18
     &             , constant = 1.0d23 * 2.0d0 * hh / cc**2, hhokB = hh / k_B )
       common / copacity / freq0, wave0, kappa0, dust2gas, iverbose
!__________________________________________________________________________________________________________________________________
!
! Construct a model function (thinbody) used in fitting.
! a(1): temperature, a(2): mass or surface density, a(3): beta, a(4): distance.
! Normally, a(3) and a(4) are kept fixed.

       arg = min ( hhokB * freqx / min ( max ( a(1), 3.0d0 ), 2.0d3 ), 1.0d2 )
       exparg = exp ( arg )
       expargm1 = exparg - 1.0d0
       bbody = constant * freqx**3 / expargm1
       dbbda1 = bbody / expargm1 * exparg * arg / a(1)
       freq = freqx / freq0
       if (a(4) .gt. 1.0001d0 * pc) then 
         taux = a(2) * dust2gas * kappa0 * freq**a(3) / a(4)**2
       else
         taux = a(2) * dust2gas * amuxmuH2 * kappa0 * freq**a(3)
       endif
       dtada2 = taux / a(2)
       dtada3 = taux * log ( freq )
       if (a(4) .gt. 1.0001d0 * pc) then 
         dtada4 = -2.0d0 * taux / a(4)
       else
         dtada4 = 0.0d0
       endif
       y = taux * bbody

! Compute partial derivatives of the model function with respect to its parameters.

       dyda(1) = taux * dbbda1
       dyda(2) = bbody * dtada2
       dyda(3) = bbody * dtada3
       dyda(4) = bbody * dtada4

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine graybody 

     &            ( freqx, a, y, dyda, na )
!__________________________________________________________________________________________________________________________________
!
! Kept here for compatibility with GETSOURCES v1.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       na, iverbose
       real*8        freqx, y, arg, k_B, hh, cc, bbody, taux, factor, freq0, wave0, kappa0, expargm1, exparg, dbbda1, dfada3, dfada4
     &             , freq, exptaux, hhokB, amuxmuH2, constant, dust2gas, a2bbody, a(na), exptau, dyda(na)
       parameter   ( k_B = 1.380658d-16, hh = 6.626d-27, cc = 2.99792458d+10, amuxmuH2 = 1.6605402D-24 * 2.8d0, hhokB = hh / k_B
     &             , constant = 1.0d23 * 2.0d0 * hh / cc**2 )
       common / copacity / freq0, wave0, kappa0, dust2gas, iverbose
!__________________________________________________________________________________________________________________________________
!
! Construct a model function (graybody) used in fitting.
! a(1): temperature, a(2): solid angle, a(3): surface density, a(4): beta.
! Normally, a(4) is kept fixed.

       arg = min ( hhokB * freqx / min ( max ( a(1), 3.0d0 ), 2.0d3 ), 1.0d2 )
       exparg = exp ( arg )
       expargm1 = exparg - 1.0d0
       bbody = constant * freqx**3 / expargm1
       dbbda1 = bbody / expargm1 * exparg * arg / a(1)
       factor = 1.0d0
       dfada3 = 0.0d0
       dfada4 = 0.0d0
       if (na .eq. 4) then 
         freq = freqx / freq0
         taux = a(3) * dust2gas * amuxmuH2 * kappa0 * freq**a(4)
         if (taux .le. 1.0d2) then
           exptau = exp ( -taux )
           factor = 1.0d0 - exptau
           exptaux = exptau * taux
           dfada3 = exptaux / a(3)
           dfada4 = exptaux * log ( freq )
         endif
       endif
       a2bbody = a(2) * bbody
       y = a2bbody * factor

! Compute partial derivatives of the model function with respect to its parameters.

       dyda(1) = a(2) * factor * dbbda1
       dyda(2) = bbody * factor
       if (na .eq. 4) then 
         dyda(3) = a2bbody * dfada3
         dyda(4) = a2bbody * dfada4
       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine covsrt 
       
     &            ( covar, npc, ma, ia, mfit )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       ma, mfit, npc, ia(ma), i, j, k
       real*8        covar(npc,npc), swap
!__________________________________________________________________________________________________________________________________
!
       do i=mfit+1,ma
         do j=1,i
           covar(i,j) = 0.0d0
           covar(j,i) = 0.0d0
         enddo
       enddo
       k = mfit
       do j=ma,1,-1
         if (ia(j) .ne. 0) then
           do i=1,ma
             swap = covar(i,k)
             covar(i,k) = covar(i,j)
             covar(i,j) = swap
           enddo
           do i=1,ma
             swap = covar(k,i)
             covar(k,i) = covar(j,i)
             covar(j,i) = swap
           enddo
           k = k - 1
         endif
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine gaussj 
       
     &            ( a, n, np, b, m, mp, crc )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*1   crc
       integer       nmax
       parameter   ( nmax = 50 )
       integer       m, mp, n, np, i, icol, irow, j, k, l, ll, indxc(nmax), indxr(nmax), ipiv(nmax)
       real*8        big, dum, pivinv, a(np,np), b(np,mp)
!__________________________________________________________________________________________________________________________________
!
       crc = ' '
       do j=1,n
         ipiv(j) = 0
       enddo
       irow = 1
       icol = 1
       do i=1,n
         big = 0.0d0
         do j=1,n
           if (ipiv(j) .ne. 1) then
             do k=1,n
               if (ipiv(k) .eq. 0) then
                 if (abs ( a(j,k) ) .ge. big) then
                   big = abs ( a(j,k) )
                   irow = j
                   icol = k
                 endif
               elseif (ipiv(k) .gt. 1) then
                 crc = '$' !<-- Singular matrix: ipiv(k) > 1'
                 return
               endif
             enddo
           endif
         enddo
         ipiv(icol) = ipiv(icol) + 1
         if (irow .ne. icol) then
           do l=1,n
             dum = a(irow,l)
             a(irow,l) = a(icol,l)
             a(icol,l) = dum
           enddo
           do l=1,m
             dum = b(irow,l)
             b(irow,l) = b(icol,l)
             b(icol,l) = dum
           enddo
         endif
         indxr(i) = irow
         indxc(i) = icol
!!         if (a(icol,icol) .eq. 0.0d0) then
         if (abs ( a(icol,icol) ) .lt. 1.0d-99) then
           crc = '*' !<-- Singular matrix: a(icol,icol) = 0'
           return
         endif
         pivinv = 1.0d0 / a(icol,icol)
         a(icol,icol) = 1.0d0
         do l=1,n
           a(icol,l) = a(icol,l) * pivinv
         enddo
         do l=1,m
           b(icol,l) = b(icol,l) * pivinv
         enddo
         do ll=1,n
           if (ll .ne. icol) then
             dum = a(ll,icol)
             a(ll,icol) = 0.0d0
             do l=1,n
               a(ll,l) = a(ll,l) - a(icol,l) * dum
             enddo
             do l=1,m
               b(ll,l) = b(ll,l) - b(icol,l) * dum
             enddo
           endif
         enddo
       enddo
       do l=n,1,-1
         if (indxr(l) .ne. indxc(l)) then
           do k=1,n
             dum = a(k,indxr(l))
             a(k,indxr(l)) = a(k,indxc(l))
             a(k,indxc(l)) = dum
           enddo
         endif
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine mrqcof 
       
     &            ( x, y, sig, ndata, a, ia, ma, alpha, beta, nalp, chisq, funcs )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN.
!
! Used by MRQMIN to evaluate the linearized fitting matrix ALPHA and vector BETA, and calculate CHI^2.
!__________________________________________________________________________________________________________________________________
!
       implicit       none
       integer        mmax
       parameter    ( mmax = 5 )
       integer        ma, nalp, ndata, ia(ma), mfit, i, j, k, l, m
       real*8         chisq, a(ma), alpha(nalp,nalp), beta(ma), sig(ndata), x(ndata), y(ndata), dy, sig2i, wt, ymod, dyda(mmax)
       external       funcs
!__________________________________________________________________________________________________________________________________
!
       mfit = 0

       do j=1,ma
         if (ia(j) .ne. 0) mfit = mfit + 1
       enddo

! Initialize (symmetric) alpha, beta.

       do j=1,mfit
         do k=1,j
           alpha(j,k) = 0.0d0
         enddo
         beta(j) = 0.0d0
       enddo

! Summation loop over all data.

       chisq = 0.0d0

       do i=1,ndata

         call funcs ( x(i), a, ymod, dyda, ma )

         sig2i = 1.0d0 / (sig(i) * sig(i))

         dy = y(i) - ymod
         j = 0
         do l=1,ma
           if (ia(l) .ne. 0) then
             j = j + 1
             wt = dyda(l) * sig2i
             k = 0

!!         write (*,*) 'aaa', wt, dyda, sig2i

             do m=1,l
               if (ia(m) .ne. 0) then
                 k = k + 1
                 alpha(j,k) = alpha(j,k) + wt * dyda(m)
               endif
             enddo

!!         write (*,*) 'bbb'

             beta(j) = beta(j) + dy * wt
           endif
         enddo

! And find chi^2.

         chisq = chisq + dy * dy * sig2i
       enddo

! Fill in the symmetric side.

       do j=2,mfit
         do k=1,j-1
           alpha(k,j) = alpha(j,k)
         enddo
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine mrqmin 
       
     &            ( x, y, sig, ndata, a, ia, ma, covar, alpha, nca, chisq, funcs, alamda, crc )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN.
!
! Levenberg-Marquardt method, attempting to reduce the value chi^2 of a fit between a set of data points x(1:ndata), y(1:ndata)
! with individual standard deviations sig(1:ndata), and a non-linear function dependent on ma coefficients a(1:ma). The input
! array ia(1:ma) indicates by nonzero entries those components of 'a' that should be fitted for and by zero entries those
! components that should be held fixed at their input values. The program returns current best fit values for the parameters
! a(1:ma) and chi^2 = 'chisq'. The arrays covar(1:nca,1:nca), alpha(1:nca,1:nca) with physical dimension 'nca' (>= the number of
! fitted parameters) are used as working space during most iterations. Supply a subroutine funcs(x,a,yfit,dyda,ma) that evaluates
! the fitting function 'yfit' and its derivatives 'dyda' with respect to the fitting parameters 'a' at 'x'. On the first call
! provide an initial guess for the parameters 'a', and set 'alamda' < 0 for initialization (which then sets alamda=0.001). If a
! step succeeds 'chisq' becomes smaller and 'alamda' decreases by a factor of 10. You must call this routine repeatedly until
! convergence is achieved. Then, make one final call with alamda=0, so that covar(1:ma,1:ma) returns the covariance matrix, and
! 'alpha' the curvature matrix. (Parameters held fixed will return zero covariances).
!
! USES covsrt, gaussj, mrqcof
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*1   crc
       integer       mmax
       parameter   ( mmax = 20 )
       integer       ma, nca, ndata, ia(ma), j, k, l, mfit
       real*8        alamda, chisq, a(ma), alpha(nca,nca), covar(nca,nca), sig(ndata), x(ndata), y(ndata), ochisq
     &             , atry(mmax), beta(mmax), da(mmax)
       save          ochisq, atry, beta, da, mfit
       external      funcs, mrqcof, gaussj, covsrt
!__________________________________________________________________________________________________________________________________
!
! Initiaization.

       if (alamda .lt. 0.0d0) then
         mfit = 0
         do j=1,ma
           if (ia(j) .ne. 0) mfit = mfit + 1
         enddo
         alamda = 0.001d0

         call mrqcof ( x, y, sig, ndata, a, ia, ma, alpha, beta, nca, chisq, funcs )

         ochisq = chisq
         do j=1,ma
           atry(j) = a(j)
         enddo
       endif

! Alter linearized fitting matrix by augmenting diagonal elements.

       do j=1,mfit
         do k=1,mfit
           covar(j,k) = alpha(j,k)
         enddo
         covar(j,j) = alpha(j,j) * (1.0d0 + alamda)
         da(j) = beta(j)
       enddo

! Matrix solution.

       call gaussj ( covar, mfit, nca, da, 1, 1, crc )

! Once converged, evaluate covariance matrix.

       if (alamda .eq. 0.0d0) then

         call covsrt ( covar, nca, ma, ia, mfit )

! Spread out 'alpha' to its full size, too.

         call covsrt ( alpha, nca, ma, ia, mfit )

         return
       endif

! Did the trial succeed?

       j=0
       do l=1,ma
         if (ia(l) .ne. 0) then
           j = j + 1
           atry(l) = a(l) + da(j)
!!           atry(l) = max ( a(l) + da(j), 1.0d-90 )
         endif
       enddo

       call mrqcof ( x, y, sig, ndata, atry, ia, ma, covar, da, nca, chisq, funcs )

       if (chisq .lt. ochisq) then

! Success, accept the new solution.

         alamda = 0.1d0 * alamda
         ochisq = chisq
         do j=1,mfit
           do k=1,mfit
             alpha(j,k) = covar(j,k)
           enddo
           beta(j) = da(j)
         enddo
         do l=1,ma
           a(l) = atry(l)
         enddo
       else

! Failure, increase 'alamda' and return.

         alamda = 10.0d0 * alamda
         chisq = ochisq
       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine spline
       
     &            ( x, y, n, yp1, ypn, y2 )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN.
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., y_i = f(x_i), with x_1 < x_2 < ... < x_n, and given values
! yp1 and ypn for the first derivative of the interpolating function at points 1 and n, respectively, this routine returns an array
! y2(1:n) of length n which contains the second derivatives of the interpolating function at the tabulated points x_i. If yp1 and/or
! ypn are equal to 10^30 or larger, the routine is signaled to set the corresponding boundary condition for a natural spline, with
! with zero second derivative at that boundary. Parameter: nmax is the largest anticipated value of n.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       n, nmax, i, k
       parameter     ( nmax = 500 )
       real*8        yp1, ypn, x(n), y(n), y2(n), p, qn, sig, un, u(nmax)
!__________________________________________________________________________________________________________________________________
!
       if (yp1 .gt. 0.99d30) then
         y2(1) = 0.0d0
         u (1) = 0.0d0
       else
         y2(1) = -0.5d0
         u (1) = (3.0d0 / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
       endif
       do i=2,n-1
         sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
         p = sig * y2(i-1) + 2.0d0
         y2(i) = (sig - 1.0d0) / p
         u (i) = (6.0d0 * ((y(i+1) - y(i)) / (x(i+1) - x(i)) - (y(i) - y(i-1)) / (x(i) - x(i-1))) 
     &         / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
       enddo
       if (ypn .gt. 0.99d30) then
         qn = 0.0d0
         un = 0.0d0
       else
         qn = 0.5d0
         un = (3.0d0 / (x(n) - x(n-1))) * (ypn - (y(n) - y(n-1)) / (x(n) - x(n-1)))
       endif
       y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + 1.0d0)
       do k=n-1,1,-1
         y2(k) = y2(k) * y2(k+1) + u(k)
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine splint
       
     &            ( xa, ya, y2a, n, x, y )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       n, k, khi, klo
       real*8        x, y, xa(n), y2a(n), ya(n), a, b, h
!__________________________________________________________________________________________________________________________________
!
       klo = 1
       khi = n

 1     if (khi - klo .gt. 1) then
         k = (khi + klo) / 2
         if (xa(k) .gt. x) then
           khi = k
         else
           klo = k
         endif
       goto 1
       endif
       h = xa(khi) - xa(klo)
       if (h .eq. 0.0d0) then
         write (*,*) 'bad xa input in splint'
         read (*,*)
       endif

       a = (xa(khi) - x) / h
       b = (x - xa(klo)) / h
       y = a * ya(klo) + b * ya(khi) + ((a**3 - a) * y2a(klo) + (b**3 - b) * y2a(khi)) * (h**2) / 6.0d0
      
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function ran1 
       
     &          ( idum )
!__________________________________________________________________________________________________________________________________
!
!    whp:zeus3d.ran1 <------------------- returns uniform random deviate
!    from numerical recipes                               december, 1992
!
!    written by: Press etal.
!    modified 1:
!
!  PURPOSE:  Returns a uniformly distrbuted random deviate between 0.0 and 1.0 (see Numerical Recipes, 1st edition for FORTRAN,
!            page 196).
!
!  INPUT VARIABLES:
!    idum        .ge. 0 => get next random number from sequence
!                .lt. 0 => initialise sequence depending on absolute value of "idum", and gets next random number
!
!  OUTPUT VARIABLES:
!    ran1        deviate with uniform distribution between 0.0 and 1.0
!
!  LOCAL VARIABLES:
!
!  EXTERNALS: [NONE]
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       imax, m1, ia1, ic1, m2, ia2, ic2, m3, ia3, ic3
       real*8        rm1, rm2
       parameter     ( imax = 97 )
       parameter     ( m1 = 259200, ia1 = 7141, ic1 = 54733, rm1 = 1.0d0 / m1 )
       parameter     ( m2 = 134456, ia2 = 8121, ic2 = 28411, rm2 = 1.0d0 / m2 )
       parameter     ( m3 = 243000, ia3 = 4561, ic3 = 51349 )
       integer       idum, iff, ix1, ix2, ix3, i
       real*8        ran1
       real*8        r(imax)
       data          iff / 0 /
       save          iff, ix1, ix2, ix3, r
!__________________________________________________________________________________________________________________________________
!
!  If sequence is to be initialised, start here ...

       if (idum .lt. 0 .or. iff .eq. 0) then
         iff = 1
         ix1 = mod ( ic1 - idum     , m1 )
         ix1 = mod ( ia1 * ix1 + ic1, m1 )
         ix2 = mod ( ix1            , m2 )
         ix1 = mod ( ia1 * ix1 + ic1, m1 )
         ix3 = mod ( ix1            , m3 )
         do i=1,imax
           ix1  = mod ( ia1 * ix1 + ic1, m1 )
           ix2  = mod ( ia2 * ix2 + ic2, m2 )
           r(i) = ( dfloat ( ix1 ) + dfloat ( ix2 ) * rm2 ) * rm1
         enddo
         idum = 1
       endif

!  Otherwise, start here ....

       ix1  = mod ( ia1 * ix1 + ic1, m1 )
       ix2  = mod ( ia2 * ix2 + ic2, m2 )
       ix3  = mod ( ia3 * ix3 + ic3, m3 )
       i    = 1 + ( imax * ix3 ) / m3
       ran1 = r(i)
       r(i) = ( dfloat ( ix1 ) + dfloat ( ix2 ) * rm2 ) * rm1

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      function gasdev
               
     &         ( idum )
!__________________________________________________________________________________________________________________________________
!
!  PURPOSE:  Returns a normally distrbuted deviate with zero mean and unit variance using RAN1 as the source of uniform deviates
!            (see Numerical Recipes, 1st edition for FORTRAN).
!__________________________________________________________________________________________________________________________________
!
      implicit none
      integer idum
      real*8 gasdev
cu    uses ran1
      integer iset
      real*8 fac,gset,rsq,v1,v2,ran1
      save iset,gset
      data iset /0/
      external ran1
!__________________________________________________________________________________________________________________________________
!
      if (iset .eq. 0) then
 1      continue
        v1 = 2.0d0 * ran1 ( idum ) - 1.0d0
        v2 = 2.0d0 * ran1 ( idum ) - 1.0d0
        rsq = v1**2 + v2**2
        if(rsq .ge. 1.0d0 .or. rsq .eq. 0.0d0) goto 1
        fac = sqrt ( -2.0d0 * log ( rsq ) / rsq )
        gset = v1 * fac
        gasdev = v2 * fac
        iset = 1
      else
        gasdev = gset
        iset = 0
      endif
      
      return
      end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine planck 
       
     &            ( t, fq, bf )
!__________________________________________________________________________________________________________________________________
!
! Planck function bf at temperature t and wave number fq = 1 / lambda.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       real*8        t, fq, bf, x1, x2, expmax, xx, ex, c2, exmx
       data          x1 / 3.972968d-16 /, x2 / 1.438833d+00 /, expmax / 5.0d+01 /
!__________________________________________________________________________________________________________________________________
!
! here: x1=2*h*c, x2=h*c/k; 1/(exp{x}-1) => exp{-x}/(1-exp{-x})

       xx = x2 * fq / t

       if (xx .le. expmax) then
         ex = exp ( -xx )
         bf = x1 * fq**3 * ex / ( 1.0d0 - ex )
       else
         c2 = expmax / x2
         exmx = exp ( -expmax )
         bf = x1 * (c2 * t)**3 * exmx / ( 1.0d0 - exmx )
       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine mess 
       
     &            ( message, crlf )
!__________________________________________________________________________________________________________________________________
!
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       iotty, iolog
       character*(*) message, crlf
       common / lmess / iotty, iolog
!__________________________________________________________________________________________________________________________________
!
       if (iotty .gt. 0) then
         if (crlf .eq. 'a' .or. crlf .eq. 'ab') write (iotty,'(a)') '         :'
         write (iotty, '(a)') message
         if (crlf .eq. 'b' .or. crlf .eq. 'ab') write (iotty,'(a)') '         :'
       end if
       if (iolog .gt. 0) then
         if (crlf .eq. 'a' .or. crlf .eq. 'ab') write (iolog,'(a)') '         :'
         write (iolog, '(a)') message
         if (crlf .eq. 'b' .or. crlf .eq. 'ab') write (iolog,'(a)') '         :'
       end if

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine fillarea 
     
     &            ( lpartition, iseedx, iseedy, nextrmax, nextr, newpixval, method, zero, boundpixval, nx, ny, n1x, n2x, n1y, n2y
     &            , nxmin, nxmax, nymin, nymax, iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )
!__________________________________________________________________________________________________________________________________
!
! Alexander Men'shchikov, DAp IRFU CEA Saclay, France, on March 12, 2009
! Based on the Tint Fill Algorithm by Alvy Ray Smith, Computer Graphics Lab, New York Institute of Technology, Old Westbury,
! NY 11568. Published in: SIGGRAPH 79 Conference Proceedings, Aug 1979, 276-282. Also Technical Memo No 6, Computer Graphics Lab,
! New York Institute of Technology, Jul 1978, and issued as tutorial notes at SIGGRAPHs 78, 80-82.
!__________________________________________________________________________________________________________________________________
!                     
       implicit      none
       logical       logic, lpartition, lnxmin, lnxmax, lnymin, lnymax, lnochange
       character*(*) method
       integer       iseedx, iseedy, ix, iy, ilx, irx, maxstack, nx, ny, nstack, iyref, ilxref, irxref, nextrmax, nextr, n1y, n2y, n
     &             , nxmin, nxmax, nymin, nymax, iter, i, j, nxmino, nxmaxo, nymino, nymaxo, irc
       parameter   ( maxstack = 10000000 )
       integer       n1x(*), n2x(*), ntouching(nextrmax,*), nsx(*)
       integer, allocatable :: ixstack(:), iystack(:)
       real*8        newpixval, oldpixval, boundpixval, zero
       real*8        iworker(nx,ny), segmn(nx,ny), segmo(nx,ny), mxco(*), myco(*), deltx(*), delty(*)
       external      fillleft, fillright, scanhi, scanlo, pop, push
!__________________________________________________________________________________________________________________________________
!
       nstack = 0
       ix = iseedx
       iy = iseedy
       nxmin = ix
       nxmax = ix
       nymin = iy
       nymax = iy
       ilx = ix
       irx = ix
       
       allocate ( ixstack(maxstack), iystack(maxstack), stat=irc )

       if (irc .ne. 0) then
         write (*,'(/a)') ' FILLAREA: ERROR: Trouble allocating memory (5).'
         stop 5
       endif

       oldpixval = iworker(ix,iy)

       n = nint ( newpixval )
       
       if (method .eq. 'positive' .or. method .eq. 'boundary' .or. 
     &     abs ( newpixval - oldpixval ) .gt. zero * abs ( newpixval + oldpixval )) then
         iyref = iy
         ilxref = ix
         irxref = ix

         call push ( ix, iy, ixstack, iystack, nstack, maxstack )
         
         do while (nstack .gt. 0)
       
           call pop ( ixstack, iystack, nstack, maxstack, ix, iy )
                                               
           if (abs ( newpixval - iworker(ix,iy) ) .gt. zero * abs ( newpixval + iworker(ix,iy) )) then
           
             call fillright ( lpartition, ix, iy, irx, nextrmax, nextr, oldpixval, newpixval, method, zero, boundpixval, nx, ny, n2x
     &                      , iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )
             
             call fillleft  ( lpartition, ix, iy, ilx, nextrmax, nextr, oldpixval, newpixval, method, zero, boundpixval, nx, ny, n1x
     &                      , iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )
             
             logic = ilx .ge. ilxref - 1 .and. irx .le. irxref + 1

             if (iy .eq. iyref + 1 .and. logic) then
                    
               call scanhi ( lpartition, iy, ilx, irx, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, ixstack
     &                     , iystack, nstack, maxstack, nx, ny, n2y, iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching
     &                     , nsx )
             
             else if (iy .eq. iyref - 1 .and. logic) then
               
               call scanlo ( lpartition, iy, ilx, irx, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, ixstack
     &                     , iystack, nstack, maxstack, nx, ny, n1y, iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching
     &                     , nsx )
             else
             
               call scanhi ( lpartition, iy, ilx, irx, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, ixstack
     &                     , iystack, nstack, maxstack, nx, ny, n2y, iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching
     &                     , nsx )

               call scanlo ( lpartition, iy, ilx, irx, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, ixstack
     &                     , iystack, nstack, maxstack, nx, ny, n1y, iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching
     &                     , nsx )
             endif
             
             iyref = iy
             ilxref = ilx
             irxref = irx
           endif
         enddo
       endif

! Determine the bounding box for the segmentation mask of a source.

       nxmin = min ( nxmin, ilx )
       nxmax = max ( nxmax, irx )
       nymin = min ( nymin, iy )
       nymax = max ( nymax, iy )

       lnxmin = .true.
       lnxmax = .true.
       lnymin = .true.
       lnymax = .true.

       do iter=1,nx*ny
       
! Expand vertical edges to find the box limits.

         if (lnxmin .or. lnymin .or. lnymax) then
           lnxmin = .false.
           nxmino = nxmin
           do j=nymin,nymax
             if (nint ( segmn(nxmin,j) ) .eq. n) then
               nxmin = max ( nxmin - 1, n1x(j) )
               if (nxmin .ne. nxmino) lnxmin = .true.
               exit
             endif
           enddo
         endif
         
         if (lnxmax .or. lnymin .or. lnymax) then
           lnxmax = .false.
           nxmaxo = nxmax
           do j=nymin,nymax
             if (nint ( segmn(nxmax,j) ) .eq. n) then
               nxmax = min ( nxmax + 1, n2x(j) )
               if (nxmax .ne. nxmaxo) lnxmax = .true.
               exit
             endif
           enddo
         endif
              
! Expand horizontal edges to find the box limits.

         if (lnymin .or. lnxmin .or. lnxmax) then
           lnymin = .false.
           nymino = nymin
           do i=nxmin,nxmax
             if (nint ( segmn(i,nymin) ) .eq. n) then
               nymin = max ( nymin - 1, n1y )
               if (nymin .ne. nymino) lnymin = .true.
               exit
             endif
           enddo
         endif
         
         if (lnymax .or. lnxmin .or. lnxmax) then
           lnymax = .false.
           nymaxo = nymax
           do i=nxmin,nxmax
             if (nint ( segmn(i,nymax) ) .eq. n) then
               nymax = min ( nymax + 1, n2y )
               if (nymax .ne. nymaxo) lnymax = .true.
               exit
             endif
           enddo
         endif

         lnochange = .not.lnxmin .and. .not.lnxmax .and. .not.lnymin .and. .not.lnymax
         if (lnochange) exit
       enddo

       deallocate ( ixstack, iystack )
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine fillright 
       
     &            ( lpartition, ix, iy, irx, nextrmax, nextr, oldpixval, newpixval, method, zero, boundpixval, nx, ny, n2x, iworker
     &            , segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )
!__________________________________________________________________________________________________________________________________
!
! Alexander Men'shchikov, DAp IRFU CEA Saclay, France, on March 12, 2009
! Based on the Tint Fill Algorithm by Alvy Ray Smith, Computer Graphics Lab, New York Institute of Technology, Old Westbury,
! NY 11568. Published in: SIGGRAPH 79 Conference Proceedings, Aug 1979, 276-282. Also Technical Memo No 6, Computer Graphics Lab,
! New York Institute of Technology, Jul 1978, and issued as tutorial notes at SIGGRAPHs 78, 80-82.
!__________________________________________________________________________________________________________________________________
!                     
       implicit      none
       logical       lpartition, lpixval
       character*(*) method
       integer       ix, iy, irx, nx, ny, mx, nextrmax, nextr, n2x(*), ntouching(nextrmax,*), nsx(*)
       real*8        oldpixval, newpixval, boundpixval, zero
       real*8        iworker(nx,ny), segmn(nx,ny), segmo(nx,ny), mxco(*), myco(*), deltx(*), delty(*)
       external      lpixval
!__________________________________________________________________________________________________________________________________
!
       mx = ix
       do while (mx .le. n2x(iy))
         if (lpixval ( lpartition, '.eq.', mx, iy, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, nx, ny, iworker
     &               , segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )) then 
           iworker(mx,iy) = newpixval
           segmn  (mx,iy) = newpixval
           mx = mx + 1
         else
           exit
         endif
       enddo
       irx = mx - 1
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine fillleft 
       
     &            ( lpartition, ix, iy, ilx, nextrmax, nextr, oldpixval, newpixval, method, zero, boundpixval, nx, ny, n1x, iworker
     &            , segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )
!__________________________________________________________________________________________________________________________________
!
! Alexander Men'shchikov, DAp IRFU CEA Saclay, France, on March 12, 2009
! Based on the Tint Fill Algorithm by Alvy Ray Smith, Computer Graphics Lab, New York Institute of Technology, Old Westbury,
! NY 11568. Published in: SIGGRAPH 79 Conference Proceedings, Aug 1979, 276-282. Also Technical Memo No 6, Computer Graphics Lab,
! New York Institute of Technology, Jul 1978, and issued as tutorial notes at SIGGRAPHs 78, 80-82.
!__________________________________________________________________________________________________________________________________
!                     
       implicit      none
       logical       lpartition, lpixval
       character*(*) method
       integer       ix, iy, ilx, nx, ny, mx, nextrmax, nextr, n1x(*), ntouching(nextrmax,*), nsx(*)
       real*8        oldpixval, newpixval, boundpixval, zero
       real*8        iworker(nx,ny), segmn(nx,ny), segmo(nx,ny), mxco(*), myco(*), deltx(*), delty(*)
       external      lpixval
!__________________________________________________________________________________________________________________________________
!
       mx = ix - 1
       if (method .eq. 'boundary') mx = ix
       do while (mx .ge. n1x(iy))
         if (lpixval ( lpartition, '.eq.', mx, iy, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, nx, ny, iworker
     &               , segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )) then
           if (method .eq. 'boundary' .and. mx .eq. ix) mx = ix - 1
           iworker(mx,iy) = newpixval
           segmn  (mx,iy) = newpixval
           mx = mx - 1
         else
           exit
         endif
       enddo
       ilx = mx + 1
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine scanhi 
       
     &            ( lpartition, iy, ilx, irx, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, ixstack, iystack
     &            , nstack, maxstack, nx, ny, n2y, iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )
!__________________________________________________________________________________________________________________________________
!
! Alexander Men'shchikov, DAp IRFU CEA Saclay, France, on March 12, 2009
! Based on the Tint Fill Algorithm by Alvy Ray Smith, Computer Graphics Lab, New York Institute of Technology, Old Westbury,
! NY 11568. Published in: SIGGRAPH 79 Conference Proceedings, Aug 1979, 276-282. Also Technical Memo No 6, Computer Graphics Lab,
! New York Institute of Technology, Jul 1978, and issued as tutorial notes at SIGGRAPHs 78, 80-82.
!__________________________________________________________________________________________________________________________________
!                     
       implicit      none
       logical       lpartition, lpixval
       character*(*) method
       integer       iy, nstack, nx, ny, ilx, irx, mx, my, maxstack, nextrmax, nextr, n2y
       integer       ixstack(maxstack), iystack(maxstack), ntouching(nextrmax,*), nsx(*)
       real*8        oldpixval, newpixval, boundpixval, zero
       real*8        iworker(nx,ny), segmn(nx,ny), segmo(nx,ny), mxco(*), myco(*), deltx(*), delty(*)
       external      lpixval, push
!__________________________________________________________________________________________________________________________________
!
       mx = ilx
       my = iy + 1
       if (my .le. n2y) then
         do while (mx .le. irx)
           do while (mx .le. irx)
             if (lpixval ( lpartition, '.ne.', mx, my, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, nx, ny
     &                   , iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )) then
               mx = mx + 1
             else
               exit
             endif
           enddo
           if (mx .gt. irx) exit

           call push ( mx, my, ixstack, iystack, nstack, maxstack )
           
           do while (mx .le. irx)
             if (lpixval ( lpartition, '.eq.', mx, my, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, nx, ny
     &                   , iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )) then
               mx = mx + 1
             else
               exit
             endif
           enddo
         enddo
       endif
       
       return
       end 
 
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine scanlo 
       
     &            ( lpartition, iy, ilx, irx, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, ixstack, iystack
     &            , nstack, maxstack, nx, ny, n1y, iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )
!__________________________________________________________________________________________________________________________________
!
! Alexander Men'shchikov, DAp IRFU CEA Saclay, France, on March 12, 2009
! Based on the Tint Fill Algorithm by Alvy Ray Smith, Computer Graphics Lab, New York Institute of Technology, Old Westbury,
! NY 11568. Published in: SIGGRAPH 79 Conference Proceedings, Aug 1979, 276-282. Also Technical Memo No 6, Computer Graphics Lab,
! New York Institute of Technology, Jul 1978, and issued as tutorial notes at SIGGRAPHs 78, 80-82.
!__________________________________________________________________________________________________________________________________
!                     
       implicit      none
       logical       lpartition, lpixval
       character*(*) method
       integer       iy, nstack, nx, ny, ilx, irx, mx, my, maxstack, nextrmax, nextr, n1y
       integer       ixstack(maxstack), iystack(maxstack), ntouching(nextrmax,*), nsx(*)
       real*8        oldpixval, newpixval, boundpixval, zero
       real*8        iworker(nx,ny), segmn(nx,ny), segmo(nx,ny), mxco(*), myco(*), deltx(*), delty(*)
       external      lpixval, push
!__________________________________________________________________________________________________________________________________
!
       mx = ilx
       my = iy - 1
       if (my .ge. n1y) then
         do while (mx .le. irx)
           do while (mx .le. irx)
             if (lpixval ( lpartition, '.ne.', mx, my, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, nx, ny
     &                   , iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )) then
               mx = mx + 1
             else
               exit
             endif
           enddo
           if (mx .gt. irx) exit
         
           call push ( mx, my, ixstack, iystack, nstack, maxstack )
         
           do while (mx .le. irx)
             if (lpixval ( lpartition, '.eq.', mx, my, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, nx, ny
     &                   , iworker, segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )) then
               mx = mx + 1
             else
               exit
             endif
           enddo
         enddo
       endif
       
       return
       end 

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function lpixval 
       
     &          ( lpartition, clogic, ix, iy, nextrmax, nextr, oldpixval, newpixval, boundpixval, method, zero, nx, ny, iworker
     &          , segmn, segmo, mxco, myco, deltx, delty, ntouching, nsx )
!__________________________________________________________________________________________________________________________________
!
! Alexander Men'shchikov, DAp IRFU CEA Saclay, France, on March 12, 2009
!__________________________________________________________________________________________________________________________________
!                     
       implicit      none
       logical       lpartition, lpixval, lntakenbyanother
       character*(*) method, clogic
       integer       ix, iy, nx, ny, nextrmax, newpixvali, n, k, l, mindist, nextr
       integer       ntouching(nextrmax,*), nsx(*)
       real*8        oldpixval, newpixval, boundpixval, zero, rxykn2min, rxyk2, dblix, dbliy
       real*8        iworker(nx,ny), segmn(nx,ny), segmo(nx,ny), mxco(*), myco(*), deltx(*), delty(*)
!__________________________________________________________________________________________________________________________________
!
       newpixvali = nint ( newpixval )                         

       k = newpixvali

       n = nint ( segmn(ix,iy) )

! If the segmentation pixel >= 1, this means that it has been assigned to a source
! and if the source number differs from the current filling number, then it's a different source: it's taken.

       if (segmn(ix,iy) .ge. 1.0d0 .and. k .ne. n) then
         lntakenbyanother = .true.
       else
         lntakenbyanother = .false.    !<-- n = 0 or newpixvali = n
       endif
              
! The next block works for existing touching sources whose pixels need to be partitioned.

       if (lpartition) then
         
         if (segmn(ix,iy) .ge. 1.0d0 .and. lntakenbyanother) then
                                                           
! If a source n is touching the source k, then consider the distances of the current pixel (ix,iy) from both sources
! to decide, whether the pixel should belong to k or to n.
                                                           
           if (nint ( iworker(ix,iy) ) .eq. n) then
             if (segmo(ix,iy) .lt. 1.0d0) then
               mindist = n
               do l=1,nextr
                 if (l .ne. n .and. ntouching(l,n) .eq. 1) then
                   dblix = dble ( ix )                           
                   dbliy = dble ( iy )                           
                   rxykn2min = (dblix - (mxco(n) + deltx(l)))**2 + (dbliy - (myco(n) + delty(l)))**2
                   rxyk2 = (dblix - mxco(l))**2 + (dbliy - myco(l))**2
                   if (rxyk2 .lt. rxykn2min) then                
                     mindist = l
                   endif
                 endif
               enddo
               if (mindist .eq. k) lntakenbyanother = .false.
             endif
           endif                                           
                                                           
! If a source n is merged with a more prominent source k, then its segmentation pixels are free to be used by the other.
                                                           
           if (nsx(n) .eq. 0) then
             if (ntouching(k,n) .eq. 2) then
               lntakenbyanother = .false.
             endif
           endif
         endif
       endif
       
       if (clogic .eq. '.ne.') then
              
! Pixel value is not good for taking the pixel in segmentation mask              
         
         if (method .eq. 'interior') then
           lpixval = abs ( iworker(ix,iy) - oldpixval   ) .gt. zero * abs ( iworker(ix,iy) + oldpixval )
         else if (method .eq. 'boundary') then
           lpixval = abs ( iworker(ix,iy) - boundpixval ) .le. zero * abs ( iworker(ix,iy) + oldpixval )
         else if (method .eq. 'positive') then
           lpixval = iworker(ix,iy) .le. zero .or. lntakenbyanother
         else if (method .eq. 'negative') then
           lpixval = iworker(ix,iy) .gt. zero
         else
           write (*,'(/a)') '   LPIXVAL: Invalid method: '//method
           stop 99
         endif

       else if (clogic .eq. '.eq.') then

! Pixel value is good for taking the pixel into segmentation mask              
         
         if (method .eq. 'interior') then
           lpixval = abs ( iworker(ix,iy) - oldpixval   ) .le. zero * abs ( iworker(ix,iy) + oldpixval )
         else if (method .eq. 'boundary') then
           lpixval = abs ( iworker(ix,iy) - boundpixval ) .gt. zero * abs ( iworker(ix,iy) + oldpixval )
         else if (method .eq. 'positive') then 
           lpixval = iworker(ix,iy) .gt. zero .and. .not.lntakenbyanother
         else if (method .eq. 'negative') then
           lpixval = iworker(ix,iy) .le. zero
         else
           write (*,'(/a)') '   LPIXVAL: Invalid method: '//method
           stop 99
         endif
       else
         write (*,'(/a)') '   LPIXVAL: Invalid clogic: '//clogic
         stop 99
       endif
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine push 
       
     &            ( ix, iy, ixstack, iystack, nstack, maxstack )
!__________________________________________________________________________________________________________________________________
!
! Alexander Men'shchikov, DAp IRFU CEA Saclay, France, on March 12, 2009
!__________________________________________________________________________________________________________________________________
!                     
       implicit      none
       integer       ix, iy, nstack, maxstack
       integer       ixstack(maxstack), iystack(maxstack)
!__________________________________________________________________________________________________________________________________
!
       nstack = nstack + 1
       if (nstack .gt. maxstack) then
         write (*,'(/a,i10)') '   PUSH: Stack overflow: maxstack:', maxstack, ix, iy
         stop 99
       endif       
       ixstack(nstack) = ix
       iystack(nstack) = iy
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine pop 

     &            ( ixstack, iystack, nstack, maxstack, ix, iy )
!__________________________________________________________________________________________________________________________________
!
! Alexander Men'shchikov, DAp IRFU CEA Saclay, France, on March 12, 2009
!__________________________________________________________________________________________________________________________________
!                     
       implicit      none
       integer       ix, iy, nstack, maxstack
       integer       ixstack(maxstack), iystack(maxstack)
!__________________________________________________________________________________________________________________________________
!
       ix = ixstack(nstack)
       iy = iystack(nstack)
       nstack = nstack - 1
       
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       SUBROUTINE fit
      
     &            ( x, y, ndata, sig, mwt, a, b, siga, sigb, chi2, q )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN.
!
! Given a set of data points x(1:ndata), y(1:ndata) with individual standard deviations sig(1:ndata), fit them to a straight line
! y = a + b * x by minimizing chi^2. Returned are a, b and their respective probable uncertainties siga and sigb, the chi-square
! chi2, and the goodness-of-fit probability q (that the fit would have chi^2 this large or larger). If mwt = 0 on input, then the
! standard deviations are assumed to be unavailable: q is returned as 1.0 and normalization of chi2 is to unit standard deviation
! on all points.
!__________________________________________________________________________________________________________________________________
!                     
       implicit      none
       INTEGER       mwt, ndata, i
       REAL*8        a, b, chi2, q, siga, sigb, sig(ndata), x(ndata), y(ndata)
CU     USES          gammq
       REAL*8        sigdat, ss, st2, sx, sxoss, sy, t, wt, gammq
       external      gammq
!__________________________________________________________________________________________________________________________________
!
       sx = 0.0d0
       sy = 0.0d0
       st2 = 0.0d0
       b = 0.0d0
       if (mwt .ne. 0) then
         ss = 0.0d0
         do 11 i=1,ndata
           wt = 1.0d0 / (sig(i)**2)
           ss = ss + wt
           sx = sx + x(i)*wt
           sy = sy + y(i)*wt
11       continue
       else
         do 12 i=1,ndata
           sx = sx + x(i)
           sy = sy + y(i)
12       continue
         ss = float ( ndata )
       endif
       sxoss = sx / ss
       if (mwt .ne. 0) then
         do 13 i=1,ndata
           t = (x(i) - sxoss) / sig(i)
           st2 = st2 + t*t
           b = b + t * y(i) / sig(i)
13       continue
       else
         do 14 i=1,ndata
           t = x(i) - sxoss
           st2 = st2 + t*t
           b = b + t * y(i)
14       continue
       endif
       b = b / st2
       a = (sy - sx * b) / ss
       siga = sqrt ( (1.0d0 + sx * sx / (ss * st2)) / ss)
       sigb = sqrt ( 1.0d0 / st2 )
       chi2 = 0.0d0
       if (mwt .eq. 0) then
         do 15 i=1,ndata
           chi2 = chi2 + (y(i) - a - b * x(i))**2
15       continue
         q = 1.0d0
         sigdat = sqrt ( chi2 / (ndata - 2) )
         siga = siga * sigdat
         sigb = sigb * sigdat
       else
         do 16 i=1,ndata
           chi2 = chi2 + ((y(i) - a - b * x(i)) / sig(i))**2
16       continue
         q = gammq ( 0.5d0 * (ndata - 2), 0.5d0 * chi2 )
       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function gammq

     &          ( a, x )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       real*8        a, gammq, x
CU     uses          gcf, gser
       real*8        gammcf, gamser, gln
       external      gser, gcf
!__________________________________________________________________________________________________________________________________
!
       if (x .lt. 0.0d0 .or. a .le. 0.0d0) then
         write (*,*) ' bad arguments in gammq'
         stop
       endif
       if (x .lt. a + 1.0d0)then
         call gser ( gamser, a, x, gln )
         gammq = 1.d0 - gamser
       else
         call gcf ( gammcf, a, x, gln )
         gammq = gammcf
       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine gcf
       
     &            ( gammcf, a, x, gln )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       itmax, i
       real*8        a, gammcf, gln, x, EPS, FPMIN
       parameter   ( itmax = 100, eps = 3.0d-7, fpmin = 1.0d-30)
CU     uses          gammln
       real*8        an, b, c, d, del, h, gammln
       external      gammln
!__________________________________________________________________________________________________________________________________
!
       gln = gammln(a)
       b = x + 1.0d0 - a
       c = 1.0d0 / FPMIN
       d = 1.0d0 / b
       h = d
       do 11 i=1,ITMAX
         an = -i * (i - a)
         b = b + 2.0d0
         d = an * d + b
         if (abs ( d ) .lt. FPMIN) d = FPMIN
         c = b + an / c
         if (abs ( c ) .lt. FPMIN) c = FPMIN
         d = 1.0d0 / d
         del = d * c
         h = h * del
         if (abs ( del - 1.0d0 ) .lt. EPS) goto 1
11     continue
       write (*,*) ' a too large, ITMAX too small in gcf'
       stop

1      gammcf = exp ( -x + a * log ( x ) - gln ) * h
     
       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       SUBROUTINE gser
      
     &            ( gamser, a, x, gln )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       itmax, n
       real*8        a, gamser, gln, x, EPS
       parameter   ( ITMAX = 100, EPS = 3.0d-7)
CU     uses gammln
       real*8        ap, del, sum, gammln
       external      gammln
!__________________________________________________________________________________________________________________________________
!
       gln = gammln ( a )
       if (x .le. 0.0d0)then
         if (x .lt. 0.0d0) then
           write (*,*) ' x < 0 in gser'
           stop
         endif
         gamser = 0.0d0
         return
       endif
       ap = a
       sum = 1.0d0 / a
       del = sum
       do 11 n=1,ITMAX
         ap = ap + 1.0d0
         del = del * x / ap
         sum = sum + del
         if (abs ( del ) .lt. abs ( sum ) * EPS) goto 1
11     continue
       write (*,*) ' a too large, ITMAX too small in gser'
       stop

1      gamser = sum * exp ( -x + a * log ( x ) - gln )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function gammln
    
     &         ( xx )
!__________________________________________________________________________________________________________________________________
!
! From Numerical Recipes for FORTRAN.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       REAL*8        gammln, xx
       INTEGER       j
       REAL*8        ser, stp, tmp, x, y, cof(6)
       SAVE          cof, stp
       DATA cof, stp / 76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     &                -.5395239384953d-5,2.5066282746310005d0 /
!__________________________________________________________________________________________________________________________________
!
       x = xx
       y = x
       tmp = x + 5.5d0
       tmp = (x + 0.5d0) * log ( tmp ) - tmp
       ser = 1.000000000190015d0
       do 11 j=1,6
         y = y + 1.0d0
         ser = ser + cof(j) / y
11     continue
       gammln = tmp + log ( stp * ser / x )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
