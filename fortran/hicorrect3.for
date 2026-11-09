!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
       program hicorrect3
!__________________________________________________________________________________________________________________________________
!
! Correct the temperature bias of a surface density image, with a factor measured on a grid of radiative transfer models.
!
! The correction of hicorrect2 derives its factor analytically, from a uniform structure whose shielding is an angular average.
! This one reads the factor from a table measured directly on the models by harvest_calib3.py, and it can use either one
! observable or three.
!
! One observable. The factor is a function of the surface density alone, exactly as in hicorrect2, but measured rather than
! derived. Validated on a second grid of models whose parameters lie between those of the first, it leaves 9.2% of the deficit at
! the center of a core against 16.3% for the analytic factor, and recovers 0.990 of the true mass integrated over the model
! against 0.947. It is smooth, it needs nothing but the image, and it applies to any surface density map.
!
! Three observables. The factor is a function of the surface density, of the temperature fitted at the same resolution, and of the
! contrast between that resolution and the coarsest, all of which the reconstruction provides. Those carry what one number cannot:
! whether the material of a pixel lies near the surface of a structure, where it is warm, or in the body of one, where it is cold.
! On the same validation it leaves 3.9% at the center of a core and flattens the radial residual by a factor of three and a half,
! which is what removes the ring of over-correction on the flanks of structures.
!
! The three-observable form can be used only where the models cover the observations. They heat isolated cores with the
! interstellar radiation field alone, so a field with internal sources is warmer than any model of the same column and falls
! outside the table. The fraction of pixels outside is reported, and where a pixel falls outside, the factor reverts smoothly to
! the one-observable value, so the result degrades gracefully rather than breaking. If that fraction is large the one-observable
! form should be used instead.
!
! The interpolation is bilinear in the surface density and the temperature on the regular grid of the table, and linear between
! the levels of the contrast. Cells holding too few pixels of the models to be reported carry no weight, and their weight passes
! to the one-observable curve, so the factor is a smooth function of a smooth image everywhere.
!
! An image that already carries a correction is refused, since applying one twice is meaningless.
!
! usage:  hicorrect3 <table> <image> [<coarse image> <temperature>] [-verb{0|1|2}]
!
!         <table>          the calibration written by harvest_calib3.py
!         <image>          the surface density image to correct
!         <coarse image>   the surface density at the coarsest resolution, which sets the contrast
!         <temperature>    the temperature image of the same resolution as the image being corrected
!
! output: <image>.corr3.fits, the corrected image, and <image>.cfactor3.fits, the factor applied to it
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       iotty, iolog, mx, my, mz
       parameter   ( iotty = 6, iolog = 0, mx = 400, my = 300, mz = 40 )
       logical       lunix, three
       character*6   cfitsversion, cverbose
       character*8   ctime
       character*10  cdate
       character*80  bunit, ctype1, ctype2, object, creator, history
       character*256 tabfile, mapfile, coarsefile, tempfile, arg, stem, fname
       character*80  chead(5)
       integer       nx, ny, nxi, nyi, i, j, k, iargc, narg, la, ls, nvalid, blank, ndate, npos, nfull, nnone
       integer       ntx, nty, ntz, ix, iy, iz
       real          fitsvers
       real*8        dx, dy, dxi, dyi, crpix1, crpix2, crval1, crval2, crota1, crota2, cd11, cd12, cd21, cd22, equinox, bzero
     &             , bscale, beamkw, funmin, funmax, datamn, datamx, wavekw, pix, rnan, rhead(17), ra, dec
     &             , tx(mx), ty(my), tz(mz), tfac(mx,my,mz), tspr(mx,my,mz), tnum(mx,my,mz), fone(mx), sone(mx)
     &             , ftwo(mx,my), stwo(mx,my)
     &             , xq, yq, zq, wx, wy, wz, wsum, fsum, ssum, ftot, stot, wtot, f, q(7)
       real*8,  allocatable :: smap(:,:), work(:,:), corr(:,:), sco(:,:), tmp(:,:)
       integer, allocatable :: valid(:,:)
       integer       have(mx,my)
       character*3   cnan
       equivalence ( cnan, rnan )
       data          cnan / 'NaN' /
       external      rfits, putima
!
       call osystem ( lunix )
       call ftvers ( fitsvers )
       write (cfitsversion,'(f6.3)') fitsvers

       cverbose   = '-verb1'
       tabfile    = ' '
       mapfile    = ' '
       coarsefile = ' '
       tempfile   = ' '

       narg = iargc ( )
       if (narg .lt. 2) then
         write (iotty,'(/a/)') '   HICORRECT3: usage: hicorrect3 <table> <image> [<coarse image> <temperature>] [-verb{0|1|2}]'
         stop
       endif

       do i=1,narg
         call getarg ( i, arg )
         if (arg(1:5) .eq. '-verb') then
           cverbose = arg(1:6)
         elseif (tabfile .eq. ' ') then
           tabfile = arg
         elseif (mapfile .eq. ' ') then
           mapfile = arg
         elseif (coarsefile .eq. ' ') then
           coarsefile = arg
         elseif (tempfile .eq. ' ') then
           tempfile = arg
         endif
       enddo
       three = coarsefile .ne. ' ' .and. tempfile .ne. ' '

       call when ( lunix, ctime, cdate, ndate, 4 )

       if (cverbose .ne. '-verb0') then
         write (iotty,'(/a)') '   HICORRECT3: correction of the single-temperature bias, calibrated on radiative transfer models'
       endif
!__________________________________________________________________________________________________________________________________
!
! The calibration. Each line of the table holds the center of a cell in the three observables, the number of pixels of the models
! that fell in it, and the 16th, 50th and 84th percentiles of the factor those pixels required. The cells form a regular grid in
! the first two coordinates and an irregular one in the third, so the grid is recovered by collecting the distinct values.
!
       call rdcal ( iotty, tabfile, mx, my, mz, ntx, nty, ntz, tx, ty, tz, tfac, tspr, tnum, q )

!
! The one-observable curve is the same table averaged over the temperature and the contrast, weighted by the pixels of the models.
!
       do ix=1,ntx
         fsum = 0.0d0
         ssum = 0.0d0
         wsum = 0.0d0
         do iy=1,nty
           do iz=1,ntz
             if (tnum(ix,iy,iz) .gt. 0.0d0) then
               fsum = fsum + tnum(ix,iy,iz) * tfac(ix,iy,iz)
               ssum = ssum + tnum(ix,iy,iz) * tspr(ix,iy,iz)
               wsum = wsum + tnum(ix,iy,iz)
             endif
           enddo
         enddo
         if (wsum .gt. 0.0d0) then
           fone(ix) = fsum / wsum
           sone(ix) = ssum / wsum
         else
           fone(ix) = 1.0d0
           sone(ix) = 0.0d0
         endif
       enddo

!
! The two-observable fallback, the factor against the surface density and the temperature, which replaces missing cells of the
! three-observable table. It is completed so that it is defined everywhere.
!
       call fall2 ( mx, my, mz, ntx, nty, ntz, tx, ty, tfac, tspr, tnum, ftwo, stwo, have )

       if (cverbose .ne. '-verb0') then
         write (iotty,'(a,i6,a,i4,a,i4,a,i3,a)') '   calibration: ', int ( q(1) ), ' cells on a grid of', ntx, ' x', nty
     &                                         , ' x', ntz, ' in surface density, temperature and contrast'
         write (iotty,'(a,1pe10.3,a)')           '   measured on ', q(2), ' pixels of the models'
         write (iotty,'(a,f6.2,a,f6.2,a,f5.1,a,f5.1,a,f7.3,a,f6.3)')
     &                                           '   covering log10 Sigma', tx(1), ' to', tx(ntx), ', T', ty(1), ' to', ty(nty)
     &                                         , ' K, contrast', tz(1), ' to', tz(ntz)
         if (three) then
           write (iotty,'(a)')                   '   using the surface density, the temperature and the contrast'
         else
           write (iotty,'(a)')                   '   using the surface density alone'
         endif
       endif
!__________________________________________________________________________________________________________________________________
!
! The image to correct.
!
       call getfitshead ( mapfile(1:len_trim(mapfile)), nx, ny, dx, dy, bunit )

       allocate ( smap(nx,ny), valid(nx,ny), corr(nx,ny), work(nx,ny), sco(nx,ny), tmp(nx,ny) )

       la = len_trim ( mapfile )
       call rfits ( nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, smap, dx, dy, object, ra, dec, mapfile(1:la)
     &            , iotty, iolog, creator, beamkw, funmin, funmax, blank, crota1, crota2, cd11, cd12, cd21, cd22, equinox, bzero
     &            , bscale, wavekw, datamn, datamx, history, cverbose )

       if (creator(1:9) .eq. 'HICORRECT') then
         write (iotty,'(/a/a/)') '   HICORRECT3: ERROR: '''//mapfile(1:la)//''' has already been corrected.'
     &                         , '   The correction is not idempotent and must be applied once only.'
         stop
       endif

       pix = dx
       ls = la
       if (mapfile(max(ls-4,1):ls) .eq. '.fits') ls = ls - 5
       stem = mapfile(1:ls)
       if (bunit .eq. ' ') bunit = 'cm-2'

       rhead( 1) = crpix1
       rhead( 2) = crpix2
       rhead( 3) = crval1
       rhead( 4) = crval2
       rhead( 5) = crota1
       rhead( 6) = crota2
       rhead( 7) = cd11
       rhead( 8) = cd12
       rhead( 9) = cd21
       rhead(10) = cd22
       rhead(11) = dx
       rhead(12) = dy
       rhead(13) = equinox
       rhead(14) = bzero
       rhead(15) = bscale
       rhead(16) = beamkw
       rhead(17) = 0.0d0
       chead( 1) = ctype1
       chead( 2) = ctype2
       chead( 3) = object
       chead( 4) = 'HICORRECT3'
       chead( 5) = history
!
! The coarse image and the temperature, if the three-observable form was asked for.
!
       if (three) then
         call getfitshead ( coarsefile(1:len_trim(coarsefile)), nxi, nyi, dxi, dyi, bunit )
         if (nxi .ne. nx .or. nyi .ne. ny) then
           write (iotty,'(/a/)') '   HICORRECT3: ERROR: the coarse image does not match the image to correct.'
           stop
         endif
         la = len_trim ( coarsefile )
         call rfits ( nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, sco, dxi, dyi, object, ra, dec
     &              , coarsefile(1:la), iotty, iolog, creator, beamkw, funmin, funmax, blank, crota1, crota2, cd11, cd12, cd21
     &              , cd22, equinox, bzero, bscale, wavekw, datamn, datamx, history, cverbose )

         call getfitshead ( tempfile(1:len_trim(tempfile)), nxi, nyi, dxi, dyi, bunit )
         if (nxi .ne. nx .or. nyi .ne. ny) then
           write (iotty,'(/a/)') '   HICORRECT3: ERROR: the temperature image does not match the image to correct.'
           stop
         endif
         la = len_trim ( tempfile )
         call rfits ( nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, tmp, dxi, dyi, object, ra, dec
     &              , tempfile(1:la), iotty, iolog, creator, beamkw, funmin, funmax, blank, crota1, crota2, cd11, cd12, cd21
     &              , cd22, equinox, bzero, bscale, wavekw, datamn, datamx, history, cverbose )
         bunit = 'cm-2'
       endif
!__________________________________________________________________________________________________________________________________
!
! The coverage. A pixel that is not positive carries no information, and in the three-observable form the coarse image and the
! temperature must be defined as well.
!
       nvalid = 0
       do j=1,ny
         do i=1,nx
           valid(i,j) = 0
           if (smap(i,j) .gt. 0.0d0 .and. .not.isnan ( smap(i,j) )) then
             if (three) then
               if (sco(i,j) .gt. 0.0d0 .and. tmp(i,j) .gt. 0.0d0 .and. .not.isnan ( sco(i,j) )
     &             .and. .not.isnan ( tmp(i,j) )) valid(i,j) = 1
             else
               valid(i,j) = 1
             endif
           endif
           if (valid(i,j) .eq. 1) nvalid = nvalid + 1
         enddo
       enddo
       if (cverbose .ne. '-verb0') write (iotty,'(a,i10,a)') '   ', nvalid, ' valid pixels'
!__________________________________________________________________________________________________________________________________
!
! The factor. In the one-observable form it is read from the curve; in the three-observable form it is interpolated bilinearly in
! the surface density and the temperature within each of the two levels of contrast that bracket the pixel, and linearly between
! them, with the weight of the cells that exist. What weight is missing goes to the two-observable fallback, which is defined
! everywhere and, through the temperature, distinguishes the radiation fields that the one-observable curve averages over.
!
       nfull = 0
       nnone = 0
       npos  = 0
       do j=1,ny
         do i=1,nx
           corr(i,j) = 1.0d0
           work(i,j) = 0.0d0
           if (valid(i,j) .eq. 0) cycle
           xq = log10 ( max ( smap(i,j), 1.0d0 ) )
           call locate ( ntx, tx, xq, ix, wx )
           if (.not.three) then
             corr(i,j) = ( 1.0d0 - wx ) * fone(ix) + wx * fone(ix+1)
             work(i,j) = ( 1.0d0 - wx ) * sone(ix) + wx * sone(ix+1)
             npos = npos + 1
             cycle
           endif
           yq = tmp(i,j)
           zq = log10 ( max ( smap(i,j), 1.0d0 ) / max ( sco(i,j), 1.0d0 ) )
           call locate ( ntz, tz, zq, iz, wz )
           call locate ( nty, ty, yq, iy, wy )
           ftot = ( 1.0d0 - wx ) * ( 1.0d0 - wy ) * ftwo(ix,iy)   + wx * ( 1.0d0 - wy ) * ftwo(ix+1,iy)
     &          + ( 1.0d0 - wx ) * wy * ftwo(ix,iy+1)             + wx * wy * ftwo(ix+1,iy+1)
           stot = ( 1.0d0 - wx ) * ( 1.0d0 - wy ) * stwo(ix,iy)   + wx * ( 1.0d0 - wy ) * stwo(ix+1,iy)
     &          + ( 1.0d0 - wx ) * wy * stwo(ix,iy+1)             + wx * wy * stwo(ix+1,iy+1)

           call cell ( mx, my, mz, tfac, tspr, tnum, ix, iy, iz, wx, wy, wz, fsum, ssum, wsum )

           wtot = min ( max ( wsum, 0.0d0 ), 1.0d0 )
           if (wtot .gt. 0.0d0) then
             f = wtot * fsum / wsum + ( 1.0d0 - wtot ) * ftot
             corr(i,j) = f
             work(i,j) = wtot * ssum / wsum + ( 1.0d0 - wtot ) * stot
           else
             corr(i,j) = ftot
             work(i,j) = stot
           endif
           if (wtot .gt. 0.999d0) nfull = nfull + 1
           if (wtot .lt. 0.001d0) nnone = nnone + 1
           npos = npos + 1
         enddo
       enddo

       if (cverbose .ne. '-verb0' .and. three) then
         write (iotty,'(a,f6.1,a,f6.1,a)') '   covered by the table:', 100.0d0 * dble ( nfull ) / dble ( max ( npos, 1 ) )
     &                                   , '% of the pixels fully,', 100.0d0 * dble ( nnone ) / dble ( max ( npos, 1 ) )
     &                                   , '% not at all'
       endif
!__________________________________________________________________________________________________________________________________
!
! The corrected image and the factor applied to it.
!
       fsum = 0.0d0
       ssum = 0.0d0
       do j=1,ny
         do i=1,nx
           if (valid(i,j) .eq. 1) then
             fsum = fsum + smap(i,j) * corr(i,j)
             ssum = ssum + smap(i,j)
             smap(i,j) = smap(i,j) * corr(i,j)
           else
             smap(i,j) = 0.0d0
           endif
         enddo
       enddo
       if (cverbose .ne. '-verb0') write (iotty,'(a,f8.4)') '   ratio of the sums, corrected to uncorrected:'
     &                                                    , fsum / max ( ssum, 1.0d-30 )

       fname = stem(1:ls)//'.corr3.fits'
       call putima ( iotty, cverbose, cfitsversion, nx, ny, smap, valid, rnan, fname, bunit, rhead, chead, blank, cdate, ctime )
       if (cverbose .ne. '-verb0') write (iotty,'(a)') '   --> '''//fname(1:len_trim(fname))//''''

       fname = stem(1:ls)//'.cfactor3.fits'
       call putima ( iotty, cverbose, cfitsversion, nx, ny, corr, valid, rnan, fname, ' ', rhead, chead, blank, cdate, ctime )
       if (cverbose .ne. '-verb0') write (iotty,'(a)') '   --> '''//fname(1:len_trim(fname))//''''

       deallocate ( smap, valid, corr, work, sco, tmp )

       if (cverbose .ne. '-verb0') write (iotty,'(/a/)') '   HICORRECT3 DONE.'
       stop
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine fall2

     &            ( mx, my, mz, ntx, nty, ntz, tx, ty, tfac, tspr, tnum, ftwo, stwo, have )
!__________________________________________________________________________________________________________________________________
!
! The two-observable fallback. Each cell is the average of the factor over the levels of contrast that exist at that surface
! density and temperature, weighted by the pixels of the models; where only filled cells (with a pixel count of zero) exist, it is
! their plain average. The table is then completed: gaps along the temperature are interpolated at fixed surface density, then
! gaps along the surface density at fixed temperature, and any cell still empty takes the value of the nearest cell, with
! distances counted in steps of the table and ties resolved in favor of the lower surface density, then the lower temperature.
! The procedure is that of hicorrect3.py, so that both give the same factor.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       mx, my, mz, ntx, nty, ntz, ix, iy, iz, j1, j2, k, cnt, ib, jb, ie, je, d2, best
       integer       have(mx,my)
       real*8        tx(mx), ty(my), tfac(mx,my,mz), tspr(mx,my,mz), tnum(mx,my,mz), ftwo(mx,my), stwo(mx,my)
       real*8        wn, fw, sw, fp, sp, t
       integer, allocatable :: done(:,:)

       do ix=1,ntx
         do iy=1,nty
           wn  = 0.0d0
           fw  = 0.0d0
           sw  = 0.0d0
           fp  = 0.0d0
           sp  = 0.0d0
           cnt = 0
           do iz=1,ntz
             if (tnum(ix,iy,iz) .lt. 0.0d0) cycle
             cnt = cnt + 1
             fp  = fp + tfac(ix,iy,iz)
             sp  = sp + tspr(ix,iy,iz)
             if (tnum(ix,iy,iz) .gt. 0.0d0) then
               wn = wn + tnum(ix,iy,iz)
               fw = fw + tnum(ix,iy,iz) * tfac(ix,iy,iz)
               sw = sw + tnum(ix,iy,iz) * tspr(ix,iy,iz)
             endif
           enddo
           have(ix,iy) = 0
           ftwo(ix,iy) = 1.0d0
           stwo(ix,iy) = 0.0d0
           if (wn .gt. 0.0d0) then
             ftwo(ix,iy) = fw / wn
             stwo(ix,iy) = sw / wn
             have(ix,iy) = 1
           elseif (cnt .gt. 0) then
             ftwo(ix,iy) = fp / dble ( cnt )
             stwo(ix,iy) = sp / dble ( cnt )
             have(ix,iy) = 1
           endif
         enddo
       enddo
!
! Gaps along the temperature, at fixed surface density.
!
       do ix=1,ntx
         j1 = 0
         do iy=1,nty
           if (have(ix,iy) .eq. 0) cycle
           if (j1 .gt. 0 .and. iy - j1 .gt. 1) then
             do k=j1+1,iy-1
               t = ( ty(k) - ty(j1) ) / ( ty(iy) - ty(j1) )
               ftwo(ix,k) = ( 1.0d0 - t ) * ftwo(ix,j1) + t * ftwo(ix,iy)
               stwo(ix,k) = ( 1.0d0 - t ) * stwo(ix,j1) + t * stwo(ix,iy)
               have(ix,k) = 1
             enddo
           endif
           j1 = iy
         enddo
       enddo
!
! Gaps along the surface density, at fixed temperature.
!
       do iy=1,nty
         j1 = 0
         do ix=1,ntx
           if (have(ix,iy) .eq. 0) cycle
           if (j1 .gt. 0 .and. ix - j1 .gt. 1) then
             do k=j1+1,ix-1
               t = ( tx(k) - tx(j1) ) / ( tx(ix) - tx(j1) )
               ftwo(k,iy) = ( 1.0d0 - t ) * ftwo(j1,iy) + t * ftwo(ix,iy)
               stwo(k,iy) = ( 1.0d0 - t ) * stwo(j1,iy) + t * stwo(ix,iy)
               have(k,iy) = 1
             enddo
           endif
           j1 = ix
         enddo
       enddo
!
! The cells still empty take the value of the nearest defined cell. Only cells defined before this step are sources, so the
! result does not depend on the order in which the empty cells are visited.
!
       allocate ( done(ntx,nty) )
       do ix=1,ntx
         do iy=1,nty
           done(ix,iy) = have(ix,iy)
         enddo
       enddo
       do ie=1,ntx
         do je=1,nty
           if (done(ie,je) .eq. 1) cycle
           best = huge ( best )
           ib = 0
           jb = 0
           do ix=1,ntx
             do iy=1,nty
               if (done(ix,iy) .eq. 0) cycle
               d2 = ( ix - ie ) ** 2 + ( iy - je ) ** 2
               if (d2 .lt. best) then
                 best = d2
                 ib = ix
                 jb = iy
               endif
             enddo
           enddo
           if (ib .gt. 0) then
             ftwo(ie,je) = ftwo(ib,jb)
             stwo(ie,je) = stwo(ib,jb)
             have(ie,je) = 1
           endif
         enddo
       enddo
       deallocate ( done )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine locate

     &            ( n, a, v, i, w )
!__________________________________________________________________________________________________________________________________
!
! The pair of values of the sorted list a that brackets v, and the fractional position of v between them.
!
! The coordinates of the calibration are not a complete regular sequence: a cell holding too few pixels of the models is absent
! from the table altogether, so the distinct values have gaps in them.  The bracketing pair is therefore found by searching the
! values that are present rather than by arithmetic on a step, and the weight uses the spacing of that pair, which is exact where
! the sequence is regular and interpolates across a gap where it is not.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       n, i, lo, hi, mid
       real*8        a(n), v, w

       if (v .le. a(1)) then
         i = 1
         w = 0.0d0
         return
       elseif (v .ge. a(n)) then
         i = n - 1
         w = 1.0d0
         return
       endif

       lo = 1
       hi = n
  10   if (hi - lo .gt. 1) then
         mid = ( lo + hi ) / 2
         if (v .ge. a(mid)) then
           lo = mid
         else
           hi = mid
         endif
         goto 10
       endif
       i = lo
       w = ( v - a(i) ) / max ( a(i+1) - a(i), 1.0d-12 )
       w = min ( max ( w, 0.0d0 ), 1.0d0 )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine cell

     &            ( mx, my, mz, tfac, tspr, tnum, ix, iy, iz, wx, wy, wz, fsum, ssum, wsum )
!__________________________________________________________________________________________________________________________________
!
! The factor and its spread, interpolated over the eight cells that surround a pixel in the three coordinates of the table, using
! only those that exist. wsum is the weight of the cells that were found, and is unity where the table is complete.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       mx, my, mz, ix, iy, iz, ki, kj, kk
       real*8        tfac(mx,my,mz), tspr(mx,my,mz), tnum(mx,my,mz), wx, wy, wz, fsum, ssum, wsum, w

       fsum = 0.0d0
       ssum = 0.0d0
       wsum = 0.0d0
       do kk=0,1
         do kj=0,1
           do ki=0,1
             if (tnum(ix+ki,iy+kj,iz+kk) .lt. 0.0d0) cycle
             w = 1.0d0
             if (ki .eq. 0) then
               w = w * ( 1.0d0 - wx )
             else
               w = w * wx
             endif
             if (kj .eq. 0) then
               w = w * ( 1.0d0 - wy )
             else
               w = w * wy
             endif
             if (kk .eq. 0) then
               w = w * ( 1.0d0 - wz )
             else
               w = w * wz
             endif
             fsum = fsum + w * tfac(ix+ki,iy+kj,iz+kk)
             ssum = ssum + w * tspr(ix+ki,iy+kj,iz+kk)
             wsum = wsum + w
           enddo
         enddo
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine rdcal

     &            ( iotty, tabfile, mx, my, mz, ntx, nty, ntz, tx, ty, tz, tfac, tspr, tnum, q )
!__________________________________________________________________________________________________________________________________
!
! Read the calibration written by harvest_calib3.py. Lines beginning with an exclamation mark are comments. Each other line holds
! the center of a cell in the three observables, the number of pixels of the models that fell in it, and the 16th, 50th and 84th
! percentiles of the factor those pixels required. The distinct values of each coordinate are collected first, so that the cells
! can be placed on the grid they came from, and the file is read twice. A cell with a pixel count of zero was interpolated along
! the temperature by harvest_calib3.py --fill: it takes part in the interpolation of the factor but not in the one-observable
! curve. Cells absent from the table are marked by a negative count.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       iotty, mx, my, mz, ntx, nty, ntz, i, j, k, io, ncell
       character*256 tabfile, line
       real*8        tx(mx), ty(my), tz(mz), tfac(mx,my,mz), tspr(mx,my,mz), tnum(mx,my,mz), q(7)
       real*8        x, y, z, n, f16, f50, f84, tot

       ntx = 0
       nty = 0
       ntz = 0
       ncell = 0
       tot = 0.0d0

       open ( 21, file=tabfile, status='old', err=901 )
  10   read ( 21, '(a)', end=20 ) line
       if (line(1:1) .eq. '!' .or. len_trim ( line ) .eq. 0) goto 10
       read ( line, *, err=10 ) x, y, z, n, f16, f50, f84
       call addval ( mx, ntx, tx, x )
       call addval ( my, nty, ty, y )
       call addval ( mz, ntz, tz, z )
       ncell = ncell + 1
       tot = tot + n
       goto 10
  20   rewind ( 21 )

       call sortval ( ntx, tx )
       call sortval ( nty, ty )
       call sortval ( ntz, tz )

       do k=1,ntz
         do j=1,nty
           do i=1,ntx
             tfac(i,j,k) = 1.0d0
             tspr(i,j,k) = 0.0d0
             tnum(i,j,k) =-1.0d0
           enddo
         enddo
       enddo

  30   read ( 21, '(a)', end=40 ) line
       if (line(1:1) .eq. '!' .or. len_trim ( line ) .eq. 0) goto 30
       read ( line, *, err=30 ) x, y, z, n, f16, f50, f84
       call findval ( ntx, tx, x, i )
       call findval ( nty, ty, y, j )
       call findval ( ntz, tz, z, k )
       if (i .gt. 0 .and. j .gt. 0 .and. k .gt. 0) then
         tfac(i,j,k) = f50
         tspr(i,j,k) = 0.5d0 * ( f84 - f16 ) / max ( f50, 1.0d-12 )
         tnum(i,j,k) = n
       endif
       goto 30
  40   close ( 21 )

       q(1) = dble ( ncell )
       q(2) = tot
       return

 901   write (iotty,'(/a/)') '   HICORRECT3: ERROR: cannot open the calibration '''//tabfile(1:len_trim(tabfile))//''''
       stop
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine addval

     &            ( m, n, a, v )
!__________________________________________________________________________________________________________________________________
!
! Add v to the list a if it is not already there, to within a tolerance that separates neighboring cells. A table with more
! distinct values than the list can hold stops the program, because the cells beyond the list would otherwise be lost silently.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       m, n, i
       real*8        a(m), v

       do i=1,n
         if (abs ( a(i) - v ) .lt. 1.0d-6) return
       enddo
       if (n .ge. m) then
         write (6,'(/a,i5,a/)') '   HICORRECT3: ERROR: the calibration has more than', m, ' distinct values of one coordinate'
         stop
       endif
       n = n + 1
       a(n) = v

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine sortval

     &            ( n, a )
!__________________________________________________________________________________________________________________________________
!
! Sort the list a into increasing order. The lists are short, so the simplest sort is the right one.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       n, i, j
       real*8        a(n), t

       do i=1,n-1
         do j=i+1,n
           if (a(j) .lt. a(i)) then
             t    = a(i)
             a(i) = a(j)
             a(j) = t
           endif
         enddo
       enddo

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine findval

     &            ( n, a, v, i )
!__________________________________________________________________________________________________________________________________
!
! The index of v in the sorted list a, or zero if it is not there.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       n, i, k
       real*8        a(n), v

       i = 0
       do k=1,n
         if (abs ( a(k) - v ) .lt. 1.0d-6) then
           i = k
           return
         endif
       enddo

       return
       end
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function planckf 
         
     &          ( tdust, wave, hplanckanck, kboltzltz, speedolight )
!__________________________________________________________________________________________________________________________________
!
! Planck function per unit frequency, in cgs units, for a wavelength in µm.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       real*8        planckf, tdust, wave, frequency, xexp, hplanckanck, kboltzltz, speedolight
!__________________________________________________________________________________________________________________________________
!
       frequency = speedolight / ( wave * 1.0d-4 )
       xexp = hplanckanck * frequency / ( kboltzltz * tdust )
       xexp = min ( max ( xexp, 1.0d-8 ), 7.0d2 )
       planckf = 2.0d0 * hplanckanck * frequency**3 / speedolight**2 / ( exp ( xexp ) - 1.0d0 )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function invplanck

     &          ( b, wave, hplanck, kboltz, speedolight )
!__________________________________________________________________________________________________________________________________
!
! The temperature whose Planck function equals b at the wavelength wave in um. The inverse of the function below.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       real*8        invplanck, b, wave, hplanck, kboltz, speedolight, nu, pref

       nu = speedolight / (wave * 1.0d-4)
       pref = 2.0d0 * hplanck * nu**3 / speedolight**2
       invplanck = hplanck * nu / (kboltz * log ( 1.0d0 + pref / max ( b, 1.0d-40 ) ))

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine interp

     &            ( n, x, y, xq, yq )
!__________________________________________________________________________________________________________________________________
!
! Linear interpolation of y against x in the logarithm of x, for a table whose x increases with the index. Beyond either end the
! value of that end is returned rather than an extrapolation, because the relation is only known over the range that was tabulated.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       n, lo, hi, mid
       real*8        x(n), y(n), xq, yq, w

       if (xq .le. x(1)) then
         yq = y(1)
         return
       endif
       if (xq .ge. x(n)) then
         yq = y(n)
         return
       endif

       lo = 1
       hi = n
  100  if (hi - lo .gt. 1) then
         mid = ( lo + hi ) / 2
         if (x(mid) .gt. xq) then
           hi = mid
         else
           lo = mid
         endif
         goto 100
       endif

       w  = ( log ( xq ) - log ( x(lo) ) ) / ( log ( x(hi) ) - log ( x(lo) ) )
       yq = y(lo) + w * ( y(hi) - y(lo) )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine rdcfg 
         
     &            ( cfgfile, mband, nband, wave, beams, offset, bandfile, maskfile, pmflag )
!__________________________________________________________________________________________________________________________________
!
! Read the hires configuration file '+hires.cfg'. Its columns are
!
!    PM  WAVE  BEAM  OFFSET  IMAGE  [MASKIMAGE]
!
! with PM the '+' or '-' flag stating whether hires used that waveband, WAVE in µm, BEAM the half-maximum beam width in arcsec,
! OFFSET the zero-point offset in MJy per steradian to be added to the image, and the two file names. Reading the same file that
! hires read removes any second, independently maintained list of wavebands. The extension may be omitted from the file names, as
! is customary in that file.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       logical       there
       character*(*) cfgfile
       character*1   pmflag(mband)
       character*256 name, bandfile(mband), maskfile(mband), fields(6)
       character*512 line, rest
       integer       mband, nband, i, j, n, ios
       real*8        wave(mband), beams(mband), offset(mband)
!__________________________________________________________________________________________________________________________________
!
       open ( 21, file=cfgfile, status='old', iostat=ios )
       if (ios .ne. 0) then
         write (6,'(a,a)') ' RDCFG: cannot open ', cfgfile
         stop
       endif

       nband = 0
  10   read (21,'(a)',end=90) line
       if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) goto 10
       nband = nband + 1
       if (nband .gt. mband) stop ' RDCFG: too many wavebands, raise mband'
       maskfile(nband) = ' '
       do n = 1, 6
         fields(n) = ' '
       enddo
       rest = trim(line)
       n = 0
       i = 1
       do
! skip whitespace
         do while (i <= len_trim(rest) .and. rest(i:i) == ' ')
           i = i + 1
         enddo
         if (i > len_trim(rest)) exit
         n = n + 1
         if (n > 6) then
           n = 6
           exit
         endif
! find end of field
         j = i
         do while (j <= len_trim(rest) .and. rest(j:j) /= ' ')
           j = j + 1
         enddo
         fields(n) = rest(i:j-1)
         i = j
       enddo
       if (n .lt. 5) stop ' RDCFG: a data line must hold at least five fields: PM WAVE BEAM OFFSET IMAGE [MASKIMAGE]'

       read(fields(1), *) pmflag(nband)
       read(fields(2), *) wave(nband)
       read(fields(3), *) beams(nband)
       read(fields(4), *) offset(nband)
       bandfile(nband) = fields(5)
       if (n .ge. 6) maskfile(nband) = fields(6)

       if (pmflag(nband) .ne. '+' .and. pmflag(nband) .ne. '-') stop
     &   ' RDCFG: the first column of a data line must be + or -; point the second argument at the +hires.cfg of the hires run'

       name = bandfile(nband)
       inquire ( file=name, exist=there )
       if (.not. there) bandfile(nband) = name(1:len_trim(name))//'.fits'

       if (len_trim(maskfile(nband)) .gt. 0) then
         name = maskfile(nband)
         inquire ( file=name, exist=there )
         if (.not. there) maskfile(nband) = name(1:len_trim(name))//'.fits'
       endif

       goto 10
  90   close ( 21 )

       if (nband .eq. 0) stop ' RDCFG: the configuration contains no waveband.'

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine putima 
         
     &            ( iotty, cverbose, cfitsversion, nx, ny, fun, valid, rnan, fname, bunit, rhead, chead, blank, cdate, ctime )
!__________________________________________________________________________________________________________________________________
!
! Write one output image. The pixels outside the coverage are made undefined first, because the values that the filling put there
! were never results. The header of the surface density image is carried in rhead and chead and is written back unchanged, so that
! every output image has the astrometry of the input:
!
!   rhead( 1) CRPIX1   rhead( 2) CRPIX2   rhead( 3) CRVAL1   rhead( 4) CRVAL2   rhead( 5) CROTA1   rhead( 6) CROTA2
!   rhead( 7) CD1_1    rhead( 8) CD1_2    rhead( 9) CD2_1    rhead(10) CD2_2    rhead(11) CDELT1   rhead(12) CDELT2
!   rhead(13) EQUINOX  rhead(14) BZERO    rhead(15) BSCALE   rhead(16) BEAM     rhead(17) WAVE
!   chead( 1) CTYPE1   chead( 2) CTYPE2   chead( 3) OBJECT   chead( 4) CREATOR  chead( 5) HISTORY
!
! The two pixel sizes are in arcsec, the wavelength in µm, and BUNIT is passed separately because it differs from one image to
! the next. The smallest and the largest value of the image are computed by wfits itself.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) cverbose, cfitsversion, fname, bunit, cdate, ctime
       character*(*) chead(5)
       integer       i, j, iotty, nx, ny, blank, valid(nx,ny)
       real*8        fun(nx,ny), rnan, rhead(17), datamn, datamx
       real*8, allocatable :: copy(:,:)
       external      wfits
!__________________________________________________________________________________________________________________________________
!
! The image is copied before the pixels outside the coverage are made undefined, so that writing an image never alters it. The
! intermediate images written with -diag2 and -diag3 are needed by the steps that follow, and would otherwise be left full of
! undefined values by their own writing.

       allocate ( copy(nx,ny) )

       do j=1,ny
         do i=1,nx
           if (valid(i,j) .eq. 0) then
             copy(i,j) = rnan
           else
             copy(i,j) = fun(i,j)
           endif
         enddo
       enddo

       datamn = 0.0d0
       datamx = 0.0d0

       if (cverbose .eq. '-verb2') write (iotty,'(a)') '   Writing output image '''//fname//''''

       call wfits ( cfitsversion, nx, ny, bunit, chead(1), chead(2), rhead(1), rhead(2), rhead(3), rhead(4), copy, rhead(11)
     &            , rhead(12), chead(3), rhead(3), rhead(4), fname, cdate, ctime, chead(4), rhead(16), blank, rhead(5), rhead(6)
     &            , rhead(7), rhead(8), rhead(9), rhead(10), rhead(13), rhead(14), rhead(15), rhead(17), datamn, datamx, chead(5) )

       deallocate ( copy )

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
       unit = 20
       group = 1
       rw = 0
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

       if (stat .ne. 0) write (iotty,'(a)') '   Error reading FITS file'

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

! An undefined pixel keeps its undefined value here, instead of being replaced by zero as it is in the other utilities that use
! this routine. HICORRECT recognizes the part of an image that the telescope did not cover by exactly those undefined values, so
! they have to survive the reading. The smallest and the largest value of the image are taken over the defined pixels only.

       do j=1,ny
         do i=1,nx
           if (.not.isnan ( fun(i,j) )) then
             funmin = min ( funmin, fun(i,j) )
             funmax = max ( funmax, fun(i,j) )
             if (int ( min ( abs ( fun(i,j) ), 1.0d9 ) ) .eq. blank .and. blank .ne. 0) then
               read (cnan,*) fun(i,j)
             endif
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

       status = 0
       blocksze = 1
       simple = .true.
       extend = .false.
       group = 1
       naxis = 2
       pcount = 0
       gcount = 1
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
cc!      bzero = 0.5d0 * (funmin + funmax)
cc!      bscale = abs ( funmax - bzero ) / 32767.0d0

       bitpix = -32     !<-- single precision floating-point values
       bzero = 0.0d0
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
