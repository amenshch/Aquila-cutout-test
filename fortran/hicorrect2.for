
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       program hicorrect
!__________________________________________________________________________________________________________________________________
!
! HIRES • CORRECTION OF THE SINGLE-TEMPERATURE BIAS OF SURFACE DENSITY IMAGES
!
! A surface density image is built by fitting one modified blackbody to the wavebands of each line of sight and dividing the
! reference waveband by it. The correction applies to any such image, whether made by hires or by any other method, at any
! resolution; it needs only the image and the list of wavebands that were fitted to make it. A real line of sight is not isothermal: its dust is warm near the surface of the cloud and cold in the
! interior, the single fitted temperature is weighted toward the warm dust that dominates the emission, and the surface density is
! therefore too low wherever the temperature varies along the line of sight. The deficit exceeds a factor of two in the densest
! material.
!
! Where the deficit is. The reconstruction of hires forms a base image at the coarsest resolution and adds to it one differential
! increment per waveband, so that the image at any resolution is
!
!   S(theta) = S(theta_coarse) + sum of the increments finer than theta.
!
! The increments are differences between an image and a smoothed copy of itself, so they very nearly cancel when integrated over
! an extended region: on the model grid they carry less than three percent of the true surface density. The whole integrated
! deficit therefore sits in the base image, and every image of the reconstruction has nearly the same fractional deficit.
!
! How the correction is applied. One surface density image is corrected per run. It is multiplied by the correction factor
! evaluated on that image itself, pixel by pixel. This works for an image of any resolution; the default of the paper is the hires
! image at 13.5 arcsec. If the temperature image fitted together with it is given, the dust temperature implied by the corrected
! surface density is also written. That temperature image must have the same resolution and the same pixels as the surface density
! image, because the temperature at each pixel is recomputed with the factor of that pixel.
!
! How the correction is computed. For a line of sight of true surface density S, the material at fractional depth u is shielded by
! the effective attenuating column
!
!   s(u) = -ln < exp ( -kapheat S l(u,psi) / 2 ) > / ( kapheat / 2 ),
!
! the angular average being taken over nang directions psi and l(u,psi) being the path from the depth u to the boundary in units of
! the half thickness. It is evaluated for a uniform sphere and for a uniform plane-parallel layer, and the two are averaged. The
! dust temperature of that material follows the relation
!
!   T(s) = tfloor + ( tsurf - tfloor ) ( 1 + s / sstar ) ** (-gamma),
!
! measured on radmc-3d models of Bonnor-Ebert spheres embedded in uniform clouds and of isolated cylinders. The emission of the line
! of sight is integrated over depth in each of the wavebands that the reconstruction fitted, one modified blackbody is fitted to the
! result exactly as the reconstruction fits it, and the surface density that the fit returns is compared with S. The ratio is the
! factor by which the fitted surface density has to be multiplied. It is tabulated once against the fitted surface density, which
! is the quantity the image holds, and interpolated onto the image. At low surface densities the analytic layer is replaced by a
! radmc-3d model of a diffuse cloud, whose surface layers are warmer than the relation predicts (see below).
!
! The geometry. The sphere and the layer are the two limits of a uniform structure. Their factors differ by up to 10%, and by less
! than 5% outside the range from 10^21.4 to 10^22.1 cm^-2, and their average is used. The arrangement of the material along the
! line of sight is not known. Where a structure lies in a weak cloud, it is corrected too little; where it lies near a face of a
! dense cloud, it is corrected too much. On radiative transfer models, these errors range from -13% to +21%.
!
! The correction must be applied once. It is not idempotent: applied to an image it has already corrected it would raise the values
! again. Every output image is stamped with CREATOR = 'HICORRECT', and the program refuses any input that carries it. The output
! images also record the radiation field of the correction (HIG0), the wavebands fitted (HIBANDS), and the four constants of the
! temperature relation (HITMIN, HIT0, HISSTAR, HIGAMMA).
!
! usage:  hicorrect2 <image> [-temp=<file>] [-bands=<list>] [-g0=<field>] [-verb{0|1|2}]
!
!         <image>      the surface density image to be corrected; the factor is evaluated on it and it is multiplied by it
!                      (the file names of <image> and <file> may be given without the extension '.fits', which is then added)
!         <file>       the dust temperature image fitted together with <image>, optional; it must have the same resolution and
!                      the same pixels as <image>
!         <list>       the wavelengths in um of the wavebands fitted to make <image>, separated by commas, for example
!                      -bands=250,350,500; 160,250,350,500 by default. The longest one is the reference waveband.
!         <field>      the radiation field of the cloud in Habing units, as estimated by g0_requirement.py; 10 by default
!
! output: written to the current directory, whatever directory the inputs come from
!         <image>.correct.fits, the corrected image
!         <image>.cfactor.fits, the image of the correction factors
!         <file>.correct.fits, the dust temperature implied by the corrected surface density, if -temp=<file> is given
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       iotty, iolog, mband, ntab, ndep, ntfit, nang
       parameter   ( iotty = 6, iolog = 0, mband = 10, ntab = 200, ndep = 512, ntfit = 3000, nang = 160 )
       integer       ngx, ngf, lt
       parameter   ( ngx = 24, ngf = 9 )
       real*8        gtab(ngx,ngf), gfld(ngf), gx0, gdx, g0, fieldfac
       parameter   ( gx0 = 20.95d0, gdx = 0.1d0 )
       logical       lunix
       character*6   cfitsversion, cverbose
       character*8   ctime
       character*10  cdate
       character*80  bunit, ctype1, ctype2, object, creator, history
       character*256 mapfile, arg, stem, fname, tempfile
       character*80  chead(5)
       integer       nx, ny, nxi, nyi, i, j, k, l, ib, iargc, narg, la, ls, nfit, nvalid, blank, ndate, igeo, ios
       real          fitsvers
       real*8        wfit(mband), waveref, bb, invplanck
     &             , dx, dy, dxi, dyi, crpix1, crpix2, crval1, crval2, crota1, crota2, cd11, cd12, cd21, cd22, equinox, bzero
     &             , bscale, beamkw, beamsd, funmin, funmax, datamn, datamx, wavekw, pix, rnan, rhead(17), ra, dec
     &             , tfloor, tsurf, sstar, gamma, kappa0, nu0, betap, dtg, muh2, amu, hplanck, kboltz, speedolight, kapheat, pi
     &             , smin, smax, sgrid(ntab), ffit(ntab), fcor(ntab), tgrid(ntfit), seff(ndep), relerr(mband)
     &             , ffits(ntab), fcors(ntab), xg, fsx, fcx, dd, d0, xb0, xb1, wb, ddx(ntab), cirx(12), cirf(12)
     &             , cnu(mband), obs(mband), model(mband), wgt(mband), amp, chi2, chibest, tbest, sig, tt, ss, f, planckf
     &             , t0, zz, psi, wps, wsum, jsph, jlay, lsph, dlay, xmu, ssph, slay
       parameter   ( tfloor = 7.049531d0, tsurf = 22.866783d0, sstar = 2.389668d+21, gamma = 1.596874d0, kappa0 = 10.0d0
     &             , nu0 = 1.0d12, betap = 2.0d0, dtg = 1.0d-2, muh2 = 2.8d0, amu = 1.66053906660d-24, hplanck = 6.62607015d-27
     &             , kboltz = 1.380649d-16, speedolight = 2.99792458d10, smin = 1.0d20, smax  = 1.0d24
     &             , kapheat = 1.0d0 / ( 0.94d21 * 1.086d0 ), pi = 3.141592653589793d0 )
       real*8,  allocatable :: smap(:,:), work(:,:), delta(:,:), corr(:,:), iref(:,:)
       integer, allocatable :: valid(:,:)
       character*3   cnan
       equivalence ( cnan, rnan )
       data          cnan / 'NaN' /
       external      planckf, rfits, wfits, interp, putima, invplanck, fieldfac, addkeys
!
! The field factor g, which adapts the correction factor to a radiation field other than G_0 = 10 in Habing units, measured on the
! grid of radiation fields (harvest_calib3_G0=*.txt.hist.npz): the pixel mean of ln f in bins of 0.05 dex of the fitted surface
! density, minus the same at G_0 = 10, smoothed by a local linear fit with Gaussian weights 0.1 dex wide times the square root of
! the number of pixels (make_field_factor.py). Above 10^22.3 cm^-2 the grid has too few pixels, and the field factor is held at its
! value at 10^22.25 cm^-2. Rows are log10 of the fitted surface density from gx0 in steps of gdx, columns are the fields gfld. It is
! interpolated linearly in log G_0 and held constant outside the table.
!
       data          gfld  / 0.98d0, 2.12d0, 4.61d0, 10.0d0, 21.7d0, 47.1d0, 102.4d0, 222.2d0, 482.5d0 /
!                   G_0 = 0.98
       data          gtab / 0.996d0, 0.987d0, 0.979d0, 0.971d0, 0.960d0, 0.945d0, 0.927d0, 0.907d0
     &             , 0.885d0, 0.856d0, 0.828d0, 0.799d0, 0.774d0, 0.741d0, 0.741d0, 0.741d0
     &             , 0.741d0, 0.741d0, 0.741d0, 0.741d0, 0.741d0, 0.741d0, 0.741d0, 0.741d0
!                   G_0 = 2.12
     &             , 0.994d0, 0.991d0, 0.988d0, 0.981d0, 0.971d0, 0.959d0, 0.944d0, 0.928d0
     &             , 0.910d0, 0.886d0, 0.861d0, 0.834d0, 0.811d0, 0.780d0, 0.780d0, 0.780d0
     &             , 0.780d0, 0.780d0, 0.780d0, 0.780d0, 0.780d0, 0.780d0, 0.780d0, 0.780d0
!                   G_0 = 4.61
     &             , 0.998d0, 0.995d0, 0.991d0, 0.988d0, 0.983d0, 0.976d0, 0.966d0, 0.959d0
     &             , 0.948d0, 0.931d0, 0.916d0, 0.904d0, 0.896d0, 0.869d0, 0.869d0, 0.869d0
     &             , 0.869d0, 0.869d0, 0.869d0, 0.869d0, 0.869d0, 0.869d0, 0.869d0, 0.869d0
!                   G_0 = 10
     &             , 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0
     &             , 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0
     &             , 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0
!                   G_0 = 21.7
     &             , 1.000d0, 1.000d0, 1.001d0, 1.003d0, 1.009d0, 1.016d0, 1.028d0, 1.046d0
     &             , 1.070d0, 1.094d0, 1.111d0, 1.123d0, 1.137d0, 1.144d0, 1.144d0, 1.144d0
     &             , 1.144d0, 1.144d0, 1.144d0, 1.144d0, 1.144d0, 1.144d0, 1.144d0, 1.144d0
!                   G_0 = 47.1
     &             , 0.999d0, 0.999d0, 1.000d0, 1.003d0, 1.011d0, 1.023d0, 1.042d0, 1.072d0
     &             , 1.117d0, 1.165d0, 1.204d0, 1.227d0, 1.240d0, 1.219d0, 1.219d0, 1.219d0
     &             , 1.219d0, 1.219d0, 1.219d0, 1.219d0, 1.219d0, 1.219d0, 1.219d0, 1.219d0
!                   G_0 = 102.4
     &             , 0.997d0, 0.995d0, 0.995d0, 0.997d0, 1.002d0, 1.013d0, 1.034d0, 1.068d0
     &             , 1.117d0, 1.168d0, 1.205d0, 1.212d0, 1.204d0, 1.172d0, 1.172d0, 1.172d0
     &             , 1.172d0, 1.172d0, 1.172d0, 1.172d0, 1.172d0, 1.172d0, 1.172d0, 1.172d0
!                   G_0 = 222.2
     &             , 0.987d0, 0.988d0, 0.988d0, 0.988d0, 0.991d0, 0.997d0, 1.010d0, 1.036d0
     &             , 1.072d0, 1.103d0, 1.117d0, 1.103d0, 1.076d0, 1.044d0, 1.044d0, 1.044d0
     &             , 1.044d0, 1.044d0, 1.044d0, 1.044d0, 1.044d0, 1.044d0, 1.044d0, 1.044d0
!                   G_0 = 482.5
     &             , 0.995d0, 0.989d0, 0.985d0, 0.983d0, 0.982d0, 0.982d0, 0.986d0, 0.999d0
     &             , 1.013d0, 1.018d0, 1.013d0, 0.987d0, 0.952d0, 0.887d0, 0.887d0, 0.887d0
     &             , 0.887d0, 0.887d0, 0.887d0, 0.887d0, 0.887d0, 0.887d0, 0.887d0, 0.887d0 /

!!! The field factor g, which adapts the correction factor to a radiation field other than G_0 = 10 in Habing units, measured on the
!!! grid of radiation fields (harvest_calib3_G0=*.txt): the pixel-weighted mean of log f in bins of 0.1 dex of the fitted surface
!!! density, smoothed over three bins, minus the same at G_0 = 10 (make_field_factor.py). Above 10^22.6 cm^-2 the grid has too few
!!! pixels, and the field factor is held at its value there. Rows are log10 of the fitted surface density from gx0 in steps of gdx,
!!! columns are the fields gfld. It is interpolated linearly in log G_0 and held constant outside the table.
!!!
!!       data          gfld  / 0.98d0, 2.12d0, 4.61d0, 10.0d0, 21.7d0, 47.1d0, 102.4d0, 222.2d0, 482.5d0 /
!!!                   G_0 = 0.98
!!       data          gtab / 0.991d0, 0.986d0, 0.978d0, 0.970d0, 0.960d0, 0.944d0, 0.920d0, 0.904d0
!!     &             , 0.904d0, 0.837d0, 0.815d0, 0.775d0, 0.756d0, 0.753d0, 0.756d0, 0.796d0
!!     &             , 0.767d0, 0.761d0, 0.761d0, 0.761d0, 0.761d0, 0.761d0, 0.761d0, 0.761d0
!!!                   G_0 = 2.12
!!     &             , 0.994d0, 0.991d0, 0.988d0, 0.981d0, 0.972d0, 0.958d0, 0.936d0, 0.924d0
!!     &             , 0.929d0, 0.891d0, 0.857d0, 0.803d0, 0.790d0, 0.786d0, 0.793d0, 0.843d0
!!     &             , 0.812d0, 0.822d0, 0.822d0, 0.822d0, 0.822d0, 0.822d0, 0.822d0, 0.822d0
!!!                   G_0 = 4.61
!!     &             , 0.999d0, 0.995d0, 0.992d0, 0.987d0, 0.981d0, 0.973d0, 0.959d0, 0.955d0
!!     &             , 0.967d0, 0.936d0, 0.889d0, 0.919d0, 0.909d0, 0.846d0, 0.851d0, 0.899d0
!!     &             , 0.916d0, 0.907d0, 0.907d0, 0.907d0, 0.907d0, 0.907d0, 0.907d0, 0.907d0
!!!                   G_0 = 10
!!     &             , 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0
!!     &             , 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0
!!     &             , 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0
!!!                   G_0 = 21.7
!!     &             , 1.001d0, 1.001d0, 1.001d0, 1.003d0, 1.007d0, 1.013d0, 1.026d0, 1.047d0
!!     &             , 1.124d0, 1.110d0, 1.084d0, 1.102d0, 1.107d0, 1.092d0, 1.099d0, 1.186d0
!!     &             , 1.168d0, 1.135d0, 1.135d0, 1.135d0, 1.135d0, 1.135d0, 1.135d0, 1.135d0
!!!                   G_0 = 47.1
!!     &             , 1.000d0, 1.000d0, 1.000d0, 1.001d0, 1.008d0, 1.017d0, 1.041d0, 1.072d0
!!     &             , 1.176d0, 1.182d0, 1.242d0, 1.179d0, 1.169d0, 1.232d0, 1.170d0, 1.220d0
!!     &             , 1.199d0, 1.165d0, 1.165d0, 1.165d0, 1.165d0, 1.165d0, 1.165d0, 1.165d0
!!!                   G_0 = 102.4
!!     &             , 0.997d0, 0.997d0, 0.997d0, 0.996d0, 1.001d0, 1.008d0, 1.034d0, 1.068d0
!!     &             , 1.178d0, 1.188d0, 1.248d0, 1.209d0, 1.162d0, 1.126d0, 1.120d0, 1.120d0
!!     &             , 1.114d0, 1.082d0, 1.082d0, 1.082d0, 1.082d0, 1.082d0, 1.082d0, 1.082d0
!!!                   G_0 = 222.2
!!     &             , 0.990d0, 0.989d0, 0.989d0, 0.986d0, 0.989d0, 0.991d0, 1.009d0, 1.035d0
!!     &             , 1.132d0, 1.126d0, 1.106d0, 1.092d0, 1.074d0, 1.032d0, 0.996d0, 0.996d0
!!     &             , 0.959d0, 0.933d0, 0.933d0, 0.933d0, 0.933d0, 0.933d0, 0.933d0, 0.933d0
!!!                   G_0 = 482.5
!!     &             , 0.994d0, 0.990d0, 0.986d0, 0.982d0, 0.982d0, 0.979d0, 0.985d0, 0.999d0
!!     &             , 1.010d0, 1.027d0, 1.029d0, 0.991d0, 0.936d0, 0.871d0, 0.838d0, 0.865d0
!!     &             , 0.834d0, 0.801d0, 0.801d0, 0.801d0, 0.801d0, 0.801d0, 0.801d0, 0.801d0
!!     &             /
!__________________________________________________________________________________________________________________________________
!
! The command line. The first positional argument is the image to be corrected; any further ones are images of the same
! reconstruction to which the same amount is added. The wavebands that were fitted to make the image are given with -bands=.
!
       call osystem ( lunix )
       call ftvers ( fitsvers )
       write (cfitsversion,'(f6.3)') fitsvers

       cverbose = '-verb1'
       mapfile  = ' '
       tempfile = ' '
       nfit     = 4
       wfit(1)  = 160.0d0
       wfit(2)  = 250.0d0
       wfit(3)  = 350.0d0
       wfit(4)  = 500.0d0
       g0       = 10.0d0

       narg = iargc ( )
       if (narg .lt. 1) then
         write (iotty,'(/a/)') '   HICORRECT: usage: hicorrect2 <image> [-temp=<file>] [-bands=<list>] [-g0=<field>]'
     &                         //' [-verb{0|1|2}]'
         stop
       endif

       do i=1,narg
         call getarg ( i, arg )
         la = len_trim ( arg )
         if (arg(1:5) .eq. '-verb') then
           cverbose = arg(1:6)
         elseif (arg(1:4) .eq. '-g0=') then
           read (arg(5:la),*) g0
         elseif (arg(1:6) .eq. '-temp=') then
           tempfile = arg(7:la)
         elseif (arg(1:7) .eq. '-bands=') then
           nfit = 1
           do k=8,la
             if (arg(k:k) .eq. ',') nfit = nfit + 1
           enddo
           if (nfit .gt. mband) then
             write (iotty,'(/a,i3,a/)') '   HICORRECT: ERROR: at most', mband, ' wavebands can be given with -bands=.'
             stop
           endif
           read (arg(8:la),*,iostat=ios) ( wfit(k), k=1,nfit )
           if (ios .ne. 0) then
             write (iotty,'(/a,a/)') '   HICORRECT: ERROR: cannot read the wavelengths from ', arg(1:la)
             stop
           endif
         elseif (mapfile .eq. ' ') then
           mapfile = arg
         else
           write (iotty,'(/a/a/)') '   HICORRECT: ERROR: only one surface density image can be corrected in one run.'
     &                         , '   Give the temperature image with -temp=<file>.'
           stop
         endif
       enddo

       if (mapfile .eq. ' ') then
         write (iotty,'(/a/)') '   HICORRECT: ERROR: no surface density image to correct was given.'
         stop
       endif

! The file names may be given without the extension '.fits', which is then added.
       la = len_trim ( mapfile )
       if (la .lt. 5 .or. mapfile(max(la-4,1):la) .ne. '.fits') mapfile = mapfile(1:la)//'.fits'
       if (tempfile .ne. ' ') then
         la = len_trim ( tempfile )
         if (la .lt. 5 .or. tempfile(max(la-4,1):la) .ne. '.fits') tempfile = tempfile(1:la)//'.fits'
       endif

       call when ( lunix, ctime, cdate, ndate, 4 )

       if (cverbose .ne. '-verb0') then
         write (iotty,'(/a)')        '   HICORRECT: correction of the single-temperature bias'
         write (iotty,'(a,i4,a)')    '   effective attenuating column from an angular average over', nang
     &                             , ' directions, sphere and layer averaged'
         write (iotty,'(a,f8.4,a,f8.4,a,f8.4,a,1pe11.4,a,0pf8.5,a)')
     &                               '   shielding law T(s) =', tfloor, ' + (', tsurf, ' -', tfloor, ') (1 + s /', sstar
     &                             , ')**(-', gamma, ') K'
         write (iotty,'(a,f9.3,a)')  '   radiation field G_0 =', g0, ' (Habing units)'
       endif
       if (g0 .lt. gfld(1) .or. g0 .gt. gfld(ngf)) then
         write (iotty,'(a,f9.3,a)')  '   HICORRECT2: WARNING: G_0 =', g0, ' lies outside the grid of radiation fields;'
         write (iotty,'(a)')         '   the field factor of the nearest field is used.'
       endif
!__________________________________________________________________________________________________________________________________
!
! The wavebands. The relation has to be built for the wavebands that were fitted to make the image, given with -bands= (160, 250,
! 350 and 500 um by default), because the bias depends on which colors were used to find the single temperature. They are sorted,
! and the longest of them is the reference waveband, the one the fit divides by.
!
       if (nfit .lt. 2) then
         write (iotty,'(/a/)') '   HICORRECT: ERROR: at least two wavebands must be given with -bands=.'
         stop
       endif
       do k=2,nfit
         do ib=nfit,k,-1
           if (wfit(ib-1) .gt. wfit(ib)) then
             bb          = wfit(ib)
             wfit(ib)    = wfit(ib-1)
             wfit(ib-1)  = bb
           endif
         enddo
       enddo
       waveref = wfit(nfit)

       do k=1,nfit
         cnu(k) = 1.0d17 * kappa0 * ( speedolight / ( wfit(k) * 1.0d-4 ) / nu0 ) ** betap * dtg * muh2 * amu
         if (wfit(k) .le. 160.5d0) then
           relerr(k) = 0.20d0
         else
           relerr(k) = 0.10d0
         endif
       enddo

       if (cverbose .ne. '-verb0') then
         write (iotty,'(a,10(1x,f5.0))') '   wavebands fitted (um):', ( wfit(k), k=1,nfit )
         write (iotty,'(a,f6.0,a)')      '   reference waveband', waveref, ' um'
       endif
!__________________________________________________________________________________________________________________________________
!
! The relation. For each of ntab true surface densities the emission is integrated over ndep layers, one modified blackbody is
! fitted to it on a grid of ntfit temperatures exactly as the reconstruction fits it, and the surface density that the fit returns
! is recorded. The table therefore runs over the fitted surface density, which is what the image holds, and gives the factor that
! restores the true one. Both arrays increase with the index, so the interpolation below is straightforward.
!
       do k=1,ntfit
         tgrid(k) = 0.98d0 * tfloor * ( 1.02d0 * tsurf / ( 0.98d0 * tfloor ) ) ** ( dble ( k - 1 ) / dble ( ntfit - 1 ) )
       enddo

! The factor that diffuse material requires, measured on a radmc-3d model of a diffuse cloud whose projection is the cloud map of
! the simulated sky: medians of the true over the fitted surface density, in bins of the fitted surface density (log10, cm^-2).
! Below 10^xb0 the factor is the average of the analytic sphere and this cloud, which replaces the analytic layer; the difference
! from the analytic average of sphere and layer is reduced smoothly to zero at 10^xb1.
!
       cirx = (/ 20.55d0, 20.65d0, 20.75d0, 20.85d0, 20.95d0, 21.05d0, 21.15d0, 21.25d0, 21.35d0, 21.45d0, 21.55d0, 21.65d0 /)
       cirf = (/ 1.049d0, 1.051d0, 1.053d0, 1.058d0, 1.066d0, 1.075d0, 1.101d0, 1.121d0, 1.146d0, 1.180d0, 1.225d0, 1.251d0 /)
       xb0 = 21.6d0
       xb1 = 22.0d0

       do igeo=1,2                                      !<~ 1: sphere and layer averaged; 2: sphere alone
       do i=1,ntab
         sig = smin * ( smax / smin ) ** ( dble ( i - 1 ) / dble ( ntab - 1 ) )

! The effective attenuating column at each depth, averaged between a uniform sphere sampled along a diameter and a uniform
! plane-parallel layer sampled along its normal. In both the average is over nang directions and t0 is the optical depth of the
! half thickness at the wavelengths that heat the dust.
!
         t0 = 0.5d0 * kapheat * sig
         do j=1,ndep
           zz   = ( dble ( j ) - 0.5d0 ) / dble ( ndep ) * 2.0d0 - 1.0d0
           jsph = 0.0d0
           jlay = 0.0d0
           wsum = 0.0d0
           do l=1,nang
             psi  = ( dble ( l ) - 0.5d0 ) / dble ( nang ) * pi
             wps  = sin ( psi )
             lsph = -abs ( zz ) * cos ( psi ) + sqrt ( max ( 1.0d0 - ( abs ( zz ) * sin ( psi ) ) ** 2, 0.0d0 ) )
             jsph = jsph + exp ( -t0 * lsph ) * wps
             wsum = wsum + wps
             xmu  = ( dble ( l ) - 0.5d0 ) / dble ( nang ) * 2.0d0 - 1.0d0
             if (xmu .gt. 0.0d0) then
               dlay = ( 1.0d0 - zz ) / abs ( xmu )
             else
               dlay = ( 1.0d0 + zz ) / abs ( xmu )
             endif
             jlay = jlay + exp ( -t0 * dlay )
           enddo
           jsph = jsph / wsum
           jlay = jlay / dble ( nang )
           ssph = -log ( max ( jsph, 1.0d-300 ) ) / max ( t0, 1.0d-12 )
           slay = -log ( max ( jlay, 1.0d-300 ) ) / max ( t0, 1.0d-12 )
           if (igeo .eq. 1) then
             seff(j) = 0.5d0 * ( ssph + slay ) * sig / 2.0d0
           else
             seff(j) = ssph * sig / 2.0d0
           endif
         enddo

         do k=1,nfit
           obs(k) = 0.0d0
           do j=1,ndep
             tt = tfloor + ( tsurf - tfloor ) * ( 1.0d0 + seff(j) / sstar ) ** (-gamma)
             obs(k) = obs(k) + planckf ( tt, wfit(k), hplanck, kboltz, speedolight )
           enddo
           obs(k) = cnu(k) * sig * obs(k) / dble ( ndep )
           wgt(k) = 1.0d0 / ( relerr(k) * obs(k) ) ** 2
         enddo

         chibest = 1.0d30
         tbest   = tgrid(1)
         do k=1,ntfit
           amp = 0.0d0
           ss  = 0.0d0
           do l=1,nfit
             model(l) = cnu(l) * planckf ( tgrid(k), wfit(l), hplanck, kboltz, speedolight )
             amp = amp + wgt(l) * obs(l) * model(l)
             ss  = ss  + wgt(l) * model(l) * model(l)
           enddo
           amp  = amp / ss
           chi2 = 0.0d0
           do l=1,nfit
             chi2 = chi2 + wgt(l) * ( obs(l) - amp * model(l) ) ** 2
           enddo
           if (chi2 .lt. chibest) then
             chibest = chi2
             tbest   = tgrid(k)
           endif
         enddo

         if (igeo .eq. 1) then
           ffit(i) = obs(nfit) / ( cnu(nfit) * planckf ( tbest, waveref, hplanck, kboltz, speedolight ) )
           fcor(i) = sig / ffit(i)
           sgrid(i) = sig
         else
           ffits(i) = obs(nfit) / ( cnu(nfit) * planckf ( tbest, waveref, hplanck, kboltz, speedolight ) )
           fcors(i) = sig / ffits(i)
         endif
       enddo
       enddo

! The calibration of the diffuse end on the radmc-3d cloud, as described above.
!
       do i=1,ntab
         xg = log10 ( ffit(i) )
         call interp ( ntab, ffits, fcors, ffit(i), fsx )
         if (xg .le. cirx(1)) then
           fcx = cirf(1)
         else if (xg .ge. cirx(12)) then
           fcx = cirf(12)
         else
           do k=1,11
             if (xg .lt. cirx(k+1)) exit
           enddo
           fcx = cirf(k) + ( cirf(k+1) - cirf(k) ) * ( xg - cirx(k) ) / ( cirx(k+1) - cirx(k) )
         endif
         ddx(i) = 0.5d0 * ( fsx + fcx ) - fcor(i)
       enddo
       d0 = ddx(1)
       do i=1,ntab
         if (log10 ( ffit(i) ) .le. xb0) d0 = ddx(i)
       enddo
       do i=1,ntab
         xg = log10 ( ffit(i) )
         if (xg .le. xb0) then
           wb = 1.0d0
           dd = ddx(i)
         else if (xg .ge. xb1) then
           wb = 0.0d0
           dd = 0.0d0
         else
           wb = 0.5d0 * ( 1.0d0 + cos ( pi * ( xg - xb0 ) / ( xb1 - xb0 ) ) )
           dd = d0
         endif
         fcor(i) = fcor(i) + wb * dd
       enddo

! The field of the cloud, and the factor made non-decreasing with the fitted surface density. Above about 3 x 10^22 cm^-2 the
! relation of uniform bodies declines, because such a body is nearly isothermal at the floor temperature. The radiative transfer
! models of concentrated cores do not decline, so the factor is held at its maximum there.
!
       do i=1,ntab
         if (g0 .ne. 10.0d0) fcor(i) = fcor(i) * fieldfac ( log10 ( ffit(i) ), g0, gtab, gfld, ngx, ngf, gx0, gdx )
         if (i .gt. 1) fcor(i) = max ( fcor(i), fcor(i-1) )
       enddo

       if (cverbose .eq. '-verb2') then
         write (iotty,'(/a)') '   the tabulated relation, fitted surface density and the factor that restores the true one:'
         do i=1,ntab,20
           write (iotty,'(4x,1pe12.4,3x,0pf9.4)') ffit(i), fcor(i)
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
! The image to be corrected.
!
       call getfitshead ( mapfile(1:len_trim(mapfile)), nx, ny, dx, dy, bunit )

       allocate ( smap(nx,ny), valid(nx,ny), delta(nx,ny), corr(nx,ny), work(nx,ny), iref(nx,ny) )

       la = len_trim ( mapfile )
       call rfits ( nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, smap, dx, dy, object, ra, dec, mapfile(1:la)
     &            , iotty, iolog, creator, beamkw, funmin, funmax, blank, crota1, crota2, cd11, cd12, cd21, cd22, equinox, bzero
     &            , bscale, wavekw, datamn, datamx, history, cverbose )

       if (creator(1:10) .eq. 'HICORRECT') then
         write (iotty,'(/a/a/)') '   HICORRECT: ERROR: '''//mapfile(1:la)//''' has already been corrected.'
     &                         , '   The correction is not idempotent and must be applied once only.'
         stop
       endif

       beamsd = beamkw
       pix = dx
! The outputs are written to the current directory, whatever directory the input comes from: the stem of their names is the name
! of the input without its directory and without the extension '.fits'.
       ls = la
       if (mapfile(max(ls-4,1):ls) .eq. '.fits') ls = ls - 5
       k  = index ( mapfile(1:ls), '/', back=.true. )
       stem = mapfile(k+1:ls)
       ls = ls - k
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
       chead( 4) = 'HICORRECT'
       chead( 5) = history
!__________________________________________________________________________________________________________________________________
!
! The coverage. A pixel that is not positive carries no information: the reconstruction did not produce a surface density there.
!
       nvalid = 0
       do j=1,ny
         do i=1,nx
           if (smap(i,j) .gt. 0.0d0 .and. .not.isnan ( smap(i,j) )) then
             valid(i,j) = 1
             nvalid = nvalid + 1
           else
             valid(i,j) = 0
           endif
         enddo
       enddo

       if (nvalid .eq. 0) then
         write (iotty,'(/a/)') '   HICORRECT: ERROR: the image has no positive pixel.'
         stop
       endif
!__________________________________________________________________________________________________________________________________
!
! The correction of the image. The factor is found by
! interpolating the table in the logarithm of the surface density, and it is held at the ends of the table rather than
! extrapolated.
!
       do j=1,ny
         do i=1,nx
           if (valid(i,j) .eq. 1) then
             call interp ( ntab, ffit, fcor, smap(i,j), f )
             corr(i,j)  = f
             delta(i,j) = ( f - 1.0d0 ) * smap(i,j)
           else
             corr(i,j)  = 1.0d0
             delta(i,j) = 0.0d0
           endif
         enddo
       enddo

       if (cverbose .ne. '-verb0') then
         funmin =  1.0d30
         funmax = -1.0d30
         ss = 0.0d0
         tt = 0.0d0
         do j=1,ny
           do i=1,nx
             if (valid(i,j) .eq. 1) then
               funmin = min ( funmin, corr(i,j) )
               funmax = max ( funmax, corr(i,j) )
               ss = ss + smap(i,j)
               tt = tt + smap(i,j) + delta(i,j)
             endif
           enddo
         enddo
         write (iotty,'(a,i10,a)')      '   valid pixels of the image:     ', nvalid, ''
         write (iotty,'(a,2(1x,f8.4))') '   correction factor, smallest and largest:', funmin, funmax
         write (iotty,'(a,f8.4)')       '   ratio of the sums, corrected to uncorrected:', tt / ss
       endif

       fname = stem(1:ls)//'.cfactor.fits'

       call putima ( iotty, cverbose, cfitsversion, nx, ny, corr, valid, rnan, fname, 'ratio', rhead, chead, blank, cdate, ctime )

       call addkeys ( fname, g0, nfit, wfit, tfloor, tsurf, sstar, gamma )

       if (cverbose .ne. '-verb0') write (iotty,'(a)') '   --> '''//fname(1:len_trim(fname))//''''

       do j=1,ny
         do i=1,nx
           work(i,j) = smap(i,j) + delta(i,j)
         enddo
       enddo
       fname = stem(1:ls)//'.correct.fits'

       call putima ( iotty, cverbose, cfitsversion, nx, ny, work, valid, rnan, fname, bunit, rhead, chead, blank, cdate, ctime )

       call addkeys ( fname, g0, nfit, wfit, tfloor, tsurf, sstar, gamma )

       if (cverbose .ne. '-verb0') write (iotty,'(a)') '   --> '''//fname(1:len_trim(fname))//''''
!__________________________________________________________________________________________________________________________________
!
! The dust temperature that goes with the corrected surface density, at the resolution of the corrected image. The reconstruction
! derived the surface density from the reference waveband with the temperature it fitted at that resolution. With f times more
! material on the line of sight, the same intensity of the reference waveband is reproduced by colder dust, whose Planck function
! is that of the fitted temperature divided by f. The fitted temperature image of the same resolution is given with -temp=<file>.
!
       if (tempfile .ne. ' ') then

         la = len_trim ( tempfile )

         call getfitshead ( tempfile(1:la), nxi, nyi, dxi, dyi, bunit )

         if (nxi .eq. nx .and. nyi .eq. ny) then

           call rfits ( nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, iref, dxi, dyi, object, ra, dec
     &                , tempfile(1:la), iotty, iolog, creator, beamkw, funmin, funmax, blank, crota1, crota2, cd11, cd12
     &                , cd21, cd22, equinox, bzero, bscale, wavekw, datamn, datamx, history, cverbose )

           if (beamsd .gt. 0.0d0 .and. beamkw .gt. 0.0d0 .and. abs ( beamkw / beamsd - 1.0d0 ) .gt. 1.0d-2) then
             write (iotty,'(a,f7.2,a,f7.2,a)') '   HICORRECT2: ERROR: the temperature image has the resolution', beamkw
     &                                       , ' arcsec, the surface density image', beamsd, ' arcsec.'
             write (iotty,'(a)') '   They must have the same resolution; no temperature map is written.'
             deallocate ( smap, valid, delta, corr, work, iref )
             stop
           endif

           do j=1,ny
             do i=1,nx
               if (valid(i,j) .eq. 1 .and. iref(i,j) .gt. 0.0d0) then
                 bb = planckf ( iref(i,j), waveref, hplanck, kboltz, speedolight ) / corr(i,j)
                 work(i,j) = invplanck ( bb, waveref, hplanck, kboltz, speedolight )
               else
                 work(i,j) = 0.0d0
               endif
             enddo
           enddo
           lt = len_trim ( tempfile )
           if (tempfile(max(lt-4,1):lt) .eq. '.fits') lt = lt - 5
           k  = index ( tempfile(1:lt), '/', back=.true. )
           fname = tempfile(k+1:lt)//'.correct.fits'

           call putima ( iotty, cverbose, cfitsversion, nx, ny, work, valid, rnan, fname, 'K', rhead, chead, blank, cdate, ctime )

           call addkeys ( fname, g0, nfit, wfit, tfloor, tsurf, sstar, gamma )

           if (cverbose .ne. '-verb0') write (iotty,'(a)') '   --> '''//fname(1:len_trim(fname))//''''
         else
           write (iotty,'(a)') '   HICORRECT2: WARNING: the temperature image does not match the corrected image,'
           write (iotty,'(a)') '   so no temperature map is written.'
         endif
       else
         if (cverbose .ne. '-verb0') write (iotty,'(a)') '   no temperature image given with -temp=<file>,'
     &                                                 //' so no temperature map is written.'
       endif
       deallocate ( smap, valid, delta, corr, work, iref )

       if (cverbose .ne. '-verb0') write (iotty,'(/a/)') '   HICORRECT DONE.'
       stop
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

       subroutine addkeys
         
     &          ( fname, g0, nfit, wfit, tfloor, tsurf, sstar, gamma )
!__________________________________________________________________________________________________________________________________
!
! Add to the header of an output image the radiation field of the correction, the wavebands fitted to make the image, and the four
! constants of the temperature relation, so that every output records how it was corrected.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       character*(*) fname
       integer       nfit, k, unit, status, block, la
       real*8        g0, wfit(nfit), tfloor, tsurf, sstar, gamma
       character*80  bands
       character*8   cw
       external      ftgiou, ftopen, ftpkyd, ftpkys, ftclos, ftfiou, printerr
!__________________________________________________________________________________________________________________________________
!
       bands = ' '
       do k=1,nfit
         write (cw,'(i0)') nint ( wfit(k) )
         la = len_trim ( bands )
         if (la .eq. 0) then
           bands = cw
         else
           bands = bands(1:la)//','//cw
         endif
       enddo

       status = 0
       call ftgiou ( unit, status )
       call ftopen ( unit, fname, 1, block, status )
       call ftpkyd ( unit, 'HIG0'   , g0    , 6, 'radiation field of the correction, Habing units', status )
       call ftpkys ( unit, 'HIBANDS', bands ,    'wavebands fitted to make the image, um', status )
       call ftpkyd ( unit, 'HITMIN' , tfloor, 8, 'T_min of the temperature relation, K', status )
       call ftpkyd ( unit, 'HIT0'   , tsurf , 8, 'T_0 of the temperature relation, K', status )
       call ftpkyd ( unit, 'HISSTAR', sstar , 7, 'Sigma_star of the temperature relation, cm^-2', status )
       call ftpkyd ( unit, 'HIGAMMA', gamma , 7, 'gamma of the temperature relation', status )
       call ftclos ( unit, status )
       call ftfiou ( unit, status )
       if (status .gt. 0) call printerr ( status )

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

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function fieldfac

     &          ( x, g0, gtab, gfld, ngx, ngf, gx0, gdx )
!__________________________________________________________________________________________________________________________________
!
! Field factor g of the correction factor for the radiation field g0, at x = log10 of the fitted surface density. The
! logarithm of the table is interpolated linearly in log g0 and in x, and held constant outside the table.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       ngx, ngf, j, i
       real*8        fieldfac, x, g0, gtab(ngx,ngf), gfld(ngf), gx0, gdx, t, a, b, xi, r1, r2

       t = log10 ( min ( max ( g0, gfld(1) ), gfld(ngf) ) )
       j = 1
       do while (j .lt. ngf - 1 .and. t .gt. log10 ( gfld(j+1) ))
         j = j + 1
       enddo
       a = ( t - log10 ( gfld(j) ) ) / ( log10 ( gfld(j+1) ) - log10 ( gfld(j) ) )

       xi = ( x - gx0 ) / gdx + 1.0d0
       xi = min ( max ( xi, 1.0d0 ), dble ( ngx ) )
       i  = min ( int ( xi ), ngx - 1 )
       b  = xi - dble ( i )

       r1 = ( 1.0d0 - b ) * log ( gtab(i,j)   ) + b * log ( gtab(i+1,j)   )
       r2 = ( 1.0d0 - b ) * log ( gtab(i,j+1) ) + b * log ( gtab(i+1,j+1) )
       fieldfac = exp ( ( 1.0d0 - a ) * r1 + a * r2 )

       return
       end
