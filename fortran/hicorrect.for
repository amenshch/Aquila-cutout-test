
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       program hicorrect
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!
! HIRES • CORRECTION OF THE SINGLE-TEMPERATURE BIAS OF SURFACE DENSITY IMAGES
!
! A surface density image is built by fitting one modified blackbody to the wavebands of each line of sight and dividing the
! reference waveband by it. The correction applies to any such image, whether made by hires or by any other method, at any
! resolution; it needs only the image and the list of wavebands that were fitted to make it. A real line of sight is not
! isothermal: its dust is warm near the surface of the cloud and cold in the interior, the single fitted temperature is weighted
! toward the warm dust that dominates the emission, and the surface density is therefore underestimated wherever the temperature
! varies along the line of sight. The underestimation exceeds a factor of two in the densest material.
!
! Where the underestimation is. The reconstruction of hires forms a base image at the coarsest resolution and adds to it one
! differential increment per waveband, so that the image at any resolution is
!
!   S(theta) = S(theta_coarse) + sum of the increments finer than theta.
!
! The increments are differences between an image and a smoothed copy of itself, so they very nearly cancel when integrated over
! an extended region: on the model grid they carry less than three percent of the true surface density. The whole integrated
! underestimation therefore sits in the base image, and every image of the reconstruction is underestimated by nearly the same
! fraction.
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
!   T(s) = tfloor + (tsurf - tfloor) (1 + s / sstar)**(-gamma),
!
! measured on radmc-3d models of Bonnor-Ebert spheres embedded in uniform clouds and of isolated cylinders. The emission of the 
! line of sight is integrated over depth in each of the wavebands that the reconstruction fitted, one modified blackbody is fitted
! to the result exactly as the reconstruction fits it, and the surface density that the fit returns is compared with S. The ratio
! is the factor by which the fitted surface density has to be multiplied. It is tabulated once against the fitted surface density,
! which is the quantity the image holds, and interpolated onto the image. At low surface densities the analytic layer is replaced
! by a radmc-3d model of a diffuse cloud, whose surface layers are warmer than the relation predicts (see below). This calibration
! is part of the correction of every image, not an adjustment for unusual wavebands.
!
! The steps, with the equations of the paper (Sect. 4): the attenuating column at every depth of a uniform sphere and a uniform
! layer, averaged (Eq. 4); the temperature at every depth (Eq. 5); the intensities of the line of sight (Eq. 6); one modified
! blackbody fitted to them and the fitted surface density from the reference waveband (Eq. 7); the factor, true over fitted (Eq.
! 8), tabulated against the fitted surface density; the calibration on the diffuse cloud below 10^22 cm^-2 (Sect. 4.3); the factor
! held at its maximum above the surface density where it is reached and multiplied by the field factor g for G_0 other than 10
! (Eq. 9); the image multiplied by the factor interpolated onto it (Eq. 10).
!
! The geometry. The sphere and the layer are the two limits of a uniform structure. Their factors differ by up to 10%, and by less
! than 5% outside the range from 10^21.4 to 10^22.1 cm^-2, and their average is used. The arrangement of the material along the
! line of sight is not known. Where a structure lies in a weak cloud, it is corrected too little; where it lies near a face of a
! dense cloud, it is corrected too much. On radiative transfer models, these errors range from -13% to +21%.
!
! The correction must be applied once. It is not idempotent: applied to an image it has already corrected it would raise the values
! again. Every output image is stamped with CREATOR = 'HICORRECT', and the program refuses any input that carries it. The output
! images also record the radiation field of the correction (HIG0), how it was obtained (HIG0SRC: given, estimated or default),
! the fields at the 16th and 84th percentiles if it was estimated (HIG0P16, HIG0P84), the wavebands fitted (HIBANDS) with their
! relative calibration uncertainties (HIERRS), and the four constants of the temperature relation (HITMIN, HIT0, HISSTAR, HIGAMMA).
!
! Usage:  hicorrect <dimage> [<timage>] [-bands=<waves>] [-errors=<relerrs>] [-g0=<field>] [-verb{0|1|2}]
!
!         <dimage>     the surface density image to be corrected; the factor is evaluated on it and it is multiplied by it
!                      (the file names of <dimage> and <timage> may be given without the extension '.fits', which is then added)
!         <timage>     the dust temperature image fitted together with <dimage>, optional; it must have the same resolution and
!                      the same pixels as <dimage>. The two images are told apart by their values: an image with more pixels
!                      between 1 and 10^4 than above 10^4 is taken as temperatures, and the run stops if <dimage> looks like
!                      temperatures or <timage> does not.
!         <waves>      the wavelengths in µm of the wavebands fitted to make <dimage>, separated by commas, for example
!                      -bands=250,350,500; by default: 160,250,350,500. The longest one is the reference waveband.
!         <relerrs>    the relative calibration uncertainties of the same wavebands, in the same order, which set their weights
!                      in the fit that made <dimage>, for example -errors=0.1,0.1,0.1 for equal weights; only their ratios
!                      matter. By default 0.2 at 160 µm and shorter and 0.1 at longer wavelengths, the weights of hires. The
!                      analytic relation is built for any wavebands and weights. The calibration on a diffuse cloud is provided
!                      for five sets: 160-500 µm with the default or equal weights, 250-500 µm with equal weights, 100-500 and
!                      70-500 µm with the default weights; for other sets that of the standard set is used, with a warning. The
!                      field factor was measured with the standard set only; for other sets it is approximate.
!         <field>      the radiation field of the cloud in Habing units. Without -g0=, the field is estimated from the images if
!                      <timage> is given, and G_0 = 10 is used otherwise. The estimate compares the fitted temperatures of the
!                      diffuse material, between 10^21.0 and 10^21.8 cm^-2, with the median fitted temperatures of the models of
!                      the grid of radiation fields at the same fitted surface density (the table ttab). For each bin of 0.2 dex
!                      with at least 100 pixels, the 16th, 50th and 84th percentiles of the temperatures of the image are
!                      converted into fields by interpolating log G_0 against the model temperatures. The field is the median
!                      of the fields of the bins weighted by their pixels, rounded to two significant digits; the 16th and 84th
!                      percentiles give its spread. The model temperatures are those of 13.5 arcsec images fitted to 160-500 µm
!                      with the standard weights; for other resolutions or wavebands the estimate is used with a warning. A
!                      warning is also printed if the healthy pixels (surface densities of at least 10^19 cm^-2 and positive
!                      temperatures) cover less than the smallest observed Herschel maps (2000 x 2000 pixels of 3 arcsec): the
!                      image is then probably a cutout, and the field should be estimated from the full map.
!
! output: written to the directory containing the input images:
!         <dimage>.correct.fits, the corrected image
!         <dimage>.cfactor.fits, the image of the correction factors
!         <timage>.correct.fits, the dust temperature implied by the corrected surface density, if <timage> is given
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       logical       lunix, estok, istemp

       character*1   cnfit
       character*3   cnan, dot
       character*6   cfitsversion, cverbose
       character*7   clibname
       character*8   ctime, cfpptime
       character*10  cdate
       character*21  compda
       character*23  g0src
       character*80  bunit, ctype1, ctype2, object, creator, history,  sbunit, chead(5)
       character*256 mapfile, stem, fname, tempfile, arg1, arg2, arg3, arg4, arg5, arg6, arg

       integer       sblank, iotty, iolog, mband, ntab, ndep, ntfit, nang, ncset, icset, nwgt, ic, lm, ngx, ngf, lt, ntb
     &             , tminpix, ndif, ngood, nx, ny, nxi, nyi, i, j, k, l, ib, la, ls, nfit, nvalid, blank, ndate, igeo, ios, lb
     &             , ia, lastc, lfn, ia1, ia2, ia3, ia4, ia5, ia6

       parameter   ( iotty = 6, iolog = 0, mband = 6, ntab = 200, ndep = 512, ntfit = 3000, nang = 160, ncset = 5 )
       parameter   ( ngx = 24, ngf = 9, ntb = 11, tminpix = 100 )
       integer       ncirn(ncset)

       real          fitsvers
       real*8        ttab(ntb,ngf), tx0, tdx, txfield, g16, g50, g84, hotfr, ee, gtab(ngx,ngf), gfld(ngf), gx0, gdx, g0, fieldfac
       parameter   ( tx0 = 21.0d0, tdx = 0.2d0, txfield = 21.8d0, gx0 = 20.95d0, gdx = 0.1d0 )

       real*8        wfit(mband), waveref, bb, invplanck
     &             , dx, dy, dxi, dyi, crpix1, crpix2, crval1, crval2, crota1, crota2, cd11, cd12, cd21, cd22, equinox, bzero
     &             , bscale, beamkw, beamsd, funmin, funmax, datamn, datamx, wavekw, pix, rnan, rhead(17), ra, dec
     &             , tfloor, tsurf, sstar, gamma, kappa0, nu0, betap, dtg, muh2, amu, hplanck, kboltz, speedolight, kapheat, pi
     &             , smin, smax, sgrid(ntab), ffit(ntab), fcor(ntab), tgrid(ntfit), seff(ndep), relerr(mband)
     &             , ffits(ntab), fcors(ntab), xg, fsx, fcx, dd, d0, xb0, xb1, wb, ddx(ntab), cirx(12), cirf(12)
     &             , relin(mband), cirtab(12,ncset), cirw(mband,ncset), cire(mband,ncset), emin
     &             , cnu(mband), obs(mband), model(mband), wgt(mband), amp, chi2, chibest, tbest, sig, tt, ss, f, planckf
     &             , t0, zz, psi, wps, wsum, jsph, jlay, lsph, dlay, xmu, ssph, slay

       parameter   ( tfloor = 7.049531d0, tsurf = 22.866783d0, sstar = 2.389668d+21, gamma = 1.596874d0, kappa0 = 10.0d0
     &             , nu0 = 1.0d12, betap = 2.0d0, dtg = 1.0d-2, muh2 = 2.8d0, amu = 1.66053906660d-24, hplanck = 6.62607015d-27
     &             , kboltz = 1.380649d-16, speedolight = 2.99792458d10, smin = 1.0d20, smax  = 1.0d24
     &             , kapheat = 1.0d0 / (0.94d21 * 1.086d0), pi = 3.141592653589793d0, dot = '•' )

! The constants. tfloor, tsurf, sstar and gamma are T_min, T_max, Sigma_star and gamma of the relation T(Sigma_A) of Eq. (5) of
! the paper, fitted by fit_attenuation_law.py to 2086 Bonnor-Ebert spheres of the main grid (the outermost shell of every sphere
! rejected) and 7 isolated cylinders computed on the full theta grid, each shell weighted by its fraction of the mass of its model,
! the spheres together and the cylinders together with half of the weight each; the scatter about the relation is 2.4%.
! kappa0, nu0 and betap give the opacity of the fit, 10 cm^2 per gram of dust at 1000 GHz (299.8 µm) with the spectral index 2,
! that is, 10 (300 µm / lambda)^2 within 0.14%; the radmc-3d models used the same law at 60 µm and longer (within 0.5%), so the
! fit and the models differ only in the temperatures along the line of sight. dtg is the dust-to-gas mass ratio and muh2 the mean
! molecular weight per H2. kapheat is the gray extinction cross section per H2 at the ultraviolet and optical wavelengths that heat
! the dust (sigma of Sect. 3.1 of the paper): A_V = Sigma / 0.94e21 cm^-2 and tau = A_V / 1.086, so 9.8e-22 cm^2. smin and smax
! bound the true surface densities of the table (ntab values), ndep is the number of depths along the line of sight, ntfit the
! number of temperatures searched by the fit, nang the number of directions of the angular average.

       real*8,  allocatable :: smap(:,:), work(:,:), delta(:,:), corr(:,:), iref(:,:)
       integer, allocatable :: valid(:,:)

       equivalence ( cnan, rnan )
       data          cnan / 'NaN' /

       external      planckf, rfits, wfits, interp, putima, invplanck, fieldfac, addkeys, estfield, istemp, lastc
!__________________________________________________________________________________________________________________________________
!
! The field factor g, which adapts the correction factor to a radiation field other than G_0 = 10 in Habing units, measured on the
! grid of radiation fields (harvest_calib3_G0=*.txt.hist.npz). The pixels are binned in the true surface density, 0.05 dex wide,
! where the populations of all fields coincide; each bin gives the mean log10 of the fitted surface density and the mean ln f, and
! ln f is smoothed against the former by a local linear fit with Gaussian weights 0.1 dex wide times the square root of the number
! of pixels; the same at G_0 = 10 is subtracted (make_field_factor.py). Above 10^22.85 cm^-2 the grid has too few pixels, and the
! field factor is held at its value there. Rows are log10 of the fitted surface density, from gx0 in steps of gdx; columns are the
! fields gfld. It is interpolated linearly in log G_0 and held constant outside the table.

       data gfld  / 0.98d0, 2.12d0, 4.61d0, 10.0d0, 21.7d0, 47.1d0, 102.4d0, 222.2d0, 482.5d0 /

       data gtab  /  0.994d0, 0.987d0, 0.979d0, 0.971d0, 0.958d0, 0.942d0, 0.924d0, 0.905d0 ! G_0 = 0.98
     &             , 0.878d0, 0.845d0, 0.815d0, 0.793d0, 0.785d0, 0.777d0, 0.768d0, 0.763d0
     &             , 0.762d0, 0.764d0, 0.770d0, 0.769d0, 0.769d0, 0.769d0, 0.769d0, 0.769d0
     
     &             , 0.995d0, 0.992d0, 0.988d0, 0.981d0, 0.970d0, 0.957d0, 0.943d0, 0.928d0 ! G_0 = 2.12
     &             , 0.906d0, 0.877d0, 0.850d0, 0.833d0, 0.824d0, 0.817d0, 0.812d0, 0.810d0
     &             , 0.811d0, 0.810d0, 0.808d0, 0.807d0, 0.807d0, 0.807d0, 0.807d0, 0.807d0
     
     &             , 0.998d0, 0.995d0, 0.991d0, 0.987d0, 0.981d0, 0.970d0, 0.962d0, 0.955d0 ! G_0 = 4.61
     &             , 0.943d0, 0.924d0, 0.907d0, 0.899d0, 0.896d0, 0.891d0, 0.888d0, 0.882d0
     &             , 0.880d0, 0.882d0, 0.889d0, 0.891d0, 0.891d0, 0.891d0, 0.891d0, 0.891d0
     
     &             , 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0 ! G_0 = 10
     &             , 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0
     &             , 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0
     
     &             , 1.000d0, 1.000d0, 1.001d0, 1.004d0, 1.009d0, 1.014d0, 1.028d0, 1.051d0 ! G_0 = 21.7
     &             , 1.082d0, 1.103d0, 1.121d0, 1.134d0, 1.140d0, 1.136d0, 1.126d0, 1.101d0
     &             , 1.080d0, 1.064d0, 1.075d0, 1.085d0, 1.085d0, 1.085d0, 1.085d0, 1.085d0
     
     &             , 0.999d0, 0.999d0, 1.000d0, 1.004d0, 1.012d0, 1.022d0, 1.044d0, 1.080d0 ! G_0 = 47.1
     &             , 1.130d0, 1.175d0, 1.215d0, 1.242d0, 1.240d0, 1.223d0, 1.205d0, 1.181d0
     &             , 1.161d0, 1.132d0, 1.111d0, 1.098d0, 1.098d0, 1.098d0, 1.098d0, 1.098d0
     
     &             , 0.996d0, 0.995d0, 0.995d0, 0.996d0, 1.001d0, 1.012d0, 1.035d0, 1.074d0 ! G_0 = 102.4
     &             , 1.130d0, 1.177d0, 1.209d0, 1.223d0, 1.206d0, 1.177d0, 1.147d0, 1.114d0
     &             , 1.083d0, 1.050d0, 1.024d0, 1.005d0, 1.005d0, 1.005d0, 1.005d0, 1.005d0
     
     &             , 0.989d0, 0.989d0, 0.988d0, 0.989d0, 0.991d0, 0.995d0, 1.010d0, 1.040d0 ! G_0 = 222.2
     &             , 1.082d0, 1.113d0, 1.120d0, 1.110d0, 1.089d0, 1.052d0, 1.018d0, 0.977d0
     &             , 0.939d0, 0.908d0, 0.884d0, 0.863d0, 0.863d0, 0.863d0, 0.863d0, 0.863d0
     
     &             , 0.995d0, 0.990d0, 0.986d0, 0.984d0, 0.981d0, 0.978d0, 0.984d0, 0.999d0 ! G_0 = 482.5
     &             , 1.014d0, 1.015d0, 1.004d0, 0.973d0, 0.943d0, 0.907d0, 0.869d0, 0.835d0
     &             , 0.798d0, 0.768d0, 0.746d0, 0.723d0, 0.723d0, 0.723d0, 0.723d0, 0.723d0 /

! The median fitted temperature (K) of the models of the grid of radiation fields, for estimating the field of a cloud when -g0= is
! not given. Rows are bins of 0.2 dex of the fitted surface density from 10^tx0 cm^-2; columns are the fields gfld; -1 where a bin
! has no models. From the images at 13.5 arcsec fitted to 160-500 µm with the standard weights (make_g0_table.py on the tables of
! harvest_calib3.py).

       data ttab /   14.900d0, 14.900d0, 13.900d0, 12.500d0, 12.300d0, 11.300d0, 10.500d0, 9.700d0  ! G_0 = 0.98
     &             , 8.897d0, 8.100d0, 7.900d0
     
     &             , 15.500d0, 15.300d0, 14.100d0, 12.900d0, 11.700d0, 11.500d0, 10.700d0, 9.900d0  ! G_0 = 2.12
     &             , 9.100d0, 8.300d0, 8.100d0
     
     &             , 16.300d0, 16.100d0, 14.700d0, 13.300d0, 12.100d0, 11.300d0, 11.100d0, 10.300d0 ! G_0 = 4.61
     &             , 9.500d0, 8.569d0, 8.500d0
     
     &             , 17.700d0, 17.300d0, 15.700d0, 14.100d0, 12.900d0, 11.900d0, 11.100d0, 11.100d0 ! G_0 = 10
     &             , 10.100d0, 9.100d0, -1.0d0
     
     &             , 19.300d0, 19.100d0, 17.100d0, 15.300d0, 14.100d0, 13.041d0, 12.100d0, 11.100d0 ! G_0 = 21.7
     &             , 11.100d0, 10.071d0, -1.0d0
     
     &             , 21.700d0, 21.300d0, 18.900d0, 16.700d0, 15.300d0, 14.500d0, 13.293d0, 12.100d0 ! G_0 = 47.1
     &             , 12.100d0, 10.900d0, -1.0d0
     
     &             , 24.300d0, 23.900d0, 20.900d0, 18.500d0, 16.900d0, 15.900d0, 14.700d0, 14.103d0 ! G_0 = 102.4
     &             , 13.300d0, 11.900d0, -1.0d0
     
     &             , 27.500d0, 27.100d0, 23.700d0, 20.500d0, 18.500d0, 17.300d0, 15.900d0, 15.300d0 ! G_0 = 222.2
     &             , 14.300d0, 13.300d0, 13.015d0
     
     &             , 33.900d0, 30.700d0, 26.700d0, 22.700d0, 20.300d0, 18.700d0, 17.900d0, 17.300d0 ! G_0 = 482.5
     &             , 15.513d0, 14.700d0, 14.100d0 /
!__________________________________________________________________________________________________________________________________
!
! The command line. The first positional argument is the surface density image to be corrected, the second one, optional, is the
! temperature image fitted together with it. The wavebands that were fitted to make the image are given with -bands=.

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

       call getarg (  1, arg1 )
       call getarg (  2, arg2 )
       call getarg (  3, arg3 )
       call getarg (  4, arg4 )
       call getarg (  5, arg5 )
       call getarg (  6, arg6 )

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
       if (ia2 .gt. 0 .and. arg2(1:ia2-1) .eq. '-verb') cverbose = arg2(1:ia2)

       mapfile  = ' '
       tempfile = ' '
       nfit     = 4
       nwgt     = 0
       wfit(1)  = 160.0d0
       wfit(2)  = 250.0d0
       wfit(3)  = 350.0d0
       wfit(4)  = 500.0d0
       g0       = 10.0d0
       g0src    = '(default)'
       g16      = 0.0d0
       g84      = 0.0d0

       if (cverbose .ge. '-verb2') then
         write (iotty,'( )')
         if (arg1(1:ia1) .eq. ':') write (iotty,'(a,$)') '  '
         write (iotty,'(a,$)') ' HICORRECT '//dot//' Correct Bias of Surface Density Images '//dot//' '//compda
       endif
       if (arg1(1:ia1) .eq. ':') stop
       if (cverbose .ge. '-verb2') write (iotty,'()')
       if (cverbose .ge. '-verb2') write (iotty,'( a)') ' Alexander Men’shchikov, DAp IRFU CEA Saclay, France.'
       if (cverbose .ge. '-verb2') write (iotty,'( a)') ' Using '//clibname(lb:7)//' library version'//cfitsversion
     &                                                //' by William D Pence.'
       if (ia1 .eq. 0) then
         write (iotty,'( )')
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USAGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')
         write (iotty,'(a)') ' hicorrect <dimage> [<timage>] [-bands=<waves>] [-errors=<relerrs>] [-g0=<field>] [-verb{0|1|2|3}]'
         write (iotty,'( )')   
         write (iotty,'(a)') ' This utility corrects images of surface densities and temperatures derived'
         write (iotty,'(a)') ' from single-temperature fitting of Herschel images for the bias caused by'
         write (iotty,'(a)') ' distribution of dust temperatures along lines of sight.'
         write (iotty,'( )')   
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')   
         write (iotty,'(a)') ' <dimage> .... surface density image to be corrected'
         write (iotty,'(a)') ' <timage> .... dust temperature image fitted together with <dimage>'
         write (iotty,'(a)') ' <waves> ..... wavelengths (µm) of the wavebands fitted to make <dimage>'
         write (iotty,'(a)') ' <relerrs> ... relative calibration uncertainties of the wavebands'
         write (iotty,'(a)') ' <field> ..... field strength; without -g0=, it is estimated from the images'
         write (iotty,'( )')   
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write (iotty,'( )')   
         write (iotty,'(a)') ' Output images are written to the directory containing the input images.'
         write (iotty,'( )')
         write (iotty,'(a)') ' <dimage>.correct.fits ... corrected surface density image'
         write (iotty,'(a)') ' <dimage>.cfactor.fits ... correction factors image'
         write (iotty,'(a)') ' <timage>.correct.fits ... dust temperature image, from <dimage>.correct.fits'
         write (iotty,'( )')   
         write (iotty,'(a)') ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         stop 11
       endif  

       if (cverbose .eq. '-verb1')
     &    write (iotty,'(a)') ' HICORRECT: '//arg1(1:ia1)//' '//arg2(1:ia2)//' '//arg3(1:ia3)//' '//arg4(1:ia4)//' '//arg5(1:ia5)
     &                                      //' '//arg6(1:ia6)
       mapfile = arg1(1:ia1)

       if (arg2(1:1) .ne. '-') then
         tempfile = arg2(1:ia2)
       endif
       arg = ''
       if (arg5(1:7) .eq. '-bands=') arg = arg5
       if (arg4(1:7) .eq. '-bands=') arg = arg4
       if (arg3(1:7) .eq. '-bands=') arg = arg3
       if (arg2(1:7) .eq. '-bands=') arg = arg2
       if (arg .ne. '') then
         ia = lastc ( arg )
         nfit = 1
         do k=8,ia
           if (arg(k:k) .eq. ',') nfit = nfit + 1
         enddo
         if (nfit .gt. mband) then
           write (iotty,'(/a,i2,a/)') '   HICORRECT: ERROR: At most', mband, ' wavebands can be given with -bands=.'
           stop
         endif
         read (arg(8:ia),*,iostat=ios) (wfit(k), k=1,nfit)
         if (ios .ne. 0) then
           write (iotty,'(/a,a/)') '   HICORRECT: ERROR: Cannot read the wavelengths from ', arg(1:ia)
           stop
         endif
         if (nfit .lt. 3) then
           write (iotty,'(/a/)') '   HICORRECT: ERROR: Too few wavebands specified (nfit < 3).'
           stop
         endif
         if (int ( wfit(nfit) ) .le. 350) then
           write (iotty,'(/a,i3,a/)') '   HICORRECT: ERROR: The longest waveband of ', int ( wfit(nfit) ), ' µm is too short.'
           stop
         endif
       endif
       arg = ''
       if (arg5(1:8) .eq. '-errors=') arg = arg5
       if (arg4(1:8) .eq. '-errors=') arg = arg4
       if (arg3(1:8) .eq. '-errors=') arg = arg3
       if (arg2(1:8) .eq. '-errors=') arg = arg2
       if (arg .ne. '') then
         ia = lastc ( arg )
         nwgt = 1
         do k=9,ia
           if (arg(k:k) .eq. ',') nwgt = nwgt + 1
         enddo
         if (nwgt .gt. mband) then
           write (iotty,'(/a,i3,a/)') '   HICORRECT: ERROR: at most', mband, ' values can be given with -errors=.'
           stop
         endif
         read (arg(9:ia),*,iostat=ios) (relin(k), k=1,nwgt)
         if (ios .ne. 0) then
           write (iotty,'(/a,a/)') '   HICORRECT: ERROR: cannot read the calibration uncertainties from ', arg(1:ia)
           stop
         endif
       endif
       arg = ''
       if (arg5(1:4) .eq. '-g0=') arg = arg5
       if (arg4(1:4) .eq. '-g0=') arg = arg4
       if (arg3(1:4) .eq. '-g0=') arg = arg3
       if (arg2(1:4) .eq. '-g0=') arg = arg2
       if (arg .ne. '') then
         ia = lastc ( arg )
         read (arg(5:ia),*) g0
         g0src = '(from command line)'
       endif

! The file names may be given without the extension '.fits', which is then added.

       la = len_trim ( mapfile )
       if (la .lt. 5 .or. mapfile(max(la-4,1):la) .ne. '.fits') mapfile = mapfile(1:la)//'.fits'
       if (tempfile .ne. ' ') then
         la = len_trim ( tempfile )
         if (la .lt. 5 .or. tempfile(max(la-4,1):la) .ne. '.fits') tempfile = tempfile(1:la)//'.fits'
       endif
       
       inquire ( file=mapfile(1:len_trim(mapfile)), exist=estok )
       
       if (.not.estok) then
         write (iotty,'(/a/)') '   HICORRECT: ERROR: '''//mapfile(1:len_trim(mapfile))//''' does not exist.'
         stop
       endif
       if (tempfile .ne. ' ') then
       
         inquire ( file=tempfile(1:len_trim(tempfile)), exist=estok )
       
         if (.not.estok) then
           write (iotty,'(/a/)') '   HICORRECT: ERROR: '''//tempfile(1:len_trim(tempfile))//''' does not exist.'
           stop
         endif
       endif

       call when ( lunix, ctime, cdate, ndate, 4 )

       if (cverbose .gt. '-verb1') then
         write (iotty,'(/a,i4,a)') '   Computing effective attenuating columns (angular averages over', nang, ' directions)'
         write (iotty,'(a,1pe10.4,a,e10.4,a,e10.4,a,e10.4,a,0pf7.5,a)')
     &                             '   Shielding law T(s) = ', tfloor, ' + (', tsurf, ' - ', tfloor, ') * (1 + s / ', sstar
     &                           , ')**(-', gamma, ') K'
       endif
!__________________________________________________________________________________________________________________________________
!
! The wavebands. The relation has to be built for the wavebands that were fitted to make the image, given with -bands= (160, 250,
! 350 and 500 µm by default), because the bias depends on which colors were used to find the single temperature. They are sorted,
! and the longest of them is the reference waveband, the one the fit divides by.

       if (nfit .lt. 2) then
         write (iotty,'(/a/)') '   HICORRECT: ERROR: at least two wavebands must be given with -bands=.'
         stop
       endif

! The relative calibration uncertainties set the weights of the wavebands in the fit, as in hires: 20% at 160 µm and shorter and
! 10% at longer wavelengths by default, or the values of -errors=, given in the order of -bands= and sorted together with them.

       if (nwgt .gt. 0) then
         if (nwgt .ne. nfit) then
           write (iotty,'(/a/)') '   HICORRECT: ERROR: -errors= must give one value for each waveband of -bands=.'
           stop
         endif
         do k=1,nfit
           if (relin(k) .le. 0.0d0) then
             write (iotty,'(/a/)') '   HICORRECT: ERROR: the values of -errors= must be positive.'
             stop
           endif
           relerr(k) = relin(k)
         enddo
       else
         do k=1,nfit
           if (wfit(k) .le. 160.5d0) then
             relerr(k) = 0.20d0
           else
             relerr(k) = 0.10d0
           endif
         enddo
       endif
       do k=2,nfit
         do ib=nfit,k,-1
           if (wfit(ib-1) .gt. wfit(ib)) then
             bb          = wfit(ib)
             wfit(ib)    = wfit(ib-1)
             wfit(ib-1)  = bb
             bb          = relerr(ib)
             relerr(ib)  = relerr(ib-1)
             relerr(ib-1)= bb
           endif
         enddo
       enddo
       waveref = wfit(nfit)

       do k=1,nfit
         cnu(k) = 1.0d17 * kappa0 * (speedolight / (wfit(k) * 1.0d-4) / nu0)**betap * dtg * muh2 * amu
       enddo

       if (cverbose .gt. '-verb1') then
         write(cnfit,'(i1)') nfit
         la = lastc ( cnfit )
         write (iotty,'(a,'//cnfit(1:la)//'(1x,i3),a,i3,a)') '   Wavebands fitted:', ( int ( wfit(k) ), k=1,nfit )
     &                                                     , ' µm  reference wavelength: ', int ( waveref ), ' µm'
         write (iotty,'(a,'//cnfit(1:la)//'(1x,f4.2))') '   Relative calibration uncertainties:', ( relerr(k), k=1,nfit )
       endif
!__________________________________________________________________________________________________________________________________
!
! The relation. For each of ntab true surface densities the emission is integrated over ndep layers, one modified blackbody is
! fitted to it on a grid of ntfit temperatures exactly as the reconstruction fits it, and the surface density that the fit returns
! is recorded. The table therefore runs over the fitted surface density, which is what the image holds, and gives the factor that
! restores the true one. Both arrays increase with the index, so the interpolation below is straightforward.

       do k=1,ntfit
         tgrid(k) = 0.98d0 * tfloor * (1.02d0 * tsurf / (0.98d0 * tfloor))**(dble ( k - 1 ) / dble ( ntfit - 1 ))
       enddo

! The calibration of the factor at low surface densities on a radmc-3d model of a diffuse cloud (Sect. 4.3 of the paper).
!
! Why. The analytic layer needs almost no correction below 10^21 cm^-2, but a radiative transfer model of an extended cloud needs
! about 5% even for thin columns and 15-25% at 2-4.5 x 10^21 cm^-2. The surface layers of the cloud are heated by unattenuated
! ultraviolet light from half of the sky and are up to 3 K warmer than the relation T(s) predicts, which a gray cross section
! cannot describe. Below 10^xb0 the analytic layer is therefore replaced by this cloud: the factor is the average of the analytic
! sphere and the cloud, and the change from the analytic average of sphere and layer is reduced smoothly to zero at 10^xb1.
!
! For which images. For every image: the calibration is part of the standard correction, not an adjustment for other wavebands.
! The factor measured on the cloud depends on the wavebands fitted and their weights, as the analytic factor does, so it was
! measured for ncset sets. The set of the image is found among them, the uncertainties compared after division by their minimum,
! because only their ratios set the weights of the fit. For any other set, the table of the standard set (the first: 160-500 µm,
! weights 20/10/10/10%) is used with a warning; the error is then several percent in diffuse material.
!
! How. The cloud is the 3D cube of make_cirrus3d.py with the central enhancement switched off: 448 x 448 columns of 12 arcsec
! (6.8 x 6.8 pc at 260 pc), 2 pc deep, whose projection is the cloud map of the simulated sky; margins of 2 pc on every side keep
! the heating through the sides away from the mapped area. The emission of each column was integrated along the line of sight
! without convolution (cloud_emission.py), so the tables are relations for single lines of sight, independent of the resolution
! of an image. Each column was fitted with one modified blackbody exactly as here, and the tables cirtab are the medians of the
! true over the fitted surface density in bins of 0.1 dex of the fitted surface density centered on cirx (log10, cm^-2), made by
! cloud_tables.py. A bin with too few columns repeats the one below it; beyond the ends of a table its end values are used.

       cirx = (/ 20.55d0, 20.65d0, 20.75d0, 20.85d0, 20.95d0, 21.05d0, 21.15d0, 21.25d0, 21.35d0, 21.45d0, 21.55d0, 21.65d0 /)
       cirw = 0.0d0
       cire = 0.0d0

! 160-500 µm, weights 20/10/10/10 percent (standard)

       ncirn(1) = 4
       cirw(1:4,1) = (/ 160.0d0, 250.0d0, 350.0d0, 500.0d0 /)
       cire(1:4,1) = (/ 2.0d0, 1.0d0, 1.0d0, 1.0d0 /)
       cirtab(:,1) = (/ 1.049d0, 1.051d0, 1.053d0, 1.058d0, 1.066d0, 1.075d0,
     &                  1.101d0, 1.121d0, 1.146d0, 1.180d0, 1.225d0, 1.251d0 /)

! 160-500 µm, equal weights

       ncirn(2) = 4
       cirw(1:4,2) = (/ 160.0d0, 250.0d0, 350.0d0, 500.0d0 /)
       cire(1:4,2) = (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)
       cirtab(:,2) = (/ 1.057d0, 1.058d0, 1.060d0, 1.067d0, 1.076d0, 1.086d0,
     &                  1.120d0, 1.139d0, 1.171d0, 1.214d0, 1.267d0, 1.294d0 /)

! 250-500 µm, equal weights

       ncirn(3) = 3
       cirw(1:3,3) = (/ 250.0d0, 350.0d0, 500.0d0 /)
       cire(1:3,3) = (/ 1.0d0, 1.0d0, 1.0d0 /)
       cirtab(:,3) = (/ 1.032d0, 1.034d0, 1.035d0, 1.039d0, 1.043d0, 1.049d0,
     &                  1.065d0, 1.076d0, 1.091d0, 1.110d0, 1.140d0, 1.160d0 /)

! 100-500 µm, weights 20/20/10/10/10 percent

       ncirn(4) = 5
       cirw(1:5,4) = (/ 100.0d0, 160.0d0, 250.0d0, 350.0d0, 500.0d0 /)
       cire(1:5,4) = (/ 2.0d0, 2.0d0, 1.0d0, 1.0d0, 1.0d0 /)
       cirtab(:,4) = (/ 1.088d0, 1.090d0, 1.094d0, 1.104d0, 1.120d0, 1.141d0,
     &                  1.194d0, 1.231d0, 1.282d0, 1.364d0, 1.430d0, 1.430d0 /)

! 70-500 µm, weights 20/20/20/10/10/10 percent

       ncirn(5) = 6
       cirw(1:6,5) = (/ 70.0d0, 100.0d0, 160.0d0, 250.0d0, 350.0d0, 500.0d0 /)
       cire(1:6,5) = (/ 2.0d0, 2.0d0, 2.0d0, 1.0d0, 1.0d0, 1.0d0 /)
       cirtab(:,5) = (/ 1.126d0, 1.128d0, 1.136d0, 1.150d0, 1.175d0, 1.216d0,
     &                  1.282d0, 1.341d0, 1.421d0, 1.536d0, 1.600d0, 1.600d0 /)

       emin = minval ( relerr(1:nfit) )
       icset = 0
       do ic=1,ncset
         if (ncirn(ic) .ne. nfit) cycle
         lm = 0
         do k=1,nfit
           if (abs ( wfit(k) - cirw(k,ic) ) .gt. 0.5d0 .or. abs ( relerr(k) / emin - cire(k,ic) ) .gt. 1.0d-3) lm = 1
         enddo
         if (lm .eq. 0) then
           icset = ic
           exit
         endif
       enddo
       if (icset .eq. 0) then
         if (cverbose .gt. '-verb1') then
           write (iotty,'(a)') '   HICORRECT: WARNING: no calibration on the diffuse cloud for these wavebands and weights;'
           write (iotty,'(a)') '   that of 160-500 µm with the standard weights is used, and diffuse material may be corrected'
           write (iotty,'(a)') '   by several percent too much or too little.'
         endif
         cirf = cirtab(:,1)
       else
         cirf = cirtab(:,icset)
       endif

       xb0 = 21.6d0
       xb1 = 22.0d0

       do igeo=1,2    !<~ 1: sphere and layer averaged; 2: sphere alone
         do i=1,ntab
           sig = smin * (smax / smin)**(dble ( i - 1 ) / dble ( ntab - 1 ))

! The effective attenuating column at each depth, averaged between a uniform sphere sampled along a diameter and a uniform
! plane-parallel layer sampled along its normal (Eq. 4 of the paper); igeo = 2 gives the sphere alone, needed for the calibration
! on the diffuse cloud. Lengths are in units of the half thickness (the radius of the sphere, half the thickness of the layer), so
! a path of length l holds the column (sig / 2) l, and t0 is the optical depth of the half thickness at the wavelengths that heat
! the dust. The depth zz runs from -1 to 1 in ndep equal steps; the fractional depth u of the paper is (zz + 1) / 2, and equal
! steps hold equal masses because the density is uniform. For every depth the transmitted fraction exp (-t0 l) is averaged over
! nang directions: for the sphere, the angle psi to the outward radial direction, weighted by sin psi (the solid angle), with the
! chord l = lsph to the boundary; for the layer, the cosine xmu of the angle to the normal, spaced uniformly (equal solid angles),
! with the path l = dlay to the front or the back face. The effective column is the one that gives the same fraction in all
! directions: -ln <exp (-t0 l)> / t0, in units of the half thickness, times sig / 2.

           t0 = 0.5d0 * kapheat * sig
           do j=1,ndep
             zz   = (dble ( j ) - 0.5d0) / dble ( ndep ) * 2.0d0 - 1.0d0
             jsph = 0.0d0
             jlay = 0.0d0
             wsum = 0.0d0
             do l=1,nang
               psi  = (dble ( l ) - 0.5d0) / dble ( nang ) * pi
               wps  = sin ( psi )
               lsph = -abs ( zz ) * cos ( psi ) + sqrt ( max ( 1.0d0 - (abs ( zz ) * sin ( psi ))**2, 0.0d0 ) )
               jsph = jsph + exp ( -t0 * lsph ) * wps
               wsum = wsum + wps
               xmu  = (dble ( l ) - 0.5d0) / dble ( nang ) * 2.0d0 - 1.0d0
               if (xmu .gt. 0.0d0) then
                 dlay = (1.0d0 - zz) / abs ( xmu )
               else
                 dlay = (1.0d0 + zz) / abs ( xmu )
               endif
               jlay = jlay + exp ( -t0 * dlay )
             enddo
             jsph = jsph / wsum
             jlay = jlay / dble ( nang )
             ssph = -log ( max ( jsph, 1.0d-300 ) ) / max ( t0, 1.0d-12 )
             slay = -log ( max ( jlay, 1.0d-300 ) ) / max ( t0, 1.0d-12 )
             if (igeo .eq. 1) then
               seff(j) = 0.5d0 * (ssph + slay) * sig / 2.0d0
             else
               seff(j) = ssph * sig / 2.0d0
             endif
           enddo

! The intensities of the line of sight (Eq. 6): optically thin emission, the temperature of each depth from the relation T(s) of
! Eq. (5); for a uniform body the integral over the fractional depth is the mean over the depths. The weights of the fit are
! inversely proportional to the squares of the calibration uncertainties, as in hires.

           do k=1,nfit
             obs(k) = 0.0d0
             do j=1,ndep
               tt = tfloor + (tsurf - tfloor) * (1.0d0 + seff(j) / sstar)**(-gamma)
               obs(k) = obs(k) + planckf ( tt, wfit(k), hplanck, kboltz, speedolight )
             enddo
             obs(k) = cnu(k) * sig * obs(k) / dble ( ndep )
             wgt(k) = 1.0d0 / (relerr(k) * obs(k))**2
           enddo

! The fit: a search over the ntfit temperatures of tgrid, which spans the whole range of the relation. For each temperature the
! best amplitude follows from linear least squares, and the temperature with the smallest chi^2 is the fitted one.

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
               chi2 = chi2 + wgt(l) * (obs(l) - amp * model(l))**2
             enddo
             if (chi2 .lt. chibest) then
               chibest = chi2
               tbest   = tgrid(k)
             endif
           enddo

! The fitted surface density (Eq. 7) is the intensity of the reference waveband divided by the emission per molecule at the fitted
! temperature (not the amplitude of the fit), as hires derives it, and the factor (Eq. 8) is the true over the fitted value.

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

! The calibration of the diffuse end on the radmc-3d cloud, as described above. ddx is the change from the analytic average of
! sphere and layer (fcor) to the calibrated factor, the average of the sphere (fsx, interpolated at the same fitted surface
! density) and the cloud (fcx, interpolated in the table cirf). Below 10^xb0 the full change is applied. Above, the change is
! frozen at its value d0 at 10^xb0 and tapered to zero at 10^xb1 by a half cosine wb, so that the factor stays continuous and
! smooth and equals the analytic average above 10^xb1.
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
           fcx = cirf(k) + (cirf(k+1) - cirf(k)) * (xg - cirx(k)) / (cirx(k+1) - cirx(k))
         endif
         ddx(i) = 0.5d0 * (fsx + fcx) - fcor(i)
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
           wb = 0.5d0 * (1.0d0 + cos ( pi * (xg - xb0) / (xb1 - xb0) ))
           dd = d0
         endif
         fcor(i) = fcor(i) + wb * dd
       enddo

! The factor made non-decreasing with the fitted surface density, and the field of the cloud. Above about 3 x 10^22 cm^-2 the
! relation of uniform bodies declines, because such a body is nearly isothermal at the floor temperature. The radiative transfer
! models of concentrated cores do not decline, so the factor is held at its maximum there. The hold applies to the factor of
! G_0 = 10 only: the field factor, measured on the radiative transfer models, declines in strong fields, and that decline is real.

       do i=2,ntab
         fcor(i) = max ( fcor(i), fcor(i-1) )
       enddo

!__________________________________________________________________________________________________________________________________
!
! The image to be corrected.

       la = len_trim ( mapfile )
       
       call getfitshead ( mapfile(1:la), nx, ny, dx, dy, bunit )

       allocate ( smap(nx,ny), valid(nx,ny), delta(nx,ny), corr(nx,ny), work(nx,ny), iref(nx,ny) )

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

! The outputs are written to the directory containing the input images.

       ls = la
       if (mapfile(max(ls-4,1):ls) .eq. '.fits') ls = ls - 5
!!       k  = index ( mapfile(1:ls), '/', back=.true. )
       k = 0
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
       if (istemp ( nx, ny, smap )) then
         write (iotty,'(/a/a/)') '   HICORRECT: ERROR: '''//mapfile(1:len_trim(mapfile))//''' looks like a temperature image.'
     &                         , '   The surface density image to be corrected must be given first.'
         stop
       endif
!__________________________________________________________________________________________________________________________________
!
! The radiation field. Without -g0=, it is estimated from the fitted temperatures of the diffuse material (Appendix C.2 of the
! paper), if the temperature image is given. Reading that image overwrites some header values of the image, which are restored.

       if (g0src .eq. '(default)' .and. tempfile .ne. ' ') then
         la = len_trim ( tempfile )
         sbunit = bunit
         sblank = blank

         call getfitshead ( tempfile(1:la), nxi, nyi, dxi, dyi, bunit )

         if (nxi .eq. nx .and. nyi .eq. ny) then

           call rfits ( nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, iref, dxi, dyi, object, ra, dec
     &                , tempfile(1:la), iotty, iolog, creator, beamkw, funmin, funmax, blank, crota1, crota2, cd11, cd12
     &                , cd21, cd22, equinox, bzero, bscale, wavekw, datamn, datamx, history, cverbose )

           if (.not.istemp ( nx, ny, iref )) then
             write (iotty,'(/a/)') '   HICORRECT: ERROR: '''//tempfile(1:la)//''' does not look like a temperature image.'
             stop
           endif

           call estfield ( nx, ny, smap, iref, valid, ttab, ntb, gfld, ngf, tx0, tdx, txfield, tminpix, cverbose, iotty
     &                   , g16, g50, g84, hotfr, ndif, estok )

           if (estok) then
             ee = 10.0d0**(floor ( log10 ( g50 ) ) - 1)
             g0 = anint ( g50 / ee ) * ee
             g0src = '(estimated from images)'
             if (cverbose .gt. '-verb1') then
               write (iotty,'(a,i0,a,f4.1,a,f4.1,a)') '   Radiation field G_0 estimated from ', ndif
     &                                  , ' pixels in the diffuse material (10^', tx0, ' to 10^', txfield, ' cm^-2):'
               write (iotty,'(a,f8.2,a,f8.2,a,f8.2,a)') '   G_0 =', g50, ' (', g16, ' to', g84
     &                                  , ' at the 16th and 84th percentiles of the temperature)'
               write (iotty,'(a,0pf6.1,a)') '   Pixels warmer than the strongest field allows:', 100.0d0 * hotfr, '%'
             endif

! The area of the healthy pixels: surface densities of at least 10^19 cm^-2 (lower values, which some methods produce around the
! usable area, are taken as bad) and positive temperatures. The smallest observed Herschel maps have about 2000 x 2000 pixels of
! 3 arcsec; a smaller area is probably a cutout, whose diffuse material may not represent the field of the cloud.

             ngood = 0
             do j=1,ny
               do i=1,nx
                 if (valid(i,j) .eq. 1 .and. smap(i,j) .ge. 1.0d19 .and. iref(i,j) .gt. 0.0d0
     &               .and. .not.isnan ( iref(i,j) )) ngood = ngood + 1
               enddo
             enddo
             if (dx .gt. 0.0d0 .and. dy .gt. 0.0d0) then
               ee = dble ( ngood ) * dx * dy / 6000.0d0**2
             else
               ee = dble ( ngood ) / 2000.0d0**2
             endif
             if (cverbose .gt. '-verb1') then
               if (ee .lt. 1.0d0) then
                 write (iotty,'(a)') '   HICORRECT: WARNING: the image appears to be a cutout of a larger map;'
     &                             //' estimate the field from the full map.'
               endif
               if (beamsd .gt. 0.0d0 .and. abs ( beamsd / 13.5d0 - 1.0d0 ) .gt. 1.0d-2) then
                 write (iotty,'(a,f7.2,a)') '   HICORRECT: WARNING: the field was estimated with the model temperatures of'
     &                                    //' 13.5 arcsec images; at', beamsd, ' arcsec it is approximate.'
               endif
               if (icset .ne. 1) then
                 write (iotty,'(a)') '   HICORRECT: WARNING: The field was estimated with the model temperatures of 160-500 µm'
     &                             //' with standard weights; for these wavebands or weights it is approximate.'
               endif
             endif
           else
             if (cverbose .gt. '-verb1') then
               write (iotty,'(a,f4.1,a,f4.1,a)') '   HICORRECT: WARNING: Too few pixels between 10^', tx0, ' and 10^', txfield
     &                                         , ' cm^-2 to estimate the field; G_0 = 10 is used'
             endif
           endif
         else
           if (cverbose .gt. '-verb1') then
             write (iotty,'(a)') '   HICORRECT: WARNING: Temperature image does not match the image;'
     &                         //' the field cannot be estimated, and G_0 = 10 is used'
           endif
         endif
         bunit = sbunit
         blank = sblank
       else
         if (g0src .eq. '(default)') then
           if (cverbose .gt. '-verb1') write (iotty,'(a)') '   No temperature image given, so the field strength'
     &                                                   //' cannot be estimated'
         endif
       endif

       if (cverbose .gt. '-verb1') then
         write (iotty,'(a,1pe8.2,a)') '   Radiation field G_0 = ', g0, ' '//g0src(1:len_trim(g0src))
         if (g0 .lt. gfld(1) .or. g0 .gt. gfld(ngf)) then
           write (iotty,'(a,1pe8.2,a)') '   HICORRECT: WARNING: G_0 =', g0, ' is outside the grid;'
     &                                //' the field factor g for the nearest G_0 is used'
         endif
         if (g0 .ne. 10.0d0 .and. icset .ne. 1) then
           write (iotty,'(a)') '   HICORRECT: WARNING: The field factor was measured for 160-500 µm and standard weights;'
     &                       //' for other wavebands or weights it is approximate'
         endif
       endif

! The field factor multiplies the factor held at its maximum (Eq. 9 of the paper).

       if (g0 .ne. 10.0d0) then
         do i=1,ntab
           fcor(i) = fcor(i) * fieldfac ( log10 ( ffit(i) ), g0, gtab, gfld, ngx, ngf, gx0, gdx )
         enddo
       endif

       if (cverbose .eq. '-verb3') then
         write (iotty,'(a)') '   Tabulated relation: fitted surface density and factor that restores the true one'
         do i=1,ntab,20
           write (iotty,'(4x,1pe11.4,0pf8.4)') ffit(i), fcor(i)
         enddo
       endif
!__________________________________________________________________________________________________________________________________
!
! The correction of the image. The factor is found by
! interpolating the table in the logarithm of the surface density, and it is held at the ends of the table rather than
! extrapolated.

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
         if (cverbose .gt. '-verb1') then
           write (iotty,'(a,1pe12.6,a,2(1x,0pf6.4))') '   Valid pixels in the image: ', dble ( nvalid )
     &                                              , '  minmax corrections:', funmin, funmax
           write (iotty,'(a,f7.4)') '   Ratio of the corrected to uncorrected sums:', tt/ss
         endif
       endif

       fname = stem(1:ls)//'.cfactor.fits'
       lfn = lastc ( fname )

       call putima ( iotty, cverbose, cfitsversion, nx, ny, corr, valid, rnan, fname(1:lfn), 'ratio', rhead, chead, blank, cdate
     &             , ctime )

       call addkeys ( fname, g0, g0src, g16, g84, nfit, wfit, relerr, tfloor, tsurf, sstar, gamma )

       do j=1,ny
         do i=1,nx
           work(i,j) = smap(i,j) + delta(i,j)
         enddo
       enddo

       fname = stem(1:ls)//'.correct.fits'
       lfn = lastc ( fname )

       call putima ( iotty, cverbose, cfitsversion, nx, ny, work, valid, rnan, fname(1:lfn), bunit, rhead, chead, blank, cdate
     &             , ctime )

       call addkeys ( fname, g0, g0src, g16, g84, nfit, wfit, relerr, tfloor, tsurf, sstar, gamma )
!__________________________________________________________________________________________________________________________________
!
! The dust temperature that goes with the corrected surface density, at the resolution of the corrected image. The reconstruction
! derived the surface density from the reference waveband with the temperature it fitted at that resolution. With f times more
! material on the line of sight, the same intensity of the reference waveband is reproduced by colder dust, whose Planck function
! is that of the fitted temperature divided by f. The fitted temperature image of the same resolution is the second image given.

       if (tempfile .ne. ' ') then

         la = len_trim ( tempfile )

         call getfitshead ( tempfile(1:la), nxi, nyi, dxi, dyi, bunit )

         if (nxi .eq. nx .and. nyi .eq. ny) then

           call rfits ( nx, ny, bunit, ctype1, ctype2, crpix1, crpix2, crval1, crval2, iref, dxi, dyi, object, ra, dec
     &                , tempfile(1:la), iotty, iolog, creator, beamkw, funmin, funmax, blank, crota1, crota2, cd11, cd12
     &                , cd21, cd22, equinox, bzero, bscale, wavekw, datamn, datamx, history, cverbose )
           
           if (.not.istemp ( nx, ny, iref )) then
             write (iotty,'(/a/)') '   HICORRECT: ERROR: '''//tempfile(1:la)//''' does not look like a temperature image.'
             stop
           endif
           if (beamsd .gt. 0.0d0 .and. beamkw .gt. 0.0d0) then
             if (abs ( beamkw / beamsd - 1.0d0 ) .gt. 1.0d-2) then
               write (iotty,'(a,f7.2,a,f7.2,a)') '   HICORRECT: ERROR: the temperature image has the resolution', beamkw
     &                                         , ' arcsec, the surface density image', beamsd, ' arcsec.'
               write (iotty,'(a)') '   They must have the same resolution; no temperature map is written.'
             endif

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
!!           k  = index ( tempfile(1:lt), '/', back=.true. )
           k = 0
           fname = tempfile(k+1:lt)//'.correct.fits'
           lfn = lastc ( fname )

           call putima ( iotty, cverbose, cfitsversion, nx, ny, work, valid, rnan, fname(1:lfn), 'K', rhead, chead, blank, cdate
     &                 , ctime )

           call addkeys ( fname, g0, g0src, g16, g84, nfit, wfit, relerr, tfloor, tsurf, sstar, gamma )

         else
           if (cverbose .gt. '-verb1') then
             write (iotty,'(a)') '   HICORRECT: WARNING: the temperature image does not match the corrected image,'
             write (iotty,'(a)') '   so no temperature map is written.'
           endif
         endif
!!       else
!!         if (cverbose .gt. '-verb1') write (iotty,'(a)') '   No temperature image given, so no temperature map is written.'
       endif

       deallocate ( smap, valid, delta, corr, work, iref )

       if (cverbose .gt. '-verb1') write (iotty,'(/a/)') ' Done.'
       stop
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function planckf 
         
     &          ( tdust, wave, hplanckanck, kboltzltz, speedolight )
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!
! Planck function per unit frequency, in cgs units, for a wavelength in µm.
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       real*8        planckf, tdust, wave, frequency, xexp, hplanckanck, kboltzltz, speedolight
!__________________________________________________________________________________________________________________________________
!
       frequency = speedolight / (wave * 1.0d-4)
       xexp = hplanckanck * frequency / (kboltzltz * tdust)
       xexp = min ( max ( xexp, 1.0d-8 ), 7.0d2 )

       planckf = 2.0d0 * hplanckanck * frequency**3 / speedolight**2 / (exp ( xexp ) - 1.0d0)

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function invplanck

     &          ( b, wave, hplanck, kboltz, speedolight )
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!
! The temperature whose Planck function equals b at the wavelength wave in µm. The inverse of the function below.
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       real*8        invplanck, b, wave, hplanck, kboltz, speedolight, nu, pref
!__________________________________________________________________________________________________________________________________
!
       nu = speedolight / (wave * 1.0d-4)
       pref = 2.0d0 * hplanck * nu**3 / speedolight**2

       invplanck = hplanck * nu / (kboltz * log ( 1.0d0 + pref/max ( b, 1.0d-40 ) ))

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine interp

     &            ( n, x, y, xq, yq )
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!
! Linear interpolation of y against x in the logarithm of x, for a table whose x increases with the index. Beyond either end the
! value of that end is returned rather than an extrapolation, because the relation is only known over the range that was tabulated.
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       integer       n, lo, hi, mid
       real*8        x(n), y(n), xq, yq, w
!__________________________________________________________________________________________________________________________________
!
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
         mid = (lo + hi) / 2
         if (x(mid) .gt. xq) then
           hi = mid
         else
           lo = mid
         endif
         goto 100
       endif

       w  = (log ( xq ) - log ( x(lo) )) / (log ( x(hi) ) - log ( x(lo) ))
       yq = y(lo) + w * (y(hi) - y(lo))

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine addkeys
         
     &          ( fname, g0, g0src, g16, g84, nfit, wfit, relerr, tfloor, tsurf, sstar, gamma )
!__________________________________________________________________________________________________________________________________
!
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!
! Add to the header of an output image the radiation field of the correction, how it was obtained, its spread if it was estimated,
! the wavebands fitted to make the image with their relative calibration uncertainties, and the four constants of the temperature
! relation, so that every output records how it was corrected.
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       character*80  bands, errs
       character*8   cw
       character*(*) fname, g0src

       integer       nfit, k, unit, status, block, la
       real*8        g0, g16, g84, wfit(nfit), relerr(nfit), tfloor, tsurf, sstar, gamma

       external      ftgiou, ftopen, ftpkyd, ftpkys, ftclos, ftfiou, printerr
!__________________________________________________________________________________________________________________________________
!
       bands = ' '
       errs  = ' '
       do k=1,nfit
         write (cw,'(f5.3)') relerr(k)
         la = len_trim ( errs )
         if (la .eq. 0) then
           errs = adjustl ( cw )
         else
           errs = errs(1:la)//','//adjustl ( cw )
         endif
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
       call ftpkyd ( unit, 'HIG0'   , g0    , 6, 'radiation field of the correction (Habing units)', status )
       call ftpkys ( unit, 'HIG0SRC', g0src ,    'radiation field given, estimated or default', status )
       if (g0src(1:9) .eq. 'estimated') then
         call ftpkyd ( unit, 'HIG0P16', g16 , 6, 'estimated field at the 16th percentile of T', status )
         call ftpkyd ( unit, 'HIG0P84', g84 , 6, 'estimated field at the 84th percentile of T', status )
       endif
       call ftpkys ( unit, 'HIBANDS', bands ,    'wavebands fitted to make the image (µm)', status )
       call ftpkys ( unit, 'HIERRS' , errs  ,    'relative calibration errors of HIBANDS', status )
       call ftpkyd ( unit, 'HITMIN' , tfloor, 8, 'T_min of the temperature relation (K)', status )
       call ftpkyd ( unit, 'HIT0'   , tsurf , 8, 'T_max of the temperature relation (K)', status )
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
! GETSF • Multi-Scale Multi-Wavelength Source & Filament Extraction • Alexander Men'shchikov, DAp IRFU CEA Saclay
!
! Write one output image. The pixels outside the coverage (not positive in the input image) are made undefined first, because
! they carry no result. The header of the surface density image is carried in rhead and chead and is written back unchanged, so that
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
! The image is copied before the pixels outside the coverage are made undefined, so that writing an image never alters it: the
! arrays written are used again afterwards (the factors corr for the temperature image).

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

       if (cverbose .gt. '-verb1') write (iotty,'(a)') '   Writing output image '''//fname//''''

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
!__________________________________________________________________________________________________________________________________
!
       t = log10 ( min ( max ( g0, gfld(1) ), gfld(ngf) ) )
       j = 1
       do while (j .lt. ngf-1 .and. t .gt. log10 ( gfld(j+1) ))
         j = j + 1
       enddo
       a = (t - log10 ( gfld(j) )) / (log10 ( gfld(j+1) ) - log10 ( gfld(j) ))

       xi = (x - gx0) / gdx + 1.0d0
       xi = min ( max ( xi, 1.0d0 ), dble ( ngx ) )
       i  = min ( int ( xi ), ngx - 1 )
       b  = xi - dble ( i )

       r1 = (1.0d0 - b) * log ( gtab(i,j)   ) + b * log ( gtab(i+1,j)   )
       r2 = (1.0d0 - b) * log ( gtab(i,j+1) ) + b * log ( gtab(i+1,j+1) )

       fieldfac = exp ( (1.0d0 - a) * r1 + a * r2 )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine estfield

     &          ( nx, ny, smap, tmap, valid, ttab, ntb, gfld, ngf, tx0, tdx, txfield, minpix, cverbose, iotty
     &          , g16, g50, g84, hotfr, ndif, estok )
!__________________________________________________________________________________________________________________________________
!
! Radiation field of a cloud from its fitted surface density image smap and temperature image tmap, as in Appendix C.2 of the
! paper. The pixels are binned in 0.2 dex of the fitted surface density. In each bin with at least minpix pixels and models in at
! least two fields, the 16th, 50th and 84th percentiles of the temperature are converted into fields by interpolating log G_0
! against the median model temperatures ttab of that bin, made non-decreasing with the field; a temperature outside them gives
! the weakest or the strongest field. The fields of the diffuse bins, below 10^txfield cm^-2, are combined into medians weighted
! by the pixels of the bins: g16, g50 and g84. hotfr is the fraction of the binned pixels warmer than the strongest field allows,
! ndif the number of diffuse pixels used; estok is false if there are none.
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       logical       estok
       character*1   fg(3)
       character*13  cnote
       character*(*) cverbose

       integer       nx, ny, ntb, ngf, minpix, iotty, ndif, valid(nx,ny)
       integer       i, j, k, n, m, npx, nall, nhot, nd, ib
       integer, allocatable :: ibin(:)

       real*8        smap(nx,ny), tmap(nx,ny), ttab(ntb,ngf), gfld(ngf), tx0, tdx, txfield, g16, g50, g84, hotfr
       real*8,  allocatable :: tb(:)
       real*8        gs(ngf), tm(ngf), tp(3), fl(3), e0, e1, wmedf, pctl, fieldof, tprev
       real*8        dv16(ntb), dv50(ntb), dv84(ntb), dw(ntb), dw2(ntb)

       external      hsort, wmedf, pctl, fieldof
!__________________________________________________________________________________________________________________________________
!
       n = nx * ny

       allocate ( ibin(n), tb(n) )

! The bin of every pixel with a positive surface density and temperature; 0 outside the bins.

       m = 0
       do j=1,ny
         do i=1,nx
           m = m + 1
           ibin(m) = 0
           if (valid(i,j) .eq. 1 .and. tmap(i,j) .gt. 0.0d0 .and. .not.isnan ( tmap(i,j) )) then
             e0 = ( log10 ( smap(i,j) ) - tx0 ) / tdx
             if (e0 .ge. 0.0d0 .and. e0 .lt. dble ( ntb )) ibin(m) = int ( e0 ) + 1
           endif
         enddo
       enddo

       if (cverbose .eq. '-verb2') then
         write (iotty,'(/a)') '   fields that the temperatures of the image require, in bins of the fitted surface density:'
         write (iotty,'(a)')  '   log10 Sigma_F       pixels   T (K) at 16, 50, 84%     G_0 at 16, 50, 84%           too warm'
       endif

       nall = 0
       nhot = 0
       nd   = 0
       ndif = 0

       do ib=1,ntb
         e0 = tx0 + tdx * dble ( ib - 1 )
         e1 = e0 + tdx

! The fields with models in this bin and their temperatures, made non-decreasing with the field.

         k = 0
         tprev = 0.0d0
         do j=1,ngf
           if (ttab(ib,j) .gt. 0.0d0) then
             k = k + 1
             gs(k) = gfld(j)
             tm(k) = max ( ttab(ib,j), tprev )
             tprev = tm(k)
           endif
         enddo
         npx = 0
         do m=1,n
           if (ibin(m) .eq. ib) then
             npx = npx + 1
             tb(npx) = tmap(mod(m-1,nx)+1,(m-1)/nx+1)
           endif
         enddo
         if (npx .lt. minpix .or. k .lt. 2) cycle
         call hsort ( npx, tb )
         tp(1) = pctl ( npx, tb, 16.0d0 )
         tp(2) = pctl ( npx, tb, 50.0d0 )
         tp(3) = pctl ( npx, tb, 84.0d0 )
         do i=1,3
           fl(i) = fieldof ( tp(i), k, gs, tm )
           fg(i) = ' '
           if (tp(i) .lt. tm(1)) fg(i) = '<'
           if (tp(i) .gt. tm(k)) fg(i) = '>'
         enddo
         m = 0
         do i=1,npx
           if (tb(i) .gt. tm(k)) m = m + 1
         enddo
         nall = nall + npx
         nhot = nhot + m
         cnote = '   (not used)'
         if (e1 .le. txfield + 1.0d-6) then
           cnote = ' '
           nd = nd + 1
           dv16(nd) = fl(1)
           dv50(nd) = fl(2)
           dv84(nd) = fl(3)
           dw(nd)   = dble ( npx )
           ndif = ndif + npx
         endif
         if (cverbose .eq. '-verb2') then
           write (iotty,'(3x,f5.1,a,f5.1,i11,3x,3f7.1,3x,3(1x,a1,1pg9.3),0pf9.1,2a)') e0, ' --', e1, npx, tp
     &          , ( fg(i), fl(i), i=1,3 ), 100.0d0 * dble ( m ) / dble ( npx ), '%', cnote
         endif
       enddo

       deallocate ( ibin, tb )

       estok = nd .gt. 0
       hotfr = 0.0d0
       if (nall .gt. 0) hotfr = dble ( nhot ) / dble ( nall )
       if (.not.estok) return

       dw2 = dw
       g16 = wmedf ( nd, dv16, dw2 )
       dw2 = dw
       g50 = wmedf ( nd, dv50, dw2 )
       dw2 = dw
       g84 = wmedf ( nd, dv84, dw2 )

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function fieldof

     &          ( t, k, gs, tm )
!__________________________________________________________________________________________________________________________________
!
! The field that reproduces the temperature t, interpolating log G_0 linearly against the non-decreasing model temperatures tm of
! the k fields gs, as numpy.interp does; the weakest or the strongest field outside them.
!__________________________________________________________________________________________________________________________________
!
       implicit      none
       integer       k, j
       real*8        fieldof, t, gs(k), tm(k)
!__________________________________________________________________________________________________________________________________
!
       if (t .lt. tm(1)) then
         fieldof = gs(1)
       elseif (t .ge. tm(k)) then
         fieldof = gs(k)
       else
         j = 1
         do while (j .lt. k - 1 .and. tm(j+1) .le. t)
           j = j + 1
         enddo

         fieldof = 10.0d0**(log10 ( gs(j) ) + (t - tm(j)) / (tm(j+1) - tm(j)) * (log10 ( gs(j+1) ) - log10 ( gs(j) )))
       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function pctl

     &          ( n, a, q )
!__________________________________________________________________________________________________________________________________
!
! The percentile q (in percent) of the n sorted values a, interpolated linearly between them, as numpy.percentile does.
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       integer       n, i
       real*8        pctl, a(n), q, p
!__________________________________________________________________________________________________________________________________
!
       if (n .eq. 1) then
         pctl = a(1)
         return
       endif
       p = q / 100.0d0 * dble ( n - 1 )
       i = min ( int ( p ), n - 2 )

       pctl = a(i+1) + (p - dble ( i )) * (a(i+2) - a(i+1))

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function wmedf

     &          ( n, v, w )
!__________________________________________________________________________________________________________________________________
!
! The median of the n values v with the weights w, interpolated between the cumulative weights at the centers of the sorted
! values, as wmedian of hicorrect2.py does. The arrays are sorted in place.
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       integer       n, i, j
       real*8        wmedf, v(n), w(n), c(n), x, y, tot
!__________________________________________________________________________________________________________________________________
!
       do i=2,n
         x = v(i)
         y = w(i)
         j = i - 1
         do while (j .ge. 1)
           if (v(j) .le. x) exit
           v(j+1) = v(j)
           w(j+1) = w(j)
           j = j - 1
         enddo
         v(j+1) = x
         w(j+1) = y
       enddo
       tot = 0.0d0
       do i=1,n
         tot  = tot + w(i)
         c(i) = tot - 0.5d0 * w(i)
       enddo
       do i=1,n
         c(i) = c(i) / tot
       enddo
       if (0.5d0 .le. c(1)) then
         wmedf = v(1)
       elseif (0.5d0 .ge. c(n)) then
         wmedf = v(n)
       else
         j = 1
         do while (c(j+1) .lt. 0.5d0)
           j = j + 1
         enddo

         wmedf = v(j) + (0.5d0 - c(j)) / (c(j+1) - c(j)) * (v(j+1) - v(j))
       endif

       return
       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       subroutine hsort

     &            ( n, a )
!__________________________________________________________________________________________________________________________________
!
! Sort the n values a in increasing order (heapsort).
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       integer       n, i, ir, j, l
       real*8        a(n), ra
!__________________________________________________________________________________________________________________________________
!
       if (n .lt. 2) return
       l  = n / 2 + 1
       ir = n
       do
         if (l .gt. 1) then
           l  = l - 1
           ra = a(l)
         else
           ra    = a(ir)
           a(ir) = a(1)
           ir    = ir - 1
           if (ir .eq. 1) then
             a(1) = ra
             return
           endif
         endif
         i = l
         j = l + l
         do while (j .le. ir)
           if (j .lt. ir) then
             if (a(j) .lt. a(j+1)) j = j + 1
           endif
           if (ra .lt. a(j)) then
             a(i) = a(j)
             i = j
             j = j + j
           else
             j = ir + 1
           endif
         enddo
         a(i) = ra
       enddo

       end

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

       function istemp

     &          ( nx, ny, a )
!__________________________________________________________________________________________________________________________________
!
! True if the image a looks like a temperature image: more of its pixels lie between 1 and 10^4 (temperatures in K) than above 10^4
! (surface densities in cm^-2). Pixels below 1, such as the wild low values around the usable area of some images, are ignored.
!__________________________________________________________________________________________________________________________________
!
       implicit      none

       logical       istemp

       integer       nx, ny, i, j, nhigh, nlow
       real*8        a(nx,ny)
!__________________________________________________________________________________________________________________________________
!
       nhigh = 0
       nlow  = 0
       do j=1,ny
         do i=1,nx
           if (a(i,j) .ge. 1.0d4) then
             nhigh = nhigh + 1
           elseif (a(i,j) .ge. 1.0d0) then
             nlow = nlow + 1
           endif
         enddo
       enddo

       istemp = nlow .gt. nhigh

       return
       end
