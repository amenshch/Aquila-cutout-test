RADMC-3D model of an isolated cylinder with the profile of the simsky4/5 filaments

FILES

userdef_module.f90                   Your module with a cylinder mode added. With cylinder = 0 (the default) it behaves exactly as
                                     before; nothing of the spherical models was changed.
radmc3d_cylinder_parameters.inp      The parameters to append to +radmc3d.inp (or radmc3d.inp) for the cylinder.
cylinder_smooth_temperature.py       Averages the Monte Carlo temperatures of the middle of the cylinder along its axis.

WHAT CHANGED IN THE MODULE

userdef_prep_model    With cylinder = 1, calls cylinder_grid, which writes a 2D grid in spherical coordinates (r, theta) with
                      theta from 0 to pi/2, so that RADMC-3D mirrors the model in the equatorial plane. The radial walls are
                      logarithmic from 1e-3 of the flat radius (the innermost cell spans out to 0.2 of the flat radius), the
                      theta walls logarithmic near the axis and linear near the equator.
userdef_setup_model   With cylinder = 1, calls cylinder_density, which writes dust_density.inp (each cell the mean over 8 x 8
                      samples equal in volume), aperture_info.inp and the file '+mc3d.out.cyl' with the true properties, then
                      continues with your code for the radiation field. The spherical optical-depth table is skipped.
userdef_dostuff       With cylinder = 1, writes mc.denstemps.tab cell by cell (R, z, densities, temperature) and
                      mc.cyl.midplane.tab (the row next to the equator), and skips the spherical optical-depth profiles.
                      The images and FITS files are written by your unchanged code.

CHECKS DONE HERE

The module compiles with the stock RADMC-3D 2.0 (with -cpp -ffree-line-length-none, a stub for your do_raytrace_tau flag, and
stubs for your library routines osystem, when, fdelete, printerr and selectk) and runs 'mctherm'. The surface density across the
axis, integrated through dust_density.inp, is 2.998e22 cm^-2 on the axis with a half-maximum width of 40.06 arcsec, against 3e22
and 40 requested; along the axis it varies by about 1 per cent. The flat radius comes out as 9.22 arcsec, the central density of
H2 as 2.26e5 cm^-3 and the mass per unit length as 82.3 Msun/pc.

The image step could not be tested here: the stock RADMC-3D refuses to write images when stdo is redirected, which your modified
version allows. It is your unchanged code.

With 3e6 photons, a test run gave the temperature of the middle part averaged along the axis: 6.9 K on the axis, 8.3 K at 20",
10.2 K at 56", 13.5 K at 158" and 18.9 K at 361". Near the axis most cells received no photon packet at that count.

RUNNING IT

1. Compile RADMC-3D with this userdef_module.f90 as usual.
2. In a new directory with your dust, wavelength and external-source files, append radmc3d_cylinder_parameters.inp to
   +radmc3d.inp, and set Rbesph = 0.0 there.
3. radmc3d mctherm
   The spherical grid has small cells near the origin and near the axis, and at 1e9 photons some of them may still be sampled
   poorly. Step 4 removes that noise in the middle part.
4. python3 cylinder_smooth_temperature.py
   It replaces the temperature of every cell of the middle part (|z| < 0.7 L/2 - R_out, 2.03 pc here) by the mass-weighted mean
   of its bin of R, keeps the original as dust_temperature.mc.dat, and writes cylinder_temperature_profile.txt. The column
   "frac T=0" shows how many cells of each bin received no packet; it should be small beyond the central few arcsec.
5. Images edge-on, with the axis vertical, 1793 x 1793 pixels of 3" at 260 pc (5379" = 1398540 AU on a side):
     radmc3d image loadlambda incl 90 npix 1793 sizeau 1398540
   and the surface density image as you make it for the spheres (lambda 999), also with incl 90 npix 1793 sizeau 1398540.
   The cylinder then runs from edge to edge of the images, as the simsky4 filaments do.
6. Your usual convolution, hires and hicorrect2 steps. The middle rows of the images (away from the ends) are the ones to
   compare with the true surface density.

SUGGESTED FIRST MODEL

cylsd = 3e22, cylfwhm = 40, cylpind = 1.75, cylrout = 0.5, cyllen = 7.234, as in radmc3d_cylinder_parameters.inp. The length
equals the side of the simsky4 images; the middle part used for the comparison is then more than 1.5 pc from both ends.

EMBEDDED FILAMENTS FOR SIMSKY6 (envelope and crater background)

Two parameters embed the cylinder in a uniform cylindrical envelope that stands for the cirrus around it:

 cylembsd = 3.5e21       ! Surface density of the envelope through the axis (cm^-2), the filament excluded; 0 for none
 cylrenv = 1.0           ! Outer radius of the envelope (pc)

Set cylembsd to the surface density of the cirrus along the axis of the filament, as printed by
python3 make_simsky6.py --cirrus. The steps for each model are then:

1. radmc3d mctherm
2. python3 cylinder_smooth_temperature.py
   (with an envelope, the middle part is |z| < 0.7 L/2 - cylrenv, 1.5 pc for the defaults)
3. radmc3d image loadlambda incl 90 npix 1793 sizeau 1398540
   and the surface density image (lambda 999) with the same options: the full model.
4. The same two image runs with the argument doback. The volume of the filament is emptied, the envelope and all temperatures
   are kept, and the images are written as mc.<wave>um.bgl.fits and mc.surfdens.bgl.fits, exactly as for the spherical models.

make_simsky6.py subtracts the .bgl images from the full ones, which leaves the emission of the filament inside the cirrus, with
its crater background removed.

ENVELOPES THAT FOLLOW THE CIRRUS ALONG EACH FILAMENT (cylenvz = 1)

1. python3 make_simsky6.py --envelopes
   writes cylenv.inp into each model directory: the cirrus beneath the filament, row by row, against the height z (pc) along the
   cylinder. It prints the median for each model.
2. In +radmc3d.inp of each model:
     cylenvz = 1             ! envelope follows cylenv.inp along the axis; the grid then covers both hemispheres
     cylembsd = <median>     ! the median printed in step 1 (a reference value; the table sets the actual envelope)
     cylrenv = 1.0           ! envelope outer radius (pc)
     cyllen = 10.3           ! longer than the sky (7.234 pc), so that the disturbed ends fall outside the images
     nthcyl = 240            ! cells per hemisphere; the grid has twice as many
3. radmc3d mctherm; python3 cylinder_smooth_temperature.py (it now averages within windows of 0.2 pc along the axis);
   the full and the doback image runs with npix 1793 sizeau 1398540, as before.
4. python3 make_simsky6.py with ENVZ = True: the full images are placed, shifted only sideways, and each one is flipped upside down
   if its envelope runs the other way than the cirrus.

The temperatures now vary along the axis, so the Monte Carlo noise is no longer averaged over the whole middle part. Use at least
1e9 photons, more if the temperature profiles look noisy.

IMAGES OF A MODEL WITHOUT MIRROR SYMMETRY (cylenvz = 1)

With the full theta grid, a ray that runs exactly in the equatorial plane stops RADMC-3D with STOP 3331 (a consistency check of
its spherical coordinates). With incl 90 and an odd npix, the middle row of the image is exactly such a ray. Shift the image
by a negligible 0.01 AU in the vertical direction with zoomau instead of sizeau, in all image runs (intensities and surface
density, full and doback):

  radmc3d image loadlambda incl 90 npix 1793 zoomau -699270 699270 -699269.99 699270.01
  radmc3d image tracecolumn lambda 999 incl 90 npix 1793 zoomau -699270 699270 -699269.99 699270.01
  (and the same two with doback)

GRID MODE: FULL THETA GRID BY DEFAULT (cylfull = 1)

A test sphere with the profile of the filaments (cyltest = 1 or 2, see below), computed with 3e6 photons, gave these mass-weighted
dust temperatures (K):

   r (arcsec)       5      20      40      80     150     300
   1D grid        8.28    9.16   10.43   12.23   14.50   18.66
   2D mirrored    7.47    8.60    9.89   11.64   13.79   17.77
   2D full        7.46    9.00   10.37   12.22   14.50   18.66

The full theta grid agrees with the 1D grid; the equatorial mirror mode is 0.5 to 0.9 K too cold. The module therefore uses the
full theta grid by default; cylfull = 0 restores the mirror mode, for comparison with earlier runs only.

To repeat the test with your RADMC-3D build: three runs with cylsd = 3e22, cylfwhm = 40, cylpind = 1.75, cylrout = 0.5,
cyllen = 1.0, cylembsd = 0, and
   (a) cyltest = 2                      1D grid
   (b) cyltest = 1, cylfull = 0         2D mirrored
   (c) cyltest = 1, cylfull = 1         2D full
then compare their mc.denstemps.tab.
