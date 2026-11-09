A 3D RADMC-3D CIRRUS WHOSE PROJECTION IS THE CIRRUS MAP OF SIMSKY6

FILES
 make_cirrus3d.py     builds the density cube from bd.fits (times 0.1235), writes amr_grid.inp, dust_density.binp, cirrus3d.json
 compare_cirrus3d.py  compares the RADMC-3D cirrus with the analytic slab of simsky6 and tests the correction on both
../cylinder/userdef_module.f90   the module, now with cirrus3d = 1 (read the grid and density instead of building a model)

STEPS (in a new directory with dustopac.inp, the opacity, wavelength_micron.inp, external_source.inp, bd.fits,
make_simsky5.py and hicorrect2.py)
1. python3 make_cirrus3d.py
2. +radmc3d.inp as for the cylinders, with Rbesph = 0.0, no cylinder line, and:
     cirrus3d = 1
     rto_style = 3          (binary dust temperatures, much smaller)
3. radmc3d mctherm          (no images are needed)
4. python3 compare_cirrus3d.py

GRID SIZE (parameters at the top of make_cirrus3d.py)
 CELL = 24, MARGIN = 1.0, DEPTH = 2.0:   about 290 x 290 x 66   = 5.6e6 cells   (a first look)
 CELL = 12, MARGIN = 2.0, DEPTH = 2.0:   about 712 x 712 x 132  = 6.7e7 cells   (the full resolution of the plan)
make_cirrus3d.py prints the size and the memory of one array before writing anything. For the large grid use at least 5e9
photons; every column sums about 130 cells, which averages the Monte Carlo noise of the cells.

TESTED HERE on a small synthetic map with the public RADMC-3D: the projection of the cube equals the map to 1e-6, the module
reads the files, mctherm runs, and compare_cirrus3d.py produces its table and images.
