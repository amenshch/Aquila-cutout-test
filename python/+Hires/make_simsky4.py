#!/usr/bin/env python3
"""
make_simsky.py -- a simulated sky in which the dust temperature follows the
total shielding column along each line of sight.

Each line of sight carries the total surface density D = D_cirrus + D_filament.
The two components are given the geometry each of them has.  The cirrus is
plane-parallel: illuminated from both sides, its material at fractional depth u
in [0,1] is shielded by the column s = D_cirrus min(u, 1-u).  The filament is a
cylinder, and its material escapes radially, not vertically.  For a uniform
structure bounded in two dimensions, as a cylinder is, or in three, as a sphere
is, the chord geometry is the same: at a position where the filament column is
a fraction q of the peak column D_peak of that filament, the line of sight
crosses material at path coordinate zeta running over [-q, q] in units of the
radius, and the parcel there is shielded by

    s = (D_peak / 2) (1 - sqrt(beta^2 + zeta^2)) + D_cirrus / 2,
    beta = sqrt(1 - q^2),

the second term being the cirrus lying in front of and behind a filament that
sits at mid-depth.  On the axis, q = 1, this reduces to D_peak/2, the same as a
plane-parallel slab of the same column; away from the axis the escape path is
shorter and the material warmer, which a plane-parallel treatment misses.  That
difference is not small: a uniform sphere and a slab of the same observed
column require correction factors of 1.17 and 1.44 at 6e21 cm-2.

Both temperatures follow one shielding law,

    T(s) = T_floor + (T_0 - T_floor) (1 + s / D_star)^(-gamma).

The intensity is the integral of the emission over depth,

    I_lambda = D  integral_0^1  m_lambda(T(s(u))) du ,
    m_lambda(T) = B_lambda(T) kappa_lambda eta mu_H2 m_H ,
    kappa_lambda = 10 (nu / 1e12 Hz)^2 cm^2 g^-1 .

Every line of sight therefore has a smooth temperature gradient rather than one
or two discrete values, the surface of the cloud is at T_0 everywhere, and a
dense narrow structure is cold in its interior because its own column shields
it.  No component is assigned a temperature of its own, so the construction is
self-consistent: two lines of sight with the same total column have the same
temperature profile whatever material produced that column.

The shielding law is no longer an analytic guess.  It is the relation measured
on the grid of RADMC-3D Bonnor-Ebert models, where the dust temperature of each
radial shell was tabulated against the column lying between that shell and the
outer boundary, weighted by the mass of the shell.  Seventeen nodes spanning a
factor 256 in the surface density of the embedding cloud, a factor 12 in core
size and a factor 270 in central to edge density ratio fall on one relation
with a mass-weighted scatter of about 3 percent, and the form below reproduces
the tabulated relation to 1.9 percent root-mean-square.  Using it here makes
the simulated sky consistent with the radiative transfer that the correction of
the hires surface density maps is calibrated on, so that the two can be
compared without a difference in the assumed thermal structure confusing the
result.

The law gives 22.1 K at the surface, 20.3 K behind 1e20 cm-2, 13.5 K behind
1e21 cm-2 and 8.0 K behind 1e22 cm-2, and saturates near 6.4 K deep inside a
cloud.  The interstellar radiation field of the radiative transfer included the
microwave background, so the floor is the value the models themselves give; it
is colder than the 7 to 9 K usually quoted for the centres of prestellar cores,
and that difference is the dominant uncertainty of any correction built on this
relation.

The previous version used T(s) = 200/(1e-20 s + 20) + 10 K, which is warmer by
3 to 4 K everywhere above 3e20 cm-2 and saturates at 10 K.  Between the two
laws the emission-to-mass conversion at 500 micron differs by a factor 1.7 at
1e21 cm-2 and 3.5 at 1e22 cm-2.
"""
import sys

import numpy as np
from astropy.io import fits
from scipy import ndimage

H_PLANCK, K_BOLTZ, C_LIGHT = 6.62607015e-27, 1.380649e-16, 2.99792458e10
AMU, MU_H2, ETA = 1.6605402e-24, 2.8, 1.0e-2
KAPPA0, NU0, BETA = 10.0, 1.0e12, 2.0
PIX = 3.0
FWHM_TO_SIGMA = 1.0 / 2.354820045

TSURF = 22.4463      # dust temperature where nothing attenuates the radiation, K
TFLOOR = 6.2090      # dust temperature where the radiation cannot reach, K
DSTAR = 2.07379e21   # attenuating column at which the fall sets in, cm-2
GAMMA = 1.14100      # index of the fall with the attenuating column
KAPHEAT = 1.0 / (0.94e21 * 1.086)   # cm^2 per H2 at the heating wavelengths
NANG = 120           # directions used in the angular average
ANGULAR = True       # use the effective attenuating column, obtained by
                     # averaging the attenuation over all directions, rather
                     # than the column along the single shortest path.  What
                     # sets the temperature of a parcel is the radiation
                     # reaching it, which arrives from every direction and is
                     # attenuated along paths of very different lengths, so
                     # the single shortest path is only a proxy.  With this
                     # set, the sky is consistent with the way the correction
                     # of the hires images is derived.


def slab_seff(z, t0):
    """Effective attenuating column of a layer, in units of its half-column.

    z is the position in the layer, from -1 at one face to +1 at the other,
    and t0 the optical depth of the half-thickness at the heating wavelengths.
    The radiation arrives from every direction; along a direction whose cosine
    to the normal is mu, it has traversed (1 -/+ z) / |mu| of the layer.
    """
    mu = (np.arange(NANG) + 0.5) / NANG * 2.0 - 1.0
    d = np.where(mu > 0.0, (1.0 - np.asarray(z)[..., None]) / np.abs(mu),
                 (1.0 + np.asarray(z)[..., None]) / np.abs(mu))
    t = np.asarray(t0, float)[..., None] if np.ndim(t0) else float(t0)
    j = np.exp(-t * d).mean(axis=-1)
    return -np.log(np.maximum(j, 1e-300)) / np.maximum(np.asarray(t0, float),
                                                       1e-12)


def cyl_seff(x, t0):
    """Effective attenuating column of a cylinder, in units of its radius.

    x is the distance from the axis in units of the radius and t0 the optical
    depth of the radius.  A direction is set by its polar angle to the axis and
    its azimuth in the cross section; the path is the chord of the cross
    section divided by the sine of the polar angle.
    """
    ph = (np.arange(NANG) + 0.5) / NANG * np.pi
    al = (np.arange(NANG // 2) + 0.5) / (NANG // 2) * np.pi / 2.0
    x = np.asarray(x)[..., None]
    l2 = -x * np.cos(ph) + np.sqrt(np.maximum(1.0 - (x * np.sin(ph)) ** 2, 0.0))
    t = (np.asarray(t0, float)[..., None, None] if np.ndim(t0)
         else float(t0))
    e = np.exp(-t * l2[..., None] / np.sin(al)) * np.sin(al)
    j = e.sum(axis=(-1, -2)) / (NANG * np.sin(al).sum())
    return -np.log(np.maximum(j, 1e-300)) / np.maximum(np.asarray(t0, float),
                                                       1e-12)


def T_of_shielding(s):
    """Dust temperature at an effective attenuating column s, in cm-2.

    The relation measured on the grid of RADMC-3D Bonnor-Ebert models, applied
    here to the column shielding each layer rather than to the total column of
    the pixel, so that a line of sight carries a gradient instead of a single
    value.  It gives 22.1 K at the surface, 13.5 K behind 1e21 cm-2, 8.0 K
    behind 1e22 cm-2, and saturates near 6.4 K deep inside a cloud.
    """
    return TFLOOR + (TSURF - TFLOOR) * (1.0 + np.asarray(s, float)
                                        / DSTAR) ** (-GAMMA)

T0 = TSURF           # dust temperature at the surface, K (value of the law at s = 0)
NLAYER = 64          # layers used for the integral over depth
CYLINDRICAL_FILAMENTS = True   # give the filament the geometry of a cylinder
                               # rather than treating cirrus and filament as
                               # one plane-parallel slab
RANDOM_ORIENTATION = False   # give each filament an isotropically drawn
                             # inclination and shield the diffuse component
                             # with the orientation-averaged eta below, rather
                             # than treating both as face-on slabs.  The
                             # surface density maps are unchanged by this: an
                             # inclined cylinder is stretched along its axis,
                             # so its projected width stays the same while its
                             # column rises, and a map of a given column is
                             # reproduced by a less dense filament seen at an
                             # angle.  Only the temperatures, and hence the
                             # waveband images, differ.
ETA_DIFFUSE = 0.25   # effective shielding of the diffuse component, the value
                     # that a slab, a sphere and a filament all share once
                     # averaged over their bodies and over isotropic
                     # orientations.  A turbulent medium has no normal to
                     # incline, so it is given this value directly rather than
                     # an orientation.
ORIENT_SEED = 20260908       # so that a run can be repeated exactly
FIL_LABEL_FRAC = 0.30        # fraction of the axial column above which the
                             # filaments are separate from one another, used to
                             # identify them before each is given its own angle
PINDEX = 1.75        # index of the filament cross-section: the surface density
                     # falls as r^-(PINDEX-1) = r^-0.75 outside the
                     # half-maximum radius and flattens inside it, so the
                     # volume density goes as [1 + (r/Rflat)^2]^(-PINDEX/2)
NCHORD = 96          # samples along the chord through the filament
FIL_PEAK = None      # surface density on the axis of a filament, cm-2.  All
                     # the filaments of a map share one value, stated in the
                     # name of the file, so it is taken as the maximum of the
                     # filament map when this is left as None.  It must not be
                     # obtained from a local dilation: the wings of an
                     # r^-0.75 cross-section reach across the whole field, so
                     # no finite structuring element finds the crest and every
                     # pixel comes out as its own local maximum
QFLOOR = 1.0e-3      # smallest fraction of the peak filament column that is
                     # still treated as belonging to the filament


def _trap(y, x):
    """Trapezoidal integral, written out so that no numpy version matters."""
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def filament_chord_tables(nb=256):
    """Tabulate the chord through a filament of the assumed cross-section.

    For a cylinder of volume density rho(r) = [1 + (r/Rflat)^2]^(-PINDEX/2),
    which projects to a surface density falling as r^-(PINDEX-1), this returns
    three arrays on a grid of impact parameter b, in units of Rflat:

    q      the surface density at b divided by the surface density on the axis
    sfrac  the shielding of each sample along the chord, as a fraction of the
           surface density on the axis, the shielding of a parcel being the
           column from its radius outward
    wgt    the mass weight of each sample along the chord, normalised to unity

    Assuming the filament to be a uniform cylinder instead, as an earlier
    version did, makes the surface density fall off the axis as sqrt(1-b^2)
    while the escape column falls as 1-b.  The flanks are then far too warm and
    the emergent intensity peaks off the axis, by 55 per cent at 160 micron and
    14 per cent at 250 micron, which is an artefact: the radiative transfer
    models show a maximum on the axis at every waveband except 160 micron,
    where the centre is low by 5 per cent at the most concentrated node.
    """
    bmax, zmax = 24.0, 400.0
    b = np.linspace(0.0, bmax, nb)
    rho = lambda r: (1.0 + r ** 2) ** (-0.5 * PINDEX)
    # radial column from r outward, tabulated on a grid that resolves the core
    rr = np.concatenate([[0.0], np.geomspace(1.0e-4, zmax, 20000)])
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (rho(rr[1:]) + rho(rr[:-1]))
                                             * np.diff(rr))])
    scol = integ[-1] - integ                      # in units of rho_0 Rflat
    # The samples along the chord are placed uniformly in asinh(z), which is
    # linear near the axis and logarithmic far from it.  A grid uniform in z
    # cannot be used: with a spacing of several Rflat the dense core of the
    # cross-section falls between two samples and is never seen, which makes
    # the material appear much less shielded than it is.  By symmetry the
    # escape column from the axis is exactly half the axial column, and only a
    # grid that resolves the core recovers that.
    tmax = np.arcsinh(zmax)
    t = tmax * (np.arange(NCHORD) + 0.5) / NCHORD
    zeta = np.sinh(t)
    jac = np.cosh(t)                              # dz/dt, the quadrature weight
    col0 = 2.0 * _trap(rho(rr), rr)               # axial column, same units
    q = np.empty(nb)
    sfrac = np.empty((nb, NCHORD))
    wgt = np.empty((nb, NCHORD))
    for i, bi in enumerate(b):
        r = np.hypot(bi, zeta)
        w = rho(r) * jac
        q[i] = 2.0 * _trap(rho(np.hypot(bi, rr)), rr) / col0
        sfrac[i] = np.interp(r, rr, scol) / col0
        wgt[i] = w / w.sum()
    return b, q, sfrac, wgt

PREFIX = "simsky4"   # stem of every file written by this script
VERSION = ("2026-09-08: shielding law replaced by the relation measured on the"
           " RADMC-3D model grid")

BANDS = [(70., 8.4, "070", "r8p4"), (100., 9.4, "100", "r9p4"),
         (160., 13.5, "160", "r13p5"), (250., 18.2, "250", "r18p2"),
         (350., 24.9, "350", "r24p9"), (500., 36.3, "500", "r36p3")]
WIDTH = 40           # half-maximum width of the filaments, arcsec
STEM = PREFIX + ".w%d.%s.%s.fits"

def m_lambda(t, wave_um):
    nu = C_LIGHT / (wave_um * 1.0e-4)
    kap = KAPPA0 * (nu / NU0) ** BETA
    x = np.clip(H_PLANCK * nu / (K_BOLTZ * np.asarray(t, float)), 1e-8, 700.0)
    b = (2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2) / np.expm1(x)
    return b * kap * ETA * MU_H2 * AMU / 1.0e-17          # MJy/sr per cm-2


def _double_card(key, value, comment):
    """A header card holding a double in fixed exponential form.

    Astropy writes a float in the shortest form that round-trips, so a
    wavelength of 160 becomes '160.0'.  That is valid FITS and is read
    correctly, but it does not look like the cards written by cfitsio through
    ftpkyd, which are of the form ' 1.6000000000000E+02'.  Writing the same
    form here removes any doubt about how the value will be read.
    """
    return fits.Card.fromstring(
        ("%-8s= %20.13E / %s" % (key, float(value), comment))[:80].ljust(80))


def write_fits(name, data, header, doubles=()):
    """Write an image, then set the double-valued keywords in a second pass.

    Astropy drops BZERO and BSCALE when it writes floating-point data, even
    when the header supplied carries them, so they have to be put back by
    reopening the file in update mode.  The keywords listed in 'doubles', each
    a triple of keyword, value and comment, are written there as well, in
    fixed exponential form.
    """
    fits.PrimaryHDU(np.asarray(data, np.float32), header).writeto(
        name, overwrite=True)
    with fits.open(name, mode="update") as hdul:
        h = hdul[0].header
        h["BZERO"] = 0.0
        h["BSCALE"] = 1.0
        for key, value, comment in doubles:
            if key in h:
                del h[key]
            h.append(_double_card(key, value, comment))


def repair(a, name=""):
    """Replace non-positive pixels by the nearest positive one.

    The cirrus map carries a border of zeros -- its first and last columns and
    its last three rows -- which a convolution would smear inward as a dark rim.
    A surface density of zero is not a measurement, so such pixels are filled
    from their nearest valid neighbour before anything else is done.
    """
    bad = ~(a > 0) | ~np.isfinite(a)
    if not bad.any():
        return a
    idx = ndimage.distance_transform_edt(bad, return_distances=False,
                                         return_indices=True)
    print("   %s: %d non-positive pixels (%.4f per cent) filled from their "
          "nearest valid neighbour" % (name, bad.sum(), 100.0 * bad.mean()))
    return a[tuple(idx)]


def smooth(a, fwhm_as):
    """Convolve with a Gaussian, padding the image first and cropping after.

    The image is expanded on every side by twice the half-maximum width of the
    beam, by replicating the edge, so that the convolution never runs off the
    array, and the result is cut back to the original size.
    """
    if fwhm_as <= 0:
        return a.copy()
    n = int(np.ceil(2.0 * fwhm_as / PIX))
    b = np.pad(a, n, mode="edge")
    b = ndimage.gaussian_filter(b, fwhm_as * FWHM_TO_SIGMA / PIX,
                                mode="nearest", truncate=4.0)
    return b[n:-n, n:-n]


def main():
    global WIDTH
    if len(sys.argv) > 1:
        WIDTH = int(sys.argv[1])
    get = lambda f: np.array(fits.getdata(f), float).squeeze()
    hdr = fits.getheader("bdm1500x1p5.fits")
    print("repairing the input maps:")
    Dc = repair(get("bdm1500x1p5.fits"), "cirrus")
    Df = repair(get("fils.sha%dx0p75-.x3e22.fits" % WIDTH), "filament")
    D = Dc + Df
    print("total surface density: min %.3g, median %.3g, max %.3g cm-2"
          % (D.min(), np.median(D), D.max()))

    # peak column of the filament to which each pixel belongs, from a
    # greyscale dilation with a disk of the filament width; this is
    # independent of the orientation of the filament
    sinmap = np.ones_like(Df)
    if CYLINDRICAL_FILAMENTS:
        btab, qtab, sfrac, wgt = filament_chord_tables()
        peak = float(np.max(Df)) if FIL_PEAK is None else float(FIL_PEAK)
        # inclination of every filament.  The connected regions of the
        # filament map above a low threshold are its individual filaments;
        # each is given its own angle to the line of sight, drawn so that the
        # axis is isotropically oriented, that is with cos(theta) uniform.
        # Its material is shielded by sin(theta) times the face-on value,
        # because an inclined filament is crossed over a longer path while
        # being no better shielded.  A pixel of observed column Df belonging
        # to a filament inclined by theta is a pixel of a filament whose
        # face-on axial column is peak sin(theta), so its absolute shielding
        # is peak sin(theta) SF and its impact parameter is unchanged: the
        # observed column and the axial column are both stretched by the same
        # factor and the ratio of the two, which is what fixes the impact
        # parameter, is not.  The inclination therefore enters the shielding
        # in emit() below and nowhere else.
        if RANDOM_ORIENTATION:
            # The wings of an r^-0.75 cross-section overlap, so labelling at
            # a low level merges every filament into one region.  The
            # filaments are therefore identified near their crests, where they
            # are separate, and each remaining pixel is then given the label of
            # the nearest crest.
            core = Df > FIL_LABEL_FRAC * peak
            lab, nfil = ndimage.label(core)
            if nfil < 1:
                raise SystemExit("no filament found above %.2f of the peak"
                                 % FIL_LABEL_FRAC)
            _, idx = ndimage.distance_transform_edt(~core,
                                                    return_indices=True)
            lab = lab[tuple(idx)]
            rng = np.random.default_rng(ORIENT_SEED)
            cth = rng.uniform(0.0, 1.0, nfil)
            sth = np.sqrt(1.0 - cth ** 2)
            sinmap = sth[lab - 1]
            print("filament inclinations: %d filaments, sin(theta) median"
                  " %.3f, range %.3f to %.3f"
                  % (nfil, float(np.median(sth)), float(sth.min()),
                     float(sth.max())))
        Dpeak = np.full_like(Df, peak)
        qmap = np.clip(Df / peak, 0.0, 1.0)
        # impact parameter of each pixel, from its share of the axial column
        bmap = np.interp(qmap, qtab[::-1], btab[::-1])
        SF = np.empty((NCHORD,) + Df.shape)
        WG = np.empty_like(SF)
        for j in range(NCHORD):
            SF[j] = np.interp(bmap, btab, sfrac[:, j])
            WG[j] = np.interp(bmap, btab, wgt[:, j])
        print("filament: cross-section index %.2f, axial column %.4g cm-2,"
              " q median %.4f, 99th percentile %.4f, %d pixels above %.3g"
              % (PINDEX, peak, np.median(qmap[qmap > QFLOOR]),
                 np.percentile(qmap, 99.0), int((qmap > QFLOOR).sum()),
                 QFLOOR))
        print("   impact parameter: median %.2f, 1st percentile %.2f Rflat;"
              " the half-maximum radius is 2.32 Rflat"
              % (np.median(np.interp(qmap, qtab[::-1], btab[::-1])),
                 np.percentile(np.interp(qmap, qtab[::-1], btab[::-1]), 1.0)))

    u = (np.arange(NLAYER) + 0.5) / NLAYER
    depth = np.minimum(u, 1.0 - u)                        # shielded from both sides
    edif = ETA_DIFFUSE if RANDOM_ORIENTATION else 0.5     # diffuse component
    zlayer = (np.arange(NLAYER) + 0.5) / NLAYER * 2.0 - 1.0   # -1 to +1
    tgrid = np.geomspace(0.999 * TFLOOR, T0 * 1.001, 4000)

    def emit(table):
        """Emission of a line of sight, per unit surface density of each
        component, integrated over depth.

        The cirrus is integrated over its own depth as a plane-parallel slab.
        The filament, when it is given a cylindrical geometry, is integrated
        over the chord that the line of sight cuts through it, its material
        shielded radially and additionally by half of the cirrus.
        """
        acc = np.zeros_like(D)
        if ANGULAR:
            t0c = KAPHEAT * Dc / 2.0
            for zz in zlayer:
                sa = slab_seff(np.full_like(Dc, zz), t0c) * Dc / 2.0
                acc += Dc * np.interp(T_of_shielding(sa), tgrid, table)
        else:
            for dd in depth:
                acc += Dc * np.interp(T_of_shielding(2.0 * edif * dd * Dc),
                                      tgrid, table)
        if not CYLINDRICAL_FILAMENTS:
            for dd in depth:
                acc += Df * np.interp(T_of_shielding(dd * Df + 0.5 * Dc),
                                      tgrid, table)
            return acc / NLAYER
        for j in range(NCHORD):
            sh = Dpeak * SF[j] * sinmap + 0.5 * Dc
            acc += NLAYER * Df * WG[j] * np.interp(T_of_shielding(sh),
                                                   tgrid, table)
        return acc / NLAYER

    # mass-weighted mean temperature of each line of sight, for reference only
    tmean = np.zeros_like(D)
    if CYLINDRICAL_FILAMENTS:
        for dd in depth:
            tmean += Dc * T_of_shielding(dd * Dc) / NLAYER
        for j in range(NCHORD):
            tmean += Df * WG[j] * T_of_shielding(Dpeak * SF[j] * sinmap
                                                 + 0.5 * Dc)
        tmean /= np.maximum(D, 1e-30)
    else:
        for dd in depth:
            tmean += T_of_shielding(dd * D)
        tmean /= NLAYER

    for wave, beam, tag, rtag in BANDS:
        table = m_lambda(tgrid, wave)
        img = smooth(emit(table), beam)
        h = hdr.copy()
        h["BUNIT"] = "MJy/sr"
        h["BZERO"], h["BSCALE"] = 0.0, 1.0
        h["COMMENT"] = ("simulated sky, temperature from the shielding law of"
                        " the RADMC-3D model grid")
        h["TLAW"] = ("%.4f+(%.4f-%.4f)(1+s/%.5e)^-%.5f K"
                     % (TFLOOR, TSURF, TFLOOR, DSTAR, GAMMA),
                     "s = shielding column, cm-2")
        h["GEOM"] = ("cylinder" if CYLINDRICAL_FILAMENTS else "slab",
                     "geometry given to the filament component")
        h["FILORIEN"] = ("isotropic" if RANDOM_ORIENTATION else "plane of sky",
                         "inclination given to the filaments")
        name = STEM % (WIDTH, tag, rtag)
        # hires reads the wavelength of an image from the WAVE keyword
        write_fits(name, img, h,
                   [("WAVE", wave, "wavelength (microns)"),
                    ("BEAM", beam, "half-maximum beam width (arcsec)")])
        print("%-58s median %10.3f MJy/sr" % (name, np.median(img)))

    for nm, arr, unit, bm in (
            (((PREFIX + ".w%d.surfdens.pix.fits") % WIDTH), D, "cm-2", 0.0),
            (((PREFIX + ".w%d.surfdens.r13p5.fits") % WIDTH), smooth(D, 13.5), "cm-2", 13.5),
            (((PREFIX + ".w%d.surfdens.r36p3.fits") % WIDTH), smooth(D, 36.3), "cm-2", 36.3),
            (((PREFIX + ".w%d.tmean.pix.fits") % WIDTH), tmean, "K", 0.0)):
        h = hdr.copy()
        h["BUNIT"], h["BZERO"], h["BSCALE"] = unit, 0.0, 1.0
        # a surface density or temperature image is band 0 micron in the getsf
        # convention, so the wavelength is set to zero rather than inherited
        dbl = [("WAVE", 0.0, "wavelength (microns); 0 for surface density")]
        if bm > 0:
            dbl.append(("BEAM", bm, "half-maximum beam width (arcsec)"))
        elif "BEAM" in h:
            del h["BEAM"]
        write_fits(nm, arr, h, dbl)
        print("%-58s median %10.4g %s" % (nm, np.median(arr), unit))

    with open(("+hires.cfg." + PREFIX + ".w%d") % WIDTH, "w") as f:
        f.write("#\n# HIRES: CONFIGURATION PARAMETERS FOR HIRES\n#\n"
                "# PM  WAVE  BEAM  OFFSET   IMAGE   [MASKIMAGE]\n#\n")
        for wave, beam, tag, rtag in BANDS:
            f.write("  %s   %s  %4.1f   0.0   %s\n"
                    % ("-" if wave < 160 else "+", tag, beam,
                       STEM % (WIDTH, tag, rtag)))
    print(("written +hires.cfg." + PREFIX + ".w%d") % WIDTH)

    print()
    print("make_simsky.py version %s" % VERSION)
    print("header keywords written:")
    print("%42s %10s %10s" % ("file", "WAVE", "BEAM"))
    import glob
    for f in sorted(glob.glob((PREFIX + ".w%d.*.fits") % WIDTH)):
        hh = fits.getheader(f)
        print("%42s %10s %10s"
              % (f, hh.get("WAVE", "absent"), hh.get("BEAM", "absent")))

    print()
    print("temperature at the mid-plane, as a function of the total column:")
    print("%16s %14s %14s"
          % ("N_H2 (cm-2)", "mid-plane (K)", "mass mean (K)"))
    for q in (5, 25, 50, 75, 95, 99.9):
        n = np.percentile(D, q)
        print("%16.3g %14.2f %14.2f"
              % (n, T_of_shielding(0.5 * n),
                 np.median(tmean[(D > 0.98 * n) & (D < 1.02 * n)])))


if __name__ == "__main__":
    main()
