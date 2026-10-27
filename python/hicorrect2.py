"""hicorrect2.py -- correct the single-temperature bias of surface density
images.

The correction applies to any surface density image obtained by fitting one
modified blackbody to the intensities of each pixel, whether made by hires or
by any other method, at any resolution. It needs only the image and the list of
wavebands that were fitted to make it.

A surface density image is obtained by fitting one modified blackbody to the
intensities of each pixel. A line of sight crossing a dense structure is not
isothermal, the fit is weighted toward the warmer material, and the surface
density comes out too low. The deficit exceeds a factor of two in the densest
material.

In the reconstruction of hires the whole integrated deficit sits in the base
image formed at the coarsest resolution, because the increments that carry the
finer resolutions are differential terms that very nearly cancel when
integrated. The cumulative images satisfy

    S(theta) = S(theta_coarse) + sum of the increments finer than theta,

so every image of the reconstruction has nearly the same fractional deficit.

How the correction is applied. One surface density image is corrected per run.
The factor f is evaluated on that image, pixel by pixel, and the image is
multiplied by it. This works for an image of any resolution; the default of the
paper is the hires image at 13.5 arcsec. If the temperature image fitted
together with it is given, the dust temperature implied by the corrected
surface density is also written. That temperature image must have the same
resolution and the same pixels as the surface density image, because the
temperature of each pixel is recomputed with the factor of that pixel.

The correction. For a line of sight of true surface density Sigma, a uniform
structure of that surface density is taken. At every depth the effective
attenuating column Sigma_a is the angular average of the column toward the
outside, for a uniform sphere along its diameter and for a uniform layer along
its normal, and the two are averaged. The dust temperature follows the relation
measured on RADMC-3D models of Bonnor-Ebert spheres embedded in uniform clouds
and of isolated cylinders,

    T(Sigma_a) = T_floor + (T_0 - T_floor) (1 + Sigma_a / S_star)^(-gamma).

Integrating the emission over depth, fitting one modified blackbody to the
result in the same wavebands that were fitted to make the image, and comparing with Sigma gives
the factor by which the fitted surface density must be multiplied. The factor
is tabulated once against the fitted surface density and interpolated onto the
image.

At low surface densities the analytic layer is replaced by a RADMC-3D model of
a diffuse cloud: below a fitted surface density of 10^21.6 cm^-2 the factor is
the average of the analytic sphere and the factor measured on that cloud, and
the difference from the analytic relation is reduced smoothly to zero at
10^22 cm^-2. The surface layers of an extended cloud are heated by unattenuated
ultraviolet light and are warmer than the relation predicts, which a gray
attenuation cannot describe.

Usage:

    python3 hicorrect2.py <image> [-temp=<temperature image>] [-bands=<list>]
                          [-weights=<list>] [-g0=<field>]

<image> is the surface density image to be corrected. The file names of the
images may be given without the extension '.fits', which is then added. The corrected image is
written as <image>.correct.fits and the image of the correction factors as
<image>.cfactor.fits.

-temp=<temperature image> gives the dust temperature image fitted together with
<image>. It must have the same resolution and the same pixels as <image>. The
dust temperature implied by the corrected surface density is then written as
<temperature image>.correct.fits.

-bands=<list> gives the wavelengths in micron of the wavebands that were fitted
to make the image, separated by commas, for example -bands=250,350,500; the
default is 160,250,350,500. The relation is built for these wavebands, and the
longest of them is the reference waveband.

-weights=<list> gives the relative calibration uncertainties of the same
wavebands, in the same order, which set their weights in the fit that made the
image, for example -weights=0.1,0.1,0.1 for equal weights. Only their ratios
matter. The default is 0.2 at 160 um and shorter and 0.1 at longer wavelengths,
the weights of hires. The analytic relation is built for any wavebands and
weights. The calibration on a diffuse cloud is provided for five sets: 160-500
um with the default or equal weights, 250-500 um with equal weights, 100-500 and
70-500 um with the default weights. For other sets, the calibration of the
standard set is used, and a warning is printed. The field factor of -g0= was
measured with the standard set only; for other sets it is approximate.

The outputs are written to the current directory, whatever directory the inputs
come from. Their names are the same as those of the Fortran version.

The relation is calibrated for the radiation field G_0 = 10 in Habing units.
For a cloud in another field, give its field with -g0=<field>, as estimated
from the diffuse material with g0_requirement.py. The factor is then multiplied
by a function of the fitted surface density and the field, measured on the
grid of radiation fields from G_0 = 0.98 to 482.5 and held constant outside it.

The correction must be applied once. It is not idempotent: applied again to an
image it has already corrected, it would raise the values a second time. Every
output carries CREATOR = 'HICORRECT' and the keyword HICORR, and any input that
carries either of them is refused, including the outputs of the Fortran version.
The outputs also record the radiation field of the correction (HIG0), the
wavebands fitted (HIBANDS) with their relative calibration uncertainties
(HIWGTS), and the four constants of the temperature relation (HITMIN, HIT0,
HISSTAR, HIGAMMA), as the Fortran version does.

Surface densities in H2 cm^-2, intensities in MJy/sr, temperatures in K.
"""

import os
import sys
import numpy as np
from astropy.io import fits

# -------------------------------- temperature against the attenuating column
# The relation T(Sigma_a), fitted by fit_attenuation_law.py 
# 2093 RADMC-3D nodes, 100693 shells, mass weighted.
# 100378 sphere shells and 315 cylinder radii; the cylinders carry 0.50 of the total weight
# outermost shell of every sphere rejected and to seven cylinders computed on the full theta grid
TFLOOR = 7.049531     # dust temperature where the radiation cannot reach, K
TSURF = 22.866783     # dust temperature where nothing attenuates it, K
SSTAR = 2.389668e+21   # attenuating column at which the fall sets in, cm^-2
GAMMA = 1.596874     # index of the fall with the attenuating column

# ------------------------------------ opacity at the heating wavelengths
KAPHEAT = 1.0 / (0.94e21 * 1.086)    # cm^2 per H2, from A_V = N/0.94e21
NANG = 160           # directions used in the angular average

# ----------------------------------------------------------------- the dust
KAPPA0 = 10.0        # cm^2 per gram of dust at NU0
NU0 = 1.0e12         # Hz
BETA = 2.0
DUST_TO_GAS = 1.0e-2
MU_H2 = 2.8

AMU = 1.66053906660e-24
H_PLANCK = 6.62607015e-27
K_BOLTZ = 1.380649e-16
C_LIGHT = 2.99792458e10

NDEPTH = 512         # layers used for the integral over depth

# The factor that diffuse material requires, measured on a RADMC-3D model of a diffuse cloud: a 3D cube of 7.2 x 7.2 x 2 pc in cells
# of 12 arcsec whose projection is the cloud map of the simulated sky, without its enhancement (make_cirrus3d.py). Each column
# of the cube was fitted with one modified blackbody, and these are the medians of the true over the fitted surface density in
# bins of the fitted surface density (log10, cm^-2), for several sets of wavebands and relative weights (cloud_tables.py). The
# surface layers of the cloud are heated by unattenuated ultraviolet light and are warmer than the relation T(Sigma_a) predicts,
# so even thin columns need about 5 percent. A bin with too few pixels repeats the one below it. The keys are the wavebands (um) and the relative calibration uncertainties.
CIRRUS_LOGX = np.array([20.55, 20.65, 20.75, 20.85, 20.95, 21.05, 21.15, 21.25, 21.35, 21.45, 21.55, 21.65])
CLOUD_TABLES = {
    # 160-500 um, weights 20/10/10/10 percent (standard)
    ((160.0, 250.0, 350.0, 500.0), (2.0, 1.0, 1.0, 1.0)):
        np.array([1.049, 1.051, 1.053, 1.058, 1.066, 1.075, 1.101, 1.121, 1.146, 1.180, 1.225, 1.251]),
    # 160-500 um, equal weights
    ((160.0, 250.0, 350.0, 500.0), (1.0, 1.0, 1.0, 1.0)):
        np.array([1.057, 1.058, 1.060, 1.067, 1.076, 1.086, 1.120, 1.139, 1.171, 1.214, 1.267, 1.294]),
    # 250-500 um, equal weights
    ((250.0, 350.0, 500.0), (1.0, 1.0, 1.0)):
        np.array([1.032, 1.034, 1.035, 1.039, 1.043, 1.049, 1.065, 1.076, 1.091, 1.110, 1.140, 1.160]),
    # 100-500 um, weights 20/20/10/10/10 percent
    ((100.0, 160.0, 250.0, 350.0, 500.0), (2.0, 2.0, 1.0, 1.0, 1.0)):
        np.array([1.088, 1.090, 1.094, 1.104, 1.120, 1.141, 1.194, 1.231, 1.282, 1.364, 1.430, 1.430]),
    # 70-500 um, weights 20/20/20/10/10/10 percent
    ((70.0, 100.0, 160.0, 250.0, 350.0, 500.0), (2.0, 2.0, 2.0, 1.0, 1.0, 1.0)):
        np.array([1.126, 1.128, 1.136, 1.150, 1.175, 1.216, 1.282, 1.341, 1.421, 1.536, 1.600, 1.600]),
}
CIRRUS_F = CLOUD_TABLES[((160.0, 250.0, 350.0, 500.0), (2.0, 1.0, 1.0, 1.0))]     # the table in use, set in main()
XBLEND = (21.6, 22.0)  # log10 of the fitted surface density over which the cloud calibration gives way to the analytic one
NTABLE = 80          # points of the tabulated relation
NTFIT = 2600         # points of the temperature grid of the fit
CAL_ERR = {70.0: 0.20, 100.0: 0.20, 160.0: 0.20, 250.0: 0.10, 350.0: 0.10, 500.0: 0.10}
                     # relative calibration uncertainty of each waveband, which sets its weight in the fit, as in hires
SMIN, SMAX = 1.0e20, 1.0e24      # range of the tabulated relation, cm^-2

# The field factor g, which adapts the correction factor to a radiation field other than G_0 = 10 (Habing units), measured on the
# grid of radiation fields (harvest_calib3_G0=*.txt.hist.npz). The pixels are binned in the true surface density, 0.05 dex wide,
# where the populations of all fields coincide; each bin gives the mean log10 of the fitted surface density and the mean ln f, and
# ln f is smoothed against the former by a local linear fit with Gaussian weights 0.1 dex wide times the square root of the number
# of pixels; the same at G_0 = 10 is subtracted (make_field_factor.py). Above 10^22.85 cm^-2 the grid has too few pixels, and the
# field factor is held at its value there. Rows are log10 of the fitted surface density, from GX0 in steps of GDX; columns are the
# fields GFIELD. It is interpolated linearly in log G_0 and held constant outside the table.
GFIELD = np.array([0.98, 2.12, 4.61, 10.0, 21.7, 47.1, 102.4, 222.2, 482.5])
GX0, GDX = 20.95, 0.1
GTAB = np.array([
    [0.994, 0.995, 0.998, 1.000, 1.000, 0.999, 0.996, 0.989, 0.995],  # 20.95
    [0.987, 0.992, 0.995, 1.000, 1.000, 0.999, 0.995, 0.989, 0.990],  # 21.05
    [0.979, 0.988, 0.991, 1.000, 1.001, 1.000, 0.995, 0.988, 0.986],  # 21.15
    [0.971, 0.981, 0.987, 1.000, 1.004, 1.004, 0.996, 0.989, 0.984],  # 21.25
    [0.958, 0.970, 0.981, 1.000, 1.009, 1.012, 1.001, 0.991, 0.981],  # 21.35
    [0.942, 0.957, 0.970, 1.000, 1.014, 1.022, 1.012, 0.995, 0.978],  # 21.45
    [0.924, 0.943, 0.962, 1.000, 1.028, 1.044, 1.035, 1.010, 0.984],  # 21.55
    [0.905, 0.928, 0.955, 1.000, 1.051, 1.080, 1.074, 1.040, 0.999],  # 21.65
    [0.878, 0.906, 0.943, 1.000, 1.082, 1.130, 1.130, 1.082, 1.014],  # 21.75
    [0.845, 0.877, 0.924, 1.000, 1.103, 1.175, 1.177, 1.113, 1.015],  # 21.85
    [0.815, 0.850, 0.907, 1.000, 1.121, 1.215, 1.209, 1.120, 1.004],  # 21.95
    [0.793, 0.833, 0.899, 1.000, 1.134, 1.242, 1.223, 1.110, 0.973],  # 22.05
    [0.785, 0.824, 0.896, 1.000, 1.140, 1.240, 1.206, 1.089, 0.943],  # 22.15
    [0.777, 0.817, 0.891, 1.000, 1.136, 1.223, 1.177, 1.052, 0.907],  # 22.25
    [0.768, 0.812, 0.888, 1.000, 1.126, 1.205, 1.147, 1.018, 0.869],  # 22.35
    [0.763, 0.810, 0.882, 1.000, 1.101, 1.181, 1.114, 0.977, 0.835],  # 22.45
    [0.762, 0.811, 0.880, 1.000, 1.080, 1.161, 1.083, 0.939, 0.798],  # 22.55
    [0.764, 0.810, 0.882, 1.000, 1.064, 1.132, 1.050, 0.908, 0.768],  # 22.65
    [0.770, 0.808, 0.889, 1.000, 1.075, 1.111, 1.024, 0.884, 0.746],  # 22.75
    [0.769, 0.807, 0.891, 1.000, 1.085, 1.098, 1.005, 0.863, 0.723],  # 22.85
    [0.769, 0.807, 0.891, 1.000, 1.085, 1.098, 1.005, 0.863, 0.723],  # 22.95
    [0.769, 0.807, 0.891, 1.000, 1.085, 1.098, 1.005, 0.863, 0.723],  # 23.05
    [0.769, 0.807, 0.891, 1.000, 1.085, 1.098, 1.005, 0.863, 0.723],  # 23.15
    [0.769, 0.807, 0.891, 1.000, 1.085, 1.098, 1.005, 0.863, 0.723],  # 23.25
])

### The field factor g, which adapts the correction factor to a radiation field other than G_0 = 10 (Habing units), measured on the
### grid of radiation fields (harvest_calib3_G0=*.txt): the pixel-weighted mean of log f in bins of 0.1 dex of the fitted surface
### density, smoothed over three bins, minus the same at G_0 = 10 (make_field_factor.py). Above 10^22.6 cm^-2 the grid has too few
### pixels, and the field factor is held at its value there. Rows are log10 of the fitted surface density, from GX0 in steps of GDX;
### columns are the fields GFIELD. It is interpolated linearly in log G_0 and held constant outside the table.
##GFIELD = np.array([0.98, 2.12, 4.61, 10.0, 21.7, 47.1, 102.4, 222.2, 482.5])
##GX0, GDX = 20.95, 0.1
##GTAB = np.array([
##    [0.991, 0.994, 0.999, 1.000, 1.001, 1.000, 0.997, 0.990, 0.994],  # 20.95
##    [0.986, 0.991, 0.995, 1.000, 1.001, 1.000, 0.997, 0.989, 0.990],  # 21.05
##    [0.978, 0.988, 0.992, 1.000, 1.001, 1.000, 0.997, 0.989, 0.986],  # 21.15
##    [0.970, 0.981, 0.987, 1.000, 1.003, 1.001, 0.996, 0.986, 0.982],  # 21.25
##    [0.960, 0.972, 0.981, 1.000, 1.007, 1.008, 1.001, 0.989, 0.982],  # 21.35
##    [0.944, 0.958, 0.973, 1.000, 1.013, 1.017, 1.008, 0.991, 0.979],  # 21.45
##    [0.920, 0.936, 0.959, 1.000, 1.026, 1.041, 1.034, 1.009, 0.985],  # 21.55
##    [0.904, 0.924, 0.955, 1.000, 1.047, 1.072, 1.068, 1.035, 0.999],  # 21.65
##    [0.904, 0.929, 0.967, 1.000, 1.124, 1.176, 1.178, 1.132, 1.010],  # 21.75
##    [0.837, 0.891, 0.936, 1.000, 1.110, 1.182, 1.188, 1.126, 1.027],  # 21.85
##    [0.815, 0.857, 0.889, 1.000, 1.084, 1.242, 1.248, 1.106, 1.029],  # 21.95
##    [0.775, 0.803, 0.919, 1.000, 1.102, 1.179, 1.209, 1.092, 0.991],  # 22.05
##    [0.756, 0.790, 0.909, 1.000, 1.107, 1.169, 1.162, 1.074, 0.936],  # 22.15
##    [0.753, 0.786, 0.846, 1.000, 1.092, 1.232, 1.126, 1.032, 0.871],  # 22.25
##    [0.756, 0.793, 0.851, 1.000, 1.099, 1.170, 1.120, 0.996, 0.838],  # 22.35
##    [0.796, 0.843, 0.899, 1.000, 1.186, 1.220, 1.120, 0.996, 0.865],  # 22.45
##    [0.767, 0.812, 0.916, 1.000, 1.168, 1.199, 1.114, 0.959, 0.834],  # 22.55
##    [0.761, 0.822, 0.907, 1.000, 1.135, 1.165, 1.082, 0.933, 0.801],  # 22.65
##    [0.761, 0.822, 0.907, 1.000, 1.135, 1.165, 1.082, 0.933, 0.801],  # 22.75
##    [0.761, 0.822, 0.907, 1.000, 1.135, 1.165, 1.082, 0.933, 0.801],  # 22.85
##    [0.761, 0.822, 0.907, 1.000, 1.135, 1.165, 1.082, 0.933, 0.801],  # 22.95
##    [0.761, 0.822, 0.907, 1.000, 1.135, 1.165, 1.082, 0.933, 0.801],  # 23.05
##    [0.761, 0.822, 0.907, 1.000, 1.135, 1.165, 1.082, 0.933, 0.801],  # 23.15
##    [0.761, 0.822, 0.907, 1.000, 1.135, 1.165, 1.082, 0.933, 0.801],  # 23.25
##])


def field_factor(logs, g0):
    """Field factor g of the correction factor for the field g0, at log10 of the fitted surface density logs."""
    x = np.asarray(logs, float)
    t = np.clip(np.log10(g0), np.log10(GFIELD[0]), np.log10(GFIELD[-1]))
    lg = np.log10(GFIELD)
    j = min(max(int(np.searchsorted(lg, t)) - 1, 0), len(lg) - 2)
    a = (t - lg[j]) / (lg[j + 1] - lg[j])
    xs = GX0 + GDX * np.arange(GTAB.shape[0])
    lr = (1.0 - a) * np.log(GTAB[:, j]) + a * np.log(GTAB[:, j + 1])
    return np.exp(np.interp(x, xs, lr))


def planck(t, wave_um):
    nu = C_LIGHT / (wave_um * 1.0e-4)
    x = np.clip(H_PLANCK * nu / (K_BOLTZ * np.asarray(t, float)), 1e-8, 500.0)
    return 2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2 / np.expm1(x)


def cnu(wave_um):
    """Intensity in MJy/sr per unit of B_nu times surface density in cm^-2."""
    nu = C_LIGHT / (wave_um * 1.0e-4)
    return (1.0e17 * KAPPA0 * (nu / NU0) ** BETA
            * DUST_TO_GAS * MU_H2 * AMU)


def temperature(s):
    """Dust temperature at an effective attenuating column s, in cm^-2."""
    return TFLOOR + (TSURF - TFLOOR) * (1.0 + np.asarray(s, float)
                                        / SSTAR) ** (-GAMMA)


def attenuating_column(sigma, geometry="mean"):
    """Effective attenuating column along a line of sight, in cm^-2.

    Returns the value at each of NDEPTH points through the middle of a uniform
    structure of total column sigma, averaged between a sphere and a
    plane-parallel layer. In both cases the angular average is taken over
    NANG directions; the sphere is sampled along a diameter and the layer
    along its normal.
    """
    t0 = KAPHEAT * sigma / 2.0                    # optical depth of the radius
    z = (np.arange(NDEPTH) + 0.5) / NDEPTH * 2.0 - 1.0
    ps = (np.arange(NANG) + 0.5) / NANG * np.pi
    wps = np.sin(ps)
    # sphere: distance to the boundary from radius |z|, in units of the radius
    lsph = (-np.abs(z)[:, None] * np.cos(ps)
            + np.sqrt(np.maximum(1.0 - (np.abs(z)[:, None]
                                        * np.sin(ps)) ** 2, 0.0)))
    jsph = (np.exp(-t0 * lsph) * wps).sum(axis=1) / wps.sum()
    # layer: the path is (1 -/+ z) / |mu|, both sides
    mu = (np.arange(NANG) + 0.5) / NANG * 2.0 - 1.0
    dlay = np.where(mu > 0.0, (1.0 - z[:, None]) / np.abs(mu),
                    (1.0 + z[:, None]) / np.abs(mu))
    jlay = np.exp(-t0 * dlay).mean(axis=1)
    ssph = -np.log(np.maximum(jsph, 1e-300)) / max(t0, 1e-12)
    slay = -np.log(np.maximum(jlay, 1e-300)) / max(t0, 1e-12)
    if geometry == "sphere":
        return ssph * sigma / 2.0
    return 0.5 * (ssph + slay) * sigma / 2.0


def parse_bands(text):
    """Wavelengths in micron of the fitted wavebands, from a list separated by commas, sorted."""
    return sorted(float(w) for w in text.split(",") if w.strip())


def cloud_key(waves, errs):
    """The key of CLOUD_TABLES for the wavebands and their relative calibration uncertainties (only the ratios matter)."""
    e = np.array(errs, float) / min(errs)
    return tuple(float(w) for w in waves), tuple(float(np.round(v, 3)) for v in e)


STANDARD_KEY = ((160.0, 250.0, 350.0, 500.0), (2.0, 1.0, 1.0, 1.0))


def build_relation(waves):
    """Fitted surface density and the factor that restores the true one.

    Both arrays are increasing in the fitted surface density, so they can be
    interpolated directly onto an observed image.
    """
    fitted, factor, _ = build_relation_t(waves)
    return fitted, factor


def build_relation_t(waves):
    """Fitted surface density, the factor that restores the true one, and the
    temperature that the fit returns.

    Two analytic relations are computed: the average of a uniform sphere and a
    uniform layer, and the sphere alone. Below the fitted surface density
    10^XBLEND[0], the factor is the average of the sphere and the RADMC-3D
    cloud, which replaces the analytic layer there. The difference between
    that and the analytic average is then reduced smoothly to zero at
    10^XBLEND[1], above which the factor is the analytic average."""
    fitted, factor, tfit = relation_geometry(waves, "mean")
    fsph, factsph, _ = relation_geometry(waves, "sphere")
    x = np.log10(fitted)
    xs = np.log10(fsph)
    fs = np.interp(x, xs, factsph)
    fc = np.interp(x, CIRRUS_LOGX, CIRRUS_F)              # held constant beyond the ends of the table
    d = 0.5 * (fs + fc) - factor
    x0, x1 = XBLEND
    d0 = np.interp(x0, x, d)
    d = np.where(x > x0, d0, d)
    w = np.where(x <= x0, 1.0, np.where(x >= x1, 0.0, 0.5 * (1.0 + np.cos(np.pi * (x - x0) / (x1 - x0)))))
    return fitted, factor + w * d, tfit


def relation_geometry(waves, geometry):
    """The analytic relation for the given geometry ("mean" or "sphere")."""
    wref = waves[-1]
    tgrid = np.geomspace(0.98 * TFLOOR, 1.02 * TSURF, NTFIT)
    model = np.array([[cnu(w) * planck(t, w) for w in waves] for t in tgrid])
    sigma = np.geomspace(SMIN, SMAX, NTABLE)
    fitted = np.empty(NTABLE)
    tfit = np.empty(NTABLE)
    for i, sig in enumerate(sigma):
        t = temperature(attenuating_column(sig, geometry))
        obs = np.array([cnu(w) * sig * np.mean(planck(t, w)) for w in waves])
        wgt = 1.0 / (np.array([CAL_ERR.get(w, 0.20) for w in waves]) * obs) ** 2
        amp = ((wgt * obs * model).sum(axis=1)
               / (wgt * model * model).sum(axis=1))
        chi2 = (wgt * (obs - amp[:, None] * model) ** 2).sum(axis=1)
        j = int(np.argmin(chi2))
        fitted[i] = (cnu(wref) * sig * np.mean(planck(t, wref))
                     / (cnu(wref) * planck(tgrid[j], wref)))
        tfit[i] = tgrid[j]
    return fitted, sigma / fitted, tfit


def invplanck(b, wave_um):
    """Temperature whose Planck function equals b at the wavelength in um."""
    nu = C_LIGHT / (wave_um * 1.0e-4)
    pref = 2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2
    return H_PLANCK * nu / (K_BOLTZ * np.log1p(pref / np.maximum(b, 1e-300)))


def stamp(header, g0, waves):
    """Record in a header how the image was corrected: the radiation field, the wavebands fitted and the temperature relation."""
    # the same CREATOR card as the Fortran version writes, which its FITS reader recognizes
    header["CREATOR"] = ("HICORRECT", "Alexander Menshchikov, DAp IRFU CEA Saclay")
    header["HICORR"] = (GAMMA, "index of the temperature relation applied")
    header["HIG0"] = (g0, "radiation field of the correction, Habing units")
    header["HIBANDS"] = (",".join("%.0f" % w for w in waves), "wavebands fitted to make the image, um")
    header["HIWGTS"] = (",".join("%g" % CAL_ERR.get(w, 0.20) for w in waves), "relative calibration errors of HIBANDS")
    header["HITMIN"] = (TFLOOR, "T_min of the temperature relation, K")
    header["HIT0"] = (TSURF, "T_max of the temperature relation, K")
    header["HISSTAR"] = (SSTAR, "Sigma_star of the temperature relation, cm^-2")
    header["HIGAMMA"] = (GAMMA, "gamma of the temperature relation")
    return header


def write_temperature(ftemp, f, good, wref, header, out, g0, waves):
    """Dust temperature that goes with the corrected surface density.

    The reconstruction derived the surface density from the reference
    waveband with the temperature it fitted at the resolution of the
    image. With f times more material on the line of sight, the
    same intensity of the reference waveband is reproduced by colder dust,
    whose Planck function is that of the fitted temperature divided by f. The
    fitted temperature image must have the resolution and the pixels of the corrected image.
    """
    t = np.nan_to_num(np.array(fits.open(ftemp)[0].data, float))
    if t.shape != f.shape:
        print(" %s skipped: its shape does not match the corrected image"
              % ftemp)
        return
    ok = good & (t > 0.0)
    tc = np.where(ok, invplanck(planck(np.where(ok, t, 1.0), wref)
                                / np.where(ok, f, 1.0), wref), 0.0)
    hh = stamp(header.copy(), g0, waves)
    hh["BUNIT"] = "K"
    fits.PrimaryHDU(np.asarray(tc, np.float32), hh).writeto(out,
                                                             overwrite=True)
    print(" corrected dust temperature: median %.2f K -> %s"
          % (np.median(tc[ok]), os.path.basename(out)))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    ftemp = [a[6:] for a in sys.argv[1:] if a.startswith("-temp=")]
    ftemp = ftemp[0] if ftemp else None
    g0 = [float(a[4:]) for a in sys.argv[1:] if a.startswith("-g0=")]
    g0 = g0[0] if g0 else 10.0
    bands = [a[7:] for a in sys.argv[1:] if a.startswith("-bands=")]
    weights = [a[9:] for a in sys.argv[1:] if a.startswith("-weights=")]
    if len(args) != 1:
        if len(args) > 1:
            print(" only one surface density image can be corrected in one run;"
                  " give the temperature image with -temp=<file>")
        else:
            print(__doc__)
        return
    # the file names may be given without the extension '.fits', which is then added
    fits_name = lambda n: n if n.endswith(".fits") else n + ".fits"
    fbase = fits_name(args[0])
    if ftemp is not None:
        ftemp = fits_name(ftemp)
    try:
        given = [float(w) for w in bands[0].split(",") if w.strip()] if bands else [160.0, 250.0, 350.0, 500.0]
    except ValueError:
        print(" cannot read the wavelengths from -bands=%s" % bands[0])
        return
    if len(given) < 2:
        print(" at least two wavebands must be given with -bands=")
        return
    if weights:
        try:
            errs = [float(e) for e in weights[0].split(",") if e.strip()]
        except ValueError:
            print(" cannot read the calibration uncertainties from -weights=%s" % weights[0])
            return
        if len(errs) != len(given) or min(errs) <= 0.0:
            print(" -weights= must give one positive value for each waveband of -bands=, in the same order")
            return
    else:
        errs = [0.20 if w <= 160.5 else 0.10 for w in given]
    order = np.argsort(given)
    waves = [given[i] for i in order]
    errs = [errs[i] for i in order]
    global CIRRUS_F, CAL_ERR
    CAL_ERR = dict(zip(waves, errs))
    key = cloud_key(waves, errs)
    if key in CLOUD_TABLES:
        CIRRUS_F = CLOUD_TABLES[key]
    else:
        CIRRUS_F = CLOUD_TABLES[STANDARD_KEY]
        print(" WARNING: no calibration on the diffuse cloud for these wavebands and weights; that of 160-500 um with the"
              " standard weights is used, and diffuse material may be corrected by several percent too much or too little")
    if g0 != 10.0 and key != STANDARD_KEY:
        print(" WARNING: the field factor was measured with 160-500 um and the standard weights; for other wavebands or"
              " weights it is approximate")
    fitted, factor = build_relation(waves)
    # The relation of uniform bodies declines above about 3e22 cm^-2, where such a body is nearly isothermal at the floor
    # temperature; the radiative transfer models of concentrated cores do not, so the factor is held at its maximum there. The hold
    # applies to the factor of G_0 = 10 only: the field factor, measured on the radiative transfer models, declines in strong
    # fields, and that decline is real.
    factor = np.maximum.accumulate(factor)
    if g0 != 10.0:
        if g0 < GFIELD[0] or g0 > GFIELD[-1]:
            print(" WARNING: G_0 = %g lies outside the grid of radiation fields (%g to %g);"
                  " the field factor of the nearest field is used" % (g0, GFIELD[0], GFIELD[-1]))
        factor = factor * field_factor(np.log10(fitted), g0)

    hdu = fits.open(fbase)[0]
    if "HICORR" in hdu.header or str(hdu.header.get("CREATOR", "")).startswith("HICORRECT"):
        print(" %s has already been corrected." % fbase)
        print(" The correction is not idempotent and must be applied once"
              " only.")
        return
    base = np.nan_to_num(np.array(hdu.data, float))
    good = base > 0.0
    f = np.where(good, np.interp(np.log10(np.maximum(base, 1.0)),
                                 np.log10(fitted), factor,
                                 left=factor[0], right=factor[-1]), 1.0)

    delta = np.where(good, (f - 1.0) * base, 0.0)

    print(" wavebands used: %s micron; reference %.0f micron; relative calibration uncertainties %s"
          % (", ".join("%.0f" % w for w in waves), waves[-1], ", ".join("%g" % CAL_ERR[w] for w in waves)))
    print(" T(Sigma_a) = %.4f + (%.4f - %.4f) (1 + Sigma_a / %.4e)^-%.5f K"
          % (TFLOOR, TSURF, TFLOOR, SSTAR, GAMMA))
    print(" effective attenuating column from an angular average over %d"
          " directions, uniform structures, sphere and layer averaged;" % NANG)
    print(" below 10^%.1f cm^-2 the layer is replaced by a RADMC-3D cloud,"
          " blended into the analytic relation up to 10^%.1f cm^-2" % XBLEND)
    print(" radiation field G_0 = %g (Habing units)%s" % (g0, "" if g0 != 10.0 else ", the field of the calibration"))
    print(" image %s: %d valid pixels" % (os.path.basename(fbase),
                                               int(good.sum())))
    p = np.percentile(f[good], [1, 25, 50, 75, 99])
    print(" correction factor  1%%: %.3f  25%%: %.3f  50%%: %.3f  75%%: %.3f"
          "  99%%: %.3f  max: %.3f" % (*p, f[good].max()))
    print(" ratio of the sums of the corrected to the uncorrected image:"
          " %.4f" % ((base + delta)[good].sum() / base[good].sum()))

    out = np.where(good, base + delta, 0.0)
    hh = stamp(hdu.header.copy(), g0, waves)
    hh["BUNIT"] = "H2/cm^2"
    # the outputs are written to the current directory, whatever directory the input comes from
    name = os.path.basename(fbase)
    stem = name[:-5] if name.endswith(".fits") else name
    fits.PrimaryHDU(np.asarray(out, np.float32), hh).writeto(stem + ".correct.fits", overwrite=True)
    print(" written %s.correct.fits" % os.path.basename(stem))
    hh2 = stamp(hdu.header.copy(), g0, waves)
    hh2["BUNIT"] = "ratio"
    fits.PrimaryHDU(np.asarray(np.where(good, f, 0.0), np.float32), hh2).writeto(stem + ".cfactor.fits", overwrite=True)
    print(" written %s.cfactor.fits" % os.path.basename(stem))
    if ftemp is None:
        print(" no temperature image given with -temp=<file>, so no temperature map is written")
        return
    th = fits.getheader(ftemp)
    b1, b2 = hdu.header.get("BEAM", 0.0), th.get("BEAM", 0.0)
    if b1 and b2 and abs(b2 / b1 - 1.0) > 0.01:
        print(" ERROR: the temperature image has the resolution %.2f arcsec, the surface density image %.2f arcsec;"
              " they must have the same resolution, so no temperature map is written" % (b2, b1))
        return
    write_temperature(ftemp, f, good, waves[-1], hdu.header,
                      (os.path.basename(ftemp)[:-5] if ftemp.endswith(".fits") else os.path.basename(ftemp)) + ".correct.fits", g0, waves)

if __name__ == "__main__":
    main()
