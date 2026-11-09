"""
hires_forward.py -- forward model of the Herschel observing chain for the
Bonnor-Ebert model grid, used to test the HIRES surface density reconstruction.

For a given node of the grid "bes_grid_final3" the module builds:

  * Sigma_true(theta)  : the true surface density of hydrogen molecules,
                         N(H2) in cm^-2, as a function of angular offset
                         theta in arcsec from the center of the core;
  * I_lambda(theta)    : the specific intensity in MJy sr^-1 in each of the
                         six Herschel wavebands, computed by optically thin
                         line-of-sight integration through the spherically
                         symmetric radiative-transfer model;
  * the same quantities sampled on a square map with a pixel size of 3 arcsec,
    after convolution to each waveband beam.

Notation follows the paper: the dust temperature is T_d(r), the gas mass
density is rho(r), the dust-to-gas mass ratio is 0.01, the mean molecular
weight per hydrogen molecule is mu = 2.8.
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------
# physical constants (cgs)
# ----------------------------------------------------------------------
H_PLANCK = 6.62607015e-27      # erg s
K_BOLTZ = 1.380649e-16         # erg K^-1
C_LIGHT = 2.99792458e10        # cm s^-1
AU_CM = 1.495978707e13         # cm
PC_CM = 3.0856775814913673e18  # cm
AMU = 1.6605402e-24            # g
MU_H2 = 2.8                    # mean molecular weight per H2
DUST_TO_GAS = 1.0e-2           # dust-to-gas mass ratio

# ----------------------------------------------------------------------
# observing chain
# ----------------------------------------------------------------------
DIST_PC = 260.0                # assumed distance, pc
PIX_ARCSEC = 3.0               # map pixel size, arcsec
WAVES = np.array([70.0, 100.0, 160.0, 250.0, 350.0, 500.0])       # micron
BEAMS = np.array([8.4, 9.4, 13.5, 18.2, 24.9, 36.3])              # arcsec FWHM
FWHM_TO_SIGMA = 1.0 / 2.354820045

# opacity assumed by HIRES: kappa_0 = 10 cm^2 per gram of dust at 1000 GHz
KAPPA0_DUST = 10.0             # cm^2 g^-1 (dust) at nu_0
NU0 = 1.0e12                   # Hz
BETA = 2.0


def kappa_abs(wave_um):
    """Dust absorption opacity in cm^2 per gram of dust at wave_um micron."""
    nu = C_LIGHT / (wave_um * 1.0e-4)
    return KAPPA0_DUST * (nu / NU0) ** BETA


def planck(temp_k, wave_um):
    """Planck function B_nu in erg s^-1 cm^-2 Hz^-1 sr^-1."""
    nu = C_LIGHT / (np.asarray(wave_um, float) * 1.0e-4)
    x = H_PLANCK * nu / (K_BOLTZ * np.asarray(temp_k, float))
    return (2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2) / np.expm1(np.clip(x, 1e-8, 700.0))


# ----------------------------------------------------------------------
# isothermal Lane-Emden solution (same construction as make_bes_grid3.py)
# ----------------------------------------------------------------------
def _lane_emden(xi_max_tab=80.0, n=200000):
    def rhs(xi, y):
        p, dp = y
        return [dp, (np.exp(-p) - 2.0 * dp / xi) if xi > 1e-9 else 0.0]
    grid = np.linspace(1e-6, xi_max_tab, n)
    sol = solve_ivp(rhs, [1e-6, xi_max_tab], [0.0, 0.0], t_eval=grid,
                    method="RK45", rtol=1e-10, atol=1e-12)
    return interp1d(sol.t, sol.y[0], kind="cubic")


_PSI = _lane_emden()


def be_density_ratio(xi):
    """rho(xi)/rho_c for the isothermal Bonnor-Ebert sphere.

    The argument is clipped to the tabulated interval; values beyond xi_max of
    a given node are never used, since the cloud pedestal applies there.
    """
    return np.exp(-_PSI(np.clip(xi, 1e-6, 79.999)))


# ----------------------------------------------------------------------
# grid node and radiative-transfer profile
# ----------------------------------------------------------------------
def load_profile_archive(path):
    """Read the concatenated profile file written by collect_profiles.py.

    Returns a dictionary mapping the node identifier to the numeric table of
    that node, with the columns of mc.denstemps.tab preserved.
    """
    import gzip
    import re
    opener = gzip.open if str(path).endswith(".gz") else open
    out, cur, rows = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("### NODE"):
                if cur is not None and rows:
                    out[cur] = np.array(rows, float)
                cur = int(re.search(r"ID=(\d+)", line).group(1))
                rows = []
            elif line.startswith("###") or line.startswith("!") or not line.strip():
                continue
            elif cur is not None:
                rows.append(line.split())
    if cur is not None and rows:
        out[cur] = np.array(rows, float)
    return out


class Node:
    """One node of the model grid, with its density and temperature profiles."""

    def __init__(self, node_id, grid_path, prof_path=None, prof_table=None):
        rows = [l.split() for l in open(grid_path) if l.strip()]
        hdr = rows[0]
        for r in rows[1:]:
            if int(r[0]) == node_id:
                rec = dict(zip(hdr, r))
                break
        else:
            raise ValueError("node %d not found" % node_id)
        self.id = node_id
        self.idx_i = int(float(rec["i"]))
        self.idx_j = int(float(rec["j"]))
        self.idx_k = int(float(rec["k"]))
        self.sd_emb = float(rec["SD_emb"])        # cm^-2, cloud surface density
        self.m_be = float(rec["M_BE"])            # Msun
        self.fwhm_pc = float(rec["FWHM_pc"])      # pc
        self.xi_max = float(rec["xi_max"])
        self.rho_c = float(rec["rho_c"])          # g cm^-3
        self.rho_emb = float(rec["rho_emb"])      # g cm^-3
        self.r0_au = float(rec["r0_AU"])
        self.r_out_au = float(rec["R_out_AU"])    # = xi_max * r0
        self.r_cloud_au = float(rec["R_cloud_AU"])
        self.peak_exc = float(rec["peak_exc"])    # cm^-2, tabulated reference
        self.contrast = float(rec["contrast"])
        self.t_be = float(rec["T_BE"])

        if prof_table is not None:
            prof = np.asarray(prof_table, float)
        else:
            prof = np.array([l.split() for l in open(prof_path)
                             if not l.startswith("!") and l.strip()], float)
        self.prof_r_as = prof[:, 5]               # arcsec, cell centered
        self.prof_rho = prof[:, 6]                # g cm^-3
        self.prof_tdust = prof[:, 11]             # K

        # angular scales, arcsec
        self.r0_as = self.r0_au / DIST_PC
        self.r_out_as = self.r_out_au / DIST_PC
        # the profile table lists cell-centered radii, so its last entry lies
        # inside the cloud boundary; the boundary itself is taken from the grid
        self.r_cloud_as = self.r_cloud_au / DIST_PC

        self._t_interp = interp1d(np.log(self.prof_r_as), self.prof_tdust,
                                  kind="linear", bounds_error=False,
                                  fill_value=(self.prof_tdust[0], self.prof_tdust[-1]))

    # ---- profiles as functions of the spherical radius in arcsec ----
    def rho_of_r(self, r_as):
        """Gas mass density in g cm^-3 at spherical radius r_as (arcsec)."""
        r_as = np.asarray(r_as, float)
        rho = np.where(r_as < self.r_out_as,
                       self.rho_c * be_density_ratio(r_as / self.r0_as),
                       self.rho_emb)
        return np.where(r_as > self.r_cloud_as, 0.0, rho)

    def tdust_of_r(self, r_as):
        """Dust temperature in K at spherical radius r_as (arcsec)."""
        r = np.clip(np.asarray(r_as, float), self.prof_r_as[0], self.prof_r_as[-1])
        return self._t_interp(np.log(r))

    # ---- line-of-sight integration ----
    def _los_nodes(self, b_as, n=2400):
        """Integration nodes z (arcsec) along the line of sight at impact
        parameter b_as, from 0 to the cloud boundary, concentrated near z = 0."""
        zmax = np.sqrt(max(self.r_cloud_as ** 2 - b_as ** 2, 0.0))
        if zmax <= 0.0:
            return np.array([0.0])
        zmin = 1.0e-4 * max(self.r0_as, 1.0e-3)
        z = np.concatenate(([0.0], np.geomspace(min(zmin, zmax * 1e-6), zmax, n)))
        # add nodes bracketing the Bonnor-Ebert boundary to resolve the jump
        if b_as < self.r_out_as:
            zb = np.sqrt(self.r_out_as ** 2 - b_as ** 2)
            z = np.concatenate((z, zb * np.array([1 - 1e-9, 1 + 1e-9])))
        return np.unique(np.clip(z, 0.0, zmax))

    def profiles(self, theta_as):
        """Return (Sigma_true, I) for the given angular offsets.

        Sigma_true : N(H2) in cm^-2, array of shape (n_theta,)
        I          : specific intensity in MJy sr^-1, shape (6, n_theta)

        The integration along the line of sight is vectorized over all impact
        parameters at once, on a single grid of the coordinate z that is
        geometrically spaced so as to resolve the center of the core.
        """
        theta_as = np.atleast_1d(np.asarray(theta_as, float))
        cm_per_as = DIST_PC * PC_CM / 206264.806
        kap = np.array([kappa_abs(w) for w in WAVES])

        zmax = self.r_cloud_as
        zmin = 1.0e-5 * max(self.r0_as, 1.0e-3)
        z = np.concatenate(([0.0], np.geomspace(zmin, zmax, 1500)))
        # nodes bracketing the Bonnor-Ebert boundary, to resolve the jump
        z = np.unique(np.concatenate((z, self.r_out_as *
                                      np.array([1 - 1e-9, 1 + 1e-9]))))

        # tabulate the density and the emissivity on a radial grid once, then
        # interpolate them onto the two-dimensional grid of radii; this avoids
        # evaluating the Lane-Emden solution and the Planck function millions
        # of times
        rt = np.geomspace(max(1e-4 * self.r0_as, 1e-6), self.r_cloud_as, 6000)
        rt = np.unique(np.concatenate((
            [0.0], rt, self.r_out_as * np.array([1 - 1e-9, 1 + 1e-9]))))
        rho_t = self.rho_of_r(rt)
        td_t = self.tdust_of_r(rt)
        emis_t = np.array([(kap[k] * DUST_TO_GAS) * rho_t * planck(td_t, w)
                           for k, w in enumerate(WAVES)])

        r = np.hypot(theta_as[:, None], z[None, :])
        rho = np.interp(r, rt, rho_t)
        sigma = 2.0 * np.trapezoid(rho, z, axis=1) * cm_per_as / (MU_H2 * AMU)
        inten = np.empty((len(WAVES), len(theta_as)))
        for k in range(len(WAVES)):
            emis = np.interp(r, rt, emis_t[k])
            inten[k] = 2.0 * np.trapezoid(emis, z, axis=1) * cm_per_as / 1.0e-17
        return sigma, inten


# ----------------------------------------------------------------------
# maps
# ----------------------------------------------------------------------
def gaussian_convolve(image, fwhm_as, pix_as):
    """Convolve a square map with a circular Gaussian of the given FWHM.

    Edge pixels are replicated outward before the convolution, reproducing the
    edge treatment of the utility 'convolve' used by the HIRES script.
    """
    if fwhm_as <= 0.0:
        return image.copy()
    sig = fwhm_as * FWHM_TO_SIGMA / pix_as
    pad = int(np.ceil(4.0 * sig))
    a = np.pad(image, pad, mode="edge")
    n = a.shape[0]
    y, x = np.mgrid[0:n, 0:n]
    c = (n - 1) / 2.0
    k = np.exp(-((x - c) ** 2 + (y - c) ** 2) / (2.0 * sig ** 2))
    k /= k.sum()
    f = np.fft.rfft2(a) * np.fft.rfft2(np.fft.ifftshift(k))
    out = np.fft.irfft2(f, a.shape)
    return out[pad:pad + image.shape[0], pad:pad + image.shape[1]]


def build_maps(node, half_as=402.0, subpix=3, verbose=False):
    """Build the true surface density map and the six intensity maps.

    The model is first sampled on a grid of pix = PIX_ARCSEC/subpix arcsec,
    convolved to each waveband beam, and then block-averaged to the 3 arcsec
    map pixel, reproducing the observing chain (convolve, then resample).

    Returns a dictionary with
      'theta'      : radius of each 3 arcsec pixel from the center, arcsec
      'sigma_true' : true N(H2) map at full resolution, cm^-2
      'sigma_13p5' : true N(H2) convolved to 13.5 arcsec, cm^-2
      'images'     : (6, n, n) intensity maps at the native beam, MJy sr^-1
    """
    fine = PIX_ARCSEC / subpix
    n3 = int(round(2 * half_as / PIX_ARCSEC))
    if n3 % 2 == 0:
        n3 += 1
    nf = n3 * subpix
    c = (nf - 1) / 2.0
    yy, xx = np.mgrid[0:nf, 0:nf]
    rad = np.hypot(xx - c, yy - c) * fine

    # tabulate the radial profiles once, then interpolate onto the fine grid
    rmax = rad.max() * 1.001
    tab = np.unique(np.concatenate((
        np.linspace(0.0, min(30.0, rmax), 900),
        np.geomspace(max(1e-3, 0.005 * node.r0_as), rmax, 1600))))
    sig_t, int_t = node.profiles(tab)

    sigma_fine = np.interp(rad, tab, sig_t)
    img_fine = np.empty((len(WAVES), nf, nf))
    for k in range(len(WAVES)):
        img_fine[k] = np.interp(rad, tab, int_t[k])

    def block(a):
        return a.reshape(n3, subpix, n3, subpix).mean(axis=(1, 3))

    out = {}
    out["n"] = n3
    out["pix"] = PIX_ARCSEC
    cc = (n3 - 1) / 2.0
    y3, x3 = np.mgrid[0:n3, 0:n3]
    out["theta"] = np.hypot(x3 - cc, y3 - cc) * PIX_ARCSEC
    out["sigma_true"] = block(sigma_fine)
    out["sigma_13p5"] = block(gaussian_convolve(sigma_fine, 13.5, fine))
    out["sigma_18p2"] = block(gaussian_convolve(sigma_fine, 18.2, fine))
    imgs = np.empty((len(WAVES), n3, n3))
    for k in range(len(WAVES)):
        imgs[k] = block(gaussian_convolve(img_fine[k], BEAMS[k], fine))
    out["images"] = imgs
    if verbose:
        print("map: %d x %d pixels of %.1f arcsec, built at %.2f arcsec sampling"
              % (n3, n3, PIX_ARCSEC, fine))
    return out
