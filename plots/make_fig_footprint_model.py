#!/usr/bin/env python3
"""
make_fig_footprint_model.py -- generate fig_footprint_model.pdf: the
convolution factor c_lambda(D) = sigma(O_lambda;D)/sigma(O_Sigma;D) that
propagates the structural-noise power spectrum from the surface-density
band to each continuum band (Sect. "the per-band footprint model" in the
paper), as a function of footprint diameter D.

FULLY SELF-CONTAINED: pure analytic computation, no cached data or model
grid needed -- just the four angular resolutions (O_lambda) and the
adopted power-spectrum slope (gamma=-2.5).

Requires: numpy, scipy, matplotlib.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate

# angular resolution O_lambda [arcsec] of each band (paper's Sect. "the
# embedding cloud"): O_160=13.5, O_250=18.2, O_350=24.9, O_500=36.3;
# O_Sigma=13.5 for the high-resolution surface-density proxy band.
O_LAMBDA = {'Sigma': 13.5, '250': 18.2, '350': 24.9, '500': 36.3}
GAMMA = -2.5                          # adopted structural-noise power-spectrum slope
A_BEAM = 4 * np.pi**2 / (8 * np.log(2))

def sigma_beam(beam, klow):
    """sigma(beam; D): structural-noise power integrated above the spatial
    frequency corresponding to footprint diameter D=1/klow, convolved to
    the given beam."""
    f = lambda k: k**(GAMMA + 1) * np.exp(-A_BEAM * beam**2 * k**2)
    v, _ = integrate.quad(f, klow, np.inf, limit=200)
    return np.sqrt(max(v, 0))

Dgrid = np.linspace(20, 200, 120)
plt.rcParams.update({'font.size': 9, 'axes.linewidth': 0.8})
fig, ax = plt.subplots(figsize=(3.4, 3.0))
colors = {'250': 'tab:blue', '350': 'tab:green', '500': 'tab:red'}

for b in ('250', '350', '500'):
    c = [sigma_beam(O_LAMBDA[b], 1.0 / D) / sigma_beam(O_LAMBDA['Sigma'], 1.0 / D)
         for D in Dgrid]
    ax.plot(Dgrid, c, color=colors[b],
            label=r'%s $\mu$m ($O_{%s}=%.1f^{\prime\prime}$)' % (b, b, O_LAMBDA[b]))
    ax.axvline(O_LAMBDA[b], color=colors[b], ls=':', lw=0.7)

ax.axhline(1.0, color='0.5', lw=0.6)
ax.set_xlabel(r'footprint diameter $D\ (^{\prime\prime})$')
ax.set_ylabel(r'$c_\lambda(D)=\sigma(O_\lambda;D)/\sigma(O_\Sigma;D)$')
ax.legend(fontsize=7, frameon=False, loc='lower right')
ax.set_ylim(0, 1.05)
plt.tight_layout()
plt.savefig('fig_footprint_model.pdf')
print("saved fig_footprint_model.pdf")

for D in (40, 57, 80):
    print("D=%d: c_250=%.2f c_350=%.2f c_500=%.2f" % (D,
        sigma_beam(O_LAMBDA['250'], 1/D) / sigma_beam(O_LAMBDA['Sigma'], 1/D),
        sigma_beam(O_LAMBDA['350'], 1/D) / sigma_beam(O_LAMBDA['Sigma'], 1/D),
        sigma_beam(O_LAMBDA['500'], 1/D) / sigma_beam(O_LAMBDA['Sigma'], 1/D)))
