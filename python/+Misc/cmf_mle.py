"""Core mass function slopes for the seven clouds, before and after correction,
by maximum likelihood.

For a power-law mass function dN/dM proportional to M^(-alpha) above a lower
limit M_min, the maximum-likelihood estimate (Clauset, Shalizi & Newman 2009) is

    alpha = 1 + n / sum_i ln(M_i / M_min),      sigma_alpha = (alpha - 1)/sqrt(n)

over the n cores with M_i >= M_min.  Results are also quoted in the logarithmic
convention dN/dlogM proportional to M^(-x), with x = alpha - 1, in which the
Salpeter value is x = 1.35.
"""
import sys, numpy as np
sys.path.insert(0,'python')
import mass_correction_pipeline_5 as mc
mc.VERBOSE = 0

def mle(M, Mmin):
    m = M[np.isfinite(M) & (M >= Mmin)]
    n = len(m)
    if n < 20:
        return np.nan, np.nan, n
    a = 1.0 + n / np.sum(np.log(m / Mmin))
    return a, (a - 1.0) / np.sqrt(n), n

CL=[('Ophiuchus',139.),('Scorpius',150.),('Aquila',260.),('OrionA',432.),
    ('California',470.),('CygnusX',1400.),('TriRegion',2000.)]
S = mc.build_recovery_samples('cats/bes_model_grid_final2_catalog',
                              'cats/bes_model_grid_final2_recovery_tables_v3')
f4 = mc.cm_size_interp(S)
lg = np.log10
DATA = {}
for name, dist in CL:
    d = 'W3W4W5-Guoyin' if name == 'TriRegion' else '%s-Guoyin' % name
    mc.load_getsf._dist = dist
    g = mc.load_getsf(['%s/%s.sw.sources.ok.cat' % (d, name),
                       '%s/%s.sw.sources.ok.add.cat' % (d, name)])
    o = mc.source_mask(g)
    u = o & (g['Nbg']>0) & (g['conc']>0) & (g['rfoot']>0) & (g['foot_radius_as']>0) \
        & (g['mass']>0)
    q = np.column_stack([lg(g['Nbg'][u]), lg(g['conc'][u]), lg(g['rfoot'][u]),
                         lg(g['foot_radius_as'][u])])
    c = 10 ** f4(q)
    DATA[name] = (dist, g['mass'][u], c, f4.last_call_bounded())

print("Maximum-likelihood slope of the core mass function, in the logarithmic")
print("convention dN/dlogM proportional to M^(-x); Salpeter is x = 1.35.")
print("M_min is the lower mass limit of the fit, in solar masses.  n is the number")
print("of cores above it.  Uncertainties are (alpha-1)/sqrt(n).\n")
for Mmin in (0.1, 0.3, 1.0):
    print("--- M_min = %.1f Msun ---" % Mmin)
    print("%-11s %6s %6s %16s %6s %16s %10s"
          % ("cloud","d(pc)","n_unc","x uncorrected","n_cor","x corrected","change"))
    for name, dist in CL:
        _, M, c, b = DATA[name]
        au, su, nu = mle(M, Mmin)
        ac, sc, nc = mle(M * c, Mmin)
        if not np.isfinite(au) or not np.isfinite(ac):
            print("%-11s %6.0f %6s %16s %6s %16s %10s"
                  % (name, dist, nu, "-", nc, "-", "-")); continue
        print("%-11s %6.0f %6d %8.2f +/- %.2f %6d %8.2f +/- %.2f %+10.2f"
              % (name, dist, nu, au-1, su, nc, ac-1, sc, (ac-au)))
    print()
