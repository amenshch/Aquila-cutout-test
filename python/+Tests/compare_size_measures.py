"""compare_size_measures.py -- test whether the model column FWHMSDbs of the
Bonnor-Ebert grid is the same quantity as the size getsf reports for the same
model after it has been injected into the map and re-extracted.

The five injection sets in tests/ are the same 24 models placed in the same
Aquila sub-field at five surface-density scalings, 0.25, 0.5, 1, 2 and 4.  Each
has a truth table giving the model tag and the pixel position of every injected
core, and an extraction catalogue of the resulting field.

For each injected core the nearest extracted source within MATCH_RADIUS is
taken as its recovery, and the comparison is between

    FWHMSDbs   major axis at half maximum of the model in the surface-density
               image, as tabulated in the grid catalogue          (arcsec)
    AFWHM03    major axis at half maximum reported by getsf in waveband 03,
               which in these extractions is the surface-density image at
               13.5 arcsec                                        (arcsec)

Note the waveband index.  In these injection extractions the wavebands are
070, 160, 161, 250, 350 and 500 micron, so the surface-density image is band
03.  In G. Zhang's whole-field Aquila extraction the wavebands are 160, 165,
250, 255, 350 and 500 micron, and the same image is band 02.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from getsf_columns import GetsfTable        # noqa: E402

PIX_ARCSEC = 3.0
MATCH_RADIUS_AS = 20.0
# name, extraction catalogue, truth table
SETS = [
    ('s2026',
     'inj/s2026_Aquila.s.sources.ok.cat',
     'inj/s2026_inj_cmf2.00_n400_s2026_truth.txt'),
    ('s2027',
     'inj/cmf2.00_n600_s2027/Aquila.s.sources.ok.cat',
     'inj/cmf2.00_n600_s2027/inj_cmf2.00_n600_s2027_truth.txt'),
    ('s2027_flat',
     'inj/cmf2.00_n600_s2027_flat/Aquila.s.sources.ok.cat',
     'inj/cmf2.00_n600_s2027_flat/inj_cmf2.00_n600_s2027_flat_truth.txt'),
]
GRID = 'cats/bes_model_grid_final2_catalog'


def load_grid_sizes(path=GRID):
    hdr = None
    for ln in open(path):
        if ln.lstrip()[:1] == '#' and 'ICSDbs' in ln.split():
            hdr = ln.lstrip('# ').split()
    rows = [l.split() for l in open(path)
            if l.strip() and l.lstrip()[:1] != '#']
    out = {}
    for r in rows:
        tag = 'i%02dj%02dk%02d' % (int(float(r[hdr.index('i')])),
                                   int(float(r[hdr.index('j')])),
                                   int(float(r[hdr.index('k')])))
        out[tag] = dict(fwhm=float(r[hdr.index('FWHMSDbs')]),
                        icsd=float(r[hdr.index('ICSDbs')]),
                        rbe=float(r[hdr.index('R_BE_as')]),
                        mbe=float(r[hdr.index('M_BE')]),
                        sde=float(r[hdr.index('SD_emb')]))
    return out


def read_truth(path):
    """id, tag, x_pix, y_pix, local_Sigma  (positions are one-based)."""
    rec = []
    for ln in open(path):
        if not ln.strip() or ln.lstrip()[0] in '#!':
            continue
        f = ln.split()
        rec.append(dict(tag=f[1], x=float(f[2]), y=float(f[3]),
                        local_sigma=float(f[4])))
    return rec


def main():
    gsz = load_grid_sizes()
    print('%-12s %6s %8s %10s %10s %12s %9s'
          % ('set', 'N_inj', 'N_match', 'med model', 'med getsf',
             'med ratio', 'scatter'))
    print('%-12s %6s %8s %10s %10s %12s %9s'
          % ('', '', '', 'FWHMSDbs', 'AFWHM03', 'getsf/model', 'dex'))
    allm, allg = [], []
    for xx, cpath, tpath in SETS:
        if not (os.path.exists(tpath) and os.path.exists(cpath)):
            print('  set %s: files missing' % xx)
            continue
        truth = read_truth(tpath)
        t = GetsfTable(cpath)

        def _col(name):
            # the injection catalogues are three getsf catalogues merged, so a
            # name that occurs in more than one of them comes back repeated,
            # with the copies of a row adjacent; take the first copy of each
            v = np.asarray(t.col(name), float)
            if len(v) == t.nrows:
                return v
            if len(v) % t.nrows == 0:
                return v.reshape(t.nrows, len(v) // t.nrows)[:, 0]
            raise ValueError('column %s has %d values for %d rows'
                             % (name, len(v), t.nrows))

        xs, ys, af = _col('XCO_P'), _col('YCO_P'), _col('AFWHM03')
        rmatch = MATCH_RADIUS_AS / PIX_ARCSEC
        mm, gm = [], []
        for r in truth:
            if r['tag'] not in gsz:
                continue
            d2 = (xs - r['x']) ** 2 + (ys - r['y']) ** 2
            k = int(np.argmin(d2))
            if d2[k] > rmatch ** 2:
                continue
            a = af[k]
            if not np.isfinite(a) or a <= 0 or a > 1e3:
                continue
            mm.append(gsz[r['tag']]['fwhm']); gm.append(a)
        mm = np.asarray(mm); gm = np.asarray(gm)
        if len(mm) < 3:
            print('%-12s %6d %8d   too few matches' % (xx, len(truth), len(mm)))
            continue
        lr = np.log10(gm / mm)
        print('%-12s %6d %8d %10.1f %10.1f %12.3f %9.3f'
              % (xx, len(truth), len(mm), np.median(mm), np.median(gm),
                 10 ** np.median(lr), np.std(lr)))
        allm.append(mm); allg.append(gm)
        # the same, split by the size of the model
        for a, b in ((0, 25), (25, 45), (45, 80), (80, 1e4)):
            k = (mm >= a) & (mm < b)
            if k.sum() >= 3:
                print('      model FWHMSDbs %4.0f - %4.0f arcsec: %3d matched, '
                      'median getsf/model = %.3f'
                      % (a, b, k.sum(), 10 ** np.median(np.log10(gm[k] / mm[k]))))
    if allm:
        mm = np.concatenate(allm); gm = np.concatenate(allg)
        lr = np.log10(gm / mm)
        print('\nover all %d matched injected cores: median getsf AFWHM03 '
              'divided by model FWHMSDbs = %.3f, scatter %.3f dex'
              % (len(lr), 10 ** np.median(lr), np.std(lr)))
        print('\n%-26s %6s %10s %10s %10s'
              % ('model FWHMSDbs (arcsec)', 'N', 'med model', 'med getsf',
                 'med ratio'))
        for a, b in ((0, 18), (18, 25), (25, 35), (35, 50), (50, 80),
                     (80, 150), (150, 1e4)):
            k = (mm >= a) & (mm < b)
            if k.sum() >= 3:
                print('%10.0f - %13.0f %6d %10.1f %10.1f %10.3f'
                      % (a, b, k.sum(), np.median(mm[k]), np.median(gm[k]),
                         10 ** np.median(np.log10(gm[k] / mm[k]))))
        # is the ratio a constant, or does it track the beam?
        import numpy.polynomial.polynomial as P
        c = np.polyfit(np.log10(mm), np.log10(gm), 1)
        print('\nlog10(AFWHM03) = %.3f * log10(FWHMSDbs) + %.3f'
              % (c[0], c[1]))
        beam = 13.5
        pred = np.sqrt(mm ** 2 + beam ** 2)
        print('if getsf reported the beam-convolved size sqrt(FWHMSDbs^2 + '
              '13.5^2), the median ratio')
        print('  of AFWHM03 to that prediction would be %.3f, scatter %.3f dex'
              % (10 ** np.median(np.log10(gm / pred)),
                 np.std(np.log10(gm / pred))))


if __name__ == '__main__':
    main()
