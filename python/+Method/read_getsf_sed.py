"""
read_getsf_sed.py — reader for the matched getsf+FITFLUXES catalog
(Aquila*.03=SED.00.cat) and injection-recovery analysis.

Merged-row layout (verified on the injected catalog, 230 tokens/row):
  getsf source table = cols 1..113  (1-based):
     11 global: NO S1 S2 SX XCO_P YCO_P WCS_ACOOR WCS_DCOOR FG SIGN GOOD
     then 6 bands x 17, band order = 070,160,161,250,350,500
     per band (b): FM SIGNM GOODM FXP_BST FXP_ERR FXT_BST FXT_ERR
                   FXT_ALT SCALE AFWHM BFWHM ASIZE BSIZE THETA FOFA FOOA FOOB
     -> col(band b, field offset o) = 11 + (b-1)*17 + o   (o: FXP_BST=4 ...
        FXT_BST=6, FXT_ALT=8, FOOA=16, FOOB=17)
  appended SED block (cols 114..):
     SED total mass  TOTL_MASS = col 121
     dust temp       DUST_TEMP = col 123
     surfdens mass (161 channel FXT_ALT, echoed) = col 53 == col 118
Verified: my own MBB fit tracks col123; col118==col53 exactly;
mean(RANDMASS, cols 133..230) ~= col121.

Key columns:
  XCO_P=5  YCO_P=6   (1-based pixel coords)
  SED_MASS=121  DUST_TEMP=123  SURFDENS_MASS=53
  161 footprint axes FOOA=61 FOOB=62 (arcsec, full axes)
"""
import numpy as np

def load_cat(path):
    rows = [l.split() for l in open(path)
            if l.strip() and not l.lstrip().startswith(('!', '#'))]
    return np.array([[ (float(x) if x not in ('ok', 'bad') else np.nan)
                       for x in r] for r in rows])

def col(A, c):            # 1-based column access
    return A[:, c-1]

# named accessors
def xco(A):  return col(A, 5)
def yco(A):  return col(A, 6)
def sed_mass(A): return col(A, 121)
def dust_temp(A): return col(A, 123)
def surfdens_mass(A): return col(A, 53)
def foo161(A): return col(A, 61), col(A, 62)   # FOOA, FOOB arcsec

def match_truth(A, tx, ty, plus_one=True, rmax=6.0):
    """Match truth (tx,ty) to nearest getsf source.  truth.csv is 0-based,
    getsf is 1-based -> set plus_one=True.  Returns index array + separations."""
    off = 1.0 if plus_one else 0.0
    gx, gy = xco(A), yco(A)
    idx, sep = [], []
    for x, y in zip(tx, ty):
        d = np.hypot(gx-(x+off), gy-(y+off)); j = int(np.argmin(d))
        idx.append(j); sep.append(d[j])
    return np.array(idx), np.array(sep)

if __name__ == '__main__':
    import csv, sys
    cat = sys.argv[1] if len(sys.argv) > 1 else 'Aquila.s.sources.ok.injected.03=SED.00.cat'
    A = load_cat(cat)
    T = list(csv.DictReader(open('truth.csv')))
    tx = np.array([float(r['x_pix']) for r in T])
    ty = np.array([float(r['y_pix']) for r in T])
    idx, sep = match_truth(A, tx, ty)
    MBE = 0.0800451
    print('matched %d/%d (max sep %.2f px)' % ((sep < 6).sum(), len(T), sep.max()))
    print('SED  mass  mean %.4f  bias %.2f' % (sed_mass(A)[idx].mean(),
                                               MBE/sed_mass(A)[idx].mean()))
    print('surfdens   mean %.4f  bias %.2f' % (surfdens_mass(A)[idx].mean(),
                                               MBE/surfdens_mass(A)[idx].mean()))
    print('dust temp  mean %.2f K' % dust_temp(A)[idx].mean())
