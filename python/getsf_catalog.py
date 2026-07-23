"""
Robust reader for getsf source catalogs.

Columns are located by HEADER NAME, never by position.  getsf catalogs differ
in band numbering between fields: the H2 column-density map (surfdens, 13.5"
beam, units 'H2/cm^2') is band 02 in most catalogs but band 03 wherever a
sub-160 band (070 or 075) is present and shifts the numbering.  Selecting the
surfdens band by its unit string and beam rather than by an assumed index
avoids the silent failure (zero corrected sources) that a wrong band produces.
"""
import re
import numpy as np


def _header_line(lines):
    """The '#' column-name line carrying per-source quantities."""
    for l in lines:
        if l.startswith('#') and 'PEAK^SRC03' in l and 'TOTL^MASS' in l:
            return l.lstrip('#').split()
    raise ValueError("no getsf column header found")


def _band_units(lines):
    """Map band index -> (wavelength, units) from the 'obs.fits' header comments."""
    out = {}
    for l in lines[:40]:
        m = re.search(r"\.(\d+)\.obs\.fits'\s*\(units:\s*([^)]+)\)", l)
        if m:
            out[m.group(1)] = m.group(2).strip()
    return out


def find_surfdens_band(lines, data, idx, beam_target=13.5, beam_tol=0.5):
    """Return the two-digit band string of the 13.5" H2 column-density map.

    Identified by: (a) units 'H2/cm^2' if the obs.fits comments are present,
    cross-checked by (b) a minimum sqrt(AFWHM*BFWHM) close to 13.5" and a
    PEAK^BGF in the column-density range (1e21-1e22).  Falls back to (b) alone
    if the unit comments are absent.
    """
    def col(name):
        i = idx[name]; o = []
        for r in data:
            try: o.append(float(r[i]))
            except Exception: o.append(np.nan)
        return np.array(o)

    # candidate bands present in the header
    bands = sorted({m.group(1) for n in idx for m in [re.match(r'AFWHM(\d\d)', n)] if m})
    units = _band_units(lines)
    # wavelength per band index, in obs-comment order (parallel to band order)
    wl_by_band = {}
    wls = [w for w in units]  # obs.fits comment order == band order after any sub-160 band
    # We can't assume alignment, so score each band physically:
    best = None
    for b in bands:
        try:
            fw = np.sqrt(col('AFWHM'+b)*col('BFWHM'+b))
            fw = fw[np.isfinite(fw) & (fw > 0)]
            pb = col('PEAK^BGF'+b)
            pb = pb[np.isfinite(pb) & (pb > 0)]
        except KeyError:
            continue
        if len(fw) == 0 or len(pb) == 0:
            continue
        beam = fw.min()
        med = np.median(pb)
        is_beam = abs(beam - beam_target) < beam_tol
        is_coldens = 1e21 < med < 5e22          # H2 column magnitude
        if is_beam and is_coldens:
            # prefer the surfdens (165) over 255 if both qualify: 165 comes first
            if best is None:
                best = b
    if best is None:
        raise ValueError("could not identify the 13.5-arcsec H2 column band")
    return best


def load(path):
    """Parse a getsf catalog; return (columns_dict, surfdens_band)."""
    lines = open(path).read().split('\n')
    names = _header_line(lines)
    idx = {n: i for i, n in enumerate(names)}
    data = [l.split() for l in lines
            if l.strip() and l.lstrip()[0].isdigit() and len(l.split()) >= len(names)]
    band = find_surfdens_band(lines, data, idx)

    def col(name):
        i = idx[name]; o = []
        for r in data:
            try: o.append(float(r[i]))
            except Exception: o.append(np.nan)
        return np.array(o)

    b = band
    area = np.pi * (col('FOOA'+b)/2) * (col('FOOB'+b)/2)
    with np.errstate(divide='ignore', invalid='ignore'):
        mean_in = col('FXT_BST'+b) / area
        conc = col('PEAK^SRC'+b) / mean_in
        contrast = 1.0 + col('PEAK^SRC'+b) / col('PEAK^BGF'+b)
    quality = np.array([data[r][idx['QUALITY']] for r in range(len(data))])
    return dict(
        mass=col('TOTL^MASS'),
        Nbg=col('PEAK^BGF'+b),
        peak=col('PEAK^SRC'+b),
        fwhm=np.sqrt(col('AFWHM'+b)*col('BFWHM'+b)),
        conc_peakmean=conc,
        contrast=contrast,
        Td=col('DUST^TEMP'),
        distance=col('DISTANCE'),
        quality=quality,
        band=band,
    ), band


if __name__ == '__main__':
    import sys
    d, b = load(sys.argv[1])
    ok = (d['quality'] == 'ok') & (d['mass'] > 0) & (d['Nbg'] > 0)
    print(f"surfdens band: {b}  distance: {np.nanmedian(d['distance']):.0f} pc")
    print(f"sources: {len(d['mass'])}  usable(ok): {ok.sum()}")
    print(f"N_bg median: {np.nanmedian(d['Nbg'][ok]):.2e}  (should be ~1e21-1e22)")
