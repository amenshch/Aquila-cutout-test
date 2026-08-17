#!/usr/bin/env python
"""
analyze_recovery.py -- match getsf extraction catalog to injection truth table,
apply RT mass correction, and report recovery statistics.

USAGE
-----
  python -B analyze_recovery.py  truth_table.txt  getsf_catalog.cat

OUTPUT
------
  recovery_SETTAG.txt  -- per-source table
  recovery_SETTAG.png  -- M_SED/M_BE and M_corr/M_BE vs M_BE figure

DEPENDENCIES
------------
  numpy, astropy, matplotlib, scipy
  fit_mass.py and bes_model_params_catalog on PYTHONPATH or same directory
"""

import sys, os, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------

# path to bes_model_params_catalog
BES_CATALOG = '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260616_RT_BES_radmc3d/bes_model_params_catalog'

# maximum center-to-center separation for a positional match (pixels)
# 1/3 of the 13.5" hires beam = 4.5" / 3"/pix = 1.5 px
MATCH_RADIUS_PIX = 1.5

# output directory ('.' = current directory)
OUT_DIR = '.'

# -----------------------------------------------------------------------
# END OF USER SETTINGS
# -----------------------------------------------------------------------


def parse_args():
    if len(sys.argv) != 3:
        print('Usage: python -B analyze_recovery.py truth_table.txt getsf_catalog.cat')
        sys.exit(1)
    return sys.argv[1], sys.argv[2]


def load_truth(path):
    truth = []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            p = line.split()
            # new format (10 cols): id model x y local_sigma SD_emb T_BE rho_BE M_BE R_BE_as
            # old format (5 cols):  id model x y local_sigma
            t = dict(
                id          = int(p[0]),
                tag         = p[1],
                x_pix       = int(p[2]),
                y_pix       = int(p[3]),
                local_sigma = float(p[4]),
            )
            if len(p) >= 10:
                t.update(dict(
                    SD_emb       = float(p[5]),
                    T_BE         = float(p[6]),
                    rho_BE       = float(p[7]),
                    M_BE         = float(p[8]),
                    R_BE_as      = float(p[9]),
                    scale_factor = float(p[10]) if len(p) > 10 else 1.0,
                ))
            else:
                # grid params will be filled from catalog lookup in main()
                t.update(dict(SD_emb=np.nan, T_BE=np.nan,
                              rho_BE=np.nan, M_BE=np.nan, R_BE_as=np.nan,
                              scale_factor=1.0))
            truth.append(t)
    return truth


def load_catalog(path):
    """Read getsf concatenated catalog.

    The combined header line starts with '#' as the first token, so each
    column name is at index (split-index - 1) relative to the data columns.
    """
    hdr_line = None
    data_lines = []
    for l in open(path):
        if l.startswith('#') and 'PEAK^BGF03' in l and 'TOTL^MASS' in l:
            hdr_line = l
        elif l.strip() and not l.lstrip().startswith(('!', '#')):
            data_lines.append(l.split())
    if hdr_line is None:
        raise ValueError('Cannot find header line with PEAK^BGF03 and TOTL^MASS')

    hdr = hdr_line.split()   # hdr[0] = '#', hdr[1] = first column name
    A = np.array([[float(x) if x not in ('ok', 'bad') else np.nan
                   for x in row] for row in data_lines])

    def col(name):
        # subtract 1 because hdr[0]='#' shifts all names by one
        if name not in hdr:
            raise KeyError(f'Column "{name}" not in header')
        return A[:, hdr.index(name) - 1]

    # XCO_P = pixel x, YCO_P = pixel y (verified against truth positions)
    return dict(
        x        = col('XCO_P'),
        y        = col('YCO_P'),
        M_SED    = col('TOTL^MASS'),
        T_SED    = col('DUST^TEMP'),
        sigma_bg = col('PEAK^BGF03'),
        fwhm250  = col('AFWHM04'),
    )


def get_model_params(tag, grid):
    """Look up grid parameters for a model tag like 'i1j2k4'."""
    ijk = tag[1:].split('j'); i = int(ijk[0])
    jk  = ijk[1].split('k'); j = int(jk[0]); k = int(jk[1])
    m = (grid['i']==i) & (grid['j']==j) & (grid['k']==k)
    if not m.any():
        return {}
    return {key: float(grid[key][m][0]) for key in
            ['M_BE','M_SED4bs','T_SED4bs','R_BE_as','FWHM250bs',
             'SD_emb','ICSDbs','T_BE','rho_BE']}


def main():
    truth_path, cat_path = parse_args()
    tag = os.path.splitext(os.path.basename(truth_path))[0]
    tag = tag.replace('inj_', '').replace('_truth', '')

    # load
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from fit_mass import load_catalog as load_grid, MassCorrector
    grid = load_grid(BES_CATALOG)
    ok   = grid['SD_emb'] != 2.4e22
    gcat = {k: v[ok] for k,v in grid.items()}
    corr = MassCorrector(gcat, sed_mass='M_SED4bs')

    truth  = load_truth(truth_path)
    # fill grid params for old-format truth tables that lack them
    for t in truth:
        if np.isnan(t['M_BE']):
            mdl = get_model_params(t['tag'], grid)
            t.update(dict(
                SD_emb  = mdl.get('SD_emb',  np.nan),
                T_BE    = mdl.get('T_BE',    np.nan) if 'T_BE' in mdl else np.nan,
                rho_BE  = np.nan,
                M_BE    = mdl.get('M_BE',    np.nan),
                R_BE_as = mdl.get('R_BE_as', np.nan),
            ))
    catdat = load_catalog(cat_path)
    gx = catdat['x']; gy = catdat['y']
    print(f'Truth entries: {len(truth)}  Catalog sources: {len(gx)}')

    # match and apply correction
    results = []
    for t in truth:
        mdl = get_model_params(t['tag'], grid)
        M_BE    = mdl.get('M_BE', np.nan)
        M_SED4  = mdl.get('M_SED4bs', np.nan)
        d = np.hypot(gx - t['x_pix'], gy - t['y_pix'])
        j = int(np.argmin(d)); sep = d[j]
        ms = catdat['M_SED'][j]
        if sep > MATCH_RADIUS_PIX:
            results.append({**t, 'status': 'NOT_DETECTED', 'sep': sep,
                            'M_BE': M_BE, 'M_SED': np.nan,
                            'T_SED': np.nan, 'M_corr': np.nan,
                            'ICSDbs': mdl.get('ICSDbs', np.nan)})
            continue
        mc = corr.correct(ms, t['local_sigma'])
        status = 'OK' if np.isfinite(mc) else 'OUT_OF_HULL'
        results.append({**t, 'status': status, 'sep': sep,
                        'M_BE': M_BE, 'M_SED': ms,
                        'T_SED': catdat['T_SED'][j], 'M_corr': mc,
                        'ICSDbs': mdl.get('ICSDbs', np.nan)})

    # print and write table
    hdr = ('# id  model         x_pix  y_pix  local_Sigma  SD_emb       '
           'M_BE    M_SED   M_corr  SED/BE  corr/BE  T_SED  status')
    lines = [
        f'# Injection recovery: {tag}',
        f'# Truth: {truth_path}',
        f'# Catalog: {cat_path}',
        hdr, '#' + '-' * (len(hdr)-1)
    ]
    ok_r = []
    for r in results:
        M_BE = r['M_BE']
        if r['status'] in ('NOT_DETECTED','CONTAM'):
            lines.append(f'  {r["id"]:2d}  {r["tag"]:12s}  {r["x_pix"]:5d}  {r["y_pix"]:5d}'
                         f'  {r["local_sigma"]:11.3e}  {r["SD_emb"]:11.3e}'
                         f'  {M_BE:6.4f}    ---     ---      ---      ---'
                         f'    ---  {r["status"]}')
        else:
            ms = r['M_SED']; mc = r['M_corr']
            sbe = ms/M_BE if M_BE>0 else np.nan
            cbe = mc/M_BE if (M_BE>0 and np.isfinite(mc)) else np.nan
            lines.append(f'  {r["id"]:2d}  {r["tag"]:12s}  {r["x_pix"]:5d}  {r["y_pix"]:5d}'
                         f'  {r["local_sigma"]:11.3e}  {r["SD_emb"]:11.3e}'
                         f'  {M_BE:6.4f}  {ms:6.4f}  '
                         f'{mc if np.isfinite(mc) else 0:6.4f}'
                         f'  {sbe:6.3f}  {cbe if np.isfinite(cbe) else 0:7.3f}'
                         f'  {r["T_SED"]:6.2f}  {r["status"]}')
            if r['status'] == 'OK':
                ok_r.append(r)

    n_det = sum(r['status'] in ('OK','OUT_OF_HULL') for r in results)
    lines += ['#',
              f'# Detected (clean): {n_det}/{len(truth)}',
              f'# Corrected (in hull): {len(ok_r)}/{n_det}']
    if ok_r:
        sbe  = np.array([r['M_SED']/r['M_BE']  for r in ok_r])
        cbe  = np.array([r['M_corr']/r['M_BE'] for r in ok_r])
        lines += [f'# Median M_SED/M_BE  = {np.median(sbe):.3f}  ({sbe.min():.3f}--{sbe.max():.3f})',
                  f'# Median M_corr/M_BE = {np.median(cbe):.3f}  ({cbe.min():.3f}--{cbe.max():.3f})']

    out_txt = os.path.join(OUT_DIR, f'recovery_{tag}.txt')
    with open(out_txt, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines[-8:]))
    print(f'\nTable: {out_txt}')

    # figure
    if not ok_r:
        print('No corrected detections -- skipping figure.')
        return

    M_BEs = np.array([r['M_BE']   for r in ok_r])
    SEDs  = np.array([r['M_SED']  for r in ok_r])
    CORRs = np.array([r['M_corr'] for r in ok_r])
    SDs   = np.array([r['SD_emb'] for r in ok_r])
    clr   = np.log10(SDs)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    cmap = 'plasma'

    ax = axes[0]
    ax.scatter(M_BEs, SEDs/M_BEs,  c=clr, cmap=cmap, s=45, marker='o',
               label='$M_\\mathrm{SED}/M_\\mathrm{BE}$',  zorder=3)
    ax.scatter(M_BEs, CORRs/M_BEs, c=clr, cmap=cmap, s=45, marker='^',
               label='$M_\\mathrm{corr}/M_\\mathrm{BE}$', zorder=3)
    ax.axhline(1.0, color='k', lw=1.2)
    ax.axhline(0.8, color='k', lw=0.8, ls='--', alpha=0.5)
    ax.set_xscale('log')
    ax.set_xlabel('$M_\\mathrm{BE}$ [M$_\\odot$]', fontsize=12)
    ax.set_ylabel('Recovered / True mass', fontsize=12)
    ax.set_title(f'Mass recovery: {tag}  (N={len(ok_r)})', fontsize=11)
    ax.legend(fontsize=9)
    sc = ax.scatter([], [], c=[], cmap=cmap, vmin=clr.min(), vmax=clr.max())
    plt.colorbar(plt.cm.ScalarMappable(
        norm=plt.Normalize(clr.min(), clr.max()), cmap=cmap),
        ax=ax, label='log$_{10}$(Σ$_\\mathrm{emb}$ [cm$^{-2}$])')

    ax = axes[1]
    sc2 = ax.scatter(M_BEs, CORRs/SEDs, c=clr, cmap=cmap, s=45, zorder=3)
    ax.axhline(1.0, color='k', lw=1.2, label='no correction')
    ax.axhline(1.15, color='k', lw=0.8, ls='--', alpha=0.5, label='+15%')
    ax.set_xscale('log')
    ax.set_xlabel('$M_\\mathrm{BE}$ [M$_\\odot$]', fontsize=12)
    ax.set_ylabel('$M_\\mathrm{corr}$ / $M_\\mathrm{SED}$', fontsize=12)
    ax.set_title('Mass correction factor', fontsize=11)
    ax.legend(fontsize=9)
    plt.colorbar(plt.cm.ScalarMappable(
        norm=plt.Normalize(clr.min(), clr.max()), cmap=cmap),
        ax=ax, label='log$_{10}$(Σ$_\\mathrm{emb}$ [cm$^{-2}$])')

    ax = axes[2]
    bins = np.logspace(np.log10(min(M_BEs.min(), SEDs.min())*0.5),
                       np.log10(max(M_BEs.max(), CORRs.max())*2), 22)
    ax.hist(M_BEs, bins=bins, histtype='step', lw=2, color='k',  label='$M_\\mathrm{BE}$ (truth)')
    ax.hist(SEDs,  bins=bins, histtype='step', lw=2, color='C0', label='$M_\\mathrm{SED}$')
    ax.hist(CORRs, bins=bins, histtype='step', lw=2, color='C3', label='$M_\\mathrm{corr}$')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('Mass [M$_\\odot$]', fontsize=12)
    ax.set_ylabel('N per bin', fontsize=12)
    ax.set_title('Recovered mass distributions', fontsize=11)
    ax.legend(fontsize=9); ax.set_ylim(0.7, None)

    plt.tight_layout()
    out_png = os.path.join(OUT_DIR, f'recovery_{tag}.png')
    plt.savefig(out_png, dpi=100)
    print(f'Figure: {out_png}')


if __name__ == '__main__':
    main()
