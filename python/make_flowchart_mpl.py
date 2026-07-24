#!/usr/bin/env python3
"""Method-overview flowchart, drawn with matplotlib -> true vector PDF.

No cairosvg / SVG dependency. Subscripts use matplotlib mathtext, which
typesets them correctly, so the tspan/baseline-shift problems do not arise.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, BoxStyle
from matplotlib.patches import FancyArrowPatch

W, H = 1120, 1240
PAL = dict(
    calib=('#edf1f7', '#3f6296'), apply=('#e9f4f0', '#2f7d6e'),
    problem=('#f9edec', '#b0524a'), match=('#f1ecf8', '#6b4f9e'),
    output=('#e7f4ea', '#3f8a54'), valid=('#f4f4f4', '#8a8a8a'))
INK, BODY, MUTE, ARR = '#1f2933', '#3a4652', '#6b7683', '#59636e'

# markup -> mathtext.  Longest tokens first so Sigma_{cloud} beats bare Sigma.
TOK = [
    ('M_{reported}', r'$M_{\mathrm{reported}}$'),
    ('\u03a3_{cloud}', r'$\Sigma_{\mathrm{cloud}}$'),
    ('\u03be_{max}', r'$\xi_{\mathrm{max}}$'),
    ('PEAK^{BGF}', r'$\mathrm{PEAK}^{\mathrm{BGF}}$'),
    ('PEAK^{SRC}', r'$\mathrm{PEAK}^{\mathrm{SRC}}$'),
    ('M_{SED}', r'$M_{\mathrm{SED}}$'), ('M_{BE}', r'$M_{\mathrm{BE}}$'),
    ('T_{SED}', r'$T_{\mathrm{SED}}$'), ('R_{foot}', r'$R_{\mathrm{foot}}$'),
    ('R_{conv}', r'$R_{\mathrm{conv}}$'), ('f_{rec}', r'$f_{\mathrm{rec}}$'),
    ('C_{M}', r'$C_{\mathrm{M}}$'),
    ('\u03a3', r'$\Sigma$'), ('\u03be', r'$\xi$'),
]

def mt(s):
    for a, b in TOK:
        s = s.replace(a, b)
    return s

fig = plt.figure(figsize=(9.7, 10.75))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis('off')

def yy(y):            # flip y-down layout coords to matplotlib y-up
    return H - y

def box(x, y, w, h, kind, dashed=False, bold=False):
    fill, stroke = PAL[kind]
    p = FancyBboxPatch((x, yy(y + h)), w, h,
                       boxstyle=BoxStyle('Round', pad=0, rounding_size=10),
                       linewidth=2.3 if bold else 1.4,
                       edgecolor=stroke, facecolor=fill,
                       linestyle=(0, (6, 4)) if dashed else '-',
                       mutation_aspect=1)
    ax.add_patch(p)

def head(cx, y, s, size=10.5, color=INK):
    ax.text(cx, yy(y), mt(s), ha='center', va='baseline',
            fontsize=size, fontweight='bold', color=color)

def body(x, y, s, size=8.8, color=BODY, anchor='left', cx=None):
    if anchor == 'center':
        ax.text(cx, yy(y), mt(s), ha='center', va='baseline',
                fontsize=size, color=color)
    else:
        ax.text(x, yy(y), mt(s), ha='left', va='baseline',
                fontsize=size, color=color)

def arrow(x1, y1, x2, y2, color=ARR, dashed=False):
    ax.add_patch(FancyArrowPatch(
        (x1, yy(y1)), (x2, yy(y2)), arrowstyle='-|>', mutation_scale=13,
        lw=1.6, color=color, shrinkA=0, shrinkB=0,
        linestyle='--' if dashed else '-'))

def elbow(x1, y1, xc, x2, y2):
    """Down from (x1,y1) to turn row, across to x2, into (x2,y2)."""
    ax.plot([x1, x1, x2], [yy(y1), yy(xc), yy(xc)], color=ARR, lw=1.6,
            solid_capstyle='round')
    ax.add_patch(FancyArrowPatch((x2, yy(xc)), (x2, yy(y2)),
                 arrowstyle='-|>', mutation_scale=13, lw=1.6, color=ARR,
                 shrinkA=0, shrinkB=0))

# ---- title
head(560, 44, 'Radiative-transfer correction of Herschel prestellar core masses', 16)
body(560, 68, 'method overview', 9.5, MUTE, 'center', 560)

# ---- problem
box(90, 96, 940, 96, 'problem')
head(560, 124, 'The problem: SED fitting systematically underestimates core masses',
     11, PAL['problem'][1])
body(120, 150, '\u2022  temperature gradient \u2014 warm outer dust dominates the '
     'emission, so the fitted T is too high and M too low', 8.8)
body(120, 174, '\u2022  flux loss \u2014 the detection-limited footprint misses the '
     'extended, low-surface-brightness outer envelope', 8.8)

# ---- column headers
head(277, 226, 'CALIBRATION  \u00b7  build the corrector once', 10.5, PAL['calib'][1])
head(843, 226, 'APPLICATION  \u00b7  per observed cloud', 10.5, PAL['apply'][1])

LX, LW, LC = 55, 444, 277
RX, RW, RC = 621, 444, 843

# ---- calibration track
box(LX, 246, LW, 94, 'calib')
head(LC, 272, 'RADMC-3D model grid', 10.5)
body(LX+22, 296, 'pressure-truncated Bonnor\u2013Ebert spheres embedded in cloud', 8.6)
body(LX+22, 318, 'axes: \u03a3_{cloud}, mass, \u03be_{max}   \u2192   true mass '
     'M_{BE} known per node', 8.6)

box(LX, 362, LW, 62, 'calib')
head(LC, 387, 'Synthetic Herschel images, 70\u2013500 \u03bcm', 10.5)
body(LC, 409, 'each convolved to its own band beam (8\u201336 arcsec)', 8.6,
     BODY, 'center', LC)

box(LX, 446, LW, 108, 'calib')
head(LC, 471, 'Mimic the extraction pipeline', 10.5)
body(LX+22, 495, 'flat interpolated background subtraction  (\u201cbs\u201d)', 8.6)
body(LX+22, 517, 'SED fit on 250 / 350 / 500 \u03bcm  \u2192  M_{SED}, T_{SED}', 8.6)
body(LX+22, 539, 'source mask  \u2192  footprint radius R_{conv}', 8.6)

box(LX, 576, LW, 100, 'calib')
head(LC, 601, 'Recovery tables', 10.5)
body(LX+22, 623, 'per node \u00d7 band \u00d7 truncation radius r :', 8.6)
body(LX+22, 645, 'f_{rec},  R_{foot} = 2r / FWHM,  concentration,  peak / rim', 8.6)
body(LX+22, 667, 'tabulated so the matching variable is chosen later', 8.2, MUTE)

box(LX, 698, LW, 74, 'calib', bold=True)
head(LC, 721, 'Corrector  (distance-invariant)', 10.5)
body(LC, 742, 'grid of ( \u03a3, concentration, R_{foot} )  \u2192  C_{M}', 8.6,
     BODY, 'center', LC)
body(LC, 762, 'no absolute mass or size enters', 8.0, MUTE, 'center', LC)

# ---- application track
box(RX, 246, RW, 60, 'apply')
head(RC, 271, 'Real getsf catalogue', 10.5)
body(RC, 292, 'any cloud, any distance', 8.6, BODY, 'center', RC)

box(RX, 328, RW, 132, 'apply')
head(RC, 353, 'Measure each source', 10.5)
body(RX+26, 377, 'M_{SED}  \u2014  reported mass', 8.6)
body(RX+26, 399, '\u03a3  =  PEAK^{BGF}  \u2014  background column', 8.6)
body(RX+26, 421, 'concentration  =  PEAK^{SRC} / footprint-mean', 8.6)
body(RX+26, 443, 'R_{foot}  =  FOOA / AFWHM', 8.6)

# ---- convergence
MX, MW, MC = 300, 520, 560
box(MX, 812, MW, 104, 'match', bold=True)
head(MC, 838, 'Match source observables to the grid', 11, PAL['match'][1])
body(MC, 862, 'measure \u03a3, concentration, R_{foot} \u2014 all ratios, so '
     'distance cancels', 8.8, BODY, 'center', MC)
body(MC, 884, 'R_{foot} \u2248 1.8\u20132.0 (set by extraction) pins the recovered '
     'fraction', 8.8, BODY, 'center', MC)
body(MC, 905, 'Bonnor\u2013Ebert profiles keep only ~50% of their flux at this '
     'R_{foot}', 8.2, MUTE, 'center', MC)

# ---- output
OX, OW, OC = 226, 668, 560
box(OX, 946, OW, 120, 'output', bold=True)
head(OC, 973, 'Corrected mass    M  =  M_{reported} \u00d7 C_{M}', 12)
body(OC, 1000, 'C_{M}  =  temperature bias (\u22481.16)  \u00d7  flux recovery '
     '(\u2248 2\u20133)', 9, INK, 'center', OC)
body(OC, 1022, 'read straight from ( \u03a3, concentration, R_{foot} ) \u2014 same '
     'C_{M} at any distance', 8.4, BODY, 'center', OC)
body(OC, 1047, 'median C_{M} flat at ~3 over 130\u20131700 pc;  carry f_{rec} per '
     'source', 8.4, MUTE, 'center', OC)

# ---- validation
VX, VW, VC = 226, 668, 560
box(VX, 1100, VW, 92, 'valid', dashed=True)
head(VC, 1125, 'Validation', 10, PAL['valid'][1])
body(VX+30, 1149, '\u2022  temperature half \u2014 frac = 1 subdivision models '
     '\u2192 2.7%', 8.6)
body(VX+30, 1171, '\u2022  flux half \u2014 injection\u2013recovery of known-flux '
     'sources (in progress)', 8.6)

# ---- connectors
arrow(277, 192, 277, 244)
arrow(843, 192, 843, 244)
arrow(277, 340, 277, 360)
arrow(277, 424, 277, 444)
arrow(277, 554, 277, 574)
arrow(277, 676, 277, 696)
arrow(843, 306, 843, 326)
elbow(277, 772, 792, 452, 810)
elbow(843, 460, 790, 668, 810)
arrow(560, 916, 560, 944)
arrow(560, 1066, 560, 1098, PAL['valid'][1], dashed=True)

fig.savefig('/home/claude/proj/fig_method_flowchart.pdf')
fig.savefig('/home/claude/proj/fig_method_flowchart.png', dpi=150)
print('wrote pdf + png')
