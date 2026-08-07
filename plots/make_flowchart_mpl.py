#!/usr/bin/env python3
"""Method-overview flowchart, drawn with matplotlib -> true vector PDF.

Revision 4: fixed a real bug in the layout checker itself. The previous
per-box registration loop only added a text object to the "check against
this box" list if it was ALREADY geometrically inside the box -- so any
text that genuinely overflowed was silently excluded from the check rather
than flagged, which is how a real overflow (Validation box bullets) passed
the "layout check passed" message. Registration is now done by recording
the ax.texts index range added by each content block (via a start-index
marker) and registering that whole range unconditionally, so overflow
cannot be invisible to the checker. Two automated checks run before the
figure is saved and abort the script if either fails: (1) no text block
exceeds its host box, (2) no two text blocks visually overlap each other.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, BoxStyle
from matplotlib.patches import FancyArrowPatch

W, H = 1160, 1358
PAL = dict(
    calib=('#edf1f7', '#3f6296'), apply=('#e9f4f0', '#2f7d6e'),
    problem=('#f9edec', '#b0524a'), key=('#fdf3e3', '#a8721f'),
    output=('#e7f4ea', '#3f8a54'), valid=('#f4f4f4', '#6b6b6b'))
INK, BODY, MUTE, ARR = '#1f2933', '#3a4652', '#6b7683', '#59636e'

TOK = [
    ('M_{reported}', r'$M_{\mathrm{reported}}$'),
    ('M_{corrected}', r'$M_{\mathrm{corrected}}$'),
    ('\u03a3_{cloud}', r'$\Sigma_{\mathrm{cloud}}$'),
    ('\u03be_{max}', r'$\xi_{\mathrm{max}}$'),
    ('PEAK^{BGF}', r'$\mathrm{PEAK}^{\mathrm{BGF}}$'),
    ('PEAK^{SRC}', r'$\mathrm{PEAK}^{\mathrm{SRC}}$'),
    ('M_{SED}', r'$M_{\mathrm{SED}}$'), ('M_{BE}', r'$M_{\mathrm{BE}}$'),
    ('T_{SED}', r'$T_{\mathrm{SED}}$'), ('R_{foot}', r'$R_{\mathrm{foot}}$'),
    ('R_{conv}', r'$R_{\mathrm{conv}}$'),
    ('R_{\u03bb}', r'$R_{\lambda}$'), ('c_{\u03bb}', r'$c_{\lambda}$'),
    ('f_{rec,\u03bb}', r'$f_{\mathrm{rec},\lambda}$'),
    ('C_{M}', r'$C_{\mathrm{M}}$'), ('F_{flux}', r'$F_{\mathrm{flux}}$'),
    ('F_{bg}', r'$F_{\mathrm{bg}}$'), ('F_{temp}', r'$F_{\mathrm{temp}}$'),
    ('\u03a3', r'$\Sigma$'), ('\u03be', r'$\xi$'), ('\u03b3', r'$\gamma$'),
    ('\u03bb', r'$\lambda$'),
]

def mt(s):
    for a, b in TOK:
        s = s.replace(a, b)
    return s

fig = plt.figure(figsize=(9.7, 11.35))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis('off')
fig.canvas.draw()
_renderer = fig.canvas.get_renderer()

def yy(y):
    return H - y

def box(x, y, w, h, kind, dashed=False, bold=False):
    fill, stroke = PAL[kind]
    p = FancyBboxPatch((x, yy(y + h)), w, h,
                       boxstyle=BoxStyle('Round', pad=0, rounding_size=10),
                       linewidth=2.4 if bold else 1.5,
                       edgecolor=stroke, facecolor=fill,
                       linestyle=(0, (6, 4)) if dashed else '-',
                       mutation_aspect=1)
    ax.add_patch(p)
    return p

def head(cx, y, s, size=11.5, color=INK):
    ax.text(cx, yy(y), mt(s), ha='center', va='baseline',
            fontsize=size, fontweight='bold', color=color)

_probe = [None]
def _text_width_px(s, size, weight='normal'):
    if _probe[0] is not None:
        _probe[0].remove()
    _probe[0] = ax.text(0, 0, s, fontsize=size, fontweight=weight, alpha=0)
    fig.canvas.draw()
    return _probe[0].get_window_extent(_renderer).width

def _wrap(words, max_width_px, size):
    lines, cur = [], []
    for w in words:
        trial = cur + [w]
        if cur and _text_width_px(mt(' '.join(trial)), size) > max_width_px:
            lines.append(' '.join(cur)); cur = [w]
        else:
            cur = trial
    if cur:
        lines.append(' '.join(cur))
    return lines

LH_DEFAULT = 23        # generous default: safe for mathtext sub/superscripts

def body_wrap(x, y, s, max_width_px, size=8.6, color=BODY, anchor='left',
              cx=None, lh=LH_DEFAULT):
    lines = _wrap(s.split(' '), max_width_px, size)
    if _probe[0] is not None:          # remove any dangling wrap-measurement
        _probe[0].remove(); _probe[0] = None   # probe before adding real text
    for i, ln in enumerate(lines):
        yi = y + i * lh
        if anchor == 'center':
            ax.text(cx, yy(yi), mt(ln), ha='center', va='baseline',
                    fontsize=size, color=color)
        else:
            ax.text(x, yy(yi), mt(ln), ha='left', va='baseline',
                    fontsize=size, color=color)
    return y + len(lines) * lh

def bullet_wrap(x, y, s, max_width_px, size=8.6, color=BODY, lh=LH_DEFAULT):
    indent_px = max(_text_width_px('\u2022  ', size),
                    _text_width_px('\u2002\u2002 ', size))
    lines = _wrap(s.split(' '), max_width_px - indent_px, size)
    if _probe[0] is not None:          # remove any dangling wrap-measurement
        _probe[0].remove(); _probe[0] = None   # probe before adding real text
    for i, ln in enumerate(lines):
        yi = y + i * lh
        prefix = '\u2022  ' if i == 0 else '\u2002\u2002 '
        ax.text(x, yy(yi), mt(prefix + ln), ha='left', va='baseline',
                fontsize=size, color=color)
    return y + len(lines) * lh

def arrow(x1, y1, x2, y2, color=ARR, dashed=False):
    ax.add_patch(FancyArrowPatch(
        (x1, yy(y1)), (x2, yy(y2)), arrowstyle='-|>', mutation_scale=14,
        lw=1.7, color=color, shrinkA=0, shrinkB=0,
        linestyle='--' if dashed else '-'))

def elbow(x1, y1, xc, x2, y2):
    ax.plot([x1, x1, x2], [yy(y1), yy(xc), yy(xc)], color=ARR, lw=1.7,
            solid_capstyle='round')
    ax.add_patch(FancyArrowPatch((x2, yy(xc)), (x2, yy(y2)),
                 arrowstyle='-|>', mutation_scale=14, lw=1.7, color=ARR,
                 shrinkA=0, shrinkB=0))

# ---- registration: record every text added since a marker, unconditionally,
# so overflow is never silently invisible to the checker.
_registry = []
def mark():
    return len(ax.texts)
def register(host, start):
    for t in ax.texts[start:]:
        _registry.append((t, host))

GAP = 30                # vertical gap between stacked boxes
PAD_BOT = 16             # padding below the last content line inside a box

# ---- problem (full width)
py0 = 20
i0 = mark()
head(580, py0+26, 'SED fitting systematically underestimates prestellar core masses',
     12.5, PAL['problem'][1])
yp = bullet_wrap(96, py0+48, 'temperature gradient: warm outer dust dominates the '
     'emission \u2192 fitted T too high, M too low', 960, 9.2)
yp = bullet_wrap(96, yp+4, 'flux loss: the detection-limited footprint misses the '
     'extended, low-surface-brightness envelope', 960, 9.2)
p_bottom = py0 + 92
pbox = box(70, py0, 1020, 92, 'problem')
register(pbox, i0)

# ---- column headers (free-standing labels, not inside any box by design)
col_y = p_bottom + GAP
_i_hdr = mark()
head(288, col_y, 'CALIBRATION  \u00b7  build the corrector once', 11.5, PAL['calib'][1])
head(872, col_y, 'APPLICATION  \u00b7  per observed cloud', 11.5, PAL['apply'][1])
_freestanding_ids = {id(t) for t in ax.texts[_i_hdr:]}

LX, LW, LC = 70, 436, 288
RX, RW, RC = 650, 440, 870
LTXT = LW - 56
RTXT = RW - 60
track_top = col_y + 18

# ================= calibration track =================
y = track_top
b1_top = y
i0 = mark()
head(LC, y+24, 'RT model grid', 11)
yk = body_wrap(LX+20, y+46, 'pressure-truncated Bonnor\u2013Ebert spheres, embedded '
     'and heated self-consistently', LTXT, 8.6)
yk = body_wrap(LX+20, yk+3, 'axes \u03a3_{cloud}, mass, \u03be_{max} \u2192 true M_{BE} '
     'known per node', LTXT, 8.6)
b1_bottom = yk + PAD_BOT
b1 = box(LX, b1_top, LW, b1_bottom - b1_top, 'calib')
register(b1, i0)

y = b1_bottom + GAP
b2_top = y
i0 = mark()
head(LC, y+24, 'Mimic the getsf extraction', 11)
yk = body_wrap(LX+20, y+46, 'background subtraction \u2192 3-band SED fit \u2192 '
     'M_{SED}, T_{SED}, at many truncation radii r', LTXT, 8.6)
b2_bottom = yk + PAD_BOT
b2 = box(LX, b2_top, LW, b2_bottom - b2_top, 'calib')
register(b2, i0)

y = b2_bottom + GAP
b3_top = y
i0 = mark()
head(LC, y+24, 'Self-consistent per-band footprint', 11, PAL['key'][1])
yk = body_wrap(LX+20, y+48, 'R_{\u03bb} solves prof_{\u03bb}(R_{\u03bb})/prof_{\u03bb}(0) '
     '= [prof_{SD}(r)/prof_{SD}(0)]\u00b7c_{\u03bb}', LTXT, 8.4, lh=LH_DEFAULT)
yk = body_wrap(LX+20, yk+3, 'c_{\u03bb}: cloud power spectrum (\u03b3\u2248-2.5) '
     'convolved to each band beam', LTXT, 8.4, lh=LH_DEFAULT)
yk = body_wrap(LX+20, yk+3, 'constant edge S/N + dust-emissivity cancellation fix '
     'the wavelength scaling', LTXT, 8.4, lh=LH_DEFAULT)
yk = body_wrap(LX+20, yk+3, '\u2192 f_{rec,\u03bb} \u2192 recovered mass M_{rec}(r)',
     LTXT, 8.4, lh=LH_DEFAULT)
b3_bottom = yk + PAD_BOT
b3 = box(LX, b3_top, LW, b3_bottom - b3_top, 'key', bold=True)
register(b3, i0)

y = b3_bottom + GAP
b4_top = y
i0 = mark()
head(LC, y+24, 'Corrector  (distance-invariant)', 11)
yk = body_wrap(LX+20, y+46, 'training set (\u03a3, concentration, R_{foot}, '
     'C_{M}=M_{BE}/M_{rec})', LTXT, 8.4, anchor='center', cx=LC)
yk = body_wrap(LX+20, yk+2, 'interpolated on log axes \u2014 no absolute mass or '
     'distance enters', LTXT, 7.8, MUTE, anchor='center', cx=LC, lh=17)
b4_bottom = yk + PAD_BOT
b4 = box(LX, b4_top, LW, b4_bottom - b4_top, 'calib')
register(b4, i0)

# ================= application track =================
y = track_top
b5_top = y
i0 = mark()
head(RC, y+24, 'Real getsf catalog', 11)
yk = body_wrap(RX+24, y+46, 'any cloud, any distance', RTXT, 8.6, anchor='center',
     cx=RC)
b5_bottom = yk + PAD_BOT
b5 = box(RX, b5_top, RW, b5_bottom - b5_top, 'apply')
register(b5, i0)

y = b5_bottom + GAP
b6_top = y
i0 = mark()
head(RC, y+24, 'Measure 3 observables per source', 11)
yr = body_wrap(RX+24, y+48, '\u03a3 = PEAK^{BGF}  (background column)', RTXT, 8.6)
yr = body_wrap(RX+24, yr+2, 'concentration = PEAK^{SRC} / footprint-mean', RTXT, 8.6)
yr = body_wrap(RX+24, yr+2, 'R_{foot} = FOOA / AFWHM', RTXT, 8.6)
yr = body_wrap(RX+24, yr+4, 'all ratios or surface brightnesses \u2014 distance '
     'cancels', RTXT, 7.8, MUTE, lh=17)
b6_bottom = yr + PAD_BOT
b6 = box(RX, b6_top, RW, b6_bottom - b6_top, 'apply')
register(b6, i0)

# ================= convergence: look up corrector =================
y = max(b4_bottom, b6_bottom) + GAP + 20
MX, MW, MC = 300, 560, 580
i0 = mark()
head(MC, y+26, 'Look up C_{M}( \u03a3, concentration, R_{foot} )', 11.5, PAL['key'][1])
body_wrap(MX+24, y+48, 'same corrector for every cloud and distance', MW-48, 8.6,
     anchor='center', cx=MC)
bM = box(MX, y, MW, 66, 'key', bold=True)
register(bM, i0)
bM_bottom = y + 66

# ================= output =================
y = bM_bottom + GAP
OX, OW, OC = 246, 668, 580
o_top = y
i0 = mark()
head(OC, y+26, 'Corrected mass    M_{corrected}  =  M_{reported} \u00d7 C_{M}', 12.5)
yo = body_wrap(OX+24, y+52, 'C_{M} = F_{flux} \u00d7 F_{bg} \u00d7 F_{temp}', OW-48,
     10.0, INK, anchor='center', cx=OC)
yo = body_wrap(OX+24, yo+4, 'flux recovery  \u00d7  background subtraction  \u00d7  '
     'temperature-gradient bias', OW-48, 8.2, MUTE, anchor='center', cx=OC)
yo = body_wrap(OX+24, yo+4, 'the same corrector applies at any distance \u2014 no '
     're-calibration per cloud', OW-48, 8.2, MUTE, anchor='center', cx=OC)
o_bottom = yo + PAD_BOT
bO = box(OX, o_top, OW, o_bottom - o_top, 'output', bold=True)
register(bO, i0)

# ================= validation =================
y = o_bottom + GAP
VX, VW, VC = 246, 668, 580
v_top = y
i0 = mark()
head(VC, y+24, 'Validation (independent of the calibration)', 10.5, PAL['valid'][1])
yv = bullet_wrap(VX+26, y+47, 'injection\u2013recovery: 108 synthetic cores, '
     'extracted with getsf \u2192 median M/M_{true} = 0.97, 82% within 2\u00d7',
     VW-140, 8.6)
yv = bullet_wrap(VX+26, yv+4, 'independent subdivision grid: 873 withheld samples '
     '\u2192 median M/M_{BE} = 1.05, 90% within 2\u00d7', VW-140, 8.6)
yv = body_wrap(VX+26, yv+6, 'neither test used to calibrate the footprint model '
     'above', VW-52, 7.8, MUTE, lh=17)
v_bottom = yv + PAD_BOT
bV = box(VX, v_top, VW, v_bottom - v_top, 'valid', dashed=True)
register(bV, i0)

# ---- connectors
arrow(288, p_bottom, 288, col_y-14)
arrow(872, p_bottom, 872, col_y-14)
arrow(288, b1_bottom, 288, b1_bottom+GAP)
arrow(288, b2_bottom, 288, b2_bottom+GAP)
arrow(288, b3_bottom, 288, b3_bottom+GAP)
arrow(872, b5_bottom, 872, b5_bottom+GAP)
mid_y = max(b4_bottom, b6_bottom) + GAP//2 + 6
elbow(288, b4_bottom, mid_y, MC-31, bM_bottom-66)
elbow(872, b6_bottom, mid_y, MC+31, bM_bottom-66)
arrow(580, bM_bottom, 580, bM_bottom+GAP)
arrow(580, o_bottom, 580, o_bottom+GAP, PAL['valid'][1], dashed=True)

# ---- shrink canvas to fit content exactly, with margin, then re-verify
H_needed = int(v_bottom + 24)
if H_needed != H:
    print('adjusting canvas H %d -> %d to match content' % (H, H_needed))

# ---- verify no text overflows its host box before saving
if _probe[0] is not None:
    _probe[0].remove()
fig.canvas.draw()
bad = []
for t, hostbox in _registry:
    tb = t.get_window_extent(_renderer)
    bb = hostbox.get_window_extent(_renderer)
    if tb.x0 < bb.x0 - 2 or tb.x1 > bb.x1 + 2 or tb.y0 < bb.y0 - 2 or tb.y1 > bb.y1 + 2:
        bad.append((t.get_text()[:50],
                    tb.x0-bb.x0, bb.x1-tb.x1, tb.y0-bb.y0, bb.y1-tb.y1))
if bad:
    print('OVERFLOW DETECTED -- fix before using this figure:')
    for txt, l, r, bot, top in bad:
        print('  %-52s  left=%.0f right=%.0f bottom=%.0f top=%.0f (negative = overflow)'
              % (txt, l, r, bot, top))
    raise SystemExit('aborting: %d overflowing text block(s)' % len(bad))

_checked_ids = {id(t) for t, _ in _registry} | _freestanding_ids
_unchecked = [t for t in ax.texts if id(t) not in _checked_ids]
assert not _unchecked, \
    'unchecked text object(s), not registered to any box: %r' \
    % [t.get_text()[:40] for t in _unchecked]

overlaps = []
all_texts = list(ax.texts)
for i in range(len(all_texts)):
    for j in range(i + 1, len(all_texts)):
        a = all_texts[i].get_window_extent(_renderer)
        b = all_texts[j].get_window_extent(_renderer)
        if a.x0 < b.x1 - 2 and b.x0 < a.x1 - 2 and a.y0 < b.y1 - 2 and b.y0 < a.y1 - 2:
            overlaps.append((all_texts[i].get_text()[:40], all_texts[j].get_text()[:40]))
if overlaps:
    print('TEXT-TEXT OVERLAP DETECTED -- fix before using this figure:')
    for s1, s2 in overlaps:
        print('  %-42s  overlaps  %-42s' % (s1, s2))
    raise SystemExit('aborting: %d overlapping text pair(s)' % len(overlaps))

print('layout check passed: no overflow, no overlap (%d/%d text blocks checked)'
      % (len(_registry), len(ax.texts)))

fig.savefig('fig_method_flowchart.pdf',
            bbox_inches='tight')
fig.savefig('fig_method_flowchart.png', dpi=150,
            bbox_inches='tight')
print('wrote pdf + png')
