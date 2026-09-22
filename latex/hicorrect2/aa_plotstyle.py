"""Common plotting style for the figures of this paper.

Ticks point inward, appear on all four sides, and are longer than the
matplotlib default.  Import and call apply() before creating any figure.
"""
import matplotlib as mpl
from matplotlib.ticker import NullFormatter


def apply():
    mpl.rcParams.update({
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
        'xtick.major.size': 5.5,
        'ytick.major.size': 5.5,
        'xtick.minor.size': 3.0,
        'ytick.minor.size': 3.0,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.7,
        'ytick.minor.width': 0.7,
        'xtick.minor.visible': True,
        'ytick.minor.visible': True,
        'axes.linewidth': 0.8,
        'font.size': 9.0,
        'legend.frameon': False,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.015,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })


def clean_log(ax, which='both'):
    """Suppress minor tick labels on log axes.

    NOTE ON AXIS LIMITS.  Because the minor labels are removed, the limits of
    every log axis must coincide with tick positions, so that a reader can read
    the tick spacing off the ends of the axis without annotation.  Use round
    values such as 0.5, 1, 2, 3 or powers of ten, never arbitrary ones such as
    0.55 or 1.8.

    On a log axis spanning less than a decade matplotlib labels the minor
    ticks as 9 x 10^-1, 7 x 10^-1 and so on.  The decade labels alone are
    unambiguous, so the minor labels are removed as clutter.  The minor ticks
    themselves are kept.
    """
    if which in ('x', 'both'):
        ax.xaxis.set_minor_formatter(NullFormatter())
    if which in ('y', 'both'):
        ax.yaxis.set_minor_formatter(NullFormatter())


def no_top_minor(ax):
    """Disable minor ticks on the primary axes where a secondary axis with a
    different scale is drawn on top, so that the two sets do not overlap."""
    ax.tick_params(axis='x', which='minor', top=False)
    ax.tick_params(axis='y', which='minor', right=False)


def panel_label(ax, text, loc='upper left', fontsize=8.0, pad=0.03):
    """Place a panel identifier and title inside the axes rather than above
    them, which saves the vertical space a title would occupy."""
    xy = {'upper left': (pad, 1 - pad), 'upper right': (1 - pad, 1 - pad),
          'lower left': (pad, pad), 'lower right': (1 - pad, pad)}[loc]
    ha = 'left' if 'left' in loc else 'right'
    va = 'top' if 'upper' in loc else 'bottom'
    return ax.text(xy[0], xy[1], text, transform=ax.transAxes, ha=ha, va=va,
                   fontsize=fontsize, linespacing=1.25)
