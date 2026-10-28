"""Shared ACS / Nano Letters figure-formatting constants and helpers.

Every number below comes from Appendix 2 "Preparing Graphics" of the ACS author
guidelines (ACS_author_guide.pdf): single-column graphics up to 240 pt
(3.33 in); double-column graphics 300-504 pt (4.167-7.0 in); maximum depth
660 pt (9.167 in) including caption; lettering no smaller than 4.5 pt in the
final published format, in Helvetica or Arial; lines no thinner than 0.5 pt.

Compliance is defined once, here, and every plotting script in the repo imports
it, so the whole submission uses one typeface, one font size, one set of line
weights and one pair of column widths.

Two deliberate choices go beyond the bare minimums:

* FONT_SIZE is 8 pt, not the 4.5 pt floor.  8 pt is the largest size that still
  lets the WAXS inset and its peak labels stay uncollided inside single-column
  Figure 2, which is what sets the ceiling.
* Figures are saved at exactly SINGLE_COL or DOUBLE_COL wide and embedded in
  the LaTeX at natural size, so 8 pt drawn is 8 pt printed.  Nothing is scaled
  at typesetting time and no bounding box is trimmed.
"""

import inspect
import os

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter

# --- geometry (in) ---
SINGLE_COL = 3.33     # 240 pt, the ACS single-column maximum
DOUBLE_COL = 7.00     # 504 pt, the ACS double-column maximum
MAX_DEPTH = 9.167     # 660 pt, including the caption

# --- type ---
FONT_SIZE = 8         # pt, every character in every figure

# --- line weights (pt); 0.5 is the ACS floor, nothing may go below it ---
LW_THIN = 0.5         # spines, ticks, grids, error bars, guide lines
LW_DATA = 1.0         # data and model curves
MEW = 0.6             # open-marker edge width

# --- markers (pt) ---
MS = 3.5              # the default: I(Q), g2
MS_DENSE = 2.4        # profiles carrying several hundred points per curve, where
                      # MS-sized markers with a 0.6 pt edge merge into a solid
                      # band.  Paired with a LW_THIN (0.5 pt) edge, the marker
                      # outlines stay open and overlapping curves stay readable.
MS_SPARSE = 4.5       # sparse points whose SHAPE encodes a variable, and which
                      # therefore have to be big enough to tell o/s/^/D/v apart


def apply_style():
    """Set every rcParam that controls typeface, size and line weight."""
    plt.rcParams.update({
        # Arial for text and for mathtext, so an axis label and the symbols
        # inside it are the same typeface.
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold', 'mathtext.default': 'it',
        # one size everywhere
        'font.size': FONT_SIZE, 'axes.labelsize': FONT_SIZE, 'axes.titlesize': FONT_SIZE,
        'xtick.labelsize': FONT_SIZE, 'ytick.labelsize': FONT_SIZE,
        'legend.fontsize': FONT_SIZE, 'legend.title_fontsize': FONT_SIZE,
        # opaque white legend boxes with a 0.5 pt black frame
        'legend.frameon': True, 'legend.facecolor': 'white', 'legend.edgecolor': 'black',
        'legend.framealpha': 1.0, 'legend.fancybox': False,
        'legend.handletextpad': 0.4, 'legend.labelspacing': 0.25,
        'legend.borderpad': 0.35, 'legend.columnspacing': 1.0,
        'patch.linewidth': LW_THIN,
        # lines and markers
        'lines.linewidth': LW_DATA, 'lines.markersize': MS,
        'lines.markeredgewidth': MEW,
        'axes.linewidth': LW_THIN,
        'xtick.major.width': LW_THIN, 'ytick.major.width': LW_THIN,
        'xtick.minor.width': LW_THIN, 'ytick.minor.width': LW_THIN,
        'grid.linewidth': LW_THIN,
        # TrueType, never Type 3 -- ACS production rejects Type 3 subsets
        'pdf.fonttype': 42, 'ps.fonttype': 42,
        'savefig.facecolor': 'white', 'savefig.transparent': False,
    })


def add_minor_grid(ax):
    """Grid on major and minor ticks, at the 0.5 pt minimum and pale enough to
    stay behind the data.

    A minor grid line must never land halfway between two labelled major ticks:
    a line at 0.015 between 0.01 and 0.02 reads as a value the axis never names.
    The linear axes here are safe, because matplotlib subdivides them by 5.  Log
    axes need care, so every hand-written locator in this repository lists only
    integer multiples of a decade, never 1.5 or 2.5 times one.
    """
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.grid(which='major', ls='-', lw=LW_THIN, color='0.85')
    ax.grid(which='minor', ls='-', lw=LW_THIN, color='0.93')


def q_log_ticks(ax, axis='x'):
    """Label 0.004, 0.01 and 0.03 on a log Q axis.

    The 8-ID Q range spans only about one decade (0.0033-0.034 Ainv), so
    matplotlib's automatic log locator finds a single decade inside it and
    labels one tick, 10^-2.  A reader then has no way to read a value off the
    axis.  The three named ticks below are integer multiples of a decade, so
    they satisfy the add_minor_grid() rule about half-decade positions, and the
    same three are used on every Q axis in the paper (Figures 3a, S4b and
    S7a-S10a) so the profiles can be compared tick for tick.
    """
    if axis == 'x':
        a = ax.xaxis
    else:
        a = ax.yaxis
    a.set_major_locator(FixedLocator([4e-3, 1e-2, 3e-2]))
    a.set_major_formatter(FixedFormatter(['0.004', '0.01', '0.03']))
    a.set_minor_locator(FixedLocator([3e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3,
                                      2e-2, 4e-2, 5e-2]))
    a.set_minor_formatter(NullFormatter())


def label_panels(axes, dx=-30, dy=5, start=0):
    """Bold (a), (b), (c)... just outside the top-left corner of each axes.

    The offset is in POINTS, not axes fractions, so every panel label in the
    paper sits the same distance out no matter how wide its panel is.  As a
    fraction, a wide panel (the 2x2 SI grids) pushed the label off the canvas
    while a narrow one (a 1-of-3 panel) tucked it against the tick labels.

    dy clears the TOPMOST Y TICK LABEL, which is centred on the top spine and
    therefore reaches about half a cap height (2.9 pt at 8 pt Arial) above it.
    Panels whose left margin is narrow -- Figure S7, where the y tick labels are
    short and tight_layout leaves only ~17 pt of margin -- put the label and
    that tick label at the same x, so the two collided at the original dy = 2.
    Raising it to 5 leaves 2.1 pt of clear space in the worst case and changes
    nothing else: no figure in the paper draws anything else in that corner.

    start skips that many letters, so a figure whose panels are created in more
    than one call can label the second group (d), (e), ... rather than starting
    over at (a).  A panel row that shares its y axis has no tick labels to clear
    on its second and later panels, and there dx = -30 would place the label
    nearer the previous panel than its own; pass a smaller dx for that row.
    """
    for ax, letter in zip(axes, 'abcdefghijklmnopqrstuvwxyz'[start:]):
        ax.annotate(f'({letter})', xy=(0, 1), xycoords='axes fraction',
                    xytext=(dx, dy), textcoords='offset points',
                    fontsize=FONT_SIZE, fontweight='bold', va='bottom', ha='left')


def save_fig(fig, filename):
    """Save as vector PDF at exactly the requested figure size, NEXT TO THE SCRIPT.

    A bare filename is resolved against the directory of the calling script, not
    the current working directory.  Every figure therefore lands beside the code
    that produced it however the script was launched -- `make`, an editor, or a
    shell sitting somewhere else.  Running a script from the wrong directory used
    to scatter stray duplicate PDFs around the tree, which then went stale while
    the real ones were regenerated.  Pass an absolute path to override.

    Deliberately no bbox_inches='tight': the media box has to stay exactly
    SINGLE_COL or DOUBLE_COL wide so the journal reproduces the graphic at
    100 % and the 8 pt lettering prints at 8 pt.  A trimmed box would be
    rescaled at typesetting and the compliance arithmetic would no longer hold.
    Call fig.tight_layout() (or set constrained layout) before this.
    """
    if not os.path.isabs(filename):
        caller = inspect.stack()[1].filename
        filename = os.path.join(os.path.dirname(os.path.abspath(caller)), filename)
    w, h = fig.get_size_inches()
    fig.savefig(filename, format='pdf')
    if abs(w - SINGLE_COL) < 0.01:
        ok_w = 'single'
    elif abs(w - DOUBLE_COL) < 0.01:
        ok_w = 'double'
    else:
        ok_w = 'OFF-SPEC'
    if h <= MAX_DEPTH:
        ok_h = 'ok'
    else:
        ok_h = 'OVER 9.167 in'
    print(f'wrote {os.path.basename(filename)}  {w:.3f} x {h:.3f} in '
          f'({72 * w:.0f} x {72 * h:.0f} pt)  width={ok_w}  depth={ok_h}')
