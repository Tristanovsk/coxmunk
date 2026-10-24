"""Draw the sun-facet-sensor geometry with the notation of algorithms.rst.

The x axis points toward the Sun azimuth (Sun azimuth = 0), the relative azimuth
phi is counted from x to the sensor azimuth, and n is the normal of the
reflecting facet (bisector of the Sun and sensor directions).

Two figures are written:

* sunglint_geometry.png: horizontal reference plane (algorithms.rst),
* ../research/tilted_reference_plane.png: horizontal plane (left) and reference
  plane tilted by beta toward the azimuth phi_p (right), with the local angles
  of the Sun and of the sensor measured from the normal n_p of the plane.

    python illustration/geometry/plot_geometry.py
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Liberation Serif', 'Times', 'DejaVu Serif'],
    'mathtext.fontset': 'stix',
})

sza, vza, azi = np.radians([40, 55, 135])
beta, phip = np.radians([25, 60])

# oblique projection: x comes toward the viewer, to the lower left
X2D = np.array([-0.55, -0.42])
O = np.zeros(3)
L = 1.25
zen = np.array([0, 0, 1.])
sun = np.array([np.sin(sza), 0, np.cos(sza)])
sen = np.array([np.sin(vza) * np.cos(azi), np.sin(vza) * np.sin(azi), np.cos(vza)])
kw_arrow = dict(arrowstyle='-|>', mutation_scale=14, lw=1.4, color='k')


def proj(p):
    p = np.asarray(p, float)
    return np.array([p[1] + p[0] * X2D[0], p[2] + p[0] * X2D[1]])


def unit(v):
    return v / np.linalg.norm(v)


def arc(a, b, r, n=60):
    """Great-circle arc of radius r from direction a to direction b."""
    a, b = unit(a), unit(b)
    ang = np.arccos(np.clip(a @ b, -1, 1))
    t = np.linspace(0, 1, n)[:, None]
    pts = (np.sin((1 - t) * ang) * a + np.sin(t * ang) * b) / np.sin(ang)
    return np.array([proj(r * p) for p in pts])


def arrow(ax, p0, p1, **kw):
    ax.annotate('', xy=proj(p1), xytext=proj(p0), arrowprops={**kw_arrow, **kw},
                annotation_clip=False)


def label(ax, p, text, dx=0, dy=0, **kw):
    x, y = proj(p)
    ax.text(x + dx, y + dy, text, ha='center', va='center', **{'fontsize': 15, **kw})


def angle(ax, a, b, r, name, off, ls='-'):
    """Arc with an arrow head from direction a to direction b, and its label."""
    pts = arc(a, b, r)
    ax.plot(*pts.T, ls, color='k', lw=1.2)
    ax.annotate('', xy=pts[-1], xytext=pts[-4], arrowprops=kw_arrow)
    mid = pts[len(pts) // 2]
    ax.text(mid[0] + off[0], mid[1] + off[1], name, fontsize=15, ha='center', va='center')


def azimuth(ax, phi, r, name, dy=-0.08):
    """Horizontal arc from the x axis (Sun azimuth) to the azimuth phi."""
    t = np.linspace(0, phi, 80)
    pts = np.array([proj(r * np.array([np.cos(a), np.sin(a), 0])) for a in t])
    ax.plot(*pts.T, color='k', lw=1.2)
    ax.annotate('', xy=pts[-1], xytext=pts[-4], arrowprops=kw_arrow)
    label(ax, r * np.array([np.cos(phi / 2), np.sin(phi / 2), 0]), name, 0.0, dy)


def plane(ax, normal, half, **kw):
    """Square of half-size `half` centred at the origin, perpendicular to normal."""
    e1 = unit(np.cross(zen, normal)) if abs(normal @ zen) < 1 else np.array([0, 1., 0])
    e2 = np.cross(normal, e1)
    corners = [half[0] * s1 * e1 + half[1] * s2 * e2 for s1, s2 in [(-1, -1), (1, -1), (1, 1), (-1, 1)]]
    ax.add_patch(Polygon([proj(c) for c in corners], closed=True, **kw))


def normal(ax, n, name):
    ax.plot(*np.array([proj(O), proj(1.05 * n)]).T, ':', color='k', lw=1.4)
    arrow(ax, 0.95 * n, 1.05 * n)
    label(ax, 1.05 * n, name, -0.05, 0.06)


def base(ax, title):
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_xlim(-0.82, 1.38)
    ax.set_ylim(-0.66, 1.54)
    ax.text(0.28, 1.45, title, fontsize=17, ha='center', va='center')
    arrow(ax, O, [1.25, 0, 0]); label(ax, [1.25, 0, 0], r'$x$', -0.06, -0.05)
    arrow(ax, O, [0, 1.15, 0]); label(ax, [0, 1.15, 0], r'$y$', 0.02, -0.08)
    arrow(ax, O, [0, 0, 1.35]); label(ax, [0, 0, 1.35], r'$z$', 0.07, 0)


def rays(ax, feet=True):
    """Incident sunlight, reflected ray toward the sensor, ground projections."""
    arrow(ax, L * sun, 0.45 * sun)
    ax.plot(*np.array([proj(0.47 * sun), proj(O)]).T, 'k', lw=1.4)
    label(ax, L * sun, 'Incident\nsunlight', -0.2, -0.02, fontsize=13)
    arrow(ax, O, L * sen)
    xs, ys = proj(L * sen)
    ax.add_patch(Polygon([[xs - 0.035, ys - 0.02], [xs + 0.035, ys - 0.02],
                          [xs + 0.035, ys + 0.045], [xs - 0.035, ys + 0.045]],
                         closed=True, fc='w', ec='k', lw=1.8, zorder=5))
    label(ax, L * sen, 'Sensor', 0.02, 0.15, fontsize=13)
    if feet:
        for d in (sun, sen):
            foot = L * np.array([d[0], d[1], 0])
            ax.plot(*np.array([proj(L * d), proj(foot)]).T, '--', color='0.4', lw=1)
            ax.plot(*np.array([proj(O), proj(foot)]).T, '--', color='0.4', lw=1)


def horizontal_panel(ax, title=None):
    """Reflecting facet over a horizontal reference plane."""
    nrm = unit(sun + sen)
    if title:
        base(ax, title)
    else:
        base(ax, '')
    plane(ax, nrm, (0.16, 0.21), fc='0.85', ec='0.3', lw=1)
    normal(ax, nrm, r'$\mathbf{n}$')
    rays(ax)
    azimuth(ax, azi, L * 0.62, r'$\phi$')
    angle(ax, zen, sun, 0.95, r'$\theta_s$', (-0.08, 0.07))
    angle(ax, zen, sen, 0.85, r'$\theta_v$', (0.16, 0.07))
    angle(ax, zen, nrm, 0.6, r'$\theta_n$', (0.02, 0.1))
    for d in (sun, sen):
        pts = arc(nrm, d, 0.38)
        ax.plot(*pts.T, ':', color='k', lw=1.4)
        ax.text(*proj(0.5 * unit(nrm + 2 * d)), r'$\omega$', fontsize=15, ha='center', va='center')


def tilted_panel(ax, title):
    """Reference plane tilted by beta, its normal leaning toward phi_p."""
    n_p = np.array([np.sin(beta) * np.cos(phip), np.sin(beta) * np.sin(phip), np.cos(beta)])
    base(ax, title)
    # horizontal square (dashed) and the same square rotated onto the tilted plane
    c, s = np.cos(phip), np.sin(phip)
    Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    Ry = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
    M = Rz @ Ry @ Rz.T
    square = 0.45 * np.array([[1, 1, 0], [-1, 1, 0], [-1, -1, 0], [1, -1, 0]])
    ax.add_patch(Polygon([proj(q) for q in square], closed=True, fc='none', ec='0.5', lw=1, ls='--'))
    ax.add_patch(Polygon([proj(M @ q) for q in square], closed=True,
                         fc='#9ec5ea', ec='#2f6ea5', lw=1.2, alpha=0.6))
    normal(ax, n_p, r'$\mathbf{n}_p$')
    rays(ax, feet=False)
    # azimuth of the normal of the plane
    foot = np.array([n_p[0], n_p[1], 0]) / np.sin(beta) * 0.7
    ax.plot(*np.array([proj(1.05 * n_p), proj(1.05 * np.array([n_p[0], n_p[1], 0]))]).T,
            '--', color='0.4', lw=1)
    ax.plot(*np.array([proj(O), proj(foot)]).T, '--', color='0.4', lw=1)
    azimuth(ax, phip, 0.6, r'$\phi_p$', dy=-0.1)
    angle(ax, zen, n_p, 0.5, r'$\beta$', (0.01, 0.08))
    angle(ax, n_p, sun, 1.0, r"$\theta_s'$", (-0.1, 0.05))
    angle(ax, n_p, sen, 0.9, r"$\theta_v'$", (0.12, 0.09))


here = Path(__file__).parent

fig, ax = plt.subplots(figsize=(6, 6))
fig.subplots_adjust(0, 0, 1, 1)
horizontal_panel(ax)
fig.savefig(here / 'sunglint_geometry.png', dpi=200, facecolor='white')

fig, axs = plt.subplots(1, 2, figsize=(12, 6))
fig.subplots_adjust(0, 0, 1, 1, wspace=0)
horizontal_panel(axs[0], 'Horizontal reference plane')
tilted_panel(axs[1], 'Tilted reference plane')
fig.savefig(here.parent / 'research' / 'tilted_reference_plane.png', dpi=200, facecolor='white')
