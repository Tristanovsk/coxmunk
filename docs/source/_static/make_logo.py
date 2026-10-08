"""Generate the coxmunk logo.

The background is a crop of *Sunrise with Sea Monsters* (J. M. W. Turner,
c. 1845, public domain) around the Sun and its glitter on the horizon, faded
into the page colour. In the wordmark, the 'o' is the Sun reflected by the sea:
a vertical gradient from sunlight gold to sea blue.

    python make_logo.py   ->  logo.png, logo_dark.png (800 px), icon.png, icon_dark.png (256 px)
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.font_manager import FontProperties
from matplotlib.ft2font import FT2Font
from matplotlib.image import imread
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib._text_helpers import layout

HERE = Path(__file__).parent
PAINTING = HERE.parents[2] / "illustration" / "kumatage" / "turner_1845_sunrise_with_sea_monsters.jpg"
FONT = "/usr/share/fonts/truetype/lato/Lato-Medium.ttf"
FP = FontProperties(fname=FONT)
SUN = "#F4B942"     # sunlight gold (top of the 'o')
GLINT = "#E8833A"   # low sun, where the light meets the sea
SEA = "#2A6F97"     # sea blue (bottom of the 'o')

# ----------------------------------------------------------------- background
_img = imread(PAINTING)[..., :3].astype(float)
_img = _img / 255 if _img.max() > 1 else _img
_h, _w = _img.shape[:2]
# position of the Sun glow in the painting (image pixels); the crop is computed
# for each canvas so that the Sun sits behind the anchor glyph
SUN_X, SUN_Y = 0.72 * _w, 0.56 * _h


# ----------------------------------------------------------------- colour
def _srgb_to_lin(c):
    c = np.asarray(c, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _lin_to_srgb(c):
    c = np.clip(c, 0, 1)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)


_M1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                [0.2119034982, 0.6806995451, 0.1073969566],
                [0.0883024619, 0.2817188376, 0.6299787005]])
_M2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                [1.9779984951, -2.4285922050, 0.4505937099],
                [0.0259040371, 0.7827717662, -0.8086757660]])


def _to_oklab(hex_):
    return _M2 @ np.cbrt(_M1 @ _srgb_to_lin(to_rgb(hex_)))


def _from_oklab(lab):
    return _lin_to_srgb(np.linalg.solve(_M1, np.linalg.solve(_M2, lab) ** 3))


def gradient(stops, n=256):
    """n colours through the colour stops, interpolated in OKLab (no muddy middle)."""
    labs = np.array([_to_oklab(c) for c in stops])
    pos = np.linspace(0, 1, len(stops))
    t = np.linspace(0, 1, n)
    lab = np.stack([np.interp(t, pos, labs[:, k]) for k in range(3)], axis=1)
    return np.array([_from_oklab(v) for v in lab])


# ----------------------------------------------------------------- glyphs
def glyphs(text, size):
    """[(char, x_offset)] using the font's advances and kerning."""
    font = FT2Font(FONT)
    font.set_size(size, 72)
    return [(it.char, it.x) for it in layout(text, font)]


def _halo(ax, path, halo, width, z):
    if halo is not None:
        ax.add_patch(PathPatch(path, fc=halo, ec=halo, lw=width, joinstyle="round",
                               zorder=z, alpha=0.9))


def draw_sun_o(ax, x, baseline, size, z=3, halo=None, halo_w=10):
    """The 'o' as the Sun over its reflection: gold at the top, sea blue at the bottom.

    Returns the centre of the 'o' (anchor of the Sun of the painting)."""
    o = TextPath((x, baseline), "o", size=size, prop=FP)
    v = np.vstack(o.to_polygons())
    _halo(ax, o, halo, halo_w, z - 1)
    ramp = gradient([SUN, GLINT, SEA], 256)[::-1]   # imshow origin='lower': bottom first
    img = np.repeat(ramp[:, None, :], 2, axis=1)
    im = ax.imshow(img, extent=[v[:, 0].min(), v[:, 0].max(), v[:, 1].min(), v[:, 1].max()],
                   origin="lower", interpolation="bilinear", zorder=z, aspect="auto")
    im.set_clip_path(PathPatch(o, transform=ax.transData))
    return (v[:, 0].min() + v[:, 0].max()) / 2, (v[:, 1].min() + v[:, 1].max()) / 2


def draw_wordmark(ax, size, cx, baseline, ink, z=3, halo=None, halo_w=10):
    """Draw 'coxmunk' centred on cx. Returns the centre of the 'o'."""
    full = np.vstack(TextPath((0, 0), "coxmunk", size=size, prop=FP).to_polygons())
    x0 = cx - (full[:, 0].min() + full[:, 0].max()) / 2
    anchor = None
    for i, (ch, dx) in enumerate(glyphs("coxmunk", size)):
        if ch == "o" and anchor is None:
            anchor = draw_sun_o(ax, x0 + dx, baseline, size, z, halo, halo_w)
        else:
            tp = TextPath((x0 + dx, baseline), ch, size=size, prop=FP)
            _halo(ax, tp, halo, halo_w, z - 1)
            ax.add_patch(PathPatch(tp, fc=ink, ec="none", zorder=z))
    return anchor


# ----------------------------------------------------------------- composition
def draw_background(ax, W, H, dark, centre, strength=None, fade_r=0.62, fade_centre=None):
    """Painting filling the canvas, its Sun at `centre`, faded radially into the page.

    The crop is the largest region of the painting that covers the canvas with
    the Sun at `centre`; the fade is centred on `fade_centre` (canvas centre)."""
    ox, oy = centre
    # image pixels per canvas pixel, limited by the distance of the Sun to the
    # borders of the painting (canvas y axis points up, image rows down)
    k = min(SUN_X / ox, (_w - SUN_X) / (W - ox), SUN_Y / (H - oy), (_h - SUN_Y) / oy)
    x0, x1 = SUN_X - ox * k, SUN_X + (W - ox) * k
    y0, y1 = SUN_Y - (H - oy) * k, SUN_Y + oy * k
    crop = _img[int(y0):int(np.ceil(y1)), int(x0):int(np.ceil(x1))]
    h, w = crop.shape[:2]
    fx, fy = fade_centre if fade_centre is not None else (W / 2, H / 2)
    yy, xx = np.mgrid[0:h, 0:w]
    gx, gy = xx / w * W, H - yy / h * H        # canvas coordinates of the pixels
    rr = np.hypot((gx - fx) / (fade_r * W), (gy - fy) / (fade_r * H))
    fade = np.clip(1.15 - rr, 0, 1) ** 1.3
    strength = strength or (0.85 if dark else 0.95)
    rgba = np.zeros((h, w, 4))
    # dark mode: dim the pale painting so that it glows on the night background
    rgba[..., :3] = crop * (0.55 if dark else 1.0)
    rgba[..., 3] = fade * strength
    ax.imshow(rgba, extent=[0, W, 0, H], origin="upper", interpolation="bilinear", zorder=1)


def canvas(W, H, bg):
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    fig.patch.set_facecolor(bg)
    return fig, ax


def _bg(dark):
    return "#0F1720" if dark else "#FFFFFF"


def render_logo(name, dark=False, W=1600, H=1600, size=330, OUT=800):
    bg = _bg(dark)
    fig, ax = canvas(W, H, bg)
    xh = TextPath((0, 0), "x", size=size, prop=FP).vertices[:, 1].max()
    # no halo: on the smooth haze of the painting it would read as an outline
    anchor = draw_wordmark(ax, size, W / 2, H / 2 - xh / 2, "#F2F2F2" if dark else "#151515")
    draw_background(ax, W, H, dark, centre=anchor)
    # PNG only: an SVG would embed the painting as a (heavy) raster anyway
    fig.savefig(HERE / f"{name}.png", facecolor=bg, dpi=100 * OUT / fig.get_figwidth() / 100)
    plt.close(fig)


def draw_reflection(ax, x, baseline, size, gap, squash=0.6, bands=10, amp=0.05, z=3):
    """Reflection of the Sun 'o' on the sea: mirrored, squashed and broken into
    horizontal bands shifted sideways by the waves, fading toward the viewer."""
    o = TextPath((x, baseline), "o", size=size, prop=FP)
    v = o.vertices.copy()
    y_bot = np.vstack(o.to_polygons())[:, 1].min()
    y_hor = y_bot - gap                                  # horizon, just below the Sun
    v[:, 1] = y_hor - (v[:, 1] - y_bot) * squash          # mirror image
    mv = np.vstack(matplotlib.path.Path(v, o.codes).to_polygons())
    x_min, x_max, y_min = mv[:, 0].min(), mv[:, 0].max(), mv[:, 1].min()
    width = x_max - x_min
    ramp = gradient([SUN, GLINT, SEA], 256)               # mirrored: blue at the top
    img = np.repeat(ramp[:, None, :], 2, axis=1)
    edges = np.linspace(y_hor, y_min, bands + 1)
    for k in range(bands):
        t = k / (bands - 1)
        dx = amp * width * (1 + 3 * t) * np.sin(2.3 * k + 0.5)
        shifted = matplotlib.path.Path(v + [dx, 0], o.codes)
        # thin gap between bands, like the troughs between wave crests
        y_hi = edges[k] - 0.15 * (edges[k] - edges[k + 1])
        im = ax.imshow(img, extent=[x_min + dx, x_max + dx, y_min, y_hor], origin="lower",
                       interpolation="bilinear", zorder=z, aspect="auto",
                       alpha=0.75 * (1 - 0.6 * t))
        im.set_clip_path(_intersect(shifted, edges[k + 1], y_hi), ax.transData)


def _intersect(glyph, ylo, yhi):
    """Glyph path clipped to the horizontal band ylo <= y <= yhi."""
    from matplotlib.path import Path as MPath
    polys = []
    for poly in glyph.to_polygons():
        clipped = _clip_poly_y(poly, ylo, yhi)
        if len(clipped) >= 3:
            polys.append(clipped)
    verts, codes = [], []
    for p in polys:
        verts += list(p) + [p[0]]
        codes += [MPath.MOVETO] + [MPath.LINETO] * (len(p) - 1) + [MPath.CLOSEPOLY]
    return MPath(verts, codes) if verts else MPath([[0, 0]], [MPath.MOVETO])


def _clip_poly_y(poly, ylo, yhi):
    """Sutherland-Hodgman clipping of a polygon to ylo <= y <= yhi."""
    def clip(points, inside, cross):
        out = []
        for i in range(len(points)):
            cur, prev = points[i], points[i - 1]
            if inside(cur):
                if not inside(prev):
                    out.append(cross(prev, cur))
                out.append(cur)
            elif inside(prev):
                out.append(cross(prev, cur))
        return out

    def at_y(y0):
        return lambda p, q: p + (q - p) * (y0 - p[1]) / (q[1] - p[1])

    pts = [np.asarray(p, float) for p in poly]
    pts = clip(pts, lambda p: p[1] >= ylo, at_y(ylo))
    if pts:
        pts = clip(pts, lambda p: p[1] <= yhi, at_y(yhi))
    return np.array(pts)


def render_icon(name, dark=False, S=1024, size=600, OUT=256):
    """Square icon: the Sun 'o' over its rippled reflection, on the painting."""
    bg = _bg(dark)
    fig, ax = canvas(S, S, bg)
    o = np.vstack(TextPath((0, 0), "o", size=size, prop=FP).to_polygons())
    x = S / 2 - (o[:, 0].min() + o[:, 0].max()) / 2
    o_h = o[:, 1].max() - o[:, 1].min()
    baseline = 0.66 * S - (o[:, 1].min() + o[:, 1].max()) / 2
    anchor = draw_sun_o(ax, x, baseline, size)
    draw_reflection(ax, x, baseline, size, gap=0.08 * o_h)
    draw_background(ax, S, S, dark, centre=anchor, fade_r=0.75,
                    strength=0.9 if dark else 0.85)
    # PNG only: an SVG would embed the painting as a (heavy) raster anyway
    fig.savefig(HERE / f"{name}.png", facecolor=bg, dpi=100 * OUT / fig.get_figwidth() / 100)
    plt.close(fig)


if __name__ == "__main__":
    render_logo("logo")
    render_logo("logo_dark", dark=True)
    render_icon("icon")
    render_icon("icon_dark", dark=True)
    print("written:", sorted(p.name for p in HERE.glob("*.png")))
