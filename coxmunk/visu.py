''' Simple sunglint computation for differentr wave slope statistics based on Cox-Munk model

Usage:
  coxmunk <sza> <wind_speed> [--stats <stats>] [--wind_azi <wind_azi>] [--figname <figname>]
  coxmunk -h | --help
  coxmunk -v | --version

Options:
  -h --help        Show this screen.
  -v --version     Show version.

  <sza>     Solar zenith angle in deg.
  <wind_speed>  Wind speed in m/s.
  --stats stats  Statistics used in Cox and Munk model:
                    - `cm_iso`, original isotropic Cox Munk statistics
                    - `cm_dir`, historical values from directional COX MUNK
                    - `bh2006`, reassessment from Breon Henriot 2006 JGR
                    [default: cm_iso]
  --wind_azi wind_azi  Azimuth of wind direction from the Principal plane
                        [default: 0]
  --figname figname   Path to save figure (ex: ./illustration/coxmunk_fig_39_14.2_bh2006_75.png)

'''

from docopt import docopt

import numpy as np
import cmocean as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.font_manager import FontProperties

import coxmunk.coxmunk as coxmunk


class plot():
    def __init__(self, font=None):
        # ----------------------------
        # set plotting styles

        plt.rcParams.update({'font.size': 18})

        if font == None:
            font = {'family': 'serif',
                    'color': 'black',
                    'weight': 'normal',
                    'size': 16,
                    }
        self.font = font

    def label_polplot(self, ax, yticks=[20., 40., 60.], ylabels=['$20^{\circ}$', '$40^{\circ}$', '$60^{\circ}$']):

        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels)
        ax.set_theta_zero_location("S")
        ax.set_theta_direction(-1)
        return

    def add_polplot(self, ax, r, theta, values, title="", scale=True, nlayers=50, cmap=cm.cm.delta, colfmt='%0.1e',
                    pad=0.1, fraction=0.034, **kwargs):

        theta = np.radians(theta)
        self.label_polplot(ax)
        ax.set_title(title, fontdict=self.font, pad=30)

        min_, max_ = values.min(), values.max()
        if np.abs(min_) < 1e-6:
            min_ = 0
        if 'vmin' in kwargs:
            min_ = kwargs['vmin']
        if 'vmax' in kwargs:
            max_ = kwargs['vmax']

        if max_ - min_ == 0:
            max_ = min_ + 1
        inc = (max_ - min_) / nlayers
        contour = np.arange(min_, max_, inc)
        cax = ax.contourf(theta, r, values, contour, extend='both', cmap=cmap, **kwargs)
        if scale:
            plt.colorbar(cax, ax=ax, format=colfmt, fraction=fraction, pad=pad)

        return cax


def main():
    from coxmunk.__init__ import __package__, __version__

    args = docopt(__doc__, version=__package__ + ' ' + __version__)
    print(args)

    sza = float(args['<sza>'])
    wind = float(args['<wind_speed>'])
    wind_azi = float(args['--wind_azi'])
    stats = args['--stats']
    figname = args['--figname']

    vza = np.linspace(0, 80, 80)
    azi = np.linspace(0, 360, 180)
    Nvza, Nazi = len(vza), len(azi)

    data = np.zeros((Nazi, Nvza))
    for i in range(Nvza):
        for j in range(Nazi):
            data[j, i] = coxmunk.sunglint(sza, vza[i], azi[j], m=1.334).sunglint(wind, wind_azi, stats=stats)

    cmap = plt.cm.gist_stern_r

    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'),figsize=(7, 6))
    plot().add_polplot(ax, vza, azi, data.T, cmap=cmap)
    plt.suptitle(r'Wind speed: {:.1f} m/s; direction: {:.1f} deg.'.format(wind,wind_azi))
    plt.tight_layout(rect=[0.0, 0.0, 0.99, 0.95])
    if figname:
        plt.savefig(figname,dpi=300)
    else:
        plt.show()

if __name__ == "__main__":
    main()
