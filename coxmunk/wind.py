# coding=utf-8
r'''
Sea surface wind speed from multidirectional and polarized sunglint.

Implementation of the forward model and of the inverse method of
[HarmelChami2012]_, and of the sunglint direction test of [HarmelChami2013]_.

The sunglint and whitecap contribution to the signal is modelled as

.. math::

    \mathbf{S}_{g+wc} = (1 - f_f)\, T\, \mathbf{S}_g + f_f\, t_d\, \mathbf{S}_{wc}

where :math:`\mathbf{S}_g` is the Cox-Munk sunglint Stokes vector
(:class:`coxmunk.coxmunk.sunglint`), :math:`T` and :math:`t_d` are the direct
and diffuse atmospheric transmittances, :math:`f_f` is the whitecap coverage
(:func:`whitecap_coverage`) and :math:`\mathbf{S}_{wc} = [R_{wc}, 0, 0]` is the
unpolarized foam reflectance.

The wind speed :math:`x` is retrieved by minimizing

.. math::

    \chi^2(x) = \left\| \mathbf{S}^*_{g+wc} - \mathbf{S}_{g+wc}(x) \right\|^2

over all the viewing directions with the Levenberg-Marquardt method
(:func:`retrieve_wind_speed`).

References
----------
.. [HarmelChami2012] Harmel, T. and Chami, M. (2012). Determination of sea surface wind
   speed using the polarimetric and multidirectional properties of satellite
   measurements in visible bands. *Geophys. Res. Lett.*, 39, L19611,
   doi:10.1029/2012GL053508.
.. [HarmelChami2013] Harmel, T. and Chami, M. (2013). Estimation of the sunglint
   radiance field from optical satellite imagery over open ocean: multidirectional
   approach and polarization aspects. *J. Geophys. Res. Oceans*, 118,
   doi:10.1029/2012JC008221.
.. [Monahan1980] Monahan, E. C. and O'Muircheartaigh, I. (1980). Optimal power-law
   description of oceanic whitecap coverage dependence on wind speed.
   *J. Phys. Oceanogr.*, 10, 2094-2099.
'''
from collections import namedtuple

import numpy as np
from scipy import optimize

from .coxmunk import sunglint

#: Stokes components that can be used in the retrieval
COMPONENTS = ('I', 'Q', 'U')

WindRetrieval = namedtuple('WindRetrieval', ['ws', 'sigma', 'flag', 'cost'])
WindRetrieval.__doc__ = 'Result of :func:`retrieve_wind_speed`.'
WindRetrieval.ws.__doc__ = 'Retrieved wind speed (m/s), NaN if the pixel is flagged.'
WindRetrieval.sigma.__doc__ = 'Uncertainty on the wind speed (m/s), see :func:`retrieve_wind_speed`.'
WindRetrieval.flag.__doc__ = 'True when the measurements are not informative on wind speed.'
WindRetrieval.cost.__doc__ = 'Value of the cost function at the solution.'


def whitecap_coverage(ws):
    r'''
    Fraction of the sea surface covered by whitecaps [Monahan1980]_.

    .. math::

        f_f = 2.95\times10^{-6}\, U^{3.52}

    Parameters
    ----------
    ws : float or array_like
        Wind speed :math:`U` (m/s).

    Returns
    -------
    float or numpy.ndarray
        Whitecap coverage :math:`f_f` (between 0 and 1).
    '''
    return np.clip(2.95e-6 * np.abs(ws) ** 3.52, 0, 1)


def _as_geometries(geometries):
    geometries = np.atleast_2d(np.asarray(geometries, dtype=float))
    if geometries.shape[-1] != 3:
        raise ValueError('geometries must be of shape (N, 3): (sza, vza, azi) in degrees')
    return geometries


def glint_whitecap(geometries, ws, wazi=0, stats='bh2006', shadow=True,
                   T=1., td=1., foam_reflectance=0.13, m=1.334):
    r'''
    Sunglint plus whitecap Stokes vector for a set of viewing geometries.

    .. math::

        \mathbf{S}_{g+wc} = (1 - f_f)\, T\, \mathbf{S}_g + f_f\, t_d\, \mathbf{S}_{wc}

    Parameters
    ----------
    geometries : array_like
        Viewing geometries, shape ``(N, 3)``: solar zenith, viewing zenith and
        relative azimuth angles in degrees (see :class:`coxmunk.coxmunk.sunglint`).
    ws : float
        Wind speed (m/s).
    wazi : float, optional
        Wind azimuth (deg), see :meth:`coxmunk.coxmunk.sunglint.sunglint`.
    stats : {'cm_iso', 'cm_dir', 'bh2006'}, optional
        Wave slope statistics.
    shadow : bool, optional
        Apply the shadowing correction.
    T, td : float or array_like, optional
        Direct and diffuse atmospheric transmittances (scalar or one value per
        geometry). Default 1 (no atmosphere).
    foam_reflectance : float, optional
        Reflectance of the foam :math:`R_{wc}`, assumed unpolarized
        (default 0.13, typical of 865 nm).
    m : float, optional
        Refractive index of water.

    Returns
    -------
    numpy.ndarray
        Stokes parameters :math:`[I, Q, U]` in reflectance units, shape ``(N, 3)``.
    '''
    geometries = _as_geometries(geometries)
    # floor avoids null slope variances for some statistics at ws = 0
    ws = max(abs(ws), 1e-2)
    Sg = np.array([sunglint(sza, vza, azi, m=m).sunglint(ws, wazi, stats=stats, shadow=shadow)[:3]
                   for sza, vza, azi in geometries])
    ff = whitecap_coverage(ws)
    Swc = np.array([foam_reflectance, 0., 0.])
    T = np.reshape(T, (-1, 1)) if np.ndim(T) else T
    td = np.reshape(td, (-1, 1)) if np.ndim(td) else td
    return (1 - ff) * T * Sg + ff * td * Swc


def retrieve_wind_speed(observed, geometries, components=COMPONENTS, first_guesses=(1., 6., 12.),
                        epsilon=0.05, ws_max=30., **kwargs):
    r'''
    Retrieve the wind speed from multidirectional (polarized) sunglint
    measurements [HarmelChami2012]_.

    The cost function

    .. math::

        \chi^2(x) = \left\| \mathbf{S}^*_{g+wc} - \mathbf{S}_{g+wc}(x) \right\|^2

    is minimized with the Levenberg-Marquardt method, using the Stokes
    parameters listed in ``components`` for all the viewing directions.

    **Initialization.** The inversion is first run from each of the
    ``first_guesses`` (1, 6 and 12 m/s by default). If none of the solutions
    departs from its first guess, the measurements are not informative on wind
    speed and the pixel is flagged. Otherwise, the final inversion starts from
    the mean of the solutions.

    **Uncertainty.** The uncertainty :math:`\sigma` is the largest departure from
    the solution :math:`x` for which the cost function stays within a fraction
    :math:`\varepsilon` of its minimum:

    .. math::

        |x - x'| \le \sigma \quad \text{for} \quad \chi^2(x') \le (1 + \varepsilon)\,\chi^2(x)

    Parameters
    ----------
    observed : array_like
        Measured sunglint plus whitecap Stokes parameters, shape ``(N, 3)``
        (:math:`I, Q, U` in reflectance units), one row per geometry.
    geometries : array_like
        Viewing geometries, shape ``(N, 3)``, see :func:`glint_whitecap`.
    components : sequence of {'I', 'Q', 'U'}, optional
        Stokes parameters used in the cost function. Use ``('I',)`` for a
        radiance-only retrieval.
    first_guesses : sequence of float, optional
        Starting wind speeds (m/s).
    epsilon : float, optional
        Tolerance :math:`\varepsilon` on the cost function for the uncertainty
        (default 5 %).
    ws_max : float, optional
        Upper limit (m/s) of the search for the uncertainty bounds.
    **kwargs
        Passed to :func:`glint_whitecap` (``wazi``, ``stats``, ``shadow``,
        ``T``, ``td``, ``foam_reflectance``, ``m``).

    Returns
    -------
    WindRetrieval
        Named tuple ``(ws, sigma, flag, cost)``.

    Notes
    -----
    The wind direction is not retrieved: it is an input of the forward model
    (``wazi``), as the sunglint pattern is ambiguous with respect to the wind
    azimuth [HarmelChami2012]_.

    For noise-free synthetic measurements, the cost function is null at the
    solution and the uncertainty is zero.
    '''
    geometries = _as_geometries(geometries)
    observed = np.atleast_2d(np.asarray(observed, dtype=float))
    idx = [COMPONENTS.index(c) for c in components]
    obs = observed[:, idx].ravel()

    def residuals(x):
        return glint_whitecap(geometries, x[0], **kwargs)[:, idx].ravel() - obs

    def cost(x):
        return np.sum(residuals([x]) ** 2)

    def solve(x0):
        res = optimize.least_squares(residuals, [x0], method='lm')
        return abs(res.x[0])

    solutions = np.array([solve(x0) for x0 in first_guesses])
    if np.all(np.abs(solutions - np.asarray(first_guesses)) < 1e-3):
        return WindRetrieval(np.nan, np.nan, True, np.nan)

    ws = solve(solutions.mean())
    chi2 = cost(ws)

    # uncertainty: distance to the points where chi2 = (1 + epsilon) * chi2_min
    level = (1 + epsilon) * chi2

    def excess(x):
        return cost(x) - level

    # excess(ws) = -epsilon * chi2 <= 0: search the crossing on each side of the solution
    sigma = 0.
    for bound in (1e-2, ws + ws_max):
        if bound == ws:
            continue
        if excess(bound) <= 0:
            # the cost function never exceeds the tolerance on this side
            sigma = np.inf
            break
        x_eps = optimize.brentq(excess, min(ws, bound), max(ws, bound), xtol=1e-4)
        sigma = max(sigma, abs(x_eps - ws))
    return WindRetrieval(ws, sigma, False, chi2)


def glint_directions(geometries, ws, dws=1., noise=4e-4, **kwargs):
    r'''
    Viewing directions that can be contaminated by sunglint [HarmelChami2013]_.

    A direction is considered free of sunglint when the Cox-Munk sunglint
    normalized radiance computed for the wind speed :math:`U - \Delta U` and
    :math:`U + \Delta U` is below the radiometric noise of the sensor. The
    margin :math:`\Delta U` (1 m/s) accounts for the typical uncertainty of
    ancillary wind speed data.

    The normalized radiance is :math:`\pi L / F_0 = \mu_s\, I`, with :math:`I`
    the glint reflectance computed by :class:`coxmunk.coxmunk.sunglint`.

    Parameters
    ----------
    geometries : array_like
        Viewing geometries, shape ``(N, 3)``, see :func:`glint_whitecap`.
    ws : float
        Ancillary wind speed (m/s).
    dws : float, optional
        Wind speed uncertainty :math:`\Delta U` (m/s).
    noise : float, optional
        Noise equivalent normalized radiance of the sensor (default 4e-4,
        PARASOL).
    **kwargs
        Passed to :meth:`coxmunk.coxmunk.sunglint.sunglint` (``wazi``,
        ``stats``, ``shadow``).

    Returns
    -------
    numpy.ndarray of bool
        True for the directions potentially contaminated by sunglint, shape ``(N,)``.
    '''
    geometries = _as_geometries(geometries)
    mu_s = np.cos(np.radians(geometries[:, 0]))
    glint = np.zeros(len(geometries), dtype=bool)
    for u in (max(ws - dws, 1e-2), ws + dws):
        I = np.array([sunglint(sza, vza, azi).sunglint(u, **kwargs)[0] for sza, vza, azi in geometries])
        glint |= mu_s * I > noise
    return glint
