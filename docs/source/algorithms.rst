Algorithm Description
=====================

This page gives the mathematical formulation of the sunglint computation, as
implemented in :py:meth:`coxmunk.coxmunk.sunglint.sunglint`. Each section links
to the method that implements it. The references are listed in
:py:mod:`coxmunk.coxmunk`.

.. mermaid::

   flowchart LR
       G["Geometry<br/>θs, θv, φ"] --> S["Scattering angle Θ<br/>facet tilt θn"]
       S --> Z["Facet slopes<br/>z_up, z_cr"]
       W["Wind<br/>U, φw"] --> Z
       Z --> P["Slope density P"]
       S --> F["Fresnel matrix R_F(ω)"]
       F --> R["Rotation to the<br/>scattering plane"]
       W --> SH["Shadowing SH"]
       G --> SH
       P --> O["Stokes vector<br/>I, Q, U, V"]
       R --> O
       SH --> O

Notation
--------

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Symbol
     - Definition
     - Unit
   * - :math:`\theta_s`, :math:`\theta_v`
     - Solar and viewing zenith angles
     - deg
   * - :math:`\mu_s`, :math:`\mu_v`
     - :math:`\cos\theta_s`, :math:`\cos\theta_v`
     - --
   * - :math:`\phi`
     - Relative azimuth, 180° when Sun and sensor are in opposition (Sun azimuth set to 0)
     - deg
   * - :math:`U`
     - Wind speed
     - m s\ :sup:`-1`
   * - :math:`\phi_w`
     - Downwind direction, counted counterclockwise from the Sun direction
     - deg
   * - :math:`\Theta`
     - Scattering angle
     - rad
   * - :math:`\omega`
     - Angle of incidence on the reflecting facet
     - rad
   * - :math:`\theta_n`
     - Tilt of the reflecting facet
     - rad
   * - :math:`z_{up}`, :math:`z_{cr}`
     - Upwind and crosswind facet slopes
     - --
   * - :math:`\sigma_{up}^2`, :math:`\sigma_{cr}^2`
     - Upwind and crosswind slope variances
     - --
   * - :math:`m`
     - Refractive index of water (default 1.334)
     - --

Reflecting facet
----------------

The glint seen in a given direction comes from the facets oriented so that the
Sun is specularly reflected toward the sensor. The scattering angle
(:py:meth:`~coxmunk.coxmunk.sunglint.scat_angle`) is

.. math::
   :label: scat_angle

   \cos\Theta = -\mu_s\mu_v - \sin\theta_s\sin\theta_v\cos\phi

from which the angle of incidence on the facet and the facet tilt are

.. math::
   :label: facet

   \omega = \frac{\pi - \Theta}{2}, \qquad
   \cos\theta_n = \frac{\mu_s + \mu_v}{2\cos\omega}

The facet slopes in the Sun frame are

.. math::
   :label: slopes

   z_x = -\frac{\sin\theta_v \cos\phi + \sin\theta_s}{\mu_s + \mu_v}, \qquad
   z_y = -\frac{\sin\theta_v \sin\phi}{\mu_s + \mu_v}

and are rotated into the wind frame:

.. math::
   :label: slopes_wind

   z_{up} = \cos\phi_w \, z_x + \sin\phi_w \, z_y, \qquad
   z_{cr} = -\sin\phi_w \, z_x + \cos\phi_w \, z_y

Wave slope statistics
---------------------

The ``stats`` argument selects the slope probability density :math:`P`.

Isotropic Cox-Munk (``cm_iso``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Original isotropic statistics of [CoxMunk1954]_:

.. math::
   :label: pdf_iso

   \sigma^2 = 0.003 + 5.12\times10^{-3}\,U, \qquad
   P = \frac{1}{\pi\sigma^2} \exp\left(-\frac{\tan^2\theta_n}{\sigma^2}\right)

The upwind and crosswind variances used by the shadowing correction are both
set to :math:`\sigma^2/2`.

Anisotropic Gram-Charlier expansion (``cm_dir``, ``bh2006``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the normalized slopes :math:`\xi = z_{up}/\sigma_{up}` and
:math:`\eta = z_{cr}/\sigma_{cr}`, the density is written with the convention of
[Munk2009]_:

.. math::
   :label: pdf_gc

   P(\xi, \eta) = \frac{e^{-(\xi^2 + \eta^2)/2}}{2\pi\,\sigma_{up}\,\sigma_{cr}}
   \Big[ 1
   + \tfrac{c_{12}}{2}\,\xi\,(1 - \eta^2)
   - \tfrac{c_{30}}{6}\,\xi\,(3 - \xi^2)
   + \tfrac{c_{40}}{24}\,(3 - 6\xi^2 + \xi^4)
   + \tfrac{c_{22}}{4}\,(1 - \xi^2)(1 - \eta^2)
   + \tfrac{c_{04}}{24}\,(3 - 6\eta^2 + \eta^4)
   \Big]

where :math:`c_{12} = -c_{21}` and :math:`c_{30} = -c_{03}`. The variances and
the skewness (:math:`c_{21}`, :math:`c_{03}`) and peakedness (:math:`c_{40}`,
:math:`c_{22}`, :math:`c_{04}`) coefficients are:

.. list-table::
   :header-rows: 1

   * - Parameter
     - ``cm_dir`` [CoxMunk1954]_
     - ``bh2006`` [BreonHenriot2006]_
   * - :math:`\sigma_{cr}^2`
     - :math:`0.003 + 1.92\times10^{-3}\,U`
     - :math:`0.003 + 1.85\times10^{-3}\,U`
   * - :math:`\sigma_{up}^2`
     - :math:`3.16\times10^{-3}\,U`
     - :math:`0.001 + 3.16\times10^{-3}\,U`
   * - :math:`c_{21}`
     - :math:`0.01 - 8.6\times10^{-3}\,U`
     - :math:`-9\times10^{-4}\,U^2`
   * - :math:`c_{03}`
     - :math:`0.04 - 0.033\,U`
     - :math:`-0.45 / (1 + e^{7 - U})`
   * - :math:`c_{40}`
     - 0.40
     - 0.30
   * - :math:`c_{22}`
     - 0.12
     - 0.12
   * - :math:`c_{04}`
     - 0.23
     - 0.40

Fresnel reflection
------------------

The amplitude reflection coefficients parallel and perpendicular to the plane
of incidence are (:py:meth:`~coxmunk.coxmunk.sunglint.fresnel`)

.. math::
   :label: fresnel_coef

   r_l = \frac{\sqrt{m^2 - \sin^2\omega} - m^2\cos\omega}
              {\sqrt{m^2 - \sin^2\omega} + m^2\cos\omega},
   \qquad
   r_r = \frac{\cos\omega - \sqrt{m^2 - \sin^2\omega}}
              {\cos\omega + \sqrt{m^2 - \sin^2\omega}}

giving the Mueller reflection matrix

.. math::
   :label: fresnel_matrix

   \mathbf{R}_F(\omega) = \frac{1}{2}
   \begin{pmatrix}
   r_l^2 + r_r^2 & r_l^2 - r_r^2 & 0 & 0 \\
   r_l^2 - r_r^2 & r_l^2 + r_r^2 & 0 & 0 \\
   0 & 0 & 2 r_l r_r & 0 \\
   0 & 0 & 0 & 2 r_l r_r
   \end{pmatrix}

Rotation to the scattering plane
--------------------------------

The Stokes vectors are defined in the meridian planes of the Sun and of the
sensor. The Fresnel matrix is rotated accordingly,
:math:`\mathbf{R} = \mathbf{L}(\pi - \sigma_2)\,\mathbf{R}_F\,\mathbf{L}(-\sigma_1)`, with

.. math::
   :label: rotation

   \mathbf{L}(\chi) =
   \begin{pmatrix}
   1 & 0 & 0 & 0 \\
   0 & \cos 2\chi & \sin 2\chi & 0 \\
   0 & -\sin 2\chi & \cos 2\chi & 0 \\
   0 & 0 & 0 & 1
   \end{pmatrix},
   \quad
   \cos\sigma_1 = \frac{-\mu_s - \mu_v\cos\Theta}{\sin\theta_v \sin\Theta},
   \quad
   \cos\sigma_2 = \frac{\mu_v + \mu_s\cos\Theta}{\sin\theta_s \sin\Theta}

The signs of :math:`\sigma_1` and :math:`\sigma_2` are reversed when
:math:`\sin\phi > 0`. No rotation is applied when :math:`\theta_v`,
:math:`\theta_s` or :math:`\Theta` is zero.

Shadowing and hiding
--------------------

When ``shadow=True``, the facets hidden from the Sun or from the sensor by the
waves are removed with the Smith function of [RossDion2005]_ and [RossDion2007]_
(:py:meth:`~coxmunk.coxmunk.sunglint.nu`,
:py:meth:`~coxmunk.coxmunk.sunglint.Lambda`):

.. math::
   :label: lambda

   \Lambda(\theta, \varphi) = \frac{e^{-\nu^2} - \nu\sqrt{\pi}\,\mathrm{erfc}(\nu)}
                                   {2\,\nu\sqrt{\pi}},
   \qquad
   \nu = \frac{1}{\sqrt{2}\,\sigma(\varphi)\tan\theta},
   \qquad
   \sigma^2(\varphi) = \sigma_{up}^2\cos^2\varphi + \sigma_{cr}^2\sin^2\varphi

where :math:`\varphi` is the azimuth of the direction relative to the upwind
axis, and :math:`\Lambda = 0` for :math:`\theta = 0`. The shadowing factor is

.. math::
   :label: shadow

   SH = \frac{1}{1 + \Lambda(\theta_s, -\phi_w) + \Lambda(\theta_v, \phi - \phi_w)}

and :math:`SH = 1` when ``shadow=False`` or :math:`\theta_v = 0`.

Sunglint Stokes vector
----------------------

Combining :eq:`facet`, :eq:`pdf_iso` or :eq:`pdf_gc`, :eq:`rotation` and
:eq:`shadow`, the sunglint in reflectance units is

.. math::
   :label: glint

   \begin{pmatrix} I \\ Q \\ U \\ V \end{pmatrix}
   = \frac{\pi \, P \, SH}{4 \, \mu_s \, \mu_v \cos^4\theta_n}
     \begin{pmatrix} R_{11} \\ R_{21} \\ R_{31} \\ R_{41} \end{pmatrix}

i.e., the first column of :math:`\mathbf{R}`, the response to an unpolarized
incident beam. The atmospheric transmittances are not applied (they are set
to 1); :py:meth:`~coxmunk.coxmunk.sunglint.atmo_trans` gives the direct
transmittance along the solar path if needed.
