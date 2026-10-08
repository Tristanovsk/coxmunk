Polarized sunglint
==================

The reflection of sunlight on the sea surface polarizes the light. This chapter
describes the polarization of the sunglint computed by ``coxmunk`` and how the
multidirectional and polarized measurements of the PARASOL satellite sensor
were used to retrieve the full Stokes vector of the sunglint without any a
priori knowledge of the sea state [HarmelChami2013]_. The
:doc:`tutorials/polarization` tutorial reproduces the main features discussed
here.

Polarization by Fresnel reflection
----------------------------------

The Stokes vector :math:`\mathbf{S} = [I, Q, U, V]^T` describes the intensity
and the polarization state of light in terms of measurable quantities. In
``coxmunk``, as in [HarmelChami2013]_, it is defined with respect to the
meridian plane of the sensor (the plane containing the viewing direction and the
local zenith). The circular component :math:`V` is negligible for the sunglint
and is not discussed further.

For an unpolarized incident sunlight, the light reflected by a facet is the
first column of the Fresnel matrix (:eq:`fresnel_matrix`). Before rotation into
the meridian plane, its degree of linear polarization is

.. math::
   :label: dolp_fresnel

   \mathrm{DoLP} = \frac{\sqrt{Q^2 + U^2}}{I} = \frac{\left| r_l^2 - r_r^2 \right|}{r_l^2 + r_r^2}

It depends on the angle of incidence :math:`\omega` on the facet only, i.e., on
the scattering angle :math:`\Theta = \pi - 2\omega`, and not on the sea state.
The reflected light is totally polarized at the Brewster angle,
:math:`\omega_B = \arctan m \approx 53°` (:math:`r_l = 0`), which corresponds to
a scattering angle of about 74°. The wind speed changes the amount of glint
(through the slope density :math:`P`), not its degree of polarization. The
angle of polarization

.. math::
   :label: aop

   \chi = \frac{1}{2} \arctan\frac{U}{Q}

is set by the rotation from the scattering plane to the meridian plane
(:eq:`rotation`). It is therefore a purely geometrical quantity.

Sunglint in the top-of-atmosphere signal
----------------------------------------

The normalized radiance measured by a satellite sensor is
:math:`\pi L / F_0`, with :math:`F_0` the extraterrestrial solar irradiance. It
equals :math:`\mu_s` times the reflectance computed by ``coxmunk``. For
incoherent light, the Stokes vectors of the different contributions add up, and
the top-of-atmosphere (TOA) signal reads

.. math::
   :label: toa

   \mathbf{S}_{TOA} = T_g(\theta_s, \theta_v) \left[ \mathbf{S}_{atm}
   + T_{down}(\theta_s) T_{up}(\theta_v) \mathbf{S}_g
   + t_{down}(\theta_s) t_{up}(\theta_v) \mathbf{S}_{wc}
   + \mathbf{S}_w^+ \right]

where :math:`T_g` is the gaseous transmittance, :math:`T_{down}` and
:math:`T_{up}` the direct transmittances (the sunglint is the direct sunlight
reflected by the surface), :math:`t_{down}` and :math:`t_{up}` the total
transmittances, and :math:`\mathbf{S}_{atm}`, :math:`\mathbf{S}_g`,
:math:`\mathbf{S}_{wc}` and :math:`\mathbf{S}_w^+` the atmospheric, sunglint,
whitecap and water-leaving contributions. Whitecaps are made of many small
bubbles that depolarize light by multiple scattering: their contribution is
unpolarized (:math:`Q_{wc} = U_{wc} = 0`).

For a Sun at 35° and wind speeds of 2, 6 and 10 m s\ :sup:`-1`, the sunglint
affects only a small range of viewing directions at low wind speed, but almost
one third of the directions at 10 m s\ :sup:`-1`. Over this range, the sunglint
Stokes parameters are generally one order of magnitude larger than the
contributions of the atmosphere and of the ocean. Some geometries are never
affected by the glint, whatever the wind speed.

The PARASOL mission
-------------------

The PARASOL sensor (CNES), the third generation of the POLDER instrument, flew
in the A-Train constellation and measured the Stokes parameters :math:`I`,
:math:`Q` and :math:`U` at 490, 670 and 865 nm (:math:`I` only in other bands
from 443 to 1020 nm). Thanks to its wide field of view (about 114°), a given
ground target is observed under up to 16 viewing directions within 4 minutes,
with a resolution of about 6 km × 7 km at nadir. Its noise equivalent normalized
radiance is about :math:`4\times10^{-4}`.

Within a multidirectional sequence, only some directions are affected by the
sunglint. They can be used to estimate the glint, while the others constrain
the atmosphere.

The POLAC-glint algorithm
-------------------------

The POLAC atmospheric correction algorithm retrieves the aerosol optical
thickness :math:`\tau_a` and the aerosol model from the multidirectional and
polarized PARASOL data, then the water-leaving radiance. In its sunglint
extension, POLAC-glint, :math:`\tau_a` is retrieved at 865 nm independently for
each viewing direction :math:`i` (:math:`\tau_{a,dir}(i)`), without any glint
in the model. The aerosol load does not depend on the viewing direction, so a
direction contaminated by sunglint stands out with an overestimated
:math:`\tau_{a,dir}` (up to five times the actual value). The directional
dispersion

.. math::
   :label: dtau

   \Delta\tau_a = \sqrt{\frac{\sum_{i=1}^{N_{dir}} \left(\tau_{a,dir}(i)
   - \mathrm{median}(\tau_{a,dir})\right)^2}{N_{dir} - 1}}

is used to remove the contaminated directions iteratively:

1. POLAC is applied to the :math:`N_{dir}` available directions;
2. :math:`\Delta\tau_a` is computed;
3. if :math:`\Delta\tau_a > 0.03 + 0.05\,\tau_a`, the direction with the highest
   :math:`\tau_{a,dir}` is removed, and POLAC is applied again to the remaining
   directions;
4. the iterations stop when :math:`\Delta\tau_a` no longer decreases, and the
   pixel is discarded when fewer than three directions remain.

The atmosphere and the water-leaving signal retrieved from the glint-free
directions are then used to compute the sunglint Stokes vector of all the
directions by inverting :eq:`toa` (whitecaps neglected up to
10 m s\ :sup:`-1`):

.. math::
   :label: sg_polac

   \mathbf{S}_g = \frac{1}{T_{down}(\theta_s) T_{up}(\theta_v)}
   \left[ \frac{\mathbf{S}_{TOA}}{T_g(\theta_s, \theta_v)} - \mathbf{S}_{atm}
   - t_{down}(\theta_s) t_{up}(\theta_v) \mathbf{S}_w^+ \right]

The major originality of the approach is that it requires no assumption on the
sea state: the sunglint is measured, not modelled.

Cloud edges
~~~~~~~~~~~

Thin clouds and cloud edges also increase the signal and can be mistaken for
sunglint. A direction flagged as contaminated is checked with the Cox-Munk
model fed with the ancillary (ECMWF) wind speed plus or minus 1 m s\ :sup:`-1`:
if both simulated glints are below the noise of the sensor, the direction is
geometrically out of the glint and is flagged as *cloud influenced*. A pixel
with more than three cloud-influenced directions is discarded. Over one day of
global acquisitions, this procedure identified 9.1 % of the ocean pixels as
cloud edges, in addition to the 69.1 % of cloudy pixels detected by the
operational processing. It is less efficient where many directions are affected
by the glint.

The geometrical test is implemented in :func:`coxmunk.wind.glint_directions`:

.. code-block:: python

   from coxmunk.wind import glint_directions

   # geometries: (sza, vza, azi) in degrees for each viewing direction
   in_glint = glint_directions(geometries, ws=ecmwf_wind_speed, dws=1., noise=4e-4)

Validation
----------

The sunglint Stokes parameters retrieved by POLAC-glint over the world ocean
(5 May 2006) were compared with the Cox-Munk model fed with ECMWF winds:

For the three Stokes parameters, the coefficient of determination
:math:`R^2` ranged from 0.84 to 0.92 and the slopes of the regression lines
from 0.88 to 0.96. The median absolute percentage differences were 22 % for
:math:`I_g`, 32 % for :math:`Q_g` and 53 % for :math:`U_g`.

The root-mean-square difference, about 0.02 for each parameter, is consistent
with the uncertainty of the Cox-Munk model propagated from the typical error of
the ECMWF wind speed (about 2 m s\ :sup:`-1`). The larger relative differences of
:math:`Q_g` and :math:`U_g` come from their small absolute values.

Spectral flatness
~~~~~~~~~~~~~~~~~

The refractive index of seawater varies only from 1.337 to 1.331 between 670 and
865 nm, which changes the Fresnel coefficients by less than 3 %. The sunglint
Stokes vector should therefore be spectrally flat. The ratios of the retrieved
sunglint parameters between 865 and 670 nm, corrected for this small Fresnel
effect, were 1.02, 1.00 and 0.98 for :math:`I_g`, :math:`Q_g` and :math:`U_g`,
with a relative dispersion below 17.3 %. For comparison, the same ratios are
centred on 0.5 for pixels out of the glint, where the signal is dominated by the
atmosphere. This was, to our knowledge, the first use of the polarized
components of the sunglint and of their spectral properties to detect the glint
and validate its retrieval.

The sensitivity to the refractive index can be checked with the ``m`` argument
of :class:`coxmunk.coxmunk.sunglint`.

Perspectives
------------

The retrieved sunglint Stokes vector can be used to determine the wind speed of
each pixel (:doc:`wind_speed`), and from it the whitecap contribution, to
reassess the wave slope statistics, to provide reliable glint data for sensor
calibration and to retrieve geophysical parameters from the glint (water vapor,
absorbing aerosols).

References
----------

See [HarmelChami2013]_ (in :py:mod:`coxmunk.wind`).
