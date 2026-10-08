Wind speed from sunglint
========================

Knowledge of the surface wind is critical for the budgets of energy transport,
oceanic primary productivity and ocean acidification. Over the ocean, wind speed
is usually measured from space by active or passive microwave sensors, with a
resolution of 25 to 50 km and a typical accuracy of 1 m s\ :sup:`-1`. Since the
sunglint pattern directly depends on wind speed (:doc:`algorithms`), it can be
inverted to estimate the wind from visible and near-infrared measurements.

Using the multidirectional and polarized measurements of PARASOL
(:doc:`polarization`), [HarmelChami2012]_ retrieved the wind speed and its
uncertainty over the full swath, at the 6 km resolution of the sensor. The
forward model and the inverse method are implemented in :py:mod:`coxmunk.wind`
and illustrated in the :doc:`tutorials/wind_speed_retrieval` tutorial.

Previous approaches
-------------------

[BreonHenriot2006]_ showed with POLDER that the unpolarized radiance measured in
the directions most sensitive to the sunglint characterizes the sea surface wind
for wind speeds below 14 m s\ :sup:`-1`. Their approach was limited to a
restricted number of viewing directions, to the radiance only (no polarization)
and to a resolution of 50 km. [HarmelChami2012]_ use all the viewing directions
of a pixel together with the polarized information.

Forward model
-------------

At the top of the atmosphere, the signal is decomposed as

.. math::
   :label: toa_2012

   \mathbf{S}_{TOA} = T_g \mathbf{S}_{atm} + T \mathbf{S}_g + t_d \mathbf{S}_{wc} + t_d \mathbf{S}_w^+

where :math:`T_g` is the gaseous transmittance, :math:`T` and :math:`t_d` the
direct and diffuse atmospheric transmittances. The POLAC atmospheric correction
algorithm [HarmelChami2013]_ retrieves the atmospheric and water-leaving terms,
so that the sunglint plus whitecap signal
:math:`\mathbf{S}^*_{g+wc} = T\mathbf{S}_g + t_d\mathbf{S}_{wc}` is derived from
the PARASOL data without any a priori assumption on the sea surface.

The sunglint :math:`\mathbf{S}_g` is computed with the Cox-Munk model and the
slope statistics of [BreonHenriot2006]_. The fraction of the surface covered by
whitecaps is given by the power law of [Monahan1980]_
(:func:`~coxmunk.wind.whitecap_coverage`):

.. math::
   :label: whitecap

   f_f = 2.95\times10^{-6}\, U^{3.52}

The foam is assumed unpolarized, with a reflectance of 0.13 at 865 nm, so that
the modelled signal is (:func:`~coxmunk.wind.glint_whitecap`)

.. math::
   :label: glint_whitecap

   \mathbf{S}_{g+wc} = (1 - f_f)\, T\, \mathbf{S}_g + f_f\, t_d\, \mathbf{S}_{wc}

Inverse method
--------------

The wind speed :math:`x` is the minimum of the cost function

.. math::
   :label: cost

   \chi^2(x) = \left\| \mathbf{S}^*_{g+wc} - \mathbf{S}_{g+wc}(x) \right\|^2

computed over the Stokes parameters :math:`I`, :math:`Q`, :math:`U` of all the
viewing directions. The problem is non-linear and is solved with the
Levenberg-Marquardt damped least-squares method.

**Uncertainty.** The uncertainty :math:`\sigma` expresses the sensitivity of the
cost function around the solution :math:`x`:

.. math::
   :label: sigma

   |x - x'| \le \sigma \quad \text{for} \quad \chi^2(x') \le (1 + \varepsilon)\,\chi^2(x)

with :math:`\varepsilon = 5\,\%`. A flat cost function gives a large
uncertainty.

**Initialization.** The inversion is first run from three wind speeds: 1, 6 and
12 m s\ :sup:`-1`. If the solutions depart from their first guesses, their mean
is the first guess of the final inversion. If they do not, the measurements of
the pixel are not informative on wind speed and the pixel is flagged.

**Wind direction.** Fairly similar sunglint patterns are obtained for different
wind azimuths, so that the wind direction is not retrieved. It is an input of
the forward model.

In ``coxmunk``:

.. code-block:: python

   from coxmunk.wind import retrieve_wind_speed

   # observed: (N, 3) array of the sunglint + whitecap I, Q, U for N viewing directions
   # geometries: (N, 3) array of (sza, vza, azi) in degrees
   res = retrieve_wind_speed(observed, geometries, wazi=0, stats='bh2006')
   print(res.ws, res.sigma, res.flag)

Pass ``components=('I',)`` to restrict the retrieval to the radiance, as in
[BreonHenriot2006]_.

Results
-------

Applied to PARASOL images, the method retrieved the wind speed for almost 80 %
of the cloud-free pixels, with uncertainties generally below 1 m s\ :sup:`-1`.

* Over the north-western Mediterranean Sea (5 May 2006), the retrieved wind
  field shows high wind speeds (about 12 m s\ :sup:`-1`) in the west and south,
  with higher uncertainties (about 0.7 m s\ :sup:`-1`), likely due to whitecaps.
  A zone of very weak wind (below 2 m s\ :sup:`-1`) appears in the Gulf of Lion,
  with uncertainties below 0.3 m s\ :sup:`-1`. The PARASOL wind field resolves
  finer structures than the ECMWF model, such as a local wind minimum and
  strong gradients.
* Against two Météo-France buoys over the year 2006 (146 match-ups), the
  correlation coefficient was greater than 0.96, the slope 0.96 and the RMSE
  1.1 m s\ :sup:`-1`.
* Against the AMSR-E microwave sensor (NASA) at the global scale, over three
  days of acquisitions and after degradation to the 25 km AMSR-E grid
  (78 044 pixels), the correlation coefficient was 0.84, the slope within 1 %
  and the RMSE 1.57 m s\ :sup:`-1`, in agreement with previous studies.

Passive sensors measuring the polarization and the multidirectionality of the
radiation at solar wavelengths are thus an alternative to estimate the surface
wind at a resolution four to eight times finer than microwave sensors. The
method was designed to prepare the multidirectional polarimetric missions PACE
(NASA) and 3MI (ESA).

References
----------

See [HarmelChami2012]_, [HarmelChami2013]_ and [Monahan1980]_ in
:py:mod:`coxmunk.wind`, and [BreonHenriot2006]_ in :py:mod:`coxmunk.coxmunk`.
