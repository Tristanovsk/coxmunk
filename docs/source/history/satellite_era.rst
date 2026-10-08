From kumatage to sunglint: the satellite era
============================================

Seen from a satellite, the kumatage becomes the **sunglint**. For ocean color
remote sensing, which monitors phytoplankton, suspended and dissolved matter
from the light leaving the water, it is again a noise: the glint is often one
order of magnitude brighter than the water signal and blinds a large part of
the swath of sensors such as MERIS or Sentinel-3/OLCI. Dedicated algorithms were
developed to remove it, e.g., from the shortwave-infrared bands of
Sentinel-2/MSI [Harmel2018]_.

Noise and signal
----------------

The same glint is also a **signal**. The multidirectional and polarized
measurements of the POLDER and PARASOL satellite sensors confirmed the Cox and
Munk statistics at the global scale [BreonHenriot2006]_ [Munk2009]_, and were
used to retrieve the sunglint Stokes vector (:doc:`../polarization`) and the sea
surface wind speed (:doc:`../wind_speed`) at a resolution of a few kilometres.

Sentinel-2: Corsica and the Ligurian currents
---------------------------------------------

At the 10 to 20 m resolution of Sentinel-2, the glint reveals the fine
structure of the sea surface roughness. In the images below (shortwave-infrared
bands), the whiter the sea, the more sunglint.

.. figure:: ../../../illustration/kumatage/s2_2019-06-21_corsica_sunglint.jpg
   :width: 100%

   Sunglint around Cape Corsica (Sentinel-2B, 21 June 2019). The dark area in the
   lee of the island is sheltered from the wind: the calm water reflects the Sun
   only near the specular direction, outside of the field of view. The stripes
   are the boundaries of the detectors of the instrument, which observe the
   surface under slightly different viewing angles.

.. figure:: ../../../illustration/kumatage/s2_2019-06-21_corsica_sunglint_zoom.jpg
   :width: 100%

   Zoom on the same image: wakes of ships and slicks modify the roughness and
   the glint.

.. list-table::
   :widths: 50 50

   * - .. figure:: ../../../illustration/kumatage/s2_2018-10-04_ligurian_sea_sunglint.jpg
          :width: 100%

          Ligurian Sea, 4 October 2018 (Sentinel-2B).
     - .. figure:: ../../../illustration/kumatage/s2_2018-10-14_ligurian_sea_sunglint.jpg
          :width: 100%

          Ligurian Sea, 14 October 2018 (Sentinel-2B): the Corsica and Ligurian
          currents and their central front modulate the surface roughness.

The Sentinel-2 images contain modified Copernicus Sentinel data (2018, 2019),
processed with Sentinel Hub.

Applications
------------

The sunglint signal is now exploited to estimate oceanic winds
[HarmelChami2012]_, to monitor oil spills (Chust and Sagarminaga, 2007; Hu et
al., 2009), to retrieve wave spectra and ocean currents [Kudryavtsev2017a]_
[Kudryavtsev2017b]_, to detect internal waves (Jackson, 2007), absorbing aerosols
(Kaufman et al., 2002) and biofilms at the surface (Neukermans et al., 2018), or
to map the bathymetry through its modulation of the waves [Shao2011]_.
Beyond the Earth, simulations show that the intensity and polarization of the
glint of a star on an exoplanet could reveal the presence of an ocean
[TreesStam2019]_: the kumatage could help assess the habitability of other
worlds.

Kumatage or sunglint?
---------------------

* The reflection of the Sun on the sea has been a major inspiration in the
  history of art, but its naming was established through scientific and
  technical usages.
* The *kumatage* denomination was abandoned for lack of scientific application.
* The satellite era woke up the lost kumatage under the name of *sunglint*,
  turning a noise into a signal on the ocean surface.
* It is now up to you to choose between kumatage and sunglint.

References
----------

.. [Harmel2018] Harmel, T., Chami, M., Tormos, T., Reynaud, N. and Danis, P.-A.
   (2018). Sunglint correction of the Multi-Spectral Instrument (MSI)-SENTINEL-2
   imagery over inland and sea waters from SWIR bands. *Remote Sens. Environ.*,
   204, 308-321.
.. [Kudryavtsev2017a] Kudryavtsev, V., et al. (2017). Sun glitter imagery of ocean
   surface waves. Part 1: Directional spectrum retrieval and validation.
   *J. Geophys. Res. Oceans*.
.. [Kudryavtsev2017b] Kudryavtsev, V., et al. (2017). Sun glitter imagery of surface
   waves. Part 2: Waves transformation on ocean currents. *J. Geophys. Res. Oceans*.
.. [Shao2011] Shao, H., et al. (2011). Sun glitter imaging of submarine sand waves
   on the Taiwan Banks: Determination of the relaxation rate of short waves.
   *J. Geophys. Res.*
.. [TreesStam2019] Trees, V. J. H. and Stam, D. M. (2019). Blue, white, and red
   ocean planets: Simulations of orbital variations in flux and polarization
   colors. *Astron. Astrophys.*, 626, A129.

* Spooner (1822). Letters to the *Correspondance astronomique, géographique,
  hydrographique et statistique du Baron de Zach*, vol. VI, p. 331, and
  *Lettre IV*, p. 65.
* Bowditch, N. (1841). *The New American Practical Navigator*, 12th new
  stereotype edition, E. & G. W. Blunt, New York.
