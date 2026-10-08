# coding=utf-8
r'''
Polarized sunglint reflectance of a wind-roughened sea surface.

The sea surface is modelled as a collection of plane facets whose slopes follow
the statistics of [CoxMunk1954]_. Light from the Sun is specularly reflected
toward the sensor by the facets tilted at the exact orientation linking the two
directions. The Stokes vector of the glint is

.. math::

    \begin{pmatrix} I \\ Q \\ U \\ V \end{pmatrix}_{glint}
    = \frac{\pi \, P(z_{up}, z_{cr}) \, SH}{4 \, \mu_s \, \mu_v \cos^4\theta_n}
      \; \mathbf{L}(\pi - \sigma_2) \, \mathbf{R}_F(\omega) \, \mathbf{L}(-\sigma_1)
      \begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}

where

* :math:`\mu_s = \cos\theta_s` and :math:`\mu_v = \cos\theta_v` are the cosines of
  the solar and viewing zenith angles,
* :math:`P` is the probability density of the facet slopes (see
  :meth:`sunglint.sunglint`),
* :math:`\theta_n` is the tilt of the reflecting facet,
* :math:`\omega` is the angle of incidence on the facet,
* :math:`\mathbf{R}_F` is the Fresnel reflection matrix (see :meth:`sunglint.fresnel`),
* :math:`\mathbf{L}` are rotation matrices from the meridian planes to the
  scattering plane,
* :math:`SH` is the shadowing/hiding factor (see :meth:`sunglint.Lambda`).

The glint is expressed in reflectance units for an unpolarized incident
irradiance, without atmospheric attenuation.

Geometry conventions
--------------------
* The Sun azimuth is 0 in the reference frame.
* The relative azimuth :math:`\phi` is 180° when the Sun and the sensor are in
  opposition (i.e. the sensor looks at the specular point).
* The wind azimuth :math:`\phi_w` is the downwind direction counted from the Sun
  direction counterclockwise.

References
----------
.. [CoxMunk1954] Cox, C. and Munk, W. (1954). Measurement of the roughness of the sea
   surface from photographs of the Sun's glitter. *J. Opt. Soc. Am.*, 44(11), 838-850.
.. [BreonHenriot2006] Bréon, F.-M. and Henriot, N. (2006). Spaceborne observations of
   ocean glint reflectance and modeling of wave slope distributions.
   *J. Geophys. Res.*, 111, C06005.
.. [Munk2009] Munk, W. (2009). An inconvenient sea truth: spread, steepness, and
   skewness of surface slopes. *Annu. Rev. Mar. Sci.*, 1, 377-415.
.. [RossDion2005] Ross, V., Dion, D. and Potvin, G. (2005). Detailed analytical approach
   to the Gaussian surface bidirectional reflectance distribution function
   specular component applied to the sea surface. *J. Opt. Soc. Am. A*, 22(11), 2442-2453.
.. [RossDion2007] Ross, V. and Dion, D. (2007). Sea surface slope statistics derived from
   Sun glint radiance measurements and their apparent dependence on sensor
   elevation. *J. Geophys. Res.*, 112, C09015.
'''
import numpy as np
from scipy import special


class sunglint:
    r'''
    Sunglint computation for a given Sun/sensor geometry.

    Parameters
    ----------
    sza : float
        Solar zenith angle :math:`\theta_s` (deg).
    vza : float
        Viewing zenith angle :math:`\theta_v` (deg).
    azi : float
        Relative azimuth :math:`\phi` (deg); 180° when Sun and sensor are in
        opposition, the Sun azimuth being 0 in the reference frame.
    m : float, optional
        Refractive index of water (default 1.334).
    tau_atm : float, optional
        Atmospheric optical thickness. Not used yet: transmittances are set to 1.

    Examples
    --------
    >>> from coxmunk import sunglint
    >>> I, Q, U, V = sunglint(sza=30, vza=30, azi=180).sunglint(ws=5, stats='bh2006')
    '''

    def __init__(self, sza, vza, azi, m=1.334, tau_atm=0):
        degrad = np.pi / 180.
        self.sza = sza * degrad
        self.mu0 = np.cos(self.sza)
        self.sin0 = np.sin(self.sza)
        self.vza = vza * degrad
        self.azi = azi * degrad
        self.m = m
        self.tau_atm = tau_atm

    def sunglint(self, ws, wazi=0, stats='cm_dir', shadow=True, slope=False):
        r'''
        Compute the Stokes vector of the sunglint for a given wind.

        **Facet orientation.** The scattering angle :math:`\Theta` (see
        :meth:`scat_angle`) gives the angle of incidence on the reflecting facet
        and its tilt :math:`\theta_n`:

        .. math::

            \omega = \frac{\pi - \Theta}{2}, \qquad
            \cos\theta_n = \frac{\mu_s + \mu_v}{2\cos\omega}

        The facet slopes in the Sun frame are

        .. math::

            z_x = -\frac{\sin\theta_v \cos\phi + \sin\theta_s}{\mu_s + \mu_v}, \qquad
            z_y = -\frac{\sin\theta_v \sin\phi}{\mu_s + \mu_v}

        and are rotated into the wind frame (upwind/crosswind) with the wind
        azimuth :math:`\phi_w`:

        .. math::

            z_{up} = \cos\phi_w \, z_x + \sin\phi_w \, z_y, \qquad
            z_{cr} = -\sin\phi_w \, z_x + \cos\phi_w \, z_y

        **Isotropic statistics** (``stats='cm_iso'``, [CoxMunk1954]_), with the
        wind speed :math:`U` in m/s:

        .. math::

            \sigma^2 = 0.003 + 5.12\times10^{-3}\,U, \qquad
            P = \frac{1}{\pi\sigma^2} \exp\left(-\frac{\tan^2\theta_n}{\sigma^2}\right)

        **Anisotropic statistics** (``stats='cm_dir'`` or ``'bh2006'``). The
        Gram-Charlier expansion is written with the convention of [Munk2009]_,
        with :math:`\xi = z_{up}/\sigma_{up}` and :math:`\eta = z_{cr}/\sigma_{cr}`:

        .. math::

            P(\xi, \eta) = \frac{e^{-(\xi^2 + \eta^2)/2}}{2\pi\,\sigma_{up}\,\sigma_{cr}}
            \Big[ 1
            + \tfrac{c_{12}}{2}\,\xi\,(1 - \eta^2)
            - \tfrac{c_{30}}{6}\,\xi\,(3 - \xi^2)
            + \tfrac{c_{40}}{24}\,(3 - 6\xi^2 + \xi^4)
            + \tfrac{c_{22}}{4}\,(1 - \xi^2)(1 - \eta^2)
            + \tfrac{c_{04}}{24}\,(3 - 6\eta^2 + \eta^4)
            \Big]

        where :math:`c_{12} = -c_{21}` and :math:`c_{30} = -c_{03}`, with the
        skewness coefficients :math:`c_{21}` and :math:`c_{03}` given in the
        convention of [CoxMunk1954]_ and [BreonHenriot2006]_:

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

        **Polarization.** The Fresnel matrix :math:`\mathbf{R}_F(\omega)` (see
        :meth:`fresnel`) is rotated from the meridian planes into the scattering
        plane, :math:`\mathbf{R} = \mathbf{L}(\pi - \sigma_2)\,\mathbf{R}_F\,\mathbf{L}(-\sigma_1)`,
        with

        .. math::

            \mathbf{L}(\chi) =
            \begin{pmatrix}
            1 & 0 & 0 & 0 \\
            0 & \cos 2\chi & \sin 2\chi & 0 \\
            0 & -\sin 2\chi & \cos 2\chi & 0 \\
            0 & 0 & 0 & 1
            \end{pmatrix},
            \qquad
            \cos\sigma_1 = \frac{-\mu_s - \mu_v\cos\Theta}{\sin\theta_v \sin\Theta},
            \qquad
            \cos\sigma_2 = \frac{\mu_v + \mu_s\cos\Theta}{\sin\theta_s \sin\Theta}

        The signs of :math:`\sigma_1` and :math:`\sigma_2` are reversed when
        :math:`\sin\phi > 0`. No rotation is applied when :math:`\theta_v`,
        :math:`\theta_s` or :math:`\Theta` is zero.

        **Glint.** The Stokes vector is the first column of :math:`\mathbf{R}`
        multiplied by the weighted slope density:

        .. math::

            \begin{pmatrix} I \\ Q \\ U \\ V \end{pmatrix}
            = \mathbf{R}_{\cdot,1} \;
              \frac{\pi\,P\,SH}{4\,\mu_s\,\mu_v\cos^4\theta_n}

        where :math:`SH = 1` if ``shadow`` is False or :math:`\theta_v = 0`;
        otherwise it is computed from :meth:`Lambda`:

        .. math::

            SH = \frac{1}{1 + \Lambda(\theta_s, -\phi_w) + \Lambda(\theta_v, \phi - \phi_w)}

        where the second argument is the azimuth of the Sun or viewing direction
        relative to the upwind axis.

        Parameters
        ----------
        ws : float
            Wind speed :math:`U` (m/s).
        wazi : float, optional
            Downwind direction :math:`\phi_w` (deg) counted counterclockwise from
            the Sun direction (e.g., ``wazi=0`` means the wind blows toward the Sun).
        stats : {'cm_iso', 'cm_dir', 'bh2006'}, optional
            Wave slope statistics:

            * ``'cm_iso'``: isotropic Cox-Munk statistics [CoxMunk1954]_,
            * ``'cm_dir'``: directional Cox-Munk statistics [CoxMunk1954]_,
            * ``'bh2006'``: reassessment of [BreonHenriot2006]_.
        shadow : bool, optional
            If True, apply the shadowing/hiding correction of
            [RossDion2005]_ and [RossDion2007]_. With ``'cm_iso'``, the upwind
            and crosswind variances are both set to :math:`\sigma^2/2`.
        slope : bool or str, optional
            Selects the output (see Returns).

        Returns
        -------
        list
            * ``slope=False`` (default): ``[I, Q, U, V]``, the glint Stokes
              vector in reflectance units.
            * ``slope=True``: ``[z_up, z_cr, R_11, P]``, the upwind and crosswind
              slopes, the first Fresnel term and the slope probability density.
            * ``slope='check'``: ``[R_11, R_21, R_31, R_41, P, cos(theta_n), SH]``,
              the first column of the rotated Fresnel matrix, the slope
              probability density, the facet tilt cosine and the shadowing factor.
        '''
        sza = self.sza
        vza = self.vza
        azi = self.azi
        wazi = (wazi) * np.pi / 180.

        muv = np.cos(vza)
        sinv = np.sin(vza)
        mu0 = self.mu0
        sin0 = self.sin0
        muazi = np.cos(azi)
        sinazi = np.sin(azi)
        scat = self.scat_angle()
        omega = (np.pi - scat) / 2  # for x in self.scat_angle()]
        thetaN = np.arccos((mu0 + muv) / (2. * np.cos(omega)))

        Rf = self.fresnel(omega)

        if stats == 'cm_iso':  # Original isotropic Cox Munk statistics
            sigma2 = 3e-3 + 512e-5 * ws
            s_cr2 = sigma2 / 2
            s_up2 = sigma2 / 2
        elif stats == 'cm_dir':  # historical values from directional COX MUNK
            s_cr2 = 0.003 + 1.92e-3 * ws
            s_up2 = 3.16e-3 * ws
            s_cr = np.sqrt(s_cr2)
            s_up = np.sqrt(s_up2)
            c21 = 0.01 - 8.6e-3 * ws
            c03 = 0.04 - 33.e-3 * ws
            c40 = 0.40
            c22 = 0.12
            c04 = 0.23
        elif stats == 'bh2006':  # from Breon Henriot 2006 JGR
            s_cr2 = 3e-3 + 1.85e-3 * ws
            s_up2 = 1e-3 + 3.16e-3 * ws
            s_cr = np.sqrt(s_cr2)
            s_up = np.sqrt(s_up2)
            c21 = -9e-4 * ws ** 2
            c03 = -0.45 / (1 + np.exp(7. - ws))
            c40 = 0.30
            c22 = 0.12
            c04 = 0.4

        zx = -1 * (sinv * muazi + sin0) / (mu0 + muv)
        zy = -1 * (sinv * sinazi) / (mu0 + muv)

        if stats == 'cm_iso':
            #TODO check consitency between sigma2 and formulation
            #Pdist_ = 1. / (np.pi *2.* sigma2) * np.exp(-1./2 * np.tan(thetaN) ** 2 / sigma2)
            Pdist_ = 1. / (np.pi * sigma2) * np.exp(- np.tan(thetaN) ** 2 / sigma2)

            z_up,z_cr=zx,zy
        else:
            sigma2 = s_cr2 + s_up2

            # TODO check why it is clockwise rotation
            z_up = np.cos(wazi) * zx + np.sin(wazi) * zy
            z_cr = -np.sin(wazi) * zx + np.cos(wazi) * zy

            if False:
                # convention from Cox & Munk, 1956, and Breon & Henriot, 2006
                eta = z_up / s_up
                xi = z_cr / s_cr

                Pdist_ = np.exp(-5e-1 * (xi ** 2 + eta ** 2)) / (2. * np.pi * s_cr * s_up) * \
                         (1. -
                          c21 * (xi ** 2 - 1.) * eta / 2. -
                          c03 * (eta ** 3 - 3. * eta) / 6. +
                          c40 * (xi ** 4 - 6. * eta ** 2 + 3.) / 24. +
                          c04 * (eta ** 4 - 6. * eta ** 2 + 3.) / 24. +
                          c22 * (xi ** 2 - 1.) * (eta ** 2 - 1.) / 4.)

            else:
                # convention from Munk 2009
                eta = z_cr / s_cr
                xi = z_up / s_up
                c12 = -c21
                c30 = -c03
                Pdist_ = np.exp(-(xi ** 2 + eta ** 2) / 2) / (2. * np.pi* s_cr * s_up) * \
                         (1. +
                          c12 * xi * (1 - eta ** 2) / 2 -
                          c30 * xi * (3 - xi ** 2) / 6. +
                          c40 * (3 - 6 * xi ** 2 + xi ** 4) / 24. +
                          c22 * (1 - xi ** 2) * (1 - eta ** 2) / 4. +
                          c04 * (3 - 6. * eta ** 2 + eta ** 4) / 24.)

        # ---------------------------------------------------------------------*
        #                    Rotation in the reference plane
        # ---------------------------------------------------------------------*
        if sinv != 0 and sin0 != 0 and scat != 0:
            L1 = np.zeros((4, 4))
            L2 = np.zeros((4, 4))
            L1[0, 0] = 1.
            L1[3, 3] = 1.
            L2[0, 0] = 1.
            L2[3, 3] = 1.

            cos_sigma1 = (-mu0 - muv * np.cos(scat)) / (sinv * np.sin(scat))
            cos_sigma2 = (muv + mu0 * np.cos(scat)) / (sin0 * np.sin(scat))

            # debug epsilon error in calcumation
            if cos_sigma1 > 1:
                cos_sigma1 = 1
            elif cos_sigma1 < -1:
                cos_sigma1 = -1
            if cos_sigma2 > 1:
                cos_sigma2 = 1
            elif cos_sigma2 < -1:
                cos_sigma2 = -1

            sigma1 = np.arccos(cos_sigma1)
            sigma2 = np.arccos(cos_sigma2)

            if (np.sin(azi) > 0):
                sigma2 = -1 * sigma2
                sigma1 = -1 * sigma1

            L1[1, 1] = np.cos(-2 * sigma1)
            L1[2, 2] = L1[1, 1]
            L2[1, 1] = np.cos(2 * (np.pi - sigma2))
            L2[2, 2] = L2[1, 1]

            L1[1, 2] = np.sin(-2 * sigma1)
            L1[2, 1] = -1. * L1[1, 2]
            L2[1, 2] = np.sin(2. * (np.pi - sigma2))
            L2[2, 1] = -1. * L2[1, 2]

            Rf = np.matmul(Rf, L1)
            Rf = np.matmul(L2, Rf)

        if shadow and vza != 0:
            # azimuths of the Sun (0) and viewing (azi) directions relative to the upwind axis,
            # consistent with the rotation used for z_up and z_cr
            Ls = self.Lambda(s_up2, s_cr2, np.cos(wazi) ** 2, sza)
            Lr = self.Lambda(s_up2, s_cr2, np.cos(azi - wazi) ** 2, vza)
            SH = 1 / (1 + Ls + Lr)
        else:
            SH = 1.

        # ---------------------------------------------------------------------*
        #        Sun glint Stokes component (in reflectance unit) at TOA
        # ---------------------------------------------------------------------*

        Pdist = SH * (np.pi * Pdist_) / (4. * mu0 * muv * np.cos(thetaN) ** 4)

        Td = 1.
        Tu = 1.
        Pdist = Td * Tu * Pdist
        Iglint = Rf[0, 0] * Pdist
        Qglint = Rf[1, 0] * Pdist
        Uglint = Rf[2, 0] * Pdist
        Vglint = Rf[3, 0] * Pdist
        # if Iglint < 0:
        #    print(vza*180/np.pi,azi*180/np.pi,Pdist_, np.exp(-5e-1 * (xi ** 2 + eta ** 2)) / (2. * np.pi * s_cr * s_up))


        if slope == "check":
            return [Rf[0, 0],Rf[1, 0],Rf[2, 0],Rf[3, 0], Pdist_,np.cos(thetaN),SH]
        elif slope:
            return [z_up, z_cr, Rf[0, 0],Pdist_]
        else:
            return [Iglint, Qglint, Uglint,Vglint]  # ,Rf

    def atmo_trans(self, tau, sza, Iglint):
        r'''
        Attenuate the glint by the direct transmittance along the solar path.

        .. math::

            I' = I \, \exp\left(-\frac{\tau}{\cos\theta_s}\right)

        Parameters
        ----------
        tau : float
            Atmospheric optical thickness :math:`\tau`.
        sza : float
            Solar zenith angle :math:`\theta_s` (rad).
        Iglint : float
            Glint reflectance :math:`I`.

        Returns
        -------
        float
            Attenuated glint reflectance :math:`I'`.
        '''
        # ---------------------------------------------------------------------*
        #                           Direct Transmittance
        # ---------------------------------------------------------------------*
        sza
        Td = np.exp(-1 * tau / np.cos(sza))
        # Tup = np.exp(-1 * (self.tau_atm) / np.cos(self.vza))

        return Iglint * Td

    def fresnel(self, angle):
        r'''
        Fresnel reflection matrix of the air-water interface for illumination
        from above.

        With the refractive index :math:`m` and the angle of incidence
        :math:`\omega`, the amplitude reflection coefficients parallel
        (:math:`r_l`) and perpendicular (:math:`r_r`) to the plane of incidence are

        .. math::

            r_l = \frac{\sqrt{m^2 - \sin^2\omega} - m^2\cos\omega}
                       {\sqrt{m^2 - \sin^2\omega} + m^2\cos\omega},
            \qquad
            r_r = \frac{\cos\omega - \sqrt{m^2 - \sin^2\omega}}
                       {\cos\omega + \sqrt{m^2 - \sin^2\omega}}

        and the Mueller reflection matrix is

        .. math::

            \mathbf{R}_F = \frac{1}{2}
            \begin{pmatrix}
            r_l^2 + r_r^2 & r_l^2 - r_r^2 & 0 & 0 \\
            r_l^2 - r_r^2 & r_l^2 + r_r^2 & 0 & 0 \\
            0 & 0 & 2 r_l r_r & 0 \\
            0 & 0 & 0 & 2 r_l r_r
            \end{pmatrix}

        Parameters
        ----------
        angle : float
            Angle of incidence :math:`\omega` on the wave facet (rad).

        Returns
        -------
        numpy.ndarray
            4x4 Fresnel reflection matrix :math:`\mathbf{R}_F`.
        '''
        m = self.m

        racine = np.sqrt(m ** 2 - (np.sin(angle)) ** 2)
        rl = (racine - m ** 2 * np.cos(angle)) / (racine + m ** 2 * np.cos(angle))
        rr = (np.cos(angle) - racine) / (np.cos(angle) + racine)

        #  Rf_pol = reflexion matrix for illumination from above
        Rf_pol = np.zeros((4, 4))
        Rf_pol[0, 0] = rl ** 2 + rr ** 2
        Rf_pol[1, 1] = Rf_pol[0, 0]
        Rf_pol[0, 1] = rl ** 2 - rr ** 2
        Rf_pol[1, 0] = Rf_pol[0, 1]
        Rf_pol[2, 2] = 2 * rl * rr
        Rf_pol[3, 3] = 2 * rl * rr

        Rf_pol = 0.5 * Rf_pol

        return Rf_pol

    def scat_angle(self):
        r'''
        Scattering angle between the incident and the reflected directions.

        .. math::

            \cos\Theta = -\cos\theta_s\cos\theta_v - \sin\theta_s\sin\theta_v\cos\phi

        where :math:`\phi` is the relative azimuth (180° when Sun and sensor are
        in opposition).

        Returns
        -------
        float
            Scattering angle :math:`\Theta` (rad).
        '''
        sza = self.sza
        vza = self.vza
        azi = self.azi
        ang = -np.cos(sza) * np.cos(vza) - np.sin(sza) * np.sin(vza) * np.cos(azi)
        # ang = np.cos(np.pi - sza) * np.cos(vza) - np.sin(np.pi - sza) * np.sin(vza) * np.cos(azi)
        ang = np.arccos(ang)

        return ang

    def nu(self, sigx2, sigy2, cosphi2, theta):
        r'''
        Shadowing parameter :math:`\nu` ([RossDion2005]_; Eq. 15 of [RossDion2007]_).

        .. math::

            \nu = \frac{1}{\sqrt{2}\,\sigma\tan\theta}, \qquad
            \sigma^2 = \sigma_x^2\cos^2\phi + \sigma_y^2\,(1 - \cos^2\phi)

        Parameters
        ----------
        sigx2 : float
            Upwind slope variance :math:`\sigma_x^2`.
        sigy2 : float
            Crosswind slope variance :math:`\sigma_y^2`.
        cosphi2 : float
            :math:`\cos^2\phi`, the squared cosine of the azimuth of the
            direction relative to the upwind axis.
        theta : float
            Zenith angle :math:`\theta` (rad).

        Returns
        -------
        float
            :math:`\nu`
        '''

        sig = np.sqrt(sigx2 * cosphi2 + sigy2 * (1 - cosphi2))

        return 1 / (np.tan(theta) * np.sqrt(2) * sig)

    def Lambda(self, sigx2, sigy2, cosphi2, theta):
        r'''
        Smith shadowing function :math:`\Lambda` (Eq. 33b of [RossDion2005]_;
        Eq. 15 of [RossDion2007]_).

        .. math::

            \Lambda(\theta, \phi) = \frac{e^{-\nu^2} - \nu\sqrt{\pi}\,\mathrm{erfc}(\nu)}
                                         {2\,\nu\sqrt{\pi}}

        with :math:`\nu` given by :meth:`nu` and :math:`\phi` the azimuth of the
        direction relative to the upwind axis. With the Sun at azimuth 0, the
        viewing direction at relative azimuth :math:`\phi_v` and the wind at
        :math:`\phi_w`, :meth:`sunglint` uses :math:`\cos^2\phi = \cos^2\phi_w`
        for the Sun and :math:`\cos^2\phi = \cos^2(\phi_v - \phi_w)` for the
        viewing direction.

        Parameters
        ----------
        sigx2 : float
            Upwind slope variance :math:`\sigma_x^2`.
        sigy2 : float
            Crosswind slope variance :math:`\sigma_y^2`.
        cosphi2 : float
            :math:`\cos^2\phi`.
        theta : float
            Zenith angle :math:`\theta` (rad).

        Returns
        -------
        float
            :math:`\Lambda(\theta)`
        '''
        # no shadowing at nadir/zenith (nu -> infinity, Lambda -> 0)
        if np.tan(theta) == 0:
            return 0.
        piroot = np.sqrt(np.pi)
        nu = self.nu(sigx2, sigy2, cosphi2, theta)
        return (np.exp(-nu ** 2) - nu * piroot * special.erfc(nu)) / (2 * nu * piroot)
