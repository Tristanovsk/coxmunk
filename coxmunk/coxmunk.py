# coding=utf-8
import numpy as np


class sunglint:

    def __init__(self, sza, vza, azi, m=1.334, tau_atm=0):
        '''

        :param sza:
        :param vza:
        :param azi: relative azimuth in deg. for convention 0° when Sun and sensor in opposition
        :param m:
        :param tau_atm:
        '''
        degrad = np.pi / 180.
        self.sza = sza * degrad
        self.vza = vza * degrad
        self.azi = (180 - azi) * degrad
        self.m = m
        self.tau_atm = tau_atm

    def sunglint(self, ws, wazi, stats='cm_dir', shadow=True):
        '''stats in `cm_iso`, `cm_dir`, `bh2006` '''
        sza = self.sza
        vza = self.vza
        azi = self.azi
        wazi = wazi* np.pi / 180.

        muv = np.cos(vza)
        sinv = np.sin(vza)
        mu0 = np.cos(sza)
        sin0 = np.sin(sza)
        muazi = np.cos(azi)
        sinazi = np.sin(azi)
        scat = self.scat_angle()
        omega = (np.pi - scat) / 2  # for x in self.scat_angle()]
        thetaN = np.arccos((mu0 + muv) / (2. * np.cos(omega)))

        Rf = self.fresnel(omega)

        if stats == 'cm_iso':  # Original isotropic Cox Munk statistics
            sigma2 = 3e-3 + 512e-5 * ws
        elif stats == 'cm_dir':  # historical values from directional COX MUNK
            sc2 = 0.003 + 1.92e-3 * ws
            su2 = 3.16e-3 * ws
            sc = np.sqrt(sc2)
            su = np.sqrt(su2)
            c21 = 0.01 - 8.6e-3 * ws
            c03 = 0.04 - 33.e-3 * ws
            c40 = 0.40
            c22 = 0.12
            c04 = 0.23
        elif stats == 'bh2006':  # from Breon Henriot 2006 JGR
            sc2 = 3e-3 + 1.85e-3 * ws
            su2 = 1e-3 + 3.16e-3 * ws
            sc = np.sqrt(sc2)
            su = np.sqrt(su2)
            c21 = -9e-4 * ws ** 2
            c03 = -0.45 / (1 + np.exp(7. - ws))
            c40 = 0.30
            c22 = 0.12
            c04 = 0.4

        if stats == 'cm_iso':
            Pdist = 1. / (np.pi * sigma2) * np.exp(-1. * np.tan(thetaN) ** 2 / sigma2)
        else:
            sigma2 = sc2 + su2

            zy = 1 * (sinv * muazi + sin0) / (mu0 + muv)
            zx = -1 * (sinv * sinazi) / (mu0 + muv)

            zx_prime = np.cos(wazi) * zx + np.sin(wazi) * zy
            zy_prime = -np.sin(wazi) * zx + np.cos(wazi) * zy

            xi = zx_prime / sc
            eta = zy_prime / su

            Pdist = np.exp(-5e-1 * (xi ** 2 + eta ** 2)) / (2. * np.pi * sc * su) * \
                    (1. -
                     c21 * (xi ** 2 - 1.) * eta / 2. -
                     c03 * (eta ** 3 - 3. * eta) / 6. +
                     c40 * (xi ** 4 - 6. * eta ** 2 + 3.) / 24. +
                     c04 * (eta ** 4 - 6. * eta ** 2 + 3.) / 24. +
                     c22 * (xi ** 2 - 1.) * (eta ** 2 - 1.) / 4.)

        # ---------------------------------------------------------------------*
        #                    Rotation in the reference plane
        # ---------------------------------------------------------------------*
        if vza != 0 and sza != 0 and scat != 0:
            L1 = np.zeros((4, 4))
            L2 = np.zeros((4, 4))
            L1[0, 0] = 1.
            L1[3, 3] = 1.
            L2[0, 0] = 1.
            L2[3, 3] = 1.

            cos_sigma1 = (np.cos(np.pi - sza) - muv * np.cos(scat)) / (sinv * np.sin(scat))
            cos_sigma2 = (muv - np.cos(np.pi - sza) * np.cos(scat)) / (np.sin(np.pi - sza) * np.sin(scat))

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

        # TODO: add the shadowing parameterization
        SH = 1.

        # ---------------------------------------------------------------------*
        #        Sun glint Stokes component (in reflectance unit) at TOA
        # ---------------------------------------------------------------------*

        Pdist = SH * (np.pi * Pdist) / (4. * mu0 * muv * np.cos(thetaN) ** 4)

        Td = 1.
        Tu = 1.
        Pdist = Td * Tu * Pdist
        Iglint = Rf[0, 0] * Pdist
        Qglint = Rf[1, 0] * Pdist
        Uglint = Rf[2, 0] * Pdist

        return Iglint

    def atmo_trans(self, tau, sza, Iglint):
        # ---------------------------------------------------------------------*
        #                           Direct Transmittance
        # ---------------------------------------------------------------------*
        sza
        Td = np.exp(-1 * tau / np.cos(sza))
        # Tup = np.exp(-1 * (self.tau_atm) / np.cos(self.vza))

        return Iglint * Td

    def fresnel(self, angle):
        ''' 
        :param angle: incident angle on the wave facets (in rad)
        :return: 
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
        sza = self.sza
        vza = self.vza
        azi = self.azi
        ang = np.cos(np.pi - sza) * np.cos(vza) - np.sin(np.pi - sza) * np.sin(vza) * np.cos(azi)
        ang = np.arccos(ang)

        return ang
