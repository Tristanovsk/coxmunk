# coding=utf-8
import numpy as np
from scipy import special

class sunglint:

    def __init__(self, sza, vza, azi, m=1.334, tau_atm=0):
        '''

        :param sza: solar zenith angle in deg.
        :param vza: viewing zenith angle in deg.
        :param azi: viewing azimuth in deg. for convention 180° when Sun and sensor in opposition;
                    Sun azimuth set to 0 in the reference frame
        :param m: refractive index of water
        :param tau_atm: TODO implement transmittances with optical thickness
        '''
        degrad = np.pi / 180.
        self.sza = sza * degrad
        self.vza = vza * degrad
        self.azi = azi * degrad
        self.m = m
        self.tau_atm = tau_atm

    def sunglint(self, ws, wazi=0, stats='cm_dir', shadow=True):
        '''

        :param ws: wind speed in m/s
        :param wazi: wind direction (in deg.) with respect to Sun direction (e.g., wazi=0° wind toward the Sun)
        :param stats: in `cm_iso`, `cm_dir`, `bh2006`
        :param shadow: if True, apply shadow correction based on
        :return:
        '''
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
            sc2 = sigma2/2
            su2 = sigma2/2
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
            # zy = muazi*np.tan(thetaN)
            # zx = sinazi*np.tan(thetaN)
            zy = -1 * (sinv * muazi + sin0) / (mu0 + muv)
            zx = -1 * (sinv * sinazi) / (mu0 + muv)

            # print('Zy',zy,zy2)
            # print('Zx',zx,zx2)
            # print('N',np.tan(thetaN),np.sqrt(zy2**2+zx2**2),np.sqrt(zy**2+zx**2))

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
        if shadow and vza!=0:
            Ls=self.Lambda(su2,sc2,1,sza)
            Lr=self.Lambda(su2,sc2,muazi**2,vza)
            print(su2,sc2,muazi**2,vza,Ls,Lr)
            SH = 1 / (1+Ls+Lr)
        else:
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
        '''
        self.azi: azimuth in rad for convention azi=180 when sun-sensenor in oppositioon
        :return: scattering angle in rad
        '''
        sza = self.sza
        vza = self.vza
        azi = self.azi
        ang = -np.cos(sza) * np.cos(vza) - np.sin(sza) * np.sin(vza) * np.cos(azi)
        #ang = np.cos(np.pi - sza) * np.cos(vza) - np.sin(np.pi - sza) * np.sin(vza) * np.cos(azi)
        ang = np.arccos(ang)

        return ang

    def nu(self,sigx2,sigy2,cosphi2,theta):
        '''
        From Ross & Dion, 2005 and Eq. 15 Ross & Dion, 2007
        :param sigx2
        :param sigy2:
        :param cosphi2:
        :param theta:
        :return:
        '''

        sig = np.sqrt(sigx2*cosphi2+sigy2*(1-cosphi2))

        return 1/(np.tan(theta)*np.sqrt(2)*sig)

    def Lambda(self,sigx2,sigy2,cosphi2,theta):
        '''
        From Eq. 33b Ross & Dion, 2005 and Eq. 15 Ross & Dion, 2007
        :param sigx2: upwind variance
        :param sigy2:crosswind variance
        :param cosphi2: square of cos phi
        :param theta: zenith angle (rad)
        :return: Lambda
        '''
        piroot = np.sqrt(np.pi)
        nu = self.nu(sigx2,sigy2,cosphi2,theta)
        return (np.exp(-nu**2)-nu*piroot*special.erfc(nu)) / (2*nu*piroot)