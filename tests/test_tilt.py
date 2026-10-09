import numpy as np
import pytest

from coxmunk import sunglint


def test_no_tilt_unchanged():
    g = sunglint(30, 40, 120)
    assert (g.sza, g.vza, g.azi, g.psi) == (g.sza_h, g.vza_h, g.azi_h, 0.)


@pytest.mark.parametrize('tilt, tilt_azi, expected', [
    (10, 0, (30, 40, 180)),  # plane facing the Sun
    (10, 180, (50, 20, 180)),  # plane facing away from the Sun
    (40, 0, (0, 70, None)),  # Sun at the zenith of the plane: azimuth undefined
])
def test_local_angles_principal_plane(tilt, tilt_azi, expected):
    g = sunglint(40, 30, 180, tilt=tilt, tilt_azi=tilt_azi)
    np.testing.assert_allclose(np.degrees([g.sza, g.vza]), expected[:2], atol=1e-9)
    if expected[2] is not None:
        assert np.degrees(g.azi) == pytest.approx(expected[2])


def test_tilted_glint_equals_horizontal_with_local_angles():
    g = sunglint(40, 30, 90, tilt=20, tilt_azi=90, pol_ref='local')
    h = sunglint(*np.degrees([g.sza, g.vza, g.azi]))
    np.testing.assert_allclose(g.sunglint(7, 30, stats='bh2006'), h.sunglint(7, 30, stats='bh2006'))


def test_hidden_surface():
    g = sunglint(80, 30, 180, tilt=20, tilt_azi=180)
    assert not g.visible
    assert g.sunglint(5) == [0., 0., 0., 0.]


def _aop_expected(sza, vza, azi):
    # light reflected by a facet is polarized perpendicular to the plane (s, v);
    # angle measured from the meridian plane of the viewing direction,
    # counterclockwise when looking from the sensor toward the surface
    sza, vza, azi = np.radians([sza, vza, azi])
    s = np.array([np.sin(sza), 0., np.cos(sza)])
    v = np.array([np.sin(vza) * np.cos(azi), np.sin(vza) * np.sin(azi), np.cos(vza)])
    e = np.cross(s, v)
    t = np.array([0., 0., 1.]) - v[2] * v
    return np.arctan2(np.dot(np.cross(t, e), -v), np.dot(t, e))


@pytest.mark.parametrize('tilt', [0, 15, 30])
def test_angle_of_polarization_global_frame(tilt):
    rng = np.random.default_rng(0)
    for _ in range(50):
        sza, vza, azi, tilt_azi = rng.uniform(5, 60), rng.uniform(5, 60), rng.uniform(0, 360), rng.uniform(0, 360)
        g = sunglint(sza, vza, azi, tilt=tilt, tilt_azi=tilt_azi)
        if not g.visible:
            continue
        I, Q, U, V = g.sunglint(5, stats='cm_iso', shadow=False)
        if np.hypot(Q, U) < 1e-8:
            continue
        diff = 0.5 * np.arctan2(U, Q) - _aop_expected(sza, vza, azi)
        assert abs((diff + np.pi / 2) % np.pi - np.pi / 2) < 1e-8
