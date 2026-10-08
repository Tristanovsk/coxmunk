import numpy as np
import pytest

from coxmunk import sunglint
from coxmunk.wind import glint_whitecap, glint_directions, retrieve_wind_speed, whitecap_coverage

# PARASOL-like multidirectional acquisition: 14 directions along a track
# 30 deg off the principal plane, Sun at 30 deg
VZA = np.linspace(-55, 55, 14)
GEOM = np.array([(30, abs(v), 150 if v >= 0 else 330) for v in VZA])


def test_whitecap_coverage():
    assert whitecap_coverage(0) == 0
    assert whitecap_coverage(10) == pytest.approx(2.95e-6 * 10 ** 3.52)


def test_glint_whitecap_reduces_to_glint_without_foam():
    S = glint_whitecap(GEOM[:3], 5, foam_reflectance=0)
    ff = whitecap_coverage(5)
    ref = [sunglint(*g).sunglint(5, 0, stats='bh2006', shadow=True)[:3] for g in GEOM[:3]]
    np.testing.assert_allclose(S, (1 - ff) * np.array(ref))


@pytest.mark.parametrize('ws', [2., 7., 12.])
def test_retrieval_noise_free(ws):
    res = retrieve_wind_speed(glint_whitecap(GEOM, ws), GEOM)
    assert not res.flag
    assert res.ws == pytest.approx(ws, abs=1e-3)


def test_retrieval_noisy():
    rng = np.random.default_rng(0)
    obs = glint_whitecap(GEOM, 8.)
    res = retrieve_wind_speed(obs + rng.normal(0, 4e-4, obs.shape), GEOM)
    assert res.ws == pytest.approx(8., abs=0.3)
    assert 0 < res.sigma < 1


def test_glint_directions():
    flags = glint_directions(GEOM, 5.)
    # directions on the Sun side of the track are free of sunglint
    assert not flags[0]
    assert flags[-1]


def test_lambda_no_nan_at_zenith():
    assert np.all(np.isfinite(sunglint(0, 30, 180).sunglint(7, 0, stats='bh2006', shadow=True)))
