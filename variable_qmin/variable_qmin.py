from gwpopulation.utils import powerlaw, truncnorm, xp
from gwpopulation.models.mass import two_component_single
import numpy as np

def two_component_primary_mass_ratio_variable_qmin(
    dataset, alpha, beta, gamma, zeta, mmin, mmax, lam, mpp, sigpp, gaussian_mass_maximum=100
):
    r"""
    Power law model for two-dimensional mass distribution, modelling primary
    mass and conditional mass ratio distribution.

    .. math::
        p(m_1, q) = p(m1) p(q | m_1)

    Parameters
    ----------
    dataset: dict
        Dictionary of numpy arrays for 'mass_1' and 'mass_ratio'.
    alpha: float
        Negative power law exponent for more massive black hole.
    mmin: float
        Minimum black hole mass.
    mmax: float
        Maximum black hole mass.
    beta: float
        Power law exponent of the mass ratio distribution.
    gamma: float
        determines the value of qmin
    lam: float
        Fraction of black holes in the Gaussian component.
    mpp: float
        Mean of the Gaussian component.
    sigpp: float
        Standard deviation of the Gaussian component.
    gaussian_mass_maximum: float, optional
        Upper truncation limit of the Gaussian component. (default: 100)
    """
    from gwpopulation.utils import xp

    params = dict(
        mmin=mmin,
        mmax=mmax,
        lam=lam,
        mpp=mpp,
        sigpp=sigpp,
        gaussian_mass_maximum=gaussian_mass_maximum,
    )
    if zeta > (1-gamma)/(mmax-mmin):
        return np.zeros_like(dataset["mass_1"])
    m2min = mmin + gamma*(dataset["mass_1"] - mmin) + zeta * (dataset["mass_1"] - mmin)**2
    p_m1 = two_component_single(dataset["mass_1"], alpha=alpha, **params)
    p_q = powerlaw(dataset["mass_ratio"], beta, 1, m2min/dataset["mass_1"])
    prob = p_m1 * p_q
    return prob