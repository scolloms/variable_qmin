from gwpopulation.utils import powerlaw, truncnorm, xp
from gwpopulation.models.mass import two_component_single
import numpy as np
from bilby.core.prior import (
    Prior, PriorDict, ConditionalPriorDict,
    Uniform, ConditionalUniform, Constraint, DeltaFunction
)
#from bilby.core.prior.analytical import DeltaFunction

def two_component_primary_mass_ratio_parabola_m2min(
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

def two_component_primary_mass_ratio_powerlaw_m2min(
    dataset, alpha, beta, gamma, mmin, mmax, lam, mpp, sigpp, gaussian_mass_maximum=100
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
    m2min = (mmax-mmin)*(((dataset["mass_1"]-mmin)/(mmax-mmin))**gamma)+mmin
    p_m1 = two_component_single(dataset["mass_1"], alpha=alpha, **params)
    p_q = powerlaw(dataset["mass_ratio"], beta, 1, m2min/dataset["mass_1"])
    prob = p_m1 * p_q
    return prob

def parabola_m2min_allgammazeta(
    dataset, alpha, beta, gamma, zeta, mmin, mmax, lam, mpp, sigpp, gaussian_mass_maximum=100
):
    r"""
    Power law model for two-dimensional mass distribution, modelling primary
    mass and conditional mass ratio distribution.

    .. math::
        p(m_1, q) = p(m1) p(q | m_1)
O4a_population/variable_qmin_inference/O4a-variablem2min.ini
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
    m2min = mmin + gamma*(dataset["mass_1"] - mmin) + zeta * (dataset["mass_1"] - mmin)**2
    p_m1 = two_component_single(dataset["mass_1"], alpha=alpha, **params)
    p_q = powerlaw(dataset["mass_ratio"], beta, 1, m2min/dataset["mass_1"])
    prob = p_m1 * p_q
    return prob

def zeta_CondPrior(reference_params,gamma,mmax,mmin):
    
    min_zeta = ((-gamma)/(mmax-mmin))-(mmin/((mmax-mmin)**2))

    mmax_condition = (1.-gamma)/(mmax-mmin)
    m2min_neg_zetamax = gamma**2/(4*mmin)
    
    return dict(minimum=m2min_neg_zetamax, maximum=mmax_condition)

def gamma_CondPrior(reference_params, mmax,mmin):
    gamma_min = -2*(mmin**0.5)/(mmax**0.5-mmin**0.5)
    gamma_max = 2*(mmin**0.5)/(mmax**0.5+mmin**0.5)
    return dict(minimum=gamma_min, maximum=gamma_max)

def zeta_with_conditional(mmax,mmin):
    """ Creates a prior dict with only gamma, mmin, mmax, and zeta with constraint"""
    zeta_prior = ConditionalPriorDict(
            dictionary=dict(
                gamma=ConditionalUniform(
                    condition_func=gamma_CondPrior,
                    minimum=-1, maximum=1
                ),
                mmax=mmax,
                mmin=mmin,
                zeta=ConditionalUniform(
                    condition_func=zeta_CondPrior, 
                    minimum=-0.017094017094, maximum=0.017094017094
                )
            )
        )
    return zeta_prior

def zeta_with_conditional(mmax,mmin):
    """ Creates a prior dict with only gamma, mmin, mmax, and zeta with constraint"""
    zeta_prior = PriorDict(
            dictionary=dict(
                gamma=Constraint(
                    condition_func=gamma_CondPrior,
                    minimum=-1, maximum=1
                ),
                mmax=mmax,
                mmin=mmin,
                zeta=ConditionalUniform(
                    condition_func=zeta_CondPrior, 
                    minimum=-0.017094017094, maximum=0.017094017094
                )
            )
        )
    return zeta_prior


