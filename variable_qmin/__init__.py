from gwpopulation.utils import powerlaw, truncnorm, xp
from gwpopulation.models.mass import two_component_single, SinglePeakSmoothedMassDistribution, BaseSmoothedMassDistribution
import numpy as np
from bilby.core.prior import (
    Prior, PriorDict, ConditionalPriorDict,
    Uniform, ConditionalUniform, Constraint, DeltaFunction
)
import scipy.special as scs
import inspect

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

def powerlaw_m2min(
    dataset, alpha, beta, gamma, mmin, mmax, delta_m, lam, mpp, sigpp, gaussian_mass_maximum=100
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
        alpha=alpha,
        mmin=mmin,
        mmax=mmax,
        lam=lam,
        mpp=mpp,
        sigpp=sigpp,
        gaussian_mass_maximum=gaussian_mass_maximum,
    )
    m2min = (mmax-mmin)*(((dataset["mass_1"]-mmin)/(mmax-mmin))**gamma)+mmin
    p_m1 = two_component_single(dataset["mass_1"],  **params)
    p_q = powerlaw(dataset["mass_ratio"], beta, 1, m2min/dataset["mass_1"])
    prob = p_m1 * p_q
    return prob

def powerlaw_m2min_smoothed(
    dataset, alpha, beta, gamma, mmin, mmax, delta_m, lam, mpp, sigpp, gaussian_mass_maximum=100
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
        alpha=alpha,
        mmin=mmin,
        mmax=mmax,
        lam=lam,
        mpp=mpp,
        sigpp=sigpp,
        gaussian_mass_maximum=gaussian_mass_maximum,
    )
    m2min = (mmax-mmin)*(((dataset["mass_1"]-mmin)/(mmax-mmin))**gamma)+mmin
    p_m1 = two_component_single(dataset["mass_1"],  **params)
    p_m1 *= smoothing(dataset["mass_1"], mmin, mmax, delta_m)
    p_q = powerlaw(dataset["mass_ratio"], beta, 1, m2min/dataset["mass_1"])
    p_q *= smoothing(
            dataset["mass_1"] * dataset["mass_ratio"],
            mmin=m2min,
            mmax=dataset["mass_1"],
            delta_m=delta_m,
        )
    prob = p_m1 * p_q
    return prob

def parabola_m2min_allgammazeta(
    dataset, alpha, beta, gamma, zeta, mmin, mmax, delta_m, lam, mpp, sigpp, gaussian_mass_maximum=100
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
        alpha=alpha,
        mmin=mmin,
        mmax=mmax,
        lam=lam,
        mpp=mpp,
        sigpp=sigpp,
        gaussian_mass_maximum=gaussian_mass_maximum,
    )
    m2min = mmin + gamma*(dataset["mass_1"] - mmin) + zeta * (dataset["mass_1"] - mmin)**2
    p_m1 = two_component_single(dataset["mass_1"],  **params)
    p_q = powerlaw(dataset["mass_ratio"], beta, 1, m2min/dataset["mass_1"])
    prob = p_m1 * p_q
    return prob

def parabola_m2min_allgammazeta_smoothed(
    dataset, alpha, beta, gamma, zeta, mmin, mmax, delta_m, lam, mpp, sigpp, gaussian_mass_maximum=100
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
        alpha=alpha,
        mmin=mmin,
        mmax=mmax,
        lam=lam,
        mpp=mpp,
        sigpp=sigpp,
        gaussian_mass_maximum=gaussian_mass_maximum,
    )
    m2min = mmin + gamma*(dataset["mass_1"] - mmin) + zeta * (dataset["mass_1"] - mmin)**2
    p_m1 = two_component_single(dataset["mass_1"],  **params)
    p_m1 *= smoothing(dataset["mass_1"], mmin, mmax, delta_m)
    p_q = powerlaw(dataset["mass_ratio"], beta, 1, m2min/dataset["mass_1"])
    p_q *= smoothing(
            dataset["mass_1"] * dataset["mass_ratio"],
            mmin=m2min,
            mmax=dataset["mass_1"],
            delta_m=delta_m,
        )
    prob = p_m1 * p_q
    return prob
    
def smoothing(masses, mmin, mmax, delta_m):
    """
    Apply a one sided window between mmin and mmin + delta_m to the
    mass pdf.

    The upper cut off is a step function,
    the lower cutoff is a logistic rise over delta_m solar masses.

    See T&T18 Eqs 7-8
    Note that there is a sign error in that paper.

    S = (f(m - mmin, delta_m) + 1)^{-1}
    f(m') = delta_m / m' + delta_m / (m' - delta_m)

    See also, https://en.wikipedia.org/wiki/Window_function#Planck-taper_window
    """
    from gwpopulation.utils import xp
    if delta_m > 0.0:
        shifted_mass = np.nan_to_num((masses - mmin) / delta_m, nan=0)
        shifted_mass = np.clip(shifted_mass, 1e-6, 1 - 1e-6)
        exponent = 1 / shifted_mass - 1 / (1 - shifted_mass)
        window = scs.expit(-exponent)
        window *= (masses >= mmin) * (masses <= mmax)
        return window
    else:
        return xp.ones(masses.shape)

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

class SmoothedPowerlawm2min(BaseSmoothedMassDistribution):

    primary_model=two_component_single

    @property
    def variable_names(self):
        vars = getattr(
            self.primary_model,
            "variable_names",
            inspect.getfullargspec(self.primary_model).args[1:],
        )
        vars += ["beta", "delta_m", "gamma"]
        vars = set(vars).difference(self.kwargs.keys())
        return vars

    def __call__(self, dataset, *args, **kwargs):
        beta = kwargs.pop("beta")
        gamma = kwargs.pop("gamma")
        mmin = kwargs.get("mmin", self.mmin)
        mmax = kwargs.get("mmax", self.mmax)
        if "jax" not in xp.__name__:
            if mmin < self.mmin:
                raise ValueError(
                    "{self.__class__}: mmin ({mmin}) < self.mmin ({self.mmin})"
                )
            if mmax > self.mmax:
                raise ValueError(
                    "{self.__class__}: mmax ({mmax}) > self.mmax ({self.mmax})"
                )
        delta_m = kwargs.get("delta_m", 0)

        pm1_kwargs = dict(
            alpha=kwargs.pop("alpha"),
            mmin=mmin,
            mmax=mmax,
            lam=kwargs.pop("lam"),
            mpp=kwargs.pop("mpp"),
            sigpp=kwargs.pop("sigpp"),
            gaussian_mass_maximum=kwargs.pop("gaussian_mass_maximum"),
        )
        p_m1 = self.p_m1(dataset, **pm1_kwargs, **self.kwargs)
        p_q = self.p_q(dataset, beta=beta, mmax=mmax, mmin=mmin, delta_m=delta_m, gamma=gamma)
        prob = p_m1 * p_q
        return prob

    def p_q(self, dataset, beta, mmax, mmin, delta_m, gamma):
        from gwpopulation.utils import xp
        m2min = (mmax-mmin)*(((dataset["mass_1"]-mmin)/(mmax-mmin))**gamma)+mmin
        p_q = powerlaw(dataset["mass_ratio"], beta, 1, m2min / dataset["mass_1"])
        p_q *= self.smoothing(
            dataset["mass_1"] * dataset["mass_ratio"],
            mmin=m2min,
            mmax=dataset["mass_1"],
            delta_m=delta_m,
        )
        
        try:
            if self.cache:
                p_q /= self.norm_p_q(beta=beta, mmin=mmin, delta_m=delta_m)
            else:
                self._cache_q_norms(dataset["mass_1"])
                p_q /= self.norm_p_q(beta=beta, mmax=mmax, mmin=mmin, delta_m=delta_m, gamma=gamma)
        except (AttributeError, TypeError, ValueError):
            self._cache_q_norms(dataset["mass_1"])
            p_q /= self.norm_p_q(beta=beta, mmax=mmax, mmin=mmin, delta_m=delta_m, gamma=gamma)

        return xp.nan_to_num(p_q)

    def norm_p_q(self, beta, mmax, mmin, delta_m, gamma):
        from gwpopulation.utils import xp
        """Calculate the mass ratio normalisation by linear interpolation"""
        m2min = (mmax-mmin)*(((self.m1s_grid-mmin)/(mmax-mmin))**gamma)+mmin
        p_q = powerlaw(self.qs_grid, beta, 1, m2min / self.m1s_grid)
        p_q *= self.smoothing(
            self.m1s_grid * self.qs_grid, mmin=m2min, mmax=self.m1s_grid, 
            delta_m=delta_m)

        norms = xp.nan_to_num(xp.trapz(p_q, self.qs, axis=0)) * (delta_m != 0) + 1 * (
            delta_m == 0
        )

        return self._q_interpolant(norms)

    def _cache_q_norms(self, masses):
        """
        Cache the information necessary for linear interpolation of the mass
        ratio normalisation
        """
        from gwpopulation.models.interped import _setup_interpolant

        self._q_interpolant = _setup_interpolant(
            self.m1s, masses, kind="cubic", backend=xp
        )


class SmoothedParabolam2min(BaseSmoothedMassDistribution):

    primary_model=two_component_single

    @property
    def variable_names(self):
        vars = getattr(
            self.primary_model,
            "variable_names",
            inspect.getfullargspec(self.primary_model).args[1:],
        )
        vars += ["beta", "delta_m", "gamma", "zeta"]
        vars = set(vars).difference(self.kwargs.keys())
        return vars

    def __call__(self, dataset, *args, **kwargs):
        beta = kwargs.pop("beta")
        gamma = kwargs.pop("gamma")
        zeta = kwargs.pop("zeta")
        mmin = kwargs.get("mmin", self.mmin)
        mmax = kwargs.get("mmax", self.mmax)
        """if "jax" not in xp.__name__:
            if mmin < self.mmin:
                raise ValueError(
                    "{self.__class__}: mmin ({mmin}) < self.mmin ({self.mmin})"
                )
            if mmax > self.mmax:
                raise ValueError(
                    "{self.__class__}: mmax ({mmax}) > self.mmax ({self.mmax})"
                )"""
        delta_m = kwargs.get("delta_m", 0)

        pm1_kwargs = dict(
            alpha=kwargs.pop("alpha"),
            mmin=mmin,
            mmax=mmax,
            lam=kwargs.pop("lam"),
            mpp=kwargs.pop("mpp"),
            sigpp=kwargs.pop("sigpp"),
            gaussian_mass_maximum=kwargs.pop("gaussian_mass_maximum"),
        )
        p_m1 = self.p_m1(dataset, **pm1_kwargs, **self.kwargs)
        p_q = self.p_q(dataset, beta=beta, mmax=mmax, mmin=mmin, delta_m=delta_m, gamma=gamma, zeta=zeta)
        prob = p_m1 * p_q
        return prob

    def p_q(self, dataset, beta, mmax, mmin, delta_m, gamma, zeta):
        from gwpopulation.utils import xp
        
        m2min = mmin + gamma*(dataset["mass_1"] - mmin) + zeta * (dataset["mass_1"] - mmin)**2
        p_q = powerlaw(dataset["mass_ratio"], beta, 1, m2min / dataset["mass_1"])
        p_q *= self.smoothing(
            dataset["mass_1"] * dataset["mass_ratio"],
            mmin=m2min,
            mmax=dataset["mass_1"],
            delta_m=delta_m,
        )
        
        try:
            if self.cache:
                p_q /= self.norm_p_q(beta=beta, mmax=mmax, mmin=mmin, delta_m=delta_m, gamma=gamma, zeta=zeta)
            else:
                self._cache_q_norms(dataset["mass_1"])
                p_q /= self.norm_p_q(beta=beta, mmax=mmax, mmin=mmin, delta_m=delta_m, gamma=gamma, zeta=zeta)
        except (AttributeError, TypeError, ValueError):
            self._cache_q_norms(dataset["mass_1"])
            p_q /= self.norm_p_q(beta=beta, mmax=mmax, mmin=mmin, delta_m=delta_m, gamma=gamma, zeta=zeta)

        return xp.nan_to_num(p_q)

    def norm_p_q(self, beta, mmax, mmin, delta_m, gamma, zeta):
        from gwpopulation.utils import xp
        
        m2min = mmin + gamma*(self.m1s_grid - mmin) + zeta * (self.m1s_grid - mmin)**2
        
        """Calculate the mass ratio normalisation by linear interpolation"""
        p_q = powerlaw(self.qs_grid, beta, 1, m2min / self.m1s_grid)
        p_q *= self.smoothing(
            self.m1s_grid * self.qs_grid, mmin=m2min, mmax=self.m1s_grid, delta_m=delta_m
        )

        norms = xp.nan_to_num(xp.trapz(p_q, self.qs, axis=0)) * (delta_m != 0) + 1 * (
            delta_m == 0
        )

        return self._q_interpolant(norms)

    def _cache_q_norms(self, masses):
        """
        Cache the information necessary for linear interpolation of the mass
        ratio normalisation
        """
        from gwpopulation.models.interped import _setup_interpolant

        self._q_interpolant = _setup_interpolant(
            self.m1s, masses, kind="cubic", backend=xp
        )
