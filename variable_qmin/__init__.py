from gwpopulation.utils import powerlaw, truncnorm, xp
from gwpopulation.models.mass import two_component_single, SinglePeakSmoothedMassDistribution, BaseSmoothedMassDistribution
import numpy as np
from bilby.core.prior import (
    Prior, PriorDict, ConditionalPriorDict,
    Uniform, ConditionalUniform, Constraint, DeltaFunction
)
import inspect

def zeta_CondUpperBound(reference_params,gamma,mmax,mmin):
    """
    Upper bound for zeta hyperparameter used in Increasing Parabola model
    """
    mmax_condition = (1.-gamma)/(mmax-mmin)
    
    return dict(minimum=0, maximum=mmax_condition)

def zeta_CondPrior(reference_params,gamma,mmax,mmin):
    """
    Upper and lower bounds for zeta hyperparameter used in Relaxed Parabola model
    """
    min_zeta = ((-gamma)/(mmax-mmin))-(mmin/((mmax-mmin)**2))

    mmax_condition = (1.-gamma)/(mmax-mmin)
    m2min_neg_zetamax = gamma**2/(4*mmin)
    
    return dict(minimum=m2min_neg_zetamax, maximum=mmax_condition)

def gamma_CondPrior(reference_params, mmax,mmin):
    """
    Upper and lower bounds for gamma (xi in paper) hyperparameter used in Relaxed Parabola model
    """
    gamma_min = -2*(mmin**0.5)/(mmax**0.5-mmin**0.5)
    gamma_max = 2*(mmin**0.5)/(mmax**0.5+mmin**0.5)
    return dict(minimum=gamma_min, maximum=gamma_max)

class SmoothedPowerlawm2min(BaseSmoothedMassDistribution):
    """
    Model class for Power law m2min model. 
    Inherits Power law plus Peak model for p(m_1) from two_component_single, and modifies
    p(q) to follow Power law m2min model.
    """
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
                p_q /= self.norm_p_q(beta=beta, mmax=mmax, mmin=mmin, delta_m=delta_m, gamma=gamma)
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
            self.m1s, masses, kind="linear", backend=xp
        )


class SmoothedParabolam2min(BaseSmoothedMassDistribution):
    """
    Model class for Parabola m2min model. 
    Inherits Power law plus Peak model for p(m_1) from two_component_single, and modifies
    p(q) to follow Parabola m2min model.

    Can be used with zeta_CondUpperBound to construct Increasing Parabola,
    or with zeta_CondPrior and gamma_CondPrior to construct Relaxed Parabola.
    """
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
            self.m1s, masses, kind="linear", backend=xp
        )
