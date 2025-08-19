# Variable minimum secondary mass models for use in `gwpopulation_pipe`

This codebase contains the model codes to allow a minimum secondary mass, and therefore mass ratio, that is variable with primary mass. This is included as a modification to a power-law mass-ratio distribution. These models are pip installable as custom models to be used with [`gwpopulation_pipe`](https://docs.ligo.org/RatesAndPopulations/gwpopulation_pipe/customizing.html).

This codebase contains two model variations, a Power law minimum secondary mass and a Parabola minimum secondary mass, are defined by `SmoothedPowerlawm2min` and `SmoothedParabolam2min` classes restpectively in `variable_qmin/__init__.py`.
This file also contains the functions defining the conditional priors on the hyperparameters used in the Parabola model.

These models were used in 'Can Big Black Holes Merge with the Smallest Black Holes?' (Colloms, Doctor, and Berry 2025) to investigate the correlation between binary black hole primary mass and mass ratio. We find that there is no support for an evolving minimum mass with GWTC-3 observations when excluding GW190814. When including GW190814, there is evidence for some structure in the mass-ratio and component-mass distributions.
