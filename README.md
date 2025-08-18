# Variable minimum secondary mass models for use in `gwpopulation_pipe`

This codebase contains the model codes to allow a minimum secondary mass, and therefore mass ratio, that is variable with primary mass. These models are pip installable as custom models to be used with `gwpopulation_pipe`.

Two model variations, a Power law minimum secondary mass and a Parabola minimum secondary mass, are defined by `SmoothedPowerlawm2min` and `SmoothedParabolam2min` classes  in `variable_qmin/__init__.py` respectively. 
This file also contains the functions defining the conditional priors on the hyperparameters used in the Parabola model.

These models were used in 'Can Big Black Holes Merge with the Smallest Black Holes?' (Colloms, Doctor, and Berry 2025) with GWTC-3 observations.
