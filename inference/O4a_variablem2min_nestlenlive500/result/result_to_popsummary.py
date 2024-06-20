from matplotlib import pyplot as plt
import numpy as np
import json
import corner
import seaborn
import pandas
import popsummary

def json_to_popsummary_hyperposterior(json_filename):
    """
    function to convert hyperposterior samples from the json format
    to popsummary h5
    """
    with open(json_filename,'r') as f:
        res = json.load(f)
    hyperparameters = list(res['posterior']['content'].keys())
    out = np.array(list(map(lambda hyperparam: res['posterior']['content'][hyperparam],hyperparameters))).T
    return out, hyperparameters

hyp_out, hyperparameters = json_to_popsummary_hyperposterior(f'{rundir}/gauss_result.json')

popfile = popsummary.popresult.PopulationResult(
    '.h5',
    hyperparameters = hyperparameters,
    events = events,
    event_parameters = event_vars,
        )
popfile.set_hyperparameter_samples(hyp_out, overwrite=True)