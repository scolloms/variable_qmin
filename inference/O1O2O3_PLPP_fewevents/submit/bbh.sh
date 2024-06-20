#! /bin/bash

echo 'Moving to /home/storm.colloms/O4a_population/variable_qmin_inference'
cd /home/storm.colloms/O4a_population/variable_qmin_inference

/home/storm.colloms/.conda/envs/O4apop-jax/bin/gwpopulation_pipe_collection O1O2O3_PLPP_fewevents/bbh_config_complete.ini --run-dir O1O2O3_PLPP_fewevents

/home/storm.colloms/.conda/envs/O4apop-jax/bin/gwpopulation_pipe_analysis O1O2O3_PLPP_fewevents/bbh_config_complete.ini --run-dir O1O2O3_PLPP_fewevents --label bbh_mass_two_component_primary_mass_ratio_redshift_powerlaw --models mass:two_component_primary_mass_ratio --models redshift:gwpopulation.models.redshift.PowerLawRedshift --vt-models mass:two_component_primary_mass_ratio --vt-models redshift:gwpopulation.models.redshift.PowerLawRedshift --vt-file O1O2O3_PLPP_fewevents/data/injections.pkl

/home/storm.colloms/.conda/envs/O4apop-jax/bin/gwpopulation_pipe_plot O1O2O3_PLPP_fewevents/bbh_config_complete.ini --run-dir O1O2O3_PLPP_fewevents --result-file O1O2O3_PLPP_fewevents/result/bbh_mass_two_component_primary_mass_ratio_redshift_powerlaw_result.hdf5 --samples O1O2O3_PLPP_fewevents/result/bbh_mass_two_component_primary_mass_ratio_redshift_powerlaw_samples.pkl 

/home/storm.colloms/.conda/envs/O4apop-jax/bin/gwpopulation_pipe_to_common_format --result-file O1O2O3_PLPP_fewevents/result/bbh_mass_two_component_primary_mass_ratio_redshift_powerlaw_result.hdf5 --n-samples 5000 --max-redshift 1.9 --minimum-mass 2.0 --maximum-mass 100.0 --injection-file O1O2O3_PLPP_fewevents/data/injections.pkl --filename O1O2O3_PLPP_fewevents/result/bbh_mass_two_component_primary_mass_ratio_redshift_powerlaw_full_posterior.hdf5 --samples-file O1O2O3_PLPP_fewevents/result/bbh_mass_two_component_primary_mass_ratio_redshift_powerlaw_samples.pkl --vt-ifar-threshold 1.0 --vt-snr-threshold 10.0 --backend jax 

cd -
