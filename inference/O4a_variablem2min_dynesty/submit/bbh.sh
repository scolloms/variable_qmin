#! /bin/bash

echo 'Moving to /home/storm.colloms/O4a_population/test_gwpopulationpipe'
cd /home/storm.colloms/O4a_population/test_gwpopulationpipe

/home/storm.colloms/.conda/envs/O4apop/bin/gwpopulation_pipe_collection O4a_variablem2min_dynesty/bbh_config_complete.ini --run-dir O4a_variablem2min_dynesty

/home/storm.colloms/.conda/envs/O4apop/bin/gwpopulation_pipe_analysis O4a_variablem2min_dynesty/bbh_config_complete.ini --run-dir O4a_variablem2min_dynesty --label bbh_mass_two_component_primary_mass_ratio_variable_qmin_redshift_powerlaw --models mass:variable_qmin.variable_qmin.two_component_primary_mass_ratio_variable_qmin --models redshift:gwpopulation.models.redshift.PowerLawRedshift --vt-models mass:variable_qmin.variable_qmin.two_component_primary_mass_ratio_variable_qmin --vt-models redshift:gwpopulation.models.redshift.PowerLawRedshift --vt-file O4a_variablem2min_dynesty/data/injections.pkl

/home/storm.colloms/.conda/envs/O4apop/bin/gwpopulation_pipe_plot O4a_variablem2min_dynesty/bbh_config_complete.ini --run-dir O4a_variablem2min_dynesty --result-file O4a_variablem2min_dynesty/result/bbh_mass_two_component_primary_mass_ratio_variable_qmin_redshift_powerlaw_result.hdf5 --samples O4a_variablem2min_dynesty/result/bbh_mass_two_component_primary_mass_ratio_variable_qmin_redshift_powerlaw_samples.pkl 

/home/storm.colloms/.conda/envs/O4apop/bin/gwpopulation_pipe_to_common_format --result-file O4a_variablem2min_dynesty/result/bbh_mass_two_component_primary_mass_ratio_variable_qmin_redshift_powerlaw_result.hdf5 --n-samples 5000 --max-redshift 1.9 --minimum-mass 2.0 --maximum-mass 100.0 --injection-file O4a_variablem2min_dynesty/data/injections.pkl --filename O4a_variablem2min_dynesty/result/bbh_mass_two_component_primary_mass_ratio_variable_qmin_redshift_powerlaw_full_posterior.hdf5 --samples-file O4a_variablem2min_dynesty/result/bbh_mass_two_component_primary_mass_ratio_variable_qmin_redshift_powerlaw_samples.pkl --vt-ifar-threshold 1.0 --vt-snr-threshold 10.0 --backend numpy 

/home/storm.colloms/.conda/envs/O4apop/bin/summarypages --config O4a_variablem2min_dynesty/bbh_config_complete.ini --webdir O4a_variablem2min_dynesty/summary --samples O4a_variablem2min_dynesty/result/bbh_mass_two_component_primary_mass_ratio_variable_qmin_redshift_powerlaw_result.hdf5 --labels bbh_mass_two_component_primary_mass_ratio_variable_qmin_redshift_powerlaw

cd -
