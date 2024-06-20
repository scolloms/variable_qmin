#! /bin/bash

echo 'Moving to /home/storm.colloms/O4a_population/variable_qmin_inference'
cd /home/storm.colloms/O4a_population/variable_qmin_inference

/home/storm.colloms/.conda/envs/O4apop/bin/gwpopulation_pipe_collection O1O2O3_parabolam2min_prod_allgammazeta/bbh_config_complete.ini --run-dir O1O2O3_parabolam2min_prod_allgammazeta

/home/storm.colloms/.conda/envs/O4apop/bin/gwpopulation_pipe_analysis O1O2O3_parabolam2min_prod_allgammazeta/bbh_config_complete.ini --run-dir O1O2O3_parabolam2min_prod_allgammazeta --label bbh_mass_parabola_m2min_allgammazeta_redshift_powerlaw --models mass:variable_qmin.parabola_m2min_allgammazeta --models redshift:gwpopulation.models.redshift.PowerLawRedshift --vt-models mass:variable_qmin.parabola_m2min_allgammazeta --vt-models redshift:gwpopulation.models.redshift.PowerLawRedshift --vt-file O1O2O3_parabolam2min_prod_allgammazeta/data/injections.pkl

/home/storm.colloms/.conda/envs/O4apop/bin/gwpopulation_pipe_plot O1O2O3_parabolam2min_prod_allgammazeta/bbh_config_complete.ini --run-dir O1O2O3_parabolam2min_prod_allgammazeta --result-file O1O2O3_parabolam2min_prod_allgammazeta/result/bbh_mass_parabola_m2min_allgammazeta_redshift_powerlaw_result.hdf5 --samples O1O2O3_parabolam2min_prod_allgammazeta/result/bbh_mass_parabola_m2min_allgammazeta_redshift_powerlaw_samples.pkl 

/home/storm.colloms/.conda/envs/O4apop/bin/gwpopulation_pipe_to_common_format --result-file O1O2O3_parabolam2min_prod_allgammazeta/result/bbh_mass_parabola_m2min_allgammazeta_redshift_powerlaw_result.hdf5 --n-samples 5000 --max-redshift 1.9 --minimum-mass 2.0 --maximum-mass 100.0 --injection-file O1O2O3_parabolam2min_prod_allgammazeta/data/injections.pkl --filename O1O2O3_parabolam2min_prod_allgammazeta/result/bbh_mass_parabola_m2min_allgammazeta_redshift_powerlaw_full_posterior.hdf5 --samples-file O1O2O3_parabolam2min_prod_allgammazeta/result/bbh_mass_parabola_m2min_allgammazeta_redshift_powerlaw_samples.pkl --vt-ifar-threshold 1.0 --vt-snr-threshold 10.0 --backend cupy 

cd -
