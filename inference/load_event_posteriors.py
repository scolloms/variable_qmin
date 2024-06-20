from gwpopulation_pipe import data_collection
import pandas as pd
import os
import h5py
import numpy as np

#these are copy pasted from the ini file, unsure how to read an ini file
sample_regex={"O3a": "/home/storm.colloms/O4a_population/test_gwpopulationpipe/o3a_samples/*.h5", "O3b": "/home/storm.colloms/O4a_population/test_gwpopulationpipe/o3b_samples/*.h5", "O4a":"/home/zoheyr.doctor/projects/O4/RatesPop/o4a-astrodist/scripts/scrape_PE/O4a_post/*.h5"}
#,"GWTC1": "/home/storm.colloms/O4a_population/test_gwpopulationpipe/gwtc1_samples/*.h5",}
preferred_labels=['PrecessingSpinIMRHM', 'C01:IMRPhenomXPHM', 'IMRPhenomXPHM']
parameters=['mass_1', 'mass_ratio', 'a_1', 'a_2', 'cos_tilt_1', 'cos_tilt_2', 'redshift']
ignore=['S230810af', 'S230830b', 'S230529ay', 'GW190425_081805', 'GW170817', 'GW200105_162426', 'GW191219_163120', 'GW200115_042309', 'S231123cg', 'S231020bw']


for label, regex in sample_regex.items():
    
    posts, meta = data_collection._load_batch_of_meta_files(
            regex=regex,
            label=label,
            labels=preferred_labels,
            keys=parameters,
            ignore=ignore,
        )
    posteriors, metadata = data_collection._load_batch_of_meta_files(regex, label, preferred_labels, parameters, ignore)
    for event_path, event in posteriors.items():
        print(event)
        print(event.keys())
        m1_q = {'mass_1':event['mass_1'],'mass_ratio':event['mass_ratio']}
        filename = os.path.basename(event_path)
        with h5py.File(f'm1_q_event_samples/{filename}', "w") as f:
            for key, value in m1_q.items():
                f.create_dataset(key, data=np.array(value))
    
    
