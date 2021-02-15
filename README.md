# Fingerprinting and Behavioural Overlap

This is the accompanying repository for the preprint "Fingerprinting and behavioural prediction rest on distinct functional systems of the human connectome". 

[![DOI](https://zenodo.org/badge/183295261.svg)](https://www.biorxiv.org/content/10.1101/2021.02.07.429922v1)



We provide example functional connectivity data in /fc, example behavioural data in /behavioural_data and example analysis input/output in /results, which should allow the reader to follow the behaviour of the scripts without real data and thus allow to adapt to their own data.

Both imaging and phenotypic data was downloaded from the HCP, where you will need to register.
https://db.humanconnectome.org/app/template/Login.vm

We extracted the phenotypic data from a file named "Unrelated_S900_Subject_multilist1_with_physio.csv" found here:
https://wiki.humanconnectome.org/display/PublicData/S900+Unrelated+Subjects+CSV

Please run the script in the following order:

1. prediction_behaviour.m
2. DP_edges.m
3. edgewise_overlap.m

All results are saved in /results. 

Any further analysis or visualization can then be based on the output in /results/overlap_results.

With any questions, feel free to contact the collaborators on this repository.

Mantwill&Gell, 15.02.2021
