This folder contains all input files used in the analyses. 

The initial, raw databases are `initial_files`. The standardized database is `DBWithMonthlyPhotoTeich_completed.csv` and the file `Information_Trait_Oiseaux_20190903_allimportantbirds_meanmass.csv` show the mass of each birds, later used to compute total biomasses. 

All other files (seasonal averages and summed timeseries when grouping per genus) are produced with `SCRIPTS/normalisation_timeseries.R` and `SCRIPTS/summed_timeseries.R`, using the previously mentioned files.
