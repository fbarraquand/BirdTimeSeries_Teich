This folder contains all input files used in the analyses. 

The initial, raw databases are in `initial_files`. The standardized database is `DBWithMonthlyPhotoTeich_completed.csv` (produced with `SCRIPTS/ExploTeich.R`). The file `Information_Trait_Oiseaux_20190903_allimportantbirds_meanmass.csv` contains the mass of each birds, later used to compute total biomasses. 

All other files (seasonal averages and summed timeseries when grouping by function or taxonomy) are produced with `SCRIPTS/normalisation_timeseries.R` and `SCRIPTS/summed_timeseries.R`, using the previously mentioned files. File names follow specific rules:
* `warm/coldseason` corresponds to the time period during which we compute the average abundance/biomass: May to August for warm season, November to February for cold season
* `abundances/biomasses` simply correspond to the quantity we are looking at: number of birds, or total biomass
* `detailed` corresponds to abundances or biomasses of each species of a given group (calidris, anas, waders, ducks, frequent birds) while `summed` corresponds to group abundances or biomasses (sums of all species within a group)
* `wtoutrarespecies/wrarespecies` corresponds to files without or with species that are too rare to be analysed with certain methods (too many missing points). 
