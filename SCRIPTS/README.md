This folder contains the scripts used to produce the results in the main text and most of the Supplementary Information.

### Data formatting

* `ExploTeich.R` turns the raw data in `./IN/Initial_files/data_ROT20160324.csv` into a clean file `./IN/DBWithMonthlyPhotoTeich_completed.csv` with monthly data for only one location (Teich)
* `analyseP_oiseaux_Teich2_synchrony.R` is the first exploratory file which was used on the data. It contains mostly time-series visualization, first uses and comparison of the Gross and Loreau synchrony indices and definition of studied species
* `analyseP_oiseaux_Teich2_synchrony_reduced.R` is based on the previous file but focuses on time-series plots (focus on waders and the group Cormoran/Heron/Egret)
* `average_abundance.R` is based on the previous file but is only used to produce Fig.1 (i.e. the time-series) in MS 
* `normalisation_timeseries.R` and `summed_timeseries.R` are used to produce input files for the annual analyses. The former differentiates between cold and warm season, with separate files for the different species in one group (anas, waders, etc.) while the later produces the files per taxonomic and functional groups (i.e., all species of one group are summed for each date). These scripts rely on the function in `outputtimeseries_seasonal.R` (sum and seasonal average). 
* `occurrence.R` produces Table 1 in the Appendices, that is the precise percentage of occurrences and abundance of each species used in the analyses.

### Data analysis (functions and figures)
* `iaaft.R` is a function taken from the function `make_surrogate_data` in the package rEDM to compute the Iterative Amplitude Adjusted Fourier Transform surrogates
* `test_synchrony_Gross.r` is based on the function `community.sync` in the package synchrony but adds a way to compute surrogates and corrects the p-value computation
* `image_mvcwt.r` is taken from the same function in the the mvcwt package but uses the BH fdr correction and adds value to the colorbars. It is now only used in `wavelets_simulation.R` and other exploratory files. It has been updated to `image_mvcwt_two_panels.r` which allows for more flexibility in positioning the wavelet plot and was improved to differentiate between high and low values of the synchrony index
* `Fig2.R` draws the second figure of the MS, with the Gross index within and between taxonomic and functional groups.
* `wavelet_wader_waterfowl.R` draws the third figure of the MS with the wavelet analysis in the wader and waterfowl communities.
* `Fig5.R` draws the fifth figure of the MS with Gross and Keitt indices for the group Cormoran/Heron/Egret
* `synchrony_community.R` plots the Gross and Keitt indices for the most frequent bird group.

### Simulations (data production and analyses)
* `simulation_data.R` produces simulated data used in the Appendices which explore the way the Gross index varies for a community reacting to different environmental signals and the effect of missing values on the wavelet analysis
* `verif_simulation.R` checks the time series produced by the previous script
* `wavelets_simulation.R` shows the effect of missing values in the wavelet analysis
* `analyse_simulation.R` computes the Gross index for the different simulations
* `plot_analyse_simulation.R` plot the values from previous file.
