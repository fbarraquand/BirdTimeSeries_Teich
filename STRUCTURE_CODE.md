# Code structure

We describe here how to obtain the main results shown in Figures in the article (main text and supplementary material). All functions are described in `SCRIPTS`, these are the most important for producing results. 

### Main text

* **Figure 1** is drawn by `SCRIPTS/average_abundance.R`.
* **Figure 2** is produced by `SCRIPTS/Fig2.R`. It calls the function `SCRIPTS/test_synchrony_Gross.r` to compute the observed synchrony value and the associated p-value. The IAAFT computation is made by `SCRIPTS\iaaft.R`.
* **Figure 3** is produced by `SCRIPTS/use_iaaft_community.r`. The top panel calls the function `SCRIPTS/test_synchrony_Gross.r` (same use as above). The bottom panel calls the function `SCRIPTS/image_mvcwt_for_colormaps.r` to draw p-values extracted from synchrony values computed with mvcwt package and randomization based on `SCRIPTS\iaaft.R`.
* **Figure 4** is produced by `SCRIPTS/use_iaaft_wavelet_wader.r`, based on `SCRIPTS/image_mvcwt_for_colormaps.r` (similarly as the bottom panel for Figure 3).
* **Figure 5** is drawn by `SCRIPTS/analyseP_oiseaux_Teich2_synchrony_reduced.R`.
* **Figure 6** is produced by `SCRIPTS/use_iaaft_triad.r` (same functioning as Figure 3).

### Supplementary material

* **Figure S1** and **Figure S2** are drawn by `SCRIPTS/analyseP_oiseaux_Teich2_synchrony_reduced.R`.
* **Figure S3** is produced by `SCRIPTS/explo/Gross_triad.R`, using `SCRIPTS/test_synchrony_Gross.r`.
* **Figure S4** is drawn by `SCRIPTS/analyseP_oiseaux_Teich2_synchrony_reduced.R`. 
* **Figure S5** is produced by `SCRIPTS/Fig2.R`. It calls the function `SCRIPTS/test_synchrony_Gross.r` to compute the observed synchrony value and the associated p-value. The IAAFT computation is made by `SCRIPTS\iaaft.R`.
* **Figure S6** is produced by `SCRIPTS/use_iaaft_wavelet_wader.r`, based on `SCRIPTS/image_mvcwt_for_colormaps.r` (same as Figure 4).
* **Figure S7** is drawn by `SCRIPTS/verif_simulation.R`. Time series themselves are produced by `SCRIPTS/simulation_data.R`
* **Figure S8**, **Figure S9** and **Figure S10** are drawn by `SCRIPTS/plot_analyse_simulation.R`, based on the simulations produced in `SCRIPTS/simulation_data.R` and analysed by `SCRIPTS/analyse_simulation.R` using `SCRIPTS/test_synchrony_Gross.r`.
* **Figure S11** is produced by `SCRIPTS/use_iaaft_simulation_with_missing_data.r`, using `SCRIPTS/image_mvcwt_for_colormaps.r` and `SCRIPTS/iaaft.r`
* **Figure S12** is drawn by `Submission_JAE/Revisions_R2/Simulations_response/time_series.R`. Time series themselves are produced by `Submission_JAE/Revisions_R2/Simulations_response/MultivariateNormalModels_SAD.R`. 
* **Figure S13**, **Figure S14** and **Figure S15** are produced by `Submission_JAE/Revisions_R2/Simulations_response/use_iaaft_simulations.r`, using `SCRIPTS/test_synchrony_Gross.r` and `SCRIPTS/image_mvcwt_for_colormaps.r`.
* **Figure S16** and **Figure S17** are produced by `Submission_EcolEvol/Simulations_response/MultivariateNormalModels_SAD_MultiTrends.R`, using `SCRIPTS/test_synchrony_Gross.r`and `SCRIPTS/image_mvcwt_for_colormaps.r`.
* **Figure S18** is produced by `Submission_EcolEvol/Simulations_response/MultivariateNormalModels_SAD.R`, relying on `SCRIPTS/image_mvcwt_for_colormaps.r`.
* **Figure S19** is produced by `ubmission_EcolEvol/Simulations_response/MultivariateNormalModels_2sp_antiphase_log_amplitude.R`, relying on `SCRIPTS/image_mvcwt_for_colormaps.r`. 

NOTE: Some functions may require to un-comment/re-comment calls to other functions, depending on possibly conflicting or repeated calls to functions in scripts. 
