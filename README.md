# Codes for "Robust Inference of Ecosystem Soil Water Stress from Eddy Covariance Data"
This repository contains the relevant codes for the following article:

Sloan, B. P., & Feng, X. (2023). Robust inference of ecosystem soil water stress from eddy covariance data. Agricultural and Forest Meteorology, 343, 109744. https://doi.org/10.1016/j.agrformet.2023.109744

The analysis is split into two main parts: 1) data pre-processing and fitting the Penman-Monteith Optimal Conductance Model under numerous assumption sets in MATLAB, and 2) extracting PMOC stress signals and performance to quantify each eddy covariances site's robustness in R.
I have tried to include all necessary codes to re-create my analysis and figures; however, I have only included to two example eddy covariance (and LAI) data sets out of the 151 sites used in Sloan and Feng (2023) given all data is freely available and the overall data sets are large. The data used in this research are freely available from the FLUXNET2015 (https://fluxnet.org/data/fluxnet2015-dataset/), AmeriFlux (https://ameriflux.lbl.gov/data/download-data/), and MODIS (https://doi.org/10.5067/MODIS/MC D15A2H.061) databases.

# Fitting the Penman-Monteith Optimal Conductance (PMOC) Model in MATLAB

## Data Pre-processsing
The codes for pre-processing of eddy covariance and LAI data are located under _01-codes/01-scripts/01-clean-data/_. The reader is referred to Sect. 2.1.1 and S1 of Sloan and Feng (2023) for full pre-processing details.

### FLUXNET2015 and AmeriFlux-FLUXNET Data
We selected a total of 151 out of 229 available eddy covariances sites from the FLUXNET2015 (https://fluxnet.org/data/fluxnet2015-dataset/) and AmeriFlux-FLUXNET (https://ameriflux.lbl.gov/data/download-data/) data sets based on adequate data coverage of variables relevant to detecting ecosystem soil water stress with the PMOC model (see Sect. S1 of Sloan and Feng, 2023). The codes included in _01-codes/01-scripts/01-clean-data/01-eddy-covariance/_ perform the following steps:

* _step_01_unpack_fluxnet2015_data_: Extracts the raw half-hourly (or hourly) FLUXNET2015 spreadsheets into a MATLAB Table. I have included a sample raw file for US-Me1 in _02-data\01-raw\01-eddy-covariance\01-fluxnet2015\FLX_US-Me1_FLUXNET2015_FULLSET_2004-2005_1-4/_. This script only runs for the sample US-Me1 data set and saves the table output as _02-data\02-processed\01-eddy-covariance\US_Me1.mat_. 
* _step_02_unpack_ameriflux_fluxnet_data_: Extracts the raw half-hourly (or hourly) AmeriFlux-FLUXNET spreadsheets into a MATLAB Table. I have included a sample raw file for US-Me2 in _02-data\01-raw\01-eddy-covariance\01-fluxnet2015\AMF_US-Me2_FLUXNET_FULLSET_2002-2020_3-5/_. This script only runs for the sample US-Me2 data set and saves the table output as _02-data\02-processed\01-eddy-covariance\US_Me2.mat_. This script is identical to the _step_01_unpack_fluxnet2015_data_, but is included as the data sets were originally stored in different locations.
* _step_03_seb_correction.m_: Re-performs the surface energy budget correction using Method 1 from the OnePipe processing pipeline (Pastorello et al., 2020). This step is necessary only for sites missing the _LE_CORR_ variable, which seems to occur when ground heat flux is not measured.
* _step_04_check_data_coverage.m_: Checks data coverage of eddy covariance variables relevent to detecting ecosystem soil water stress with the PMOC model (e.g., soil moisture, net radiation). This script outputs the table _01-codes\01-scripts\00_setting_files\ec_site_coverage.mat_ that summarizes the data coverage and is used to filter the final 151 eddy covariance sites.

At the time of analysis, there were 229 potential eddy covariance sites (_01-codes\01-scripts\00_setting_files\potential-ec-sites.mat_). Site selection was an iterative process of checking data coverage (_ec_site_coverage.mat_) as well as other relevant metadata (vegetation height, crop rotations, etc). The final eddy covariance sites are provided in _01-codes\01-scripts\00_setting_files\final-ec-site-properties.mat_ and shown in Table S1 of Sloan and Feng (2023).

### MODIS Leaf Area Index (LAI) Data
The MODIS LAI data was one of the selected treatment levels to control for the influence of vegetation dynamics on soil water stress inference (see Fig. 2 of paper). The LAI data was downloaded for each site with the ORNL DAAC's Terrestrial Ecology Subsetting & Visualization Services (TESViS, https://modis.ornl.gov/). The codes included in _01-codes/01-scripts/01-clean-data/02-modis-lai/_ perform the following steps:

* _step_01_download_modis_lai.m_: Automatically downloads the raw MODIS LAI and fPAR (MOD15A3H) as well as EVI and NDVI (MOD13Q1) data products. I have included sample outputs from this code under _02-data\01-raw\01-modis-lai\_ for both US-Me1 and US-Me1.
* _step_02_filter_modis_lai.m_: Smooths the noisy LAI data using a short- and long-term filter similar to the methods in Ukkola et al. (2022). I have included sample outputs and figures for US-Me1 and US-Me2 in _02-data\02-processed\01-modis-lai\_. The _site_id_LAI_FLX15.mat_ files are used to normalize the relevant fluxes in the PMOC model. See Sect. S6 of the paper for details on LAI processing.

## Fit PMOC Model by Soil Moisture Percentiles
The PMOC parameters are estimated for all 151 sites under the 2,304 unique assumption sets for at most 10 soil moisture bins. The full details on the PMOC model, assumption sets, and parameter estimation are located in Sect. 2 and S5 of Sloan and Feng (2023). Here, I will briefly cover the relevant codes for running the parameter estimation and exporting results for use in R located in  _01-codes/01-scripts/02-fit-pmoc-model/_. Note, the codes can take awhile to execute even with parallelization given that at each site nearly 2,304 fits need to be run for each soil moisture bin (at most 10). As a rough estimate, US-Me2 required ~10 minutes to fit all parameter sets when using 12 cores in parallel on my HP Zbook Firefly with 32 GB or RAM. 

### Estimate PMOC Parameters
The main parameter estimation script is _01-codes/01-scripts/02-fit-pmoc-model/step_01_fit_pmoc_parameters.m_. This script is a wrapper that sets up parallelization (if wanted), settings, and loops through the parameter and model performance estimation at each eddy covariance site. The wrapper script leverages several other nested functions that actually perform the fit for each of the 2,304 assumption sets in each soil moisture bin. These functions are primarily located in _01-codes/02-functions/02-pmoc-model/_:

* _fit_pmoc_factorial_: The function called by the wrapper script that pre-allocates all matrices, loops through the fit algorithm for each of the 2,304 assumption sets, and saves the resulting parameter and performance estimate for a site. We will discuss the outputs in more detail below. Note, the user can set a parfor loop on Line 68 of this file (i.e., parallelize the loop for assumption sets) to speed up simulation.
* _fit_pmoc_sm_bin_: The function called by _fit_pmoc_factorial_ to loop through the fit algorithm for each soil moisture bin for a fixed assumption set and eddy covariance site. There are multiple options for defining soil moisture bins, but we used soil moisture percentiles in the paper (see Sect. 4.4 for discussion).
* _fit_pmoc_train_test_: The function called by _fit_pmoc_sm_bin_ to perform the training and testing of the parameter estimation for a fixed site, assumption set, and soil moisture bin. Here we use sampling with replacement to extract a training set and test the fit on out-of-bag samples. This essentially results in a 2/3 train-test split. The random sampling means that results will not equal exactly those in the paper as the random seed is different. However, the results are the same with respect to statistics.
* _fit_pmoc_: The function called by _fit_pmoc_train_test_ to estimate the parameters using various MATLAB nonlinear solvers for a fixed site, assumption set, moisture bin, and training data set. The function specifies the PMOC model based on the assumption set and formats the inputs for the selected MATLAB solver. The function leverages multiple other functions in the _01-penman-monteith_ and _02-optimal-conductance_ folders. We will not cover these supporting functions, but details are given in the header of each script. 
* _train_fit_: A untility function used by _fit_pmoc_train_test_ to train and test _fit_pmoc_. There are options to perform bootstrapping of fits as well, but they are not used in Sloan and Feng (2023) and may not function properly given multiple iterations of this code.






