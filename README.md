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


