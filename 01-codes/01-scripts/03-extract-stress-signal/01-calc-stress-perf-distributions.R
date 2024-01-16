################################################################################ 
# Name: 01-calc-stress-perf-distributions
# By: Brandon Sloan
# Last Updated: 1/16/24
#
# Description: Calculates the stress signal and mean performance for the 
# 2,304 total assumption sets (AS) for the 150 FLUXNET sites in 
# Sloan and Feng (2023). This code is set up to be parallelizable with furrr. 
#
# Outputs:
#   Slp - Stress signal for each AS derived from segmented or linear regression
#         for the three ecosystem parameters
#   Prf - Mean performance metrics based on stress signal segmentation.
################################################################################

# Paths and files
#================
gc(reset = TRUE)  
rm(list = ls())
source("./01-codes/02-functions/04-robustness/robustness-helper-fxns.R")
loadPacks(c("segmented","dplyr","tictoc","furrr"))
select <- dplyr::select
load("./01-codes/01-scripts/00-setting-files/final_ec_site_properties_for_r.RData")
datpath = "./03-outputs/01-estimated-pmoc-parameters/02-excel/"
svpath = "./03-outputs/01-stress-perf-dist/"

# Settings for parallelization
#=================================
nc <- unname(availableCores())
plan("multicore",workers = nc - 1)
#plan("sequential")

# Analysis settings
#==================
# Contrast coding for regression
options(contrasts = rep("contr.sum", 2))

# Independent and dependent variable names
x <-  c("SWC_nep","Fct1","Fct2","Fct3","Fct4","Fct5","Fct6","Fct7","Fct8")
y.slp <-  c("G_1","G_1_VPD_m","G_v")
y.prf <-  c("NSE","LCE","AfT1","AfT2","AfT3","Ap")

# Regression equations for breakpoint, slope, and performance analysis
bp.frm <-  as.formula('pvVal ~ SWC_nep + Fct1 + Fct2 + Fct3 + Fct4 + Fct5 + Fct6 + Fct7 + Fct8')
slp.frm <-  as.formula('pvVal ~ SWC_nep')
prf.frm <-  as.formula('mtcVal ~ 1')

# Percent performance increase in R^2 required to have a breakpoint
bp.prfThresh <- 0.05


# Load and process flux tower data
#=================================
# Create site list
fn = list.files(datpath,pattern = "*csv")

# Pre-allocate
bp.All = data.frame()
badsites = c()

# Loop through sites
tic()
for(ii in 1:length(fn)){

# Load data
fnii = fn[ii]
SiteID = substr(fnii, start = 1, stop = 6)
D <- na.omit(getFitResults(datpath, SiteID)) %>% 
  mutate(ExpID = as.factor(ExpID))

# Separate and lengthen a slope and performance dataset
D.slp <-  D %>% select(Fct1:ExpID, all_of(y.slp)) %>%
  pivot_longer(cols = all_of(y.slp),
               names_to = "pvID",
               values_to = "pvVal")
D.prf <-  D %>% select(Fct1:ExpID, SWC_nep, all_of(y.prf)) %>%
  pivot_longer(cols = all_of(y.prf),
               names_to = "mtcID",
               values_to = "mtcVal")

tryCatch({

  # Define Wet-Dry breakpoint
  #==========================
  SWC.all <- sort(unique(D$SWC_nep))
  nbins <- length(SWC.all)
  
  # If adequate bins, test which breakpoint gives best performance increase
  if (nbins > 5) {
    # Ensure at least three bins for wet and dry
    SWC.bins <- SWC.all[seq(3, nbins - 3, 1)]
    
    # Fit breakpoint model starting with Davies test
    fit.bp <- D.slp %>%
      group_by(pvID) %>%
      nest() %>% ungroup() %>%
      mutate(fit  = map(data,
                        ~ lm(bp.frm,
                             data = .x)),
             bp = map2(data, fit,  ~ bpfxn(.y, .x, SWC.bins)))
    
    # Extract breakpoints and see if any provide significant benefit
    bps <-
      fit.bp %>% unnest(bp) %>% select(pvID, bp, relBIC, absR2, relDev, relR2) %>% ungroup()
    bp.prf <- as.matrix(bps[, c("absR2")])
    bp.idx <-
      bp.prf > bp.prfThresh
    bp.wts <- t(bp.prf / sum(bp.prf))
    if (sum(bp.idx > 0)) {
      bp <- as.numeric(bp.wts %*% bps$bp)
    } else{
      bp <- NaN
    }
  } else{
    # If insufficient points, do not have a breakpoint
    bp <- NaN
  }

# Assign Wet, Dry, or All flag to data points
if(is.nan(bp)) {
  D.slp$WetID = factor("All")
  D.prf$WetID = factor("All")
} else{
  D.slp$WetID = with(D.slp,factor(ifelse(SWC_nep<=bp, "Dry", "Wet")))
  D.prf$WetID = with(D.prf,factor(ifelse(SWC_nep<=bp, "Dry", "Wet")))
}


# Fit and extract plant variable response to soil moisture
#=========================================================

# Helper functions for extracting slope or intercept
exSlp <- function(data,slp.frm,trm){
  slp <- lm(slp.frm, data = data) %>%
    tidy() %>%
    filter(term == trm)
}
posexSlp <- possibly(exSlp,otherwise = NULL)
  
# Extract plant variable beta and SE for each run
tic()
fit.slp <- D.slp %>%
  group_by(pvID,WetID,ExpID) %>%
  nest() %>% 
  ungroup() %>% 
  mutate(Slp  = future_map(
    data,
    ~posexSlp(.x,
      slp.frm,
      "SWC_nep"
    )
  )
  )
toc()

Slp <- fit.slp %>% 
  unnest(Slp,keep_empty = TRUE) %>% 
  mutate(beta = estimate,se = std.error) %>% 
  select(pvID,WetID,ExpID,beta,se) %>% 
  ungroup()

# Calculate maximum median and sd values of plant variables
MaxMed <-
  D.slp %>% 
  group_by(pvID, SWC_nep) %>% 
  summarize(MaxMed = median(pvVal, na.rm = TRUE)) %>% 
  ungroup(SWC_nep) %>% 
  summarize(MaxMed = max(MaxMed)) %>% 
  ungroup()

RefSD <-
  D.slp %>% 
  group_by(pvID,WetID) %>% 
  summarize(RefSD = sd(pvVal, na.rm = TRUE)) %>% 
  ungroup() 

# Final slope table
Slp <- Slp %>% 
  left_join(MaxMed,by = "pvID") %>% 
  left_join(RefSD,by = c("pvID","WetID")) %>% 
  mutate(beta.nrm = beta/MaxMed)

# Fit and extract average functional and predictive performance
#==============================================================

# Extract performance beta and SE for each run
tic()
fit.prf <- D.prf %>%
  group_by(mtcID, WetID, ExpID) %>%
  nest() %>%
  ungroup() %>%
  mutate(Prf  = future_map(data,
                           ~ posexSlp(.x,
                                      prf.frm,
                                      "(Intercept)")))
toc()

Prf <- fit.prf %>% 
  unnest(Prf,keep_empty = TRUE) %>% 
  mutate(beta = estimate,se = std.error) %>% 
  select(mtcID,WetID,ExpID,beta,se) %>% 
  ungroup()

# Save as R data - Commented out as not to overwrite the given files
# save(Slp,Prf,file = paste0(svpath,SiteID,"_SlpPrfDist.RData"))

}
,error=function(e){badsites = rbind(badsites,SiteID)})

}
toc()
# Show which sites failed to run
print(badsites)
