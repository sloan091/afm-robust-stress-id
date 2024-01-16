################################################################################ 
# Name: 01-extract-dominant-signals.R
# By: Brandon Sloan
# Last Updated: 1/16/24
#
# Description: This script loads the individual performance and stress signals  
# for the 2304 assumptions sets at the 151 FLUXNET site, calculates key 
# statistics (e.g., median, IQR), and extracts the dominant stress signal 
# (slope) and corresponding performance metric. These results are used to create
# the figures in Sloan and Feng (2023).
#
# Outputs:
#   Slp.Stat - Stress signal stats for 150 FLUXNET sites
#   Prf.Stat - Corresponding performance stats for Slp.Stat
#   Slp.Cat  - Stress categorization based on locations and signs of stress
#   Slp.Sel  - Dominant stress signals and corresponding performances extracted
#              from Slp.Stat and Prf.Stat
################################################################################


# Clear workspace, add functions, and load libraries
gc(reset = TRUE)  
rm(list = ls())
source("./01-codes/02-functions/04-robustness/robustness-helper-fxns.R")
loadPacks(c("lemon","broom","ggh4x","tictoc","entropy"))
load("./01-codes/01-scripts/00-setting-files/final_ec_site_properties_for_r.RData")
select <- dplyr::select


# Thresholds for reliability
#===========================

# Practically significant slope threshold [% pv drop/SWC_nep]
th.Slp <- 0.2 / 0.5

# Performance thresholds
Prf.vars <- c("LCE", "NSE", "AfT1", "AfT2", "AfT3")
th.Prf <-
  tibble(mtcID = as.factor(Prf.vars), 
         th = c(0.4, 0, 0.7, 0.7, 0.7))


# Load file and setpaths
datpath = "./03-outputs/02-robust-framework/01-stress-perf-dist/"
fn = list.files(datpath,"*.RData")
SiteID = unique(substr(fn,start = 1,stop = 6))

# Pre-allocate breakpoints
Slp.Stat = data.frame()
Prf.Stat = data.frame()

#===============================================================================
# 1. Calculate and combine slope and performance stats for each site
#===============================================================================
tic()
# Loop through results for different sites
for (ii in 1:length(fn)){
  
  # Read in files
  load(paste0(datpath,fn[ii]))
  
  # y2d = discretize(Slp$beta, numBins=10)
  # entropy(y2d, method = "shrink")
  # Calculate slope statistics
  Slp.ii <- Slp  %>%
    group_by(WetID, pvID) %>%
    summarize(
      Slp_mn = median(beta,na.rm = TRUE),
      Slp_nrm = median(beta.nrm,na.rm = TRUE),
      Slp_sd = IQR(beta,na.rm = TRUE),
      Slp_kt = (quantile(beta,0.9,na.rm = TRUE) - quantile(beta,0.1,na.rm = TRUE))/Slp_sd,
      Slp_ent = discretize(na.omit(beta),numBins = 20) %>% entropy(),
      Slp_cv = min(Slp_sd / abs(Slp_mn), 10),
      MaxMed = max(MaxMed),
      SigID = abs(Slp_mn) >= th.Slp * abs(MaxMed),
      RobSlp = case_when(SigID & Slp_mn > 0 ~ 1 - ecdf(beta)(th.Slp * abs(MaxMed)),
                         SigID & Slp_mn < 0 ~ ecdf(beta)(-th.Slp * abs(MaxMed)),
                         TRUE ~ ecdf(beta)(th.Slp * abs(MaxMed)) - ecdf(beta)(-th.Slp * abs(MaxMed)))
    ) %>%
    ungroup() 
  
  
  # Update information metrics Afp and AfT to have 1 be maximum good value
  Prf <-
    Prf %>% filter(mtcID %in% Prf.vars) %>% 
    mutate(beta = case_when(mtcID %in% c("AfT1", "AfT2", "AfT3") 
                                         ~ 1 - abs(beta), 
                                         TRUE ~ beta)) %>% 
    left_join(th.Prf,by = "mtcID")
  
  Prf.ii <- Prf %>% select(-se) %>%
    pivot_wider(
      names_from = mtcID,
      values_from = c(beta, th),
      names_glue = "{mtcID}_{.value}"
    ) %>% 
    mutate(Aft.mn = (AfT1_beta + AfT2_beta + AfT3_beta) / 3) 
  
  # Calculate joint probabilities
  Jnt.Prf <- Prf.ii %>%
    group_by(WetID) %>%
    summarize(RobPrf.jnt = sum(Aft.mn > AfT1_th &
                                 LCE_beta > LCE_th, na.rm = TRUE) / sum(!is.na(LCE_beta))) %>%
    ungroup()
  
  Jnt.Slp <-
    Slp %>% left_join(Slp.ii %>% select(WetID, pvID, Slp_mn, SigID),
                      by = c("WetID", "pvID")) %>%
    select(pvID, WetID, ExpID, beta.nrm, Slp_mn, SigID) %>%
    left_join(Prf.ii, by = c("WetID", "ExpID")) %>%
    group_by(pvID, WetID) %>%
    summarize(
      SigID = unique(SigID),
      Slp_mn = unique(Slp_mn),
      RobSlpPrf.jnt = case_when(
        SigID & Slp_mn > 0 ~ sum(Aft.mn > AfT1_th &
                                   LCE_beta > LCE_th &
                                   beta.nrm >= th.Slp,
                                 na.rm = TRUE) / sum(!is.na(LCE_beta)),
        SigID &
          Slp_mn < 0 ~ sum(Aft.mn > AfT1_th &
                             LCE_beta > LCE_th &
                             beta.nrm <= -th.Slp,
                           na.rm = TRUE) / sum(!is.na(LCE_beta)),
        TRUE ~ sum(
          Aft.mn > AfT1_th &
            LCE_beta > LCE_th &
            beta.nrm >= -th.Slp &
            beta.nrm <= th.Slp,
          na.rm = TRUE
        ) / sum(!is.na(LCE_beta))
      )
    ) %>% ungroup()
  

  # Average performance metrics
  Prf.ii <- Prf %>%
    group_by(WetID, mtcID) %>%
    summarize(
      mn = median(beta,na.rm = TRUE),
      sd = IQR(beta,na.rm = TRUE),
      kt = (quantile(beta,0.9,na.rm = TRUE) - quantile(beta,0.1,na.rm = TRUE))/sd,
      ent = discretize(na.omit(beta),numBins = 20) %>% entropy(),
      cv = min(sd / abs(mn), 10),
      Rob = 1 - ecdf(beta)(max(th))
    ) %>%
    ungroup() %>% 
    pivot_wider(
      names_from = mtcID,
      values_from = c(mn, sd, kt, ent, cv, Rob),
      names_glue = "{mtcID}_{.value}"
    ) %>% 
    left_join(Jnt.Prf,by = "WetID")
  
  # Add site ID as a factor and concatenate
  Prf.ii$SiteID = as.factor(SiteID[ii])
  Slp.ii$SiteID = as.factor(SiteID[ii])
  Slp.ii <- Slp.ii %>%
    left_join(Jnt.Slp %>% select(-SigID, -Slp_mn), by = c("pvID", "WetID"))
  Prf.Stat = rbind(Prf.Stat, Prf.ii)
  Slp.Stat = rbind(Slp.Stat, Slp.ii)
}
toc()


#===============================================================================
# 2. Categorize average stress slopes
#===============================================================================

# For piecewise slope fits (ie, two slopes)
tmpDF <- left_join(Slp.Stat, Prf.Stat, by = c("SiteID", "WetID"))
SC.pw <- tmpDF %>% 
  filter(WetID != "All") %>%
  group_by(SiteID, pvID) %>%
  summarize(
    Sigct = sum(SigID),
    dSigID = SigID[WetID == "Dry"],
    wSigID = SigID[WetID == "Wet"],
    dsgn = sign(Slp_mn[WetID == "Dry"]),
    wsgn = sign(Slp_mn[WetID == "Wet"]),
    MaxMed = max(MaxMed)
  ) %>%
  mutate(dsgn = ifelse(dSigID, dsgn, 0),
         wsgn = ifelse(wSigID, wsgn, 0)) %>%
  mutate(SlpCatID = as.character((
    case_when(
      Sigct == 2 & wsgn > 0 & dsgn > 0 ~ "1.c",
      Sigct == 2 & wsgn < 0 & dsgn < 0 ~ "2.c",
      Sigct == 2 & wsgn < 0 & dsgn > 0 ~ "3.a",
      Sigct == 2 & wsgn > 0 & dsgn < 0 ~ "3.b",
      Sigct == 1 & wsgn == 0 & dsgn < 0 ~ "2.a",
      Sigct == 1 & wsgn == 0 & dsgn > 0 ~ "1.a",
      Sigct == 1 & wsgn < 0 & dsgn == 0 ~ "2.b",
      Sigct == 1 & wsgn > 0 & dsgn == 0 ~ "1.b",
      TRUE ~ "4.a"
    )
  ))) %>% 
  ungroup()

# For single slope fits
SC.sng <- Slp.Stat %>% 
  filter(WetID == "All") %>%
  group_by(SiteID, pvID) %>%
  summarize(
    Sigct = sum(SigID),
    sgn = sign(Slp_mn),
    SigID = SigID,
    MaxMed = max(MaxMed)
  ) %>%
  mutate(sgn = ifelse(SigID, sgn, 0)) %>%
  mutate(SlpCatID = as.character((
    case_when(Sigct == 1 & sgn > 0 ~ "1.d",
              Sigct == 1 & sgn < 0 ~ "2.d",
              TRUE ~ "4.d")
  ))) %>% 
  ungroup()

# Combine and clean up stress categorization 
SC <- rbind(SC.pw[c("SiteID", "pvID", "MaxMed", "SlpCatID")],
            SC.sng[c("SiteID", "pvID", "MaxMed", "SlpCatID")])
SC <- separate(SC, SlpCatID, c("SlpID", "LocID"), remove = FALSE)
SC[c("SlpCatID", "SlpID", "LocID")] <-
  lapply(SC[c("SlpCatID", "SlpID", "LocID")], factor)
SC$pvID <-
  factor(SC$pvID,
         levels = c(
           "G_1",
           "G_1_VPD_m",
           "G_v"
         ))

# Add on site characteristics, relevant slope thresholds, and create final df
Slp.Cat <- left_join(SC, Site, by = c("SiteID")) %>%
  group_by(SiteID) %>%
  mutate(conID = factor(length(unique(SlpID)))) %>%
  ungroup()

#===============================================================================
# 3. Extract practically significant stress slopes and performance metrics
#===============================================================================


# Extraction function based on stress category
extractDom <- function(x1, x2, dom, loc, WetID) {
  dom = unique(dom)
  loc = unique(loc)
  val = if (dom %in% c("1", "2") & loc == "a") {
    x1[WetID == "Dry"]
  }
  else if (dom %in% c("1", "2") & loc == "b") {
    x1[WetID == "Wet"]
  }
  else if (dom %in% c("1", "2") &
           loc == "c") {
    x1[which.max(abs(x2))]
  }
  else if (dom %in% c("1", "2") & loc == "d") {
    x1
  }
  else if (dom %in% c("4")) {
    x1[which.max(abs(x2))]
  }
  else{
    NA
  }
  return(val)
}

# Temporary df for processing
tmpDF <- Slp.Stat %>%
  left_join(
    Prf.Stat,
    by = c("SiteID", "WetID")
  ) %>%
  left_join(
    select(Slp.Cat, contains("ID")),
    by = c("SiteID", "pvID"),
    suffix = c("", ".rep")
  ) %>%
  select(!ends_with(".rep")) 


# Select significant slopes for situations with either a single dominant slope 
# signs or no dominants slopes (SlpID = 1,2,4)
Slp.1ds <- tmpDF %>% group_by(SiteID, pvID) %>%
  summarize(across(contains(c(
    "_mn", "_sd", "_cv", "_kt", "_ent", "Rob"
  )),  ~ extractDom(.x, Slp_mn, SlpID, LocID, WetID))) %>%
  ungroup() %>% 
  na.omit()

# Select significant slopes for situations with two dominants slopes (SlpID = 3)
Slp.2ds <- tmpDF %>%
  filter(SlpID == "3") %>% na.omit() %>%
  select(SiteID,
         pvID,contains(c(
           "_mn", "_sd", "_cv", "_kt", "_ent", "Rob"
         )))

# Combine, add site characteristics and calculate robustness metric
Slp.Sel <- rbind(Slp.1ds, Slp.2ds) %>%
  left_join(Slp.Cat,
            by = c("SiteID", "pvID"),
            suffix = c("", ".rep")) %>%
  mutate(
    RobPrf = ((AfT1_Rob + AfT2_Rob + AfT3_Rob)/3 + LCE_Rob)/2,
    Slp_nrm = Slp_mn / MaxMed,
    Slp_sd_nrm = Slp_sd / MaxMed,
  ) %>%
  select(!ends_with(".rep"))


# Save files as RData - Commented out as not to overwrite the given files
# save(
#   Prf.Stat,
#   Slp.Stat,
#   Slp.Cat,
#   Slp.Sel,
#   file = "./03-outputs/02-robust-framework/02-robust-summary/SlpPrfDist_Stats.RData")

# Plot illustrating that the stress and performance metric robustness are nearly
# independent.
ggplot(data = Slp.Sel, aes(x = RobSlp*RobPrf.jnt,y = RobSlpPrf.jnt,color = pvID)) +
  geom_point() + labs(x = "Pr(Stress > th)*Pr( Perf > th)", y = "Pr(Stress > th, Perf > th)") 
