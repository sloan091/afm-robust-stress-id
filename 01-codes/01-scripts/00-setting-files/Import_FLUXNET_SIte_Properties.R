################################################################################ 
# Name: Import_FLUXNET_Site_Properties
# By: Brandon Sloan
# Last Updated: 10/17/22
# Description: This script loads in FLUXNET site characteristics from MATLAB and
# creates factors to be used in plotting anlaysis.
################################################################################
gc(reset = TRUE)  
rm(list = ls())
library("readxl");library("dplyr"); library("forcats")
Site = read_excel("./Data/R_FLUXNET2015_SiteProp.xlsx")
Site[c(1:3,13:15)] = lapply(Site[c(1:3,13:15)],factor)
summary(Site)

# Rename SiteID to match existing R data
Site = rename(Site,SiteID = Site_ID)

# Recode climate keys
Site$ClimID = fct_recode(
  Site$ClimID,
  "Tropical" = "1",
  "Arid" = "2",
  "Temperate" = "3",
  "Cold" = "4",
  "Polar" = "5"
)

# Format the Site data into categorical variables
Site$DI.f = cut(Site$DI,breaks = c(0,0.5,1,1.5,Inf),
                 labels = c("0<DI<0.5","0.5<DI<1","1<DI<1.5","DI>1.5"))
Site$DI.f2 = cut(Site$DI,breaks = c(0,0.8,1.2,Inf),
                labels = c("0<DI<0.8","0.8<DI<1.2","DI>1.2"))
Site$DI.f3 = cut(Site$DI,breaks = c(0,1,2,Inf),
                 labels = c("0<DI<1","1<DI<2","DI>2"))
Site$DI.flx = cut(Site$DIFlx,breaks = c(0,1,2,Inf),
                 labels = c("0<DI<1","1<DI<2","DI>2"))
Site$DI.hyb = cut(Site$DIhyb,breaks = c(0,1,2,Inf),
                  labels = c("0<DI<1","1<DI<2","DI>2"))
Site$EI.f = cut(Site$EI,breaks = c(0,0.25,0.5,0.75,1,Inf),
              labels = c("0<EI<0.25","0.25<EI<0.5","0.5<EI<0.75","0.75<EI<1","EI>1"))
Site$ET_ASI.f = cut(Site$ET_ASI,breaks = c(-Inf,0.001,0.1,0.3,Inf),
                  labels = c("ET_ASI=0","0<ET_ASI<0.1","0.1<ET_ASI<0.3","ET_ASI>0.3"))
Site$PET_ASI.f = cut(Site$PET_ASI,breaks = c(-Inf,0.001,0.1,0.3,Inf),
                   labels = c("PET_ASI=0","0<PET_ASI<0.1","0.1<PET_ASI<0.3","PET_ASI>0.3"))
Site$GPP_ASI.f = cut(Site$GPP_ASI,breaks = c(-Inf,0.001,0.1,0.3,Inf),
                   labels = c("GPP_ASI=0","0<GPP_ASI<0.1","0.1<GPP_ASI<0.3","GPP_ASI>0.3"))
Site$nyrs.f = cut(Site$nyrs,breaks = c(0,2,5,10,Inf),
                labels = c("<2yr","2-5yr","5-10yr",">10yr"))
Site$h_v.f = cut(Site$h_v,breaks = c(0,1,5,15,Inf),
                  labels = c("<1m","1-5m","5-15m",">15m"))

# Calculate maximum possible correlation coefficient
Site <- Site %>% mutate(r.max = 1 / (sqrt(1 + (1 - expVar))),
                        LCE.max = 1 - sqrt((expVar * r.max - 1) ^ 2 +
                        (r.max / expVar - 1) ^ 2 + (EBslope - 1) ^ 2))

# Save file
save(Site,file = "./Data/R_FLUXNET2015_SiteProp.RData")
