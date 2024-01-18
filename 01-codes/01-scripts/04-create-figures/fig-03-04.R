################################################################################ 
# Name: fig-03-04
# By: Brandon Sloan
# Last Updated: 1/17/24
#
# Description: Creates the slope (Fig. 3) and performance distributions (Fig.4) 
# for the case study sites in Sect. 3.3 and Sloan and Feng (2023).
################################################################################

# Load libraries and add datasets
gc(reset = TRUE)  
rm(list = ls())
source("./01-codes/02-functions/04-robustness/robustness-helper-fxns.R")
load("./01-codes/01-scripts/00-setting-files/final_ec_site_properties_for_r.RData")
load("./03-outputs/02-robust-framework/02-robust-summary/SlpPrfDist_Stats.RData")
svpath <- "./03-outputs/03-figures/"
loadPacks(c("ggridges","egg"))
select <- dplyr::select
prf.vars = c("LCE","AfT1","AfT2","AfT3")
svpath <- "./03-outputs/03-figures/"

# Plot settings
# Font size
fs <- 8
lw <- 0.3
lt30 <- "dotted"
mk <- 21
mkst <- 0.25
scl <- 1

# Select stress periods and load data
#====================================
# Select Sites
SiteSel <- c("US_Me1","US_Me2","US_Me3","US_Me5","US_Me6")
Slp.Sel2 <- Slp.Sel %>% filter(SiteID %in% SiteSel,!(SlpID == 3 &
                                                       Slp_nrm < 0)) %>%
  select(SiteID, pvID,Slp_nrm, Slp_kt,Slp_ent,Slp_sd, SlpID, LocID,
         all_of(paste0(prf.vars,'_mn')),
         all_of(paste0(prf.vars,'_kt')),
         all_of(paste0(prf.vars,'_ent'))) %>%
  mutate(WetIDSel = case_when(
    LocID == "a" ~ "Dry",
    LocID == "b" ~ "Wet",
    LocID == "d" ~ "All",
    TRUE ~ "Dry"
  ))

# Use 2 out of 3 rule to determine which category each belongs to
Prf.Sel <- Slp.Sel2 %>% group_by(SiteID) %>% 
  summarize(Unq = length(unique(WetIDSel))) %>% 
  left_join(Slp.Sel2,by = "SiteID") %>% ungroup() %>% 
  group_by(SiteID,Unq,WetIDSel) %>% 
  summarize(n = n()) %>% 
  ungroup() %>%
  left_join(Slp.Sel2,by = c("SiteID","WetIDSel")) %>% 
  group_by(SiteID, Unq) %>%
  summarize(WetIDSel = WetIDSel[which.max(n)],
            across(contains(prf.vars),~mean(.[which.max(n)],na.rm = TRUE))) %>%
  ungroup() 

# Helper functions to load slope files
Slp.ldfxn <- function(ID,flt){
  load(paste0("./03-outputs/02-robust-framework/01-stress-perf-dist/",ID,"_SlpPrfDist.RData"))
  flt <- flt %>% filter(SiteID == ID)
  Slp.flt <- Slp %>% left_join(flt,by = c("pvID")) %>% 
    filter(WetID == WetIDSel) %>% mutate(SiteID = as.factor(ID))
  return(Slp.flt)
}

# Helper functions to load performance files
Prf.ldfxn <- function(ID,flt){
  load(paste0("./03-outputs/02-robust-framework/01-stress-perf-dist/",ID,"_SlpPrfDist.RData"))
  flt <- flt %>% filter(SiteID == ID)
  Prf.flt <- Prf %>% filter(WetID == flt$WetIDSel) %>% 
    mutate(SiteID = as.factor(ID))
  return(Prf.flt)
}

Slp <- map_dfr(SiteSel,~Slp.ldfxn(.x,Slp.Sel2)) %>% 
  mutate(pvID = recode_factor(
    pvID,
    "G_1" = "G[1]",
    "G_1_VPD_m" = "G[1] / VPD^{m}",
    "G_v" = "G[v]"
  ))

Slp.Selplt <- Slp.Sel2 %>% 
  mutate(pvID = recode_factor(
    pvID,
    "G_1" = "G[1]",
    "G_1_VPD_m" = "G[1] / VPD^{m}",
    "G_v" = "G[v]"
  ))


# Plot stress signal variability
#===============================

Slp.TrtLvl = Slp %>% filter(ExpID ==594)

ef <- function(x,probs = probs){
  out = ecdf(x)(probs)
  out = 0.4
}
g.slp <- ggplot(Slp,aes(x = beta.nrm,y = SiteID,fill = fct_rev(SiteID)))  +
  stat_density_ridges(
    scale = scl,
    calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = c(0.5),
    alpha = 0.7,
    linetype = "solid",
    lwd = lw,
    from = -1,
    to = 3
  ) +
  geom_point(data = Slp.TrtLvl,pch = mk, stroke = mkst, aes(fill = fct_rev(SiteID))) + 
  coord_cartesian(xlim = c(-1,3),ylim = c(1.5,6.2)) +
  scale_x_continuous(breaks = c(0,2),expand = c(0,0)) +
  labs(x = "Stress Signal", y = "Kernel Density") +
  scale_fill_brewer(type = "qual",palette = "Paired", name = "Site ID") +
  geom_vline(xintercept = 0.4,color = "red",linetype = "dotted",lwd = lw) +
  scale_linetype_manual(name = NULL,
                        breaks=c('30th', 'Median'),
                        values=c('Median'='solid', '30th'=lt30)) + 
  facet_wrap(~pvID,dir = "v",labeller = label_parsed) +
  theme_cowplot(font_size = fs,line_size = lw) + theme(legend.position = 'right',legend.justification = "center",
                            axis.text.y = element_blank()) 
tag_facet(g.slp, open = "",size = fs/2.75) +   
  theme_cowplot(font_size = fs, line_size = lw) + theme(
    legend.position = 'right',
    legend.justification = "center",
    axis.text.y = element_blank()
  ) 
# Save figure
ggsave2(paste0(svpath, "fig-03.pdf"),
        width = 3, height = 5, units = "in")


# Plot Performance variability
#=============================
Prf <- map_dfr(SiteSel,~Prf.ldfxn(.x,Prf.Sel)) %>% 
  mutate(across(contains("AfT"),~1 - abs(.)))
Prf[Prf$mtcID %in% c("AfT1","AfT2","AfT3"),"beta"] = 
  1 - abs(Prf[Prf$mtcID %in% c("AfT1","AfT2","AfT3"),"beta"])
Prf.plt <- Prf %>% group_by(SiteID,mtcID) %>% 
  summarize(mtcVal = median(beta,na.rm = TRUE)) %>% 
  filter(mtcID %in% prf.vars)

# Performance thresholds
th.Prf <-
  tibble(mtcID = as.factor(c("LCE",
                             "AfT1" = "A * minute[f * {list(,)}*T * {list(,)}*1]",
                             "AfT2" = "A * minute[f * {list(,)}*T * {list(,)}*2]",
                             "AfT3" = "A * minute[f * {list(,)}*T * {list(,)}*3]")), 
         th = c(0.4, 0.7, 0.7, 0.7))

Prf.TrtLvl = Prf %>% filter(ExpID ==594,mtcID %in% prf.vars) %>% 
  mutate(mtcID = recode_factor(
    mtcID,
    "LCE" = "LCE",
    "AfT1" = "A * minute[f * {list(,)}*T * {list(,)}*1]",
    "AfT2" = "A * minute[f * {list(,)}*T * {list(,)}*2]",
    "AfT3" = "A * minute[f * {list(,)}*T * {list(,)}*3]"
  ))

Prf.plt <- Prf  %>% filter(mtcID %in% prf.vars) %>% 
  mutate(mtcID = recode_factor(
    mtcID,
    "LCE" = "LCE",
    "AfT1" = "A * minute[f * {list(,)}*T * {list(,)}*1]",
    "AfT2" = "A * minute[f * {list(,)}*T * {list(,)}*2]",
    "AfT3" = "A * minute[f * {list(,)}*T * {list(,)}*3]"
  ))

g.prf <- Prf.plt %>%
ggplot(aes(x = beta,y = SiteID,fill = fct_rev(SiteID)))  +
  stat_density_ridges(
    scale = scl,
    calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = c(0.5),
    alpha = 0.7,
    linetype = "solid",
    lwd = lw,
    from = -1.1,
    to = 1
  ) +
  geom_point(data = Prf.TrtLvl,pch = mk, stroke = mkst, aes(fill = fct_rev(SiteID))) + 
  coord_cartesian(xlim = c(-1.1,1),ylim = c(1.5,6.2)) +
  scale_x_continuous(breaks = c(-1,0,1),expand = c(0,0)) +
  labs(x = "Performance", y = "Kernel Density") +
  scale_fill_brewer(type = "qual",palette = "Paired", name = "Site ID") +
  geom_vline(data = th.Prf,aes(xintercept = th), color = "red",
  linetype = "dotted",lwd = lw) +
  facet_wrap(~mtcID,dir = "v",labeller = label_parsed) +
  theme_cowplot(font_size = fs,line_size = lw) + 
  theme(legend.position = 'right',legend.justification = "center",
                                        axis.text.y = element_blank())
tag_facet(g.prf, open = "", size = fs/2.75) +   
  theme_cowplot(font_size = fs, line_size = lw) + theme(
    legend.position = 'right',
    legend.justification = "center",
    axis.text.y = element_blank()
  ) 
# Save figure
ggsave2(paste0(svpath,"fig-04.pdf"),
        width = 6, height = 4, units = "in")
