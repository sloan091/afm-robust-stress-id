################################################################################ 
# Name: Fig8_Example_Distributions
# By: Brandon Sloan
# Last Updated: 1/19/23
# Description: Plots the distributions of slopes and performance for selected
# sites.
################################################################################

# Load libraries and add datasets
gc(reset = TRUE)  
rm(list = ls())
source("./Fxns/FLUXNET_Treatment_Avg_Plot_Functions.R")
loadPacks("ggridges","ggdensity")
select <- dplyr::select
load("./Setting_Files/R_FLUXNET2015_SiteProp.RData")
load("./Outputs/Robust_AFM/FLX_Factorial_SlpPrf_Stats.RData")
svpath = "./Figs/AFM/Fig1/"

prf.vars = c("LCE","AfT1","AfT2","AfT3")
v <- "_v3"

# Plot settings
# Font size
fs <- 6
lw <- 0.25
lt30 <- "dotted"
mk <- 21
mkst <- 0.25
scl <- 1
f.wd = 1
f.ar = 4 / 3
f.buff = 0.25

# Select stress periods and load data
#====================================
# Select Sites
SiteSel <- c("US_Me2")
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
  load(paste0("./Outputs/AS_Fits/",ID,"_AllFits.RData"))
  flt <- flt %>% filter(SiteID == ID)
  Slp.flt <- Slp %>% left_join(flt,by = c("pvID")) %>% 
    filter(WetID == WetIDSel) %>% mutate(SiteID = as.factor(ID))
  return(Slp.flt)
}

# Helper functions to load performance files
Prf.ldfxn <- function(ID,flt){
  load(paste0("./Outputs/AS_Fits/",ID,"_AllFits.RData"))
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

Slp.TrtLvl = Slp %>% filter(ExpID %in% c(594,604))

ef <- function(x,probs = probs){
  out = ecdf(x)(probs)
  out = 0.4
}
gps1 <- ggplot(Slp,aes(x = beta.nrm,y = pvID,fill = pvID))  +
stat_density_ridges(
  scale = scl,
  calc_ecdf = TRUE,
  quantile_lines = TRUE,
  quantiles = c(0.3),
  alpha = 0.7,
  linetype = lt30,
  lwd = 0.99*lw
) +
  stat_density_ridges(
    scale = scl,
    calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = c(0.5),
    alpha = 0.7,
    linetype = "solid",
    lwd = lw,
    fill = NA
  ) +
  geom_point(data = Slp.TrtLvl,pch = mk, stroke = mkst) + 
  coord_cartesian(xlim = c(-1,2.5),ylim = c(1.5,3.5)) +
  scale_x_continuous(breaks = c(0,2)) +
  labs(x = "Stress Signal", y = "Kernel Density") +
  scale_fill_brewer(type = "qual",palette = "Paired", name = "Site ID") +
  scale_linetype_manual(name = NULL,
                        breaks=c('30th', 'Median'),
                        values=c('Median'='solid', '30th'=lt30)) + 
  #facet_wrap(~pvID,dir = "v",labeller = label_parsed) +
  theme_cowplot(font_size = fs,line_size = lw) + 
  theme(legend.position = 'right',legend.justification = "center",
                            axis.text.y = element_blank()) +
geom_vline(xintercept = 0.4, lwd = lw, linetype = "dashed") + 
  force_panelsizes(cols = unit(f.wd, "in"), rows = unit(f.wd/ f.ar * 3, "in"))
plot(gps1)
ggsave2(paste0(svpath,"Fig1_Slp_Robust.pdf"))


# Plot Performance variability
#=============================
Prf <- map_dfr(SiteSel,~Prf.ldfxn(.x,Prf.Sel))
Prf[Prf$mtcID %in% c("AfT1","AfT2","AfT3"),"beta"] = 
  1 - abs(Prf[Prf$mtcID %in% c("AfT1","AfT2","AfT3"),"beta"])

# Calculate average info metric
Prf.dis <- Prf %>% select(-se) %>% 
  pivot_wider(values_from = beta,names_from = mtcID) %>% 
  mutate(AfT = (AfT1 + AfT2 + AfT3)/3)

Prf.TrtLvl = Prf.dis %>% filter(ExpID %in% c(594,604))


gxy <- ggplot(Prf.dis, aes(x = LCE, y = AfT)) +
  geom_hdr(xlim = c(0, 1), ylim = c(-1, 1),
                 probs = c(0.95, 0.75, 0.5, 0.25, 0.05)) + 
  # geom_hdr_lines(xlim = c(0, 1), ylim = c(-1, 1),
  #                probs = c(0.99,0.95, 0.75, 0.5, 0.25, 0.05, 0.01), size = lw) + 
  coord_cartesian(ylim = c(-1,1), xlim = c(0,1)) + 
  geom_vline(xintercept = 0.4, lwd = lw, linetype = "dashed") +
  geom_hline(yintercept = 0.7, lwd = lw, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  labs(x = "LCE", y = TeX(r"($\bar{A'}_{f,T}$)")) +
  theme_cowplot(fs, line_size = lw) + 
  # theme(legend.position = "bottom",
  #       legend.justification = "center") + 
  force_panelsizes(cols = unit(f.wd, "in"), rows = unit(f.wd/ f.ar, "in"))
plot(gxy)
ggsave2(paste0(svpath,"Fig1_Prf_Rob_Jnt.pdf"))
# A plotting option
gx <- ggplot(Prf.dis,aes(x = LCE))  + 
  geom_density(alpha = 0.2, fill = "red") + 
  geom_point(data = Prf.TrtLvl,aes(y = 0)) + 
  coord_cartesian(xlim = c(0,1)) + 
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "LCE", y = "Density") +
  theme_cowplot(fs, line_size = lw) + 
  theme(legend.position = "top",
        legend.justification = "center",
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
geom_vline(xintercept = 0.4, lwd = lw, linetype = "dashed") + 
  force_panelsizes(cols = unit(f.wd, "in"), rows = unit(f.wd/ f.ar, "in"))
plot(gx)
ggsave2(paste0(svpath,"Fig1_Prf_Rob_LCE.pdf"))

# A plotting option
gy <- ggplot(Prf.dis,aes(x = AfT))  + 
  geom_density(alpha = 0.2, fill = "red") + 
  geom_point(data = Prf.TrtLvl,aes(y = 0)) + 
  coord_flip(xlim = c(-1,1)) + 
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  labs(y = "Density", x = "AfT") +
  theme_cowplot(fs, line_size = lw) + 
  theme(legend.position = "none",
        legend.justification = "center",
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
geom_vline(xintercept = 0.7, lwd = lw, linetype = "dashed") + 
#force_panelsizes(cols = unit(f.wd, "in"), rows = unit(f.wd/ f.ar, "in"))
force_panelsizes(cols = unit(f.wd/ f.ar, "in"), rows = unit(f.wd/ f.ar, "in"))
plot(gy)
ggsave2(paste0(svpath,"Fig1_Prf_Rob_AfT.pdf"))

# wrap_plots(
#   gx, 
#   plot_spacer(), 
#   gxy,
#   gy,
#   nrow = 2
# )
# ggsave2("./Figs/Paper/Fig1/Fig1_Prf_Robust.pdf")
# ggsave2("./Figs/Paper/Fig1/Fig1_Joint.pdf",width = 7,height = 5.25,units = "in")
