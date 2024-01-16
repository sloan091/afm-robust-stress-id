################################################################################ 
# Name: FigS2andS6
# By: Brandon Sloan
# Last Updated: 3/2/23
# Description: Creates the observation period and maximum achievable LCE for the
# 150 FLUXNET sites used in Sloan and Feng (2023).
################################################################################


# Clear workspace, add functions, and load libraries
gc(reset = TRUE)  
rm(list = ls())
source("./Fxns/FLUXNET_Factorial_Analysis_Functions.R")
loadPacks(c("lemon","broom","ggh4x","segmented","dplyr"))
load("./Data/R_FLUXNET2015_SiteProp.RData")
svpath <- "./Figs/Paper/Fig2/"
v <- "_v2"

# Axis properties
fs <- 8
#ft <- c("CMU Classical Serif")
lw.ax <- 0.25

# Jitter point properties
stk = 0.05
mrk = 20
sz = 0.5

# Boxplot properties
lw = 0.1
ft = 1

# Fig. S2: Observation Years
#===========================
g.hst <- Site %>%
  ggplot(aes(x = nyrs)) +
  geom_histogram(bins = 10, fill = 'white', color = 'black') +
  theme_cowplot(fs, line_size = lw.ax) +
  labs(x = "Observed Years", y = "Number of Sites") +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())

g.cdf <-
  Site %>% ggplot(aes(x = nyrs)) + stat_ecdf(geom = "step") +
  theme_cowplot(fs, line_size = lw.ax) +
  labs(x = "Observed Years", y = "Empirical CDF")
g.hst / g.cdf + plot_annotation(tag_levels = 'a', tag_suffix = ')')
ggsave2(
  paste0(svpath, "FigS_nyrsObserved", v, ".png"),
  width = 3,
  height = 4,
  units = "in"
)


# Fig. S6: Max LCE
#==================
Site %>% 
  ggplot(aes(x = IGBP, y = LCE.max,fill = DI.f3)) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge2(preserve = "total"),
    linewidth = lw,
    fatten = ft
  ) +
  geom_point(
    position = position_jitterdodge(),
    pch = mrk,
    stroke = stk,
    size = sz
  ) +
  labs(x = "IGBP", y = "LCE Upper Bound")  +
  scale_fill_brewer(type= "div",direction = -1, name = "DI = PET/P") + 
  coord_cartesian(ylim = c(0,1))  +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme_cowplot(fs, line_size = lw.ax) +
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5))
ggsave2(paste0(svpath,"FigS_LCE_mx_DIvIGBP",v,".png"),width = 5,height = 2,units = "in")


