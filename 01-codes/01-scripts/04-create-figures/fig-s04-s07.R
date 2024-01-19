################################################################################ 
# Name: fig-s06-s09
# By: Brandon Sloan
# Last Updated: 1/16/24
#
# Description: Creates supplementary Figure S6 and S9 from Sloan and Feng (2023).
################################################################################


# Clear workspace, add functions, and load libraries
gc(reset = TRUE)  
rm(list = ls())
source("./01-codes/02-functions/04-robustness/robustness-helper-fxns.R")
loadPacks(c("lemon","broom","ggh4x","segmented","dplyr"))
load("./01-codes/01-scripts/00-setting-files/final_ec_site_properties_for_r.RData")
svpath <- "./03-outputs/03-figures/"

# Axis properties
fs <- 8
lw.ax <- 0.25

# Jitter point properties
stk = 0.05
mrk = 20
sz = 0.5

# Boxplot properties
lw = 0.1
ft = 1

# Fig. S6: Observed Years
#========================
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
  paste0(svpath, "fig-s06.png"),
  width = 3,
  height = 4,
  units = "in"
)

# Fig. S7: Max LCE
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

ggsave2(paste0(svpath,"fig-s09.png"),width = 5,height = 2,units = "in")

