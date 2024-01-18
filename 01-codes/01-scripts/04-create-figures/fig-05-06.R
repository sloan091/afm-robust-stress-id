#===============================================================================
# Name: fig-05-06
# By: Brandon Sloan
# Last Updated: 1/17/24
#
# Description: This script creates the ecosystem specific boxplots of stress signal
# (Fig. 5) and performance (Fig. 6) from Sloan and Feng (2023)
#===============================================================================

# Load libraries and add datasets
gc(reset = TRUE)  
rm(list = ls())
source("./01-codes/02-functions/04-robustness/robustness-helper-fxns.R")
load("./01-codes/01-scripts/00-setting-files/final_ec_site_properties_for_r.RData")
load("./03-outputs/02-robust-framework/02-robust-summary/SlpPrfDist_Stats.RData")
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

# Create condensed factors for plotting
#======================================
fD.Slp = Slp.Sel %>%
  filter(IGBP %in% c("ENF", "DBF", "GRA", 
                     "CRO"),!(SlpID == 3 &
                                          Slp_nrm < 0),!(ClimID == "Polar")) 

# Select only on plant variable for performance plots
fD.Prf <- fD.Slp %>% filter(pvID == "G_1")

# Add new dryness index
fD.Slp <-
  fD.Slp %>% left_join(Site[, c("SiteID", "DI.f3")], by = "SiteID", suffix = c("", ".y"))
fD.Prf <-
  fD.Prf %>% left_join(Site[, c("SiteID", "DI.f3")], by = "SiteID", suffix = c("", ".y"))


# Fig. 4: Slope and SNR by biomes
#================================
g1 <- fD.Slp %>%
  ggplot(aes(x = pvID, y = Slp_nrm, fill = DI.f3))  +
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
  labs(x = "Plant Parameter", y = "Median Stress Signal")  +
  scale_fill_brewer(type = "div",
                    direction = -1,
                    name = "DI = PET/P") +
  theme_cowplot(fs, line_size = lw.ax) +
  facet_grid( ~ IGBP, scales = "fixed") +
  coord_cartesian(ylim = c(-1, 3))  +
  scale_y_continuous(breaks = c(0, 2)) +
  scale_x_discrete(labels = c(
    "G_1" = TeX(r"($G_1$)"),
    "G_1_VPD_m" = TeX(r"($G_1/VPD^m$)"),
    "G_v" = TeX(r"($G_v$)")
  )) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom"
  ) +
  geom_hline(aes(yintercept = 0.4, linetype = "Dry Stress"),
             colour = 'red',
             lwd = lw.ax) +
  geom_hline(aes(yintercept = -0.4, linetype = "Wet Stress"),
             colour = 'blue',
             lwd = lw.ax) +
  scale_linetype_manual(
    name = "Critical Threshold",
    values = c(2, 2),
    guide = guide_legend(override.aes = list(color = c("red", "blue")))
  )

g2 <- fD.Slp %>%
  ggplot(aes(x = pvID, y = 1 / Slp_cv, fill = DI.f3)) +
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
  labs(x = "Ecosystem Parameter", y = "Stress SNR")  +
  scale_fill_brewer(type = "div",
                    direction = -1,
                    name = "DI = PET/P") +
  theme_cowplot(fs, line_size = lw.ax) +
  coord_cartesian(ylim = c(0, 1.5))  +
  scale_y_continuous(breaks = c(0, 1)) +
  facet_grid( ~ IGBP, scales = "fixed") +
  scale_x_discrete(labels = c(
    "G_1" = TeX(r"($G_1$)"),
    "G_1_VPD_m" = TeX(r"($G_1/VPD^m$)"),
    "G_v" = TeX(r"($G_v$)")
  )) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
fig4 <- g1 / g2 + plot_layout(guides = "collect") &
  theme(legend.position = 'bottom') &
  plot_annotation(tag_levels = 'a',tag_suffix = ')')
ggsave(
  paste0(svpath, "fig-05.pdf"),
  plot = fig4,
  width = 7,
  height = 4,
  units = "in",
  #device = cairo_pdf
)


# Fig. 5: Prf and SNR by biomes
#================================

# LCE
gpL <- fD.Prf %>% 
  ggplot(aes(x = IGBP, y = LCE_mn,fill = DI.f3)) +
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
  geom_hline(aes(yintercept = 0.4, linetype = "Acceptable Performance"),
             colour = 'red',
             lwd = lw.ax) +
  scale_linetype_manual(
    name = "Critical Threshold",
    values = 2,
   guide = guide_legend(override.aes = list(color = c("red")))
  ) +
  labs(x = "IGBP", y = "Median LCE")  +
  scale_fill_brewer(type= "div",direction = -1, name = "DI = PET/P") + 
  coord_cartesian(ylim = c(0.2,1))  +
  scale_y_continuous(breaks = c(0.2,0.6,1)) +
  theme_cowplot(fs, line_size = lw.ax) +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank()) 


# d. LCE IQR
gpLsnr <- fD.Prf %>% 
  ggplot(aes(x = IGBP, y = 1/LCE_cv,fill = DI.f3)) +
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
  labs(x = "IGBP", y = "LCE SNR")  +
  scale_fill_brewer(type= "div",direction = -1, name = "DI = PET/P") + 
  coord_cartesian(ylim = c(0,3))  +
  scale_y_continuous(breaks = c(0,1.5,3)) +
  theme_cowplot(fs, line_size = lw.ax) +
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5))


# c. AfT Average
gpAf <- fD.Prf %>% 
  ggplot(aes(x = IGBP, y = (AfT1_mn + AfT2_mn + AfT3_mn)/3,fill = DI.f3)) +
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
  geom_hline(aes(yintercept = 0.7, linetype = "Acceptable Performance"),
             colour = 'red',
             lwd = lw.ax) +
  scale_linetype_manual(
    name = "Critical Threshold",
    values = 2,
    guide = guide_legend(override.aes = list(color = c("red")))
  ) +
  labs(x = "IGBP", y = TeX(r"(Median $\bar{A'}_{f,T}$)"))  +
  scale_fill_brewer(type= "div",direction = -1, name = "DI = PET/P") + 
  coord_cartesian(ylim = c(0.2,1))  +
  scale_y_continuous(breaks = c(0.2,0.6,1)) +
  theme_cowplot(fs, line_size = lw.ax) +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())

# c. AfT Average SNR
gpAfsnr <- fD.Prf %>% 
  ggplot(aes(x = IGBP, y = (1/AfT1_cv + 1/AfT2_cv + 1/AfT3_cv)/3,fill = DI.f3)) +
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
  labs(x = "IGBP", y = TeX(r"($\bar{A'}_{f,T}$ SNR)"))  +
  scale_fill_brewer(type= "div",direction = -1, name = "DI = PET/P") + 
  coord_cartesian(ylim = c(0,12))  +
  scale_y_continuous(breaks = c(0,6,12)) +
  theme_cowplot(fs, line_size = lw.ax) +
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5))
(gpL + gpAf)/(gpLsnr + gpAfsnr)+ plot_layout(guides = "collect") & theme(legend.position = 'bottom') &
  plot_annotation(tag_levels = 'a',tag_suffix = ')')

ggsave2(paste0(svpath,"fig-06.pdf"),width = 6,height = 4,units = "in")



