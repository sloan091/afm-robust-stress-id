################################################################################ 
# Name: fig-01
# By: Brandon Sloan
# Last Updated: 1/16/24
#
# Description: Creates versions of the FLUXNET 2015 sites in terms of ecosystem
# and climatic categories for Sloan and Feng (2023)
################################################################################


# Clear workspace, add functions, and load libraries
gc(reset = TRUE)  
rm(list = ls())
source("./01-codes/02-functions/04-robustness/robustness-helper-fxns.R")
load("./01-codes/01-scripts/00-setting-files/final_ec_site_properties_for_r.RData")
loadPacks(c("lemon","broom","ggh4x","segmented","dplyr"))
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

# Site plot - IGBP vs. DI
#========================
ggplot(data = Site, aes(x = IGBP, fill = DI.f3)) +
  geom_bar(color = 'black',lwd = lw) +
  labs(x = "IGBP", y = "No. of Sites") +
  theme_cowplot(fs, line_size = lw.ax) +
  theme(
    legend.text.align = 0,
    legend.position = c(0.6, 0.75),
    axis.text.x = element_text(angle = 30, vjust = 0.5)
  ) +
  scale_fill_brewer(type = "div",
                    direction = -1,
                    name = "DI = PET/P")  +
  scale_y_continuous(lim = c(0, 40), breaks = c(0, 20, 40))

ggsave2(
  paste0(svpath, "fig-01.pdf"),
  width = 2.75,
  height = 2,
  units = "in"
)

