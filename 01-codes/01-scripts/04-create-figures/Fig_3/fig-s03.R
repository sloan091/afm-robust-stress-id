#===============================================================================
# Name: fig-s03
# By: Brandon Sloan
# Last Updated: 1/16/23
#
# Description: This script creates Figure S3 of Sloan and Feng (2023) 
# summarizing the count and agreement of soil water stress signals
#===============================================================================

# Load libraries and add datasets
gc(reset = TRUE)  
rm(list = ls())
source("./01-codes/02-functions/04-robustness/robustness-helper-fxns.R")
load("./03-outputs/02-robust-framework/02-robust-summary/SlpPrfDist_Stats.RData")


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

# Colors for slope categories
bc2 <- scale_fill_manual(
  name = "Stress Class",
  labels = c('Dry', 'Wet', "Both", "Negligible", "Unsure"),
  values = c("#e41a1c", "#377eb8", "#984ea3", "#deebf7", "white"),
  drop = FALSE
)

# Figure 2a: Show count of stress response by plant variable
Slp.Cat <- Slp.Cat %>% mutate(SlpID = fct_expand(SlpID,"5"))
g2a <-
  ggplot(Slp.Cat, aes(x = pvID, fill = SlpID)) + 
  geom_bar(color = "black", position = "stack", lwd = lw) +
  scale_x_discrete(labels = c(
    "G_1" = TeX(r"($G_1$)"),
    "G_1_VPD_m" = TeX(r"($G_1/VPD^m$)"),
    "G_v" = TeX(r"($G_v$)")
  )) + bc2 +
  scale_y_continuous(breaks = c(0,75,150)) +
  labs(x = "Plant Parameter", y = "Number of Sites") +
  theme_cowplot(fs, line_size = lw.ax) +
  theme(legend.position = "none")

counts <- Slp.Cat %>% group_by(pvID,SlpID) %>% summarize(n = n(),prc = n()/151)

# Figure 2b: Count consistent plant variable response
Slp.Cat2 <- Slp.Cat %>% filter(pvID %in% "G_1") %>% 
  mutate(SlpID2 = ifelse(conID %in% "3","5",SlpID))
g2b <-
  ggplot(Slp.Cat2, aes(x = conID, fill = SlpID2)) +
  geom_bar(color = "black", position = "stack", lwd = lw) +
  scale_x_discrete(labels = c("1" = "3 of 3", "2" = "2 of 3", "3" = "No Match")) +
  scale_y_continuous(breaks = c(0,40,80))+
  bc2 +
  labs(x = "Matching Stress Response", y = "") +
  theme_cowplot(fs, line_size = lw.ax) +
  theme(legend.position = "top")

# Combine Figure 2
(g2a + g2b) + plot_layout(guides = "collect")&theme(legend.position = "right") &
  plot_annotation(tag_levels = 'a',tag_suffix = ')')
ggsave2(paste0("./03-outputs/03-figures/fig-s03.pdf"),
        width = 4,height = 2,units = "in" )

