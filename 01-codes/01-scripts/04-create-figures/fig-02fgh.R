################################################################################
# Name: fig-02-fgh
# By: Brandon Sloan
# Last Updated: 1/19/2024
#
# Description: Creates Figures 2f, 2g, and 2h in Sloan and Feng (2023)
################################################################################


# Clear workspace, add functions, and load libraries
gc(reset = TRUE)
rm(list = ls())
source("./01-codes/02-functions/04-robustness/robustness-helper-fxns.R")
loadPacks(c("lemon", "broom", "ggh4x", "segmented", "dplyr","ggridges"))
load("./01-codes/01-scripts/00-setting-files/final_ec_site_properties_for_r.RData")
load("./03-outputs/02-robust-framework/02-robust-summary/SlpPrfDist_Stats.RData")
svpath <- "./03-outputs/03-figures/"
select <- dplyr::select

# Plot settings
fs = 6
ls = 0.25
ms = 1
f.ht = 0.75
f.ar = 4 / 3
f.buff = 0.25
mk <- 21
mkst <- 0.25
scl <- 1
lt30 <- "dotted"

# Plot options
xlm = c(-1,3)
xlmp = c(0,1)


# Load file and setpaths
datpath = "./03-outputs/01-estimated-pmoc-parameters/02-excel/"
fn = list.files(datpath, pattern = "*csv")

# Analysis settings
yvars = c("G_1", "G_1_VPD_m", "G_v")

# Names for treatments
TrtStr = c("Gc", "Flx", "Prm", "VPD", "Fit", "SEB", "GPP", "LAI")

# Select US-Me2
ii = 124
fnii = fn[ii]
Site.ID = substr(fnii, start = 1, stop = 6)
D <- na.omit(getFitResults(datpath, Site.ID))

# Select the site
Slp.Stat.Sel <- Slp.Stat %>% filter(SiteID %in% Site.ID)

# Restructure data for plotting
Dplt <- D %>%
  as_tibble()  %>%
  pivot_longer(
    cols = all_of(yvars),
    names_to = "pvID",
    values_to = "pvVal",
    names_transform = list(VarID = as.factor)
  ) %>%
  left_join(Slp.Stat.Sel[, c("pvID", "MaxMed")], by = c("pvID"))

# Median by SM bins
D.med <- Dplt %>%
  group_by(pvID, SWC_nep) %>%
  summarize(Med = median(pvVal, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pvID %in% "G_1")



# Select G1 only w/ Lin
Dplt1 <- Dplt %>%
  filter(
    pvID %in% c("G_1"),
    Fct1 == 1,
    Fct2 == 1,
    Fct3 == 4,
    Fct4 == 2,
    Fct5 == 1,
    Fct6 == 1,
    Fct7 == 2,
    Fct8 == 1
  ) %>%
  mutate(pvID = recode_factor(pvID,
                              "G_1" = "G[1]"))
Dplt1$WetID = as.factor(ifelse(Dplt1$SWC_nep < 0.5, "Wet1", "Dry1"))

# Select all parameters fit w/ Medlyn
Dplt2 <- Dplt %>%
  filter(
    pvID %in% c("G_1"),
    Fct1 == 2,
    Fct2 == 1,
    Fct3 == 1,
    Fct4 == 2,
    Fct5 == 1,
    Fct6 == 1,
    Fct7 == 2,
    Fct8 == 1
  ) %>%
  mutate(pvID = recode_factor(pvID,
                              "G_1" = "G[1]"))
Dplt2$WetID = as.factor(ifelse(Dplt2$SWC_nep < 0.5, "Wet2", "Dry2"))


# Figure 2f: Slope extraction for two example assumption sets
#===========================================================
# Plot all assumption set fits
g1 <- Dplt %>%
  filter(pvID %in% c("G_1")) %>%
  mutate(pvID = recode_factor(pvID,
                              "G_1" = "G[1]")) %>%
  ggplot(aes(x = SWC_nep, y = pvVal)) +
  geom_jitter(
    width = 0.03,
    pch = 4,
    color = "gray",
    size = ms * 0.2,
    stroke = 0.1
  ) +
  labs(x = TeX(r"($\theta_p$)"),
       y = TeX(r"($G_1 \, [kPa^m]$)")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 10)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 5, 10)) +
  theme_cowplot(fs, line_size = ls) +
  theme(legend.text.align = 0, legend.position = 'none')
plot(g1)

# Create separate red fit
g2 <-  g1 +
  geom_point(
    data = Dplt1,
    pch = 20,
    color = "red",
    size = ms
  ) +
  geom_smooth(
    data = Dplt1,
    method = lm,
    formula = y ~ splines::bs(
      x,
      df = 2,
      degree = 1,
      knots = 0.5
    ),
    size = ls,
    se = FALSE,
    color = "red"
  )


g3 <-  g2 +
  geom_point(
    data = Dplt2,
    pch = 20,
    color = "blue",
    size = ms
  ) +
  geom_smooth(
    data = Dplt2,
    method = lm,
    formula = y ~ splines::bs(
      x,
      df = 2,
      degree = 1,
      knots = 0.5
    ),
    size = ls,
    se = FALSE,
    color = "blue"
  )
plot(g3)

g4 <-  g3 +
  geom_point(
    data = D.med,
    aes(y = Med),
    pch = 18,
    color = "black",
    size = ms
  ) +
  force_panelsizes(cols = unit(f.ar * f.ht, "in"), rows = unit(f.ht, "in"))
plot(g4)

ggsave2(
  paste0(svpath,"fig-02-f.jpg"),
  width = f.ar * f.ht + f.buff,
  height = f.ht + f.buff,
  units = "in",
  dpi = 1200
)


# Figure 2h: AfT extraction for two example assumption sets
#=========================================================

# Plot all treatment assumptions
gp1 <- Dplt %>% 
  filter(pvID %in% c("G_1")) %>%
  mutate(pvID = recode_factor(pvID,
                              "G_1" = "G[1]")) %>%
  ggplot(aes(x = SWC_nep, y =  1 - (AfT1 + AfT2 + AfT3) / 3)) +
  geom_jitter(
    width = 0.03,
    pch = 4,
    color = "gray",
    size = ms * 0.2,
    stroke = 0.1
  ) +
  labs(x = TeX(r"($\theta_p$)"),
       y = TeX(r"($\bar{A'}_{f,T}$)")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.5, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_cowplot(fs, line_size = ls) +
  theme(legend.text.align = 0, legend.position = 'none')
plot(gp1)

# Create separate red fit
gp2 <-  gp1 +
  geom_point(
    data = Dplt1,
    pch = 20,
    color = "red",
    size = ms
  ) +
  geom_smooth(
    aes(color = WetID),
    data = Dplt1,
    method = lm,
    formula = y ~ 1,
    size = ls,
    se = FALSE
  )


gp3 <-  gp2 +
  geom_point(
    data = Dplt2,
    pch = 20,
    color = "blue",
    size = ms
  ) +
  geom_smooth(
    aes(color = WetID),
    data = Dplt2,
    method = lm,
    formula = y ~ 1,
    size = ls,
    se = FALSE
  ) +
  scale_color_manual(values = c(
    'Wet2' = 'blue',
    'Dry2' = 'blue',
    'Wet1' = 'red',
    'Dry1' = 'red'
  )) +
  geom_hline(yintercept = 0.7, size = ls, linetype = "dashed") + 
  force_panelsizes(cols = unit(f.ar * f.ht, "in"), rows = unit(f.ht, "in"))
plot(gp3)

ggsave2(
  paste0(svpath,"fig-02h.jpg"),
  width = f.ar * f.ht + f.buff + 0.05,
  height = f.ht + f.buff,
  units = "in",
  dpi = 1200
)


# Figure 2g: LCE extraction for two example assumption sets
#=========================================================

# Plot all treatment assumptions
gp1 <- Dplt %>%
  filter(pvID %in% c("G_1")) %>%
  mutate(pvID = recode_factor(pvID,
                              "G_1" = "G[1]")) %>%
  ggplot(aes(x = SWC_nep, y =  LCE)) +
  geom_jitter(
    width = 0.03,
    pch = 4,
    color = "gray",
    size = ms * 0.2,
    stroke = 0.1
  ) + 
  labs(x = TeX(r"($\theta_p$)"), y = TeX(r"($LCE$)")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_cowplot(fs, line_size = ls) +
  theme(legend.text.align = 0, legend.position = 'none')
plot(gp1)

# Create separate red fit
gp2 <-  gp1 +
  geom_point(
    data = Dplt1,
    pch = 20,
    color = "red",
    size = ms
  ) +
  geom_smooth(
    aes(color = WetID),
    data = Dplt1,
    method = lm,
    formula = y ~ 1,
    size = ls,
    se = FALSE
  )

# Create separate blue fit
gp3 <-  gp2 +
  geom_point(
    data = Dplt2,
    pch = 20,
    color = "blue",
    size = ms
  ) +
  geom_smooth(
    aes(color = WetID),
    data = Dplt2,
    method = lm,
    formula = y ~ 1,
    size = ls,
    se = FALSE
  ) +
  scale_color_manual(values = c(
    'Wet2' = 'blue',
    'Dry2' = 'blue',
    'Wet1' = 'red',
    'Dry1' = 'red'
  )) +
  geom_hline(yintercept = 0.4, size = ls, linetype = "dashed") + 
  force_panelsizes(cols = unit(f.ar * f.ht, "in"), rows = unit(f.ht, "in"))
plot(gp3)
ggsave2(
  paste0(svpath,"fig-02g.jpg"),
  width = f.ar * f.ht + f.buff,
  height = f.ht + f.buff,
  units = "in",
  dpi = 1200
)

