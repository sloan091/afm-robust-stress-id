#===============================================================================
# Name: Fig4_Stress_IGBPvDI
# By: Brandon Sloan
# Last Updated: 2/5/23
# Description: This script tries to explain the presence of soil water stress in
# terms of site characteristics IGBP and DI. I am trimming to ENF, DBF, GRA, and
# CRO land use types as they are the most dominant
#===============================================================================

# Load libraries and add datasets
gc(reset = TRUE)  
rm(list = ls())
source("./Fxns/FLUXNET_Treatment_Avg_Plot_Functions.R")
loadPacks("glmnet","leaps","caret","broom", "furrr","tidymodels","easystats",
          "patchwork")
load(paste0("./Outputs/Robust_AFM/SlpPrfDist_Stats.RData"))
load("./Setting_Files/R_FLUXNET2015_SiteProp.RData")
svpath <- "./Figs/AFM/minor-revision/Supp/"

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

# # Create condensed factors for plotting
# #======================================
# fD.Slp = Slp.Sel %>%
#   filter(!(SlpID == 3 &
#              Slp_nrm < 0), !(ClimID == "Polar")) 
# 
# # Select only on plant variable for performance plots
# fD.Prf <- fD.Slp %>% filter(pvID == "G_1")
# fD.Prf2 <-
#   fD.Prf %>%  pivot_longer(
#     cols = contains(all_of(c(
#       "LCE", "AfT1", "AfT2", "AfT3"
#     ))) & !contains(".max"),
#     names_to = c("mtcID", "mtcType"),
#     names_sep = "_",
#     values_to = "Score"
#   )
# 
# # Add new dryness index
# fD.Slp <-
#   fD.Slp %>% left_join(Site[, c("SiteID", "DI.f3")], by = "SiteID", suffix = c("", ".y"))
# fD.Prf <-
#   fD.Prf %>% left_join(Site[, c("SiteID", "DI.f3")], by = "SiteID", suffix = c("", ".y"))

df_slp <- Slp.Sel |> filter(pvID %in% "G_v") |> 
  mutate(Slp_nrm = Slp_mn/MaxMed)

# Discrete model testing

mdl_dis <- lm(Slp_nrm ~ IGBP * DI.f3,data = df_slp)
summary(mdl_dis)
anova(mdl_dis)
plot(estimate_means(mdl_dis)) + 
  ylab("Effect on Soil Water Stress Signal") +
  scale_color_brewer(type = "div",
                    direction = -1,
                    name = "DI = PET/P") +
  theme_cowplot(11)


check_model(mdl_dis)


# continuous model
mdl_cnt <- lm(Slp_nrm ~ IGBP + DI,data = df_slp)
summary(mdl_cnt)
anova(mdl_cnt)
check_model(mdl_cnt)

boots <- bootstraps(df_slp,times = 5e2)
lm_mod <- linear_reg()
slp_rec <- recipe(Slp_nrm ~ IGBP + DI, data = df_slp) |> 
  step_zv()

slp_wf <- workflow() |> 
  add_recipe(slp_rec) |> 
  add_model(lm_mod)

get_lm_coefs <- function(x) {
  mdl <- x |> 
    # get the lm model object
    extract_fit_engine()
  prms <- mdl |>  
    # transform its format
    tidy()
  prf <- glance(mdl) |> select(adj.r.squared,p.value)
  prms <- prms |> mutate(R2_adj = prf$adj.r.squared,
                         F_p = prf$p.value)
  
  
}
tidy_ctrl <- control_grid(extract = get_lm_coefs)

slp_fit <- slp_wf |> 
  fit_resamples(resamples = boots, control = tidy_ctrl)
slp_fit
collect_metrics(slp_fit)

slps <- slp_fit |> unnest(.extracts) |> unnest(.extracts)

g1 <- slps %>%
  filter(term == "DI") %>% 
  ggplot(aes(x = estimate,y = after_stat(count/500))) +
  geom_histogram(bins = 10) +
  labs(y = "Frequency", x = TeX(r"(DI Effect on $G_1$ Soil Water Stress Signal)")) +
  theme_cowplot(9)

g2 <- slp_fit |> unnest(.metrics) |> filter(.metric %in% "rsq") |> 
  ggplot(aes(x = .estimate, y = after_stat(count/500))) +
  geom_histogram(bins = 10) +
  labs(y = NULL, x = TeX(r"($R^2$ on Assessment Set)")) +
  theme_cowplot(9)

g3 <- slps %>%
  filter(term == "DI") %>% 
  mutate(F_sig = ifelse(F_p <= 0.05, "Significant", "Insignificant")) |> 
  ggplot(aes(x = F_sig)) +
  geom_bar(aes(y = after_stat(count/500))) +
  labs(y = NULL, x = "F-test Results") +
  theme_cowplot(9)

fig <- g1 + g2 + g3 +
plot_annotation(tag_levels = 'a',tag_suffix = ')')
ggsave(
  paste0(svpath, "ecosystem-sig-DIcnt.pdf"),
  plot = fig,
  width = 6.5,
  height = 3,
  units = "in",
  #device = cairo_pdf
)

# Analysis with all assumption sets

files <- list.files("./Outputs/SlpPrfDist_Fits", pattern = "*.RData",
                    full.names = TRUE)

rd <- function(file){
  load(file)
  Prf <- Prf |> filter(mtcID %in% c("LCE"))
  Slp <- Slp |> filter(pvID %in% c("G_1"))
}

prf <- map(files, \(x) rd(x)) |> list_rbind(names_to = "rowid") |> 
  left_join(Site |> rowid_to_column() |> select(IGBP,DI.f3,DI,rowid), by = "rowid")

# continuous model all sims
mdl_cnt <- lm(beta.nrm ~ IGBP + DI,data = prf)
summary(mdl_cnt)
anova(mdl_cnt)
check_heteroscedasticity(mdl_cnt)




test <- aov(beta.nrm~IGBP*DI.f3,data = prf)
TukeyHSD(test)
test <- tidy(mdl)
rem <- test$term[test$estimate  %in% NA]
x <- model.matrix(mdl)
drop <- which(colnames(x) %in% rem)
x_mod <-  x[,-drop]
mdl2 <- Slp.Sel %>% filter(pvID %in% "G_v") %>%
  lm(Slp_mn ~ x_mod,data = .)
# Anova(mdl, type = "III")
# # Extract the residuals
# aov_residuals <- residuals(object = mdl)
# # Run Shapiro-Wilk test
# shapiro.test(x = aov_residuals )

