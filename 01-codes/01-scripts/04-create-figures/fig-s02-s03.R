#===============================================================================
# Name: fig-s02-s03
# By: Brandon Sloan
# Last Updated: 2/5/23
#
# Description: This script creates Figure S2-S3 in the supplement of Sloan and 
# Feng (2023)
#===============================================================================

# Load libraries and add datasets
gc(reset = TRUE)  
rm(list = ls())
source("./01-codes/02-functions/04-robustness/robustness-helper-fxns.R")
load("./01-codes/01-scripts/00-setting-files/final_ec_site_properties_for_r.RData")
load("./03-outputs/02-robust-framework/02-robust-summary/SlpPrfDist_Stats.RData")
loadPacks("glmnet","leaps","caret","broom", "furrr","tidymodels","easystats",
          "patchwork")
select <- dplyr::select
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

df_slp <- Slp.Sel |> filter(pvID %in% "G_1") |> 
  mutate(Slp_nrm = Slp_mn/MaxMed)

# Discrete model testing
mdl_dis <- lm(Slp_nrm ~ IGBP + DI.f3,data = df_slp)
# These two outputs are combined for Figure S2
summary(mdl_dis)
anova(mdl_dis)

# Helpful Visualization
plot(estimate_means(mdl_dis)) + 
  ylab("Effect on Soil Water Stress Signal") +
  scale_color_brewer(type = "div",
                    direction = -1,
                    name = "DI = PET/P") +
  theme_cowplot(11)

check_model(mdl_dis)


# Continuous model analysis (Figure S3)
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
  paste0(svpath, "fig-s03.pdf"),
  plot = fig,
  width = 6.5,
  height = 3,
  units = "in",
  #device = cairo_pdf
)

