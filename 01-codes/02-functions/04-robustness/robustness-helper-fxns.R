library("latex2exp")

# Load required libraries
#========================
loadPacks = function(...){
  # Loading the essential libraries
  baselibs =  c(
    "hrbrthemes",
    "RColorBrewer",
    "latex2exp",
    "ggplot2",
    "dplyr",
    "forcats",
    "readxl",
    "tibble",
    "tidyverse",
    "ggh4x",
    "ggpubr",
    "tidyr",
    "cowplot",
    "patchwork",
    "broom",
    ...
  )
  # Additional packages
  lapply(baselibs, require, character.only = TRUE) 
}


# Load MATLAB fit results from Excel
#===================================
getFitResults = function(fp,site){
  fn = list.files(fp)
  fnsel = fn[grep(site,fn)]
  D = read_csv(paste0(fp,'/',fnsel))
  D[1:8] = lapply(D[1:8],factor)
  return(D)}


# Find best breakpoint by checking relevant SWC bins
#===================================================
bpfxn <- function(mdl,.x,bins){
  bp <- unname(tidy(davies.test(
    mdl,
    seg.Z = ~ SWC_nep,
    values = bins
  ))$statistic)
  
  # Use breakpoint from Davies test
  mdl.sg <-
    segmented(
      mdl,
      seg.Z = ~ SWC_nep,
      psi = bp,
      control = seg.control(it.max = 0)
    )
  
  # Store summary statistics of each model
  relBIC <- (glance(mdl)$BIC - glance(mdl.sg)$BIC)/abs(glance(mdl)$BIC)
  absR2 <-
    (glance(mdl.sg)$adj.r.squared - glance(mdl)$adj.r.squared)
  relR2 <- absR2/glance(mdl)$adj.r.squared
  relDev <-
    (glance(mdl)$deviance - glance(mdl.sg)$deviance) / abs(glance(mdl)$deviance)
  bpStats <- tibble(bp,relBIC,absR2,relDev,relR2)
  return(bpStats)
}

# Get beta main effects from lm using emmeans or emtrends
#========================================================
getBetaME = function(mdl,Fct,type,iv,df,FctStr){
  emm_options(msg.interaction = FALSE)
  beta <-
    if (type == "Slope") {
      map_dfr(
        Fct,
        ~ tidy(emtrends(mdl, ., var = iv, data = df)) %>%
          rename(
            TrtLvlID = contains("Fct"),
            beta = contains(".trend"),
            se = std.error
          ),
        .id = "TrtID"
      )
    } else{
      map_dfr(Fct,  ~ tidy(emmeans(mdl, ., data = df))%>%
                rename(
                  TrtLvlID = contains("Fct"),
                  beta = contains("estimate"),
                  se = std.error
                ),
              .id = "TrtID"
      )
    }
  # Correct factor levels to make unique
  beta <- beta %>% mutate(TrtID = as.factor(TrtID))
  levels(beta$TrtID) <-  FctStr
  beta$TrtLvlID <- fct_inorder(fct_cross(beta$TrtID,beta$TrtLvlID,sep = ''),
                               ordered = NA)
  return(beta)
}


# Get beta main effects from lm by uncontrasting coefficients
#============================================================
# Helper fxn uncontrast regression coefficients
unconstrastSlp <- function(Fct, tdyprms, relBeta) {
  tdyprms <- tdyprms %>% filter(str_detect(term,relBeta))
  ac <- tdyprms %>% filter(term == relBeta)
  ex <-
    tdyprms %>% filter(str_detect(term, Fct)) %>% select(term, estimate, std.error) %>% 
    rename(se = std.error)
  miss <- 0 - sum(ex$estimate)
  miss.se <- mean(ex$se)
  nn <- paste0("SWC_nep:", Fct, as.character(nrow(ex) + 1))
  out <- add_row(ex, term = nn, estimate = miss, se = miss.se) %>%
    mutate(beta = estimate + ac$estimate, TrtLvlID = str_sub(term,-1,-1))
  return(out)
}
unconstrastInt <- function(Fct, tdyprms, relBeta) {
  ac <- tdyprms %>% filter(term == "(Intercept)")
  tdyprms <- tdyprms %>% filter(str_detect(term,relBeta))
  ex <-
    tdyprms %>% filter(str_detect(term, Fct)) %>% select(term, estimate, std.error) %>% 
    rename(se = std.error)
  miss <- 0 - sum(ex$estimate)
  miss.se <- mean(ex$se)
  nn <- paste0(Fct, as.character(nrow(ex) + 1))
  out <- add_row(ex, term = nn, estimate = miss, se = miss.se) %>%
    mutate(beta = estimate + ac$estimate, TrtLvlID = str_sub(term,-1,-1))
  return(out)
}

# Wrapper function to iterate over multiple factors
getBetaMEFast <- function(Fct, tdyprms, relBeta, FctStr,type) {
  if (type == "Slope") {
    beta <- map_dfr(Fct,  ~ unconstrastSlp(., tdyprms, relBeta),
                    .id = "TrtID")
  } else{
    beta <- map_dfr(Fct,  ~ unconstrastInt(., tdyprms, relBeta),
                    .id = "TrtID")
  }
  beta <- beta %>% mutate(TrtID = as.factor(TrtID))
  levels(beta$TrtID) <-  FctStr
  beta$TrtLvlID <-
    fct_inorder(fct_cross(beta$TrtID, beta$TrtLvlID, sep = ''),
                ordered = NA)
  return(beta)
}


# Plotting script for breakpoints
#================================
plot.bp <- function(D,x,y,c,bp){
g.bp <- D %>%
  ggplot(aes_string(x = x, y = y, color = c)) +
  geom_point() +
  geom_smooth(
    method = "loess",
    se = FALSE,
    size = 1,
    linetype = 1
  ) +
  geom_smooth(method = "lm",
              se = FALSE,
              size = 1,
              linetype = 2) +
  scale_color_manual(
    values = c("#ca0020",
               "#f4a582",
               "#a6cee3",
               "#1f78b4"),
    name = "Prms",
    labels = c(
      '1' = TeX(r"(All)"),
      '2' = TeX(r"($G_1, m$)"),
      '3' = TeX(r"($G_o, G_1$)"),
      '4' = TeX(r"($G_1$)")
    )
  ) +
  geom_vline(
    xintercept = bp,
    linetype = "dashed",
    color = "black",
    size = 2
  ) +
  theme_cowplot() 

return(g.bp)
}

# Re-ordering and renaming treatments and levels
#===============================================
reorder_treatments <- function(In) {
  # Re-ordering treatments
  TrtOrd <- c("Fct5",
              "Fct6",
              "Fct7",
              "Fct8",
              "Fct1",
              "Fct3",
              "Fct2",
              "Fct4")
  

  
  # Re-ordering treatment levels
  TrtLvlOrd <- c(
    "Fct5:1",
    "Fct5:2",
    "Fct5:3",
    "Fct5:4",
    "Fct6:1",
    "Fct6:2",
    "Fct6:3",
    "Fct7:1",
    "Fct7:2",
    "Fct8:1",
    "Fct8:2",
    "Fct1:1",
    "Fct1:2",
    "Fct3:1",
    "Fct3:2",
    "Fct3:3",
    "Fct3:4",
    "Fct2:1",
    "Fct2:2",
    "Fct4:1",
    "Fct4:2",
    "Fct4:3"
  )
  
  Out <-
    In %>% mutate(
      TrtID = fct_relevel(TrtID, TrtOrd),
      TrtLvlID = fct_relevel(TrtLvlID, TrtLvlOrd)
    )

  return(Out)
  
}

reord.TrtOnly.fxn <- function(In) {
  # Re-ordering treatments
  TrtOrd <- c("Fct5",
              "Fct6",
              "Fct7",
              "Fct8",
              "Fct1",
              "Fct3",
              "Fct2",
              "Fct4")
  
  Out <-
    In %>% mutate(
      TrtID = fct_relevel(TrtID, TrtOrd)
    )
  
  return(Out)
  
}

# Load treatment level names
# Rename treatments
TrtNames <- c(
  Fct1 = TeX(r"(5. $G_c$ Equation)"),
  Fct2 = "7. Response Variable",
  Fct3 = "6. Fitting Parameters",
  Fct4 = "8. VPD",
  Fct5 = "1. Fit Algorithm",
  Fct6 = "2. SEB Closure",
  Fct7 = "3. Growing Season",
  Fct8 = "4. LAI"
)

# Renaming Treatment Levels
TrtLvlNames <-  c(
  "Fct1:1" = "Lin18",
  "Fct1:2" = "Med17",
  "Fct2:1" = TeX(r"($G_c$)"),
  "Fct2:2" = "ET",
  "Fct3:1" = "All",
  "Fct3:2" = TeX(r"($G_1, m$)"),
  "Fct3:3" = TeX(r"($G_o, G_1$)"),
  "Fct3:4" = TeX(r"($G_1$)"),
  "Fct4:1" = TeX(r"($VPD_a$)"),
  "Fct4:2" = TeX(r"($VPD_l$)"),
  "Fct4:3" = TeX(r"($Iter \, VPD_l$)"),
  "Fct5:1" = "NLS",
  "Fct5:2" = "Robust NLS",
  "Fct5:3" = "L1",
  "Fct5:4" = "LCE",
  "Fct6:1" = "None",
  "Fct6:2" = "Bowen",
  "Fct6:3" = "All to LE",
  "Fct7:1" = "No GPP filter",
  "Fct7:2" = "50% GPP filter",
  "Fct8:1" = "No LAI",
  "Fct8:2" = "LAI"
)

# Heatmap helper fxn - Continuous color
htHelp = function(htMat, xnm, ynm, htnm, ttl, xlab, ylab) {
  gh = ggplot(htMat,
              aes_string(x = xnm, y = ynm, fill = htnm)) +
    geom_tile(color = "gray", aes(height = 1, width = 1)) +
    coord_fixed() +
    labs(x = xlab, y = ylab) +
    guides(fill = guide_colorbar(
      title = ttl,
      title.position = "top",
      title.hjust = 0.5,
      title.theme = element_text(size  = 8)
    )) +
    theme_cowplot(8) +
    theme(
      legend.position = "top",
      legend.key.width = unit(0.1, "in"),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      legend.justification = "center"
    )
  return(gh)
  
}

# Heatmap helper fxn - Discrete color
htdisHelp = function(htMat, xnm, ynm, htnm, ttl, xlab, ylab) {
  gh = ggplot(htMat,
              aes_string(x = xnm, y = ynm, fill = htnm)) +
    geom_tile(color = "gray", aes(height = 1, width = 1)) +
    coord_fixed() +
    labs(x = xlab, y = ylab) +
    theme_cowplot(8) +
    theme(
      legend.position = "top",
      legend.key.width = unit(0.1, "in"),
      legend.justification = "center",
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )
  return(gh)
  
}
