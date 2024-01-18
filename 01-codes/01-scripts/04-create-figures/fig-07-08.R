#===============================================================================
# Name: fig-07-08
# By: Brandon Sloan
# Last Updated: 1/17/24
#
# Description: This script creates the robust ecosystem soil water stress
# heat maps shown in Figures 7-8 in Sloan and Feng (2023). Part of the script
# creates the overall robustness rating used in the paper.
#===============================================================================

# Load libraries and add datasets
gc(reset = TRUE)
rm(list = ls())
source("./01-codes/02-functions/04-robustness/robustness-helper-fxns.R")
load("./01-codes/01-scripts/00-setting-files/final_ec_site_properties_for_r.RData")
load("./03-outputs/02-robust-framework/02-robust-summary/SlpPrfDist_Stats.RData")
loadPacks("DescTools","ggpattern")
select <- dplyr::select
svpath <- "./03-outputs/03-figures/"

# Inner point stroke size
stk <- 0.8

# Multiplier for marker size
mlt <- 0.5

# Thresholds
#===========

# Practically significant slope threshold [% pv drop/SWC_nep]
th.Slp <- 0.2 / 0.5

# Performance thresholds
Prf.vars <- c("LCE_mn", "AfT1_mn", "AfT2_mn", "AfT3_mn")
th.Prf <-
  tibble(mtcID = as.factor(Prf.vars), 
         th = c(0.4, 0.7, 0.7, 0.7))

# Robust thresholds
RobPrf.th <- 0.7
RobSlp.th <- 0.7


#===============================================================================
# 1. Create overall Robust Stress Signal Score
#===============================================================================

# Specify robustness for each plant variable 
# Rb1: Both th met, Rb2: Prf th met, Rb3: Slp th met

RbSlp <- Slp.Sel  %>% filter(!(SlpID == 3 & Slp_nrm < 0)) %>%
  select(SiteID,pvID,IGBP,RobSlp,RobPrf.jnt,SlpID,Slp_nrm) %>% 
  mutate(RbPrf.ID = RobPrf.jnt >= RobPrf.th,
         RbSlp.ID = RobSlp >= RobSlp.th) %>%
  mutate(RbID = as.character(
    case_when(
      RbPrf.ID & RbSlp.ID ~ "Rb1",RbPrf.ID & !RbSlp.ID ~ "Rb2",
      !RbPrf.ID &
        RbSlp.ID ~ "Rb3",
      TRUE ~ ""
    )
  ),
  StrID = as.character(
    case_when(
      SlpID %in% c("1", "3") ~ "Dry Stress",
      SlpID %in% c("2") ~ "Wet Stress",
      SlpID %in% c("4") ~ "No Stress",
      TRUE ~ "Undetermined"
    )
  ),
Rbrnk =
    case_when(
      RbPrf.ID & RbSlp.ID ~ 3,RbPrf.ID & !RbSlp.ID ~ 2,
      !RbPrf.ID &
        RbSlp.ID ~ 1,
      TRUE ~ 0
    )
  ) %>% mutate(RbStrID = as.factor(paste(StrID, RbID)), 
                across(c(StrID, RbID), as.factor)) %>% 
  mutate(RbStrID = fct_relevel(
    RbStrID,
    c(
      "Dry Stress Rb1",
      "Dry Stress Rb2",
      "Dry Stress Rb3",
      "Dry Stress ",
      "No Stress Rb1",
      "No Stress Rb2",
      "No Stress Rb3",
      "No Stress ",
      "Undetermined",
      "Wet Stress ",
      "Wet Stress Rb3",
      "Wet Stress Rb2",
      "Wet Stress Rb1"
    )
  ))

grbfxn <- function(x){
  
  # Specify number of unique stress responses
  n1 <- unique(x$RbStr.n) 
  n2 <- unique(x$Str.n) 
    
  # If two or more Robust Stress categories match
  x1 <- x %>% summarize(id = which.max(n),
                        RbStrID = RbStrID[id],
                        StrID = StrID[id],
                        Slp_nrm = Slp_nrm[id],
                        Rbrnk = Rbrnk[id])
  
  # If only two or more Stress categories match, grab minimum robustness category
  x2 <- x %>% group_by(StrID) %>% 
    summarize(n = n(),
                  StrID = unique(StrID),
                  Slp_nrm = mean(Slp_nrm),
                  Rbrnktmp = min(Rbrnk),
                  RbStrID = RbStrID[which.min(Rbrnk)]) %>% ungroup() %>% 
    summarize(id = which.max(n),
              RbStrID = RbStrID[id],
              StrID = StrID[id],
              Slp_nrm = Slp_nrm[id],
              Rbrnk = Rbrnktmp[id])
  
  # If no matches in stress categories, we are unsure
  x3 <- x %>% summarize(id = 0,
                        RbStrID = "Unsure",
                        StrID = "Unsure",
                        Slp_nrm = mean(Slp_nrm,na.rm = TRUE),
                        Rbrnk = 0)
  if(n1 < 3) {
    out <- x1
  }else if(n1 == 3 & n2 < 3) {
    out <- x2
  } else{
    out <- x3
  }
  
  return(out)
}
  
# Use 2 out of 3 rule to determine which category each belongs to
RbSlpCat.1 <- RbSlp %>% group_by(SiteID) %>%
  summarize(RbStr.n = n_distinct(RbStrID),
            Str.n = n_distinct(StrID)) %>% ungroup() %>%
  left_join(RbSlp, by = "SiteID") %>% ungroup() %>%
  group_by(SiteID, IGBP, RbStrID) %>%
  summarize(
    n = n(),
    RbStr.n = unique(RbStr.n),
    Str.n = unique(Str.n),
    StrID = unique(StrID),
    Slp_nrm = mean(Slp_nrm),
    Rbrnk = min(Rbrnk)
  ) %>% ungroup()

RbSlpCat <- RbSlpCat.1 %>%
  group_by(SiteID,IGBP) %>% nest() %>% 
  ungroup() %>%
  mutate(RbStr = map(data,  ~ grbfxn(.x)))  %>%
  unnest(RbStr) %>% select(-data,-id) %>%
  mutate(StrID = as.factor(StrID),
    RbID =     as.factor(
      case_when(Rbrnk == 3 ~ "Both",
                Rbrnk == 2 ~ "Prf",
                Rbrnk == 1 ~ "Slp",
                Rbrnk == 0 ~ "",
                TRUE ~ "")
    ),
    RbStrID = fct_recode(
      fct_relevel(
        fct_cross(RbID, StrID, sep = " "),
        c(
          "Both Dry Stress",
          "Prf Dry Stress",
          "Slp Dry Stress",
          " Dry Stress",
          "Both No Stress",
          "Prf No Stress",
          "Slp No Stress",
          " No Stress",
          " Unsure",
          " Wet Stress",
          "Slp Wet Stress",
          "Prf Wet Stress",
          "Both Wet Stress"
        )
      ),
      "Dry Stress" = " Dry Stress",
      "Wet Stress" = " Wet Stress",
      "No Stress" = " No Stress",
      "Unsure" = " Unsure"
    )
  )

# Commented out so the given file is not overwritten
#save(RbSlpCat,file = "./03-outputs/02-robust-framework/02-robust-summary/AFM_RC_0_7.RData")
ggplot(data = RbSlpCat,aes(x = StrID,fill = RbStrID)) + geom_bar()

#===============================================================================
# 2. Select and standardize performance metrics
#===============================================================================

perfVars = c(
  "LCE_mn",
  "AfT1_mn",
  "AfT2_mn",
  "AfT3_mn",
  "RobPrf.jnt",
  "LCE_Rob",
  "AfT1_Rob",
  "AfT2_Rob",
  "AfT3_Rob"
)

grbfxn.prf <- function(x,prfVars){
  
  # Specify number of unique stress responses
  n1 <- unique(x$RbStr.n) 
  n2 <- unique(x$Str.n) 
  
  # If two or more Robust Stress categories match
  x1 <- x %>% summarize(id = which.max(n),
                        RbStrID = RbStrID[id],
                        StrID = StrID[id],
                        Rbrnk = Rbrnk[id],
                        across(all_of(perfVars),~.x[id]))
  
  # If only two or more Stress categories match, grab minimum robustness category
  x2 <- x %>% group_by(StrID) %>% 
    summarize(n = n(),
              StrID = unique(StrID),
              Rbrnktmp = min(Rbrnk),
              RbStrID = RbStrID[which.min(Rbrnk)],
              across(all_of(perfVars),~.x[which.min(Rbrnk)])) %>% 
    ungroup() %>% 
    summarize(id = which.max(n),
              RbStrID = RbStrID[id],
              StrID = StrID[id],
              Rbrnk = Rbrnktmp[id],
              across(all_of(perfVars),~.x[id]))
  
  # If no matches in stress categories, we are unsure
  x3 <- x %>% summarize(id = 0,
                        RbStrID = "Unsure",
                        StrID = "Unsure",
                        across(all_of(perfVars),~.x[which.min(Rbrnk)]),
                        Rbrnk = 0)
  if(n1 < 3) {
    out <- x1
  }else if(n1 == 3 & n2 < 3) {
    out <- x2
  } else{
    out <- x3
  }
  
  return(out)
}

htPerf.1 <-  Slp.Sel %>% filter(!(SlpID == 3 &
                                  Slp_nrm < 0)) %>% 
  select(SiteID,IGBP,pvID,all_of(perfVars)) %>% 
  left_join(RbSlp %>% select(SiteID,pvID,RbStrID,StrID,Rbrnk),by = c("SiteID","pvID"))

htPerf.2 <- htPerf.1 %>% 
  group_by(SiteID) %>%
  summarize(RbStr.n = n_distinct(RbStrID),
            Str.n = n_distinct(StrID)) %>% 
  ungroup() %>%
  left_join(htPerf.1, by = "SiteID") %>% 
  ungroup() %>%
  group_by(SiteID, IGBP, RbStrID) %>%
  summarize(
    n = n(),
    RbStr.n = unique(RbStr.n),
    Str.n = unique(Str.n),
    StrID = unique(StrID),
    Rbrnk = min(Rbrnk),
    across(all_of(perfVars),~mean(.,na.rm = TRUE))
  ) %>% ungroup()

htPerf <- htPerf.2 %>% 
  group_by(SiteID,IGBP) %>% 
  nest() %>% 
  ungroup() %>%
  mutate(RbStr = map(data,  ~ grbfxn.prf(.x,perfVars)))  %>%
  unnest(RbStr) %>% 
  select(-data,-id) %>%
  mutate(
    across(contains("_sd"),~1 - (.))) %>%
  pivot_longer(cols = all_of(perfVars),
               names_to = "mtcID",
               values_to = "MtcVal") %>%
  mutate(mtcID = fct_relevel(mtcID, perfVars)) %>%
  left_join(th.Prf,by = "mtcID")



#===============================================================================
# 3. Create heat maps
#===============================================================================

# Unique IGBP
IGBP = unique(RbSlp$IGBP)

for (ii in IGBP){
# Select IGBP for analysis
IGBPsel = ii

# Plot panel sizes inches
boxdim = 0.117


# A. Overall Slope Robustness Score
#===================================
pltDat <- RbSlpCat %>% filter(IGBP == IGBPsel) %>%
  arrange(desc(RbStrID),Slp_nrm)
SlpOrd = as.character(pltDat$SiteID)
pltDat$SiteID = fct_relevel(pltDat$SiteID,SlpOrd)
pltDat$dummyID = 1
g1 <-  pltDat %>%
  htdisHelp("dummyID",
            "SiteID",
            "StrID",
            "Slope Score",
            "",
            "Site ID") +
  scale_fill_manual(
    'Stress Type',
    values = c("Dry Stress" = "#d7191c",
               "Unsure" = "#f0f0f0",
               "No Stress" = "#ffffbf",
               "Wet Stress" = "#2c7bb6"),
    labels = c("Dry Stress" = "Dry",
               "Unsure" = "Unsure",
               "No Stress" = "Minor",
               "Wet Stress" = "Wet"),
    drop = FALSE
  ) +
  force_panelsizes(rows = unit(length(unique(pltDat$SiteID)) * boxdim, "in"),
                   cols = unit(boxdim, "in")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "left"
  ) 
g1 <- g1 + geom_point(aes(shape = RbID,size = RbID),color = "black",stroke = stk) + 
  scale_shape_manual(name = "Robust Type",values = c(Both = 10, Prf = 3, Slp = 1),
                     labels = c(Both = "Full",Prf = "Perf Only",Slp = "Slope Only"),drop = FALSE) + 
  scale_size_manual(name = "Robust Type",values = c(Both = 5*mlt, Prf = 3*mlt, Slp = 5*mlt),
                    labels = c(Both = "Full",Prf = "Perf Only",Slp = "Slope Only"),drop = FALSE) + 
  guides(fill = guide_legend(override.aes = list(shape = NA)))


# B1. Normalized slope plot
#=========================
pltDat <-
  RbSlp %>% filter(IGBP == IGBPsel & !(SlpID == 3 &
                                        Slp_nrm < 0)) %>%
  mutate(SiteID = fct_relevel(SiteID,SlpOrd))
g2 <-  pltDat %>%
  htHelp("pvID",
         "SiteID",
         "Slp_nrm",
         "Slope",
         "Plant Variable",
         "") +
  scale_fill_gradientn(
    colors = rev(hcl.colors(20, "RdYlBu")),
    breaks = c(-2, 2),
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_x_discrete(name = "Plant Variable",
                   labels = c(
                     "G_1" = TeX(r"($G_1$)"),
                     "G_1_VPD_m" = TeX(r"($G_1/VPD^m$)"),
                     "G_v" = TeX(r"($G_v$)")
                   )) +
  force_panelsizes(rows = unit(length(unique(pltDat$SiteID)) * boxdim, "in"),
                   cols = unit(3 * boxdim, "in")) +
  theme(axis.text.y = element_blank())

g2 <- g2 + geom_point(aes(shape = RbSlp.ID,size = RbSlp.ID),color = "black",stroke = stk) + 
  scale_shape_manual(values = c("TRUE" = 1)) + 
  scale_size_manual(values = c("TRUE" = 5*mlt))+ guides(shape = "none",size = "none")



# C1. Performance metrics
#================================

selVars = c("RobPrf.jnt","LCE_Rob","AfT1_Rob","AfT2_Rob","AfT3_Rob")
pltDat <- htPerf %>% filter(IGBP == IGBPsel & mtcID %in% selVars) %>%
  mutate(SiteID = fct_relevel(SiteID,SlpOrd),PRFlag = ifelse(MtcVal >= RobPrf.th,"Rl",""))
# Save markers for perf plot
prfmrk <- pltDat %>% filter(!(mtcID %in% "RobPrf.jnt")) %>% select(PRFlag)

selVars = c("LCE_mn","AfT1_mn","AfT2_mn","AfT3_mn")
pltDat <- htPerf %>% filter(IGBP == IGBPsel & mtcID %in% selVars) %>%
  mutate(SiteID = fct_relevel(SiteID,SlpOrd),
         PRFlag = prfmrk$PRFlag) 

g3 <- pltDat %>%
  htHelp("mtcID",
         "SiteID",
         "MtcVal",
         "Performance",
         "Metric",
         "") +
  scale_fill_gradientn(
    colors = rev(hcl.colors(5, "Greens3")),
    breaks = c(0, 1),
    limits = c(0, 1),
    oob = scales::squish
  ) +
  scale_x_discrete(
    labels = c(
      "LCE_mn" = TeX(r"($LCE$)"),
      "AfT1_mn" = TeX(r"($A'_{f,T,1}$)"),
      "AfT2_mn" = TeX(r"($A'_{f,T,2}$)"),
      "AfT3_mn" = TeX(r"($A'_{f,T,3}$)")
    )
  ) +
  force_panelsizes(rows = unit(length(unique(pltDat$SiteID)) * boxdim, "in"),
                   cols = unit(4 * boxdim, "in")) +
  theme(axis.text.y = element_blank(),axis.title.x = element_text(margin = margin(t = 25)))
g3 <- g3 + geom_point(aes(shape = PRFlag,size = PRFlag),color = "black",stroke = stk) + 
  scale_shape_manual(values = c(Rl = 3)) + 
  scale_size_manual(values = c(Rl = 3*mlt))+ guides(shape = "none",size = "none")


# Combine all heatmaps
g.print <- (g1 + g2 + g3) + plot_layout(ncol = 3) + plot_annotation(
  title = as.character(IGBPsel)) 

# Save figure
ggsave2(paste0(svpath,"/fig-07-08-",IGBPsel,".pdf"),
        width = 8, height = 11, units = "in")

}
