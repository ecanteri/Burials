library(data.table)
library(ggplot2)
library(pbapply)
library(lme4)
library(lattice)
library(coefplot2)
library(cowplot)
library(patchwork)
setwd("~/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Projects/Graves")

## ---- ##
## DATA ##
## ---- ##

## Burials
graves <- fread("./Data/BurialCulture.csv", na.strings = "")
graves$DepositionType <- as.factor(graves$DepositionType)
graves$BodyPositioning <- as.factor(graves$BodyPositioning)
graves$BurialSide <- as.factor(graves$BurialSide)

# Remove outlier
graves <- graves[!is.na(YearBP)]
graves <- graves[YearBP < 15000]
summary(graves$YearBP)

# 10 most frequent cultures (for burial orientation data)
culture_keep <- readRDS("./Data/culture_keep.RDS")

## Final dataset
# subset data to keep only individuals belonging to those 10 cultures and that have ancestry information
DT <- graves[Culture %in% culture_keep][, .(Num, Longitude, Latitude, YearBP,
                                            Left = BurialSide_left,
                                            Right = BurialSide_right,
                                            Back = BurialSide_back,
                                            mobility = mobility_MDS_250y_retrospective,
                                            EHG = ANCE2_v2,
                                            WHG = ANCE4_v2,
                                            LVN = ANCE6_v2,
                                            CHG = ANCE8_v2,
                                            Culture)]
DT <- DT[!is.na(WHG)]
DT$Culture <- factor(DT$Culture, 
                     labels = c("Linearbandkeramik",
                                "Lažňany",
                                "Sopot",
                                "Bell Beaker",
                                "Únětice",
                                "Corded Ware",
                                "Baden",
                                "Dnieper-Donets",
                                "Trichterbecher",
                                "Baltic Meso-Neolithic"),
                     levels = unique(DT$Culture))

DT <- DT[!is.na(Left)]
DT <- DT[!is.na(mobility)]
DT$YearBin <- plyr::round_any(DT$YearBP, 250)

## Define covariates and scale
namescov <- c("Longitude", "Latitude", "YearBP", "mobility", "EHG", "WHG", "LVN", "CHG")
newnames <- paste0(namescov, "_scaled")
covs <- DT[, ..namescov]
covs <- scale(covs)
colnames(covs) <- newnames
DT <- cbind(DT, covs)

## Exploring variance
# Boxplots of the log of right side (response variable) versus LVN and Culture.
dat.tf <- within(DT,
                 {
                   # culture x lvn
                   gna <- interaction(Culture,LVN_scaled)
                   gna <- reorder(gna, Back, mean)
                   # culture x lvn x latitude x longitude
                   pna <- interaction(Culture,LVN_scaled,Longitude_scaled, Latitude_scaled)
                   pna <- reorder(pna, Back, mean)
                 })

ggplot(data = dat.tf, aes(factor(x = gna), y = log(Back + 1))) +
  geom_boxplot(colour = "skyblue2", outlier.shape = 21,
               outlier.colour = "skyblue2") +
  ylab("log (Back)\n") + # \n creates a space after the title
  xlab("\nCulture x LVN") + # space before the title
  theme_bw() + theme(axis.text.x = element_blank()) +
  stat_summary(fun = mean, geom = "point", colour = "red")

ggplot(data = dat.tf, aes(factor(x = pna), y = log(Back + 1))) +
  geom_boxplot(colour = "skyblue2", outlier.shape = 21,
               outlier.colour = "skyblue2") +
  ylab("log (Back)\n") + # \n creates a space after the title
  xlab("\nCulture x LVN x Longitude x Latitude") + # space before the title
  theme_bw() + theme(axis.text.x = element_blank()) +
  stat_summary(fun = mean, geom = "point", colour = "red")


## Modelling

## Easiest example
m0 <- glm(Back ~ 
            LVN_scaled +
            EHG_scaled +
            CHG_scaled +
            mobility_scaled +
            Longitude_scaled*Latitude_scaled +
            YearBP_scaled,
          data = DT,
          family = binomial)
summary(m0)

## Adding random effects
m1 <- glmer(Back ~ 
              LVN_scaled +
              EHG_scaled +
              CHG_scaled +
              mobility_scaled +
              Longitude_scaled*Latitude_scaled +
              YearBP_scaled +
              (1 | Culture), # +
              # (1 | Longitude_scaled:Latitude_scaled) +
              # (1 | YearBP_scaled),
            data = DT,
            family = binomial)
summary(m1)

# Variance terms
coefplot2(m1, ptype = "vcov", intercept = TRUE, main = "Random effect variance")
# Fixed effects
coefplot2(m1, intercept = TRUE, main = "Fixed effect coefficient")

# dotplot code
pp <- list(layout.widths = list(left.padding = 0, right.padding = 0),
           layout.heights = list(top.padding = 0, bottom.padding = 0))
r <- ranef(m1, condVar = TRUE)
d <- dotplot(r, par.settings = pp)
d

## Testing other models
m2 <- update(m1, . ~ . - LVN_scaled) 
m3 <- update(m1, . ~ . - EHG_scaled)
m4 <- update(m1, . ~ . - CHG_scaled) 
m5 <- update(m1, . ~ . - mobility_scaled)
m6 <- update(m1, . ~ . - Longitude_scaled)
m7 <- update(m1, . ~ . - Latitude_scaled)
m8 <- update(m1, . ~ . - YearBP_scaled)
m9 <- update(m1, . ~ . - Longitude_scaled:Latitude_scaled)
m10 <- glm(Back ~ 
             LVN_scaled +
             EHG_scaled +
             CHG_scaled +
             mobility_scaled +
             Longitude_scaled*Latitude_scaled +
             YearBP_scaled,
           data = DT,
           family = binomial)
# m10 <- update(m1, . ~ . - (1|Culture))
# m11 <- update(m1, . ~ . - (1|Longitude_scaled:Latitude_scaled))
# m12 <- update(m1, . ~ . - (1|YearBP_scaled))
aic_tab  <- MuMIn::model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10) #, m11, m12)
aic_table <- (round(aic_tab[ , c("AICc", "delta", "df")], digits = 2))
aic_table$Model <- rownames(aic_table)
aic_table$Formula <- as.character(lapply(list(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10), formula))[gtools::mixedorder(aic_table$Model)]
aic_table$Round <- "Round 1"
aic_table$Model <- factor(aic_table$Model, levels = c("m0", "m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10")) #, "m11", "m12"))

# With drop function
dd_LRT <- drop1(m1, test = "Chisq")
dd_LRT

# Plot
p1 <- ggplot(data = aic_table[aic_table$Model != "m0",]) +
  geom_line(aes(x = Model, y = AICc, group = 1)) +
  geom_point(aes(x = Model, y = AICc),
             size = 3) +
  geom_hline(data = aic_table[aic_table$Model == "m1",], 
             aes(yintercept = AICc),
             linetype = 2,
             alpha = 0.5) +
  facet_grid(cols = vars(Round)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))

## Based on the AICc scores and the drop test, it seems that Culture and ancestry do not have an effect, given that removing them improves the fit.
## Instead, we see an effect of Longitude, Latitude (less extent) and their interaction.

## Build another model
# mn1 <- glmer(Back ~
#                LVN_scaled +
#                Longitude_scaled*Latitude_scaled +
#                (1|Culture),
#              data = DT,
#              family = binomial)
# summary(mn1)

mn1 <- glm(Back ~
             EHG_scaled +
             CHG_scaled +
             Longitude_scaled +
             YearBP_scaled +
             Longitude_scaled:Latitude_scaled,
           data = DT,
           family = binomial)
summary(mn1)

# Variance terms
coefplot2(mn1, ptype = "vcov", intercept = TRUE, main = "Random effect variance")
# Fixed effects
coefplot2(mn1, intercept = TRUE, main = "Fixed effect coefficient")

# dotplot code
pp <- list(layout.widths = list(left.padding = 0, right.padding = 0),
           layout.heights = list(top.padding = 0, bottom.padding = 0))
r2 <- ranef(mn1, condVar = TRUE)
d2 <- dotplot(r2, par.settings = pp)
d2

## Testing other models
mn2 <- update(mn1, . ~ . - EHG_scaled)
mn3 <- update(mn1, . ~ . - CHG_scaled) 
mn4 <- update(mn1, . ~ . - Longitude_scaled) 
mn5 <- update(mn1, . ~ . - YearBP_scaled) 
mn6 <- update(mn1, . ~ . - Longitude_scaled:Latitude_scaled)

# mn6 <- glm(Back ~
#              LVN_scaled +
#              Longitude_scaled*Latitude_scaled,
#            data = DT,
#            family = binomial)
aic_tab2  <- MuMIn::model.sel(mn1, mn2, mn3, mn4, mn5, mn6)
aic_table2 <- (round(aic_tab2[ , c("AICc", "delta", "df")], digits = 2))
aic_table2$Model <- rownames(aic_table2)
aic_table2$Formula <- as.character(lapply(list(mn1, mn2, mn3, mn4, mn5, mn6), formula))
aic_table2$Round <- "Round 2"
aic_table2$Model <- factor(aic_table2$Model, levels = c("mn1", "mn2", "mn3", "mn4", "mn5", "mn6"))

# With drop function
dd_LRT2 <- drop1(mn1, test = "Chisq")
dd_LRT2

# Plot
p2 <- ggplot(data = aic_table2) +
  geom_point(aes(x = Model, y = AICc),
             size = 3) +
  geom_line(aes(x = Model, y = AICc, group = 1)) +
  geom_hline(data = aic_table2[aic_table2$Model == "mn1",], 
             aes(yintercept = AICc),
             linetype = 2,
             alpha = 0.5) +
  facet_grid(cols = vars(Round)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))


## Final model
fm <- mn1
summary(fm)
# # Variance terms
# coefplot2(fm, ptype = "vcov", intercept = TRUE, main = "Random effect variance")
# # Fixed effects
# coefplot2(fm, intercept = TRUE, main = "Fixed effect coefficient")
# 
# # dotplot code
# pp <- list(layout.widths = list(left.padding = 0, right.padding = 0),
#            layout.heights = list(top.padding = 0, bottom.padding = 0))
# r3 <- ranef(fm, condVar = TRUE)
# d3 <- dotplot(r3, par.settings = pp)
# d3
# 
# # Random effects
# random.dt <- as.data.frame(r3)
# rtab <- coeftab(fm, ptype="ranef",ctype="quad", clevel=c(0.95))
# random.dt <- cbind(random.dt, rtab[, 3:4])

# p3 <- ggplot(data = random.dt) +
#   geom_linerange(aes(x = condval, y = grp, xmin = `2.5%`, xmax = `97.5%`)) +
#   geom_point(aes(x = condval, y = grp),
#              size = 3,
#              shape = 21,
#              fill = "dodgerblue") +
#   geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
#   # scale_x_continuous(limits = c(-2.25, 2.25)) +
#   facet_grid(cols = vars(grpvar)) +
#   labs(x = "Conditional mode (Intercept)") +
#   theme_bw() +
#   theme(axis.title.y = element_blank(),
#         strip.text = element_text(face = "bold", size = 14),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14))

# Fixed effects
# fixed.dt <- as.data.frame(fixef(fm, condVar = T))
fixed.dt <- coeftab(fm, ptype="fixef",ctype="quad", clevel=c(0.95))
# fixed.dt <- cbind(fixed.dt, ftab[, 3:4])
# colnames(fixed.dt)[1] <- "condval"
fixed.dt$Variable <- rownames(fixed.dt)
fixed.dt$Variable[1] <- stringr::str_extract(fixed.dt$Variable[1], "Intercept")
fixed.dt$Parameter <- "Fixed-effect"

p4 <- ggplot(data = fixed.dt) +
  geom_linerange(aes(x = Estimate, y = Variable, xmin = `2.5%`, xmax = `97.5%`)) +
  geom_point(aes(x = Estimate, y = Variable),
             size = 3,
             shape = 21,
             fill = "dodgerblue") +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  scale_x_continuous(limits = c(-6, 6)) +
  facet_grid(cols = vars(Parameter), scales = "free") +
  labs(x = "Estimate") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

(p1 + p2 + plot_layout(axis_titles = "collect"))/(p4 + p3) + plot_annotation(tag_levels = "A", title = "Back burial side")
