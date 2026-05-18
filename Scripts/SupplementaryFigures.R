library(BIADconnect)
library(data.table)
library(ggplot2)
library(ggpattern)
library(paletteer)
library(ggrepel)
library(segmented)
library(terra)
library(sf)
library(spatstat)
library(inlabru)
library(INLA)
library(fmesher)
library(pbapply)
library(mobest)
library(knitr)
library(factoextra)
library(FactoMineR)
library(UpSetR)
library(webshot)

setwd("~/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Projects/Graves")
conn <- init.conn()
proj <- "+proj=aea +lon_0=44.296875 +lat_1=43.7864128 +lat_2=69.9512657 +lat_0=56.8688392 +datum=WGS84 +units=km +no_defs"

## ---- ##
## DATA ##
## ---- ##

## Coastlines and continents poly
conts <- rnaturalearth::ne_download(scale = 10, type = "coastline", "physical", returnclass = "sf")
land <- rnaturalearth::ne_download(scale = 10, type = "land", "physical", returnclass = "sf")
grat <- rnaturalearth::ne_download(scale = 10, type = "graticules_10", "physical", returnclass = "sf")

## Graves
graves <- fread("./Data/Burial.csv", na.strings = "")
graves$DepositionType <- as.factor(graves$DepositionType)
graves$BodyPositioning <- as.factor(graves$BodyPositioning)
graves$BurialSide <- as.factor(graves$BurialSide)

# Sites table (BIAD)
sites <- setDT(query.database("SELECT * FROM `Sites`", conn))

# Culture and Period
culture <- fread("./Data/burial rites_culture_period.csv", na.strings = "\\N")
culture <- culture[!duplicated(IndividualID)]

## Add info to graves table
graves <- merge(graves, sites[, .(SiteID, Country)], by = "SiteID", all.x = T)
rm(sites)
graves <- merge(graves, culture, by = "IndividualID", all.x = T)

## Remove outlier
graves <- graves[!is.na(YearBP)]
graves <- graves[YearBP < 15000]
summary(graves$YearBP)

## ---------------------------
## Figure S1 ##

# Check proportion of NAs
n <- names(which(sapply(graves, is.factor)))
cols <- unlist(sapply(n, function(n) grep(paste0(n, "_"), colnames(graves))))
colids <- colnames(graves)[26:53]

p2 <- melt(graves[, ..colids])
p2[, Type := sapply(strsplit(as.character(p2$variable), "_"), "[", 1)]
p2[, variable := sapply(strsplit(as.character(p2$variable), "_"), "[", 2)]
p2$Type <- factor(p2$Type, levels = unique(p2$Type), labels = c("Deposition Type",
                                                                "Body Positioning",
                                                                "Burial Side" ))
p2$variable <- factor(p2$variable, levels = unique(p2$variable))

pS1 <- ggplot(data = p2) + 
  geom_bar(aes(x = variable,
               fill = as.factor(value)),
           stat = "count") +
  scale_fill_manual("", 
                    values = paletteer::paletteer_d("beyonce::X6", 4)[3:4],
                    na.value = as.character(paletteer::paletteer_d("beyonce::X6", 4)[2]),
                    labels = c("No", "Yes", "NA")) +
  facet_wrap(~Type, scales = "free_x", space = "free_x") +
  ylab("Number of individuals") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

## ---------------------------
## Figure S2 ##

coll <- graves[, .(EHG = ANCE2_v2, WHG = ANCE4_v2,
                   NEOL = ANCE6_v2, CHG = ANCE8_v2)]

pS2.a <- ggplot(data = coll) +
  geom_point(aes(x = WHG, y = EHG)) +
  theme_light() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        aspect.ratio = 1)

pS2.b <- ggplot(data = coll) +
  geom_point(aes(x = WHG, y = CHG)) +
  theme_light() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        aspect.ratio = 1)

pS2.c <- ggplot(data = coll) +
  geom_point(aes(x = WHG, y = NEOL)) +
  theme_light() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        aspect.ratio = 1)

png("./Figures/Supplementary/Figure_S2.png", width = 16, height = 5, res = 330, units = 'in')
egg::ggarrange(plots = list(pS2.a, pS2.b, pS2.c), ncol = 3, labels = c("A", "B", "C"),
               label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
dev.off()

## ---------------------------
## Figure S3 ##

performance.dt <- fread("./Results/ModelPerformance.csv")
performance.dt$Variable <- factor(performance.dt$Variable,
                                  levels = c("Left", "Right", "Back"))

pS3 <- ggplot(data = performance.dt) +
  geom_violin(aes(y = Performance, x = Model, fill = Variable),
              width = 0.8) +
  geom_boxplot(aes(y = Performance, x = Model, fill = Variable),
               width = 0.25,
               linewidth = 1,
               color = "grey20"
  ) +
  scale_fill_manual(values = c("#aedbf1", "#1a7bbb", "#b5e48c")) +
  facet_wrap(~Variable, ncol = 3, scales = "fixed") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18, face = "bold")
  )

png("./Figures/Supplementary/Figure_S3.png", width = 14, height = 10, res = 330, units = "in")
pS3
dev.off()

## ---------------------------
## Figure S4 ##

# 10 most frequent cultures (for burial orientation data)
culture_keep <- readRDS("./Data/culture_keep.RDS")

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
                                "Sopot",
                                "Lažňany",
                                "Baltic Meso-Neolithic",
                                "Dnieper-Donets",
                                "Trichterbecher",
                                "Baden",
                                "Corded Ware",
                                "Bell Beaker",
                                "Únětice"),
                     levels = c("Linearbandkeramik",
                                "Sopot",
                                "Lažňany",
                                "Baltic Mesolithic-Neolithic",
                                "Dnieper-Donets",
                                "Trichterbecher",
                                "Baden",
                                "Corded Ware/Single Grave/Battle Axe/Fatyanovo",
                                "Bell Beaker complex",
                                "Únětice"))

DT <- DT[!is.na(Left)]
DT <- DT[!is.na(mobility)]
DT$YearBin <- plyr::round_any(DT$YearBP, 250)

## Spatial points
pnts <- project(
  vect(
    DT,
    geom = c("Longitude", "Latitude"),
    crs = crs(rast())
  ),
  proj
)

pnts <- intersect(pnts, poly)
pnts <- st_as_sf(pnts)

pnts.dt <- setDT(as.data.frame(pnts))
pnts.dt <- melt(pnts.dt, measure.vars = c("Left", "Right", "Back"))
pnts.dt <- pnts.dt[!is.na(value)]

p4.1 <- ggplot(data = pnts.dt) +
  # geom_sf(data = land_prj, inherit.aes = T, fill = "grey90") +
  geom_sf(data = countries_prj, fill = "grey90") +
  geom_sf(
    data = pnts.dt$geometry,
    inherit.aes = T,
    shape = 21,
    fill = pals::ocean.thermal(6)[4],
    size = 2.5,
    alpha = 0.8
  ) +
  coord_sf(
    xlim = st_bbox(pnts.dt$geometry)[c(1,3)],
    ylim = st_bbox(pnts.dt$geometry)[c(2,4)]
  ) +
  labs(tag = "A") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        # plot.margin = margin(t = 0, b = 0, l = 0, r = 0),
        # aspect.ratio = 1,
        plot.tag = element_text(face = 'bold', size = 18)
  )

p4.2 <- ggplot(data = pnts.dt) +
  geom_density(
    data = pnts.dt[value == 1], aes(x = YearBP),
    fill = "grey80",
    color = "transparent"
  ) +
  geom_point(aes(
    x = YearBP, y = 0,
    fill = as.factor(value)
  ), shape = 21, size = 1.5, alpha = 0.8, color = "grey30") +
  scale_fill_manual(values = pals::ocean.thermal(6)[c(2,5)]) +
  scale_x_reverse() +
  facet_grid(rows = vars(variable), switch = "both") +
  xlab("Years (cal BP)") +
  labs(tag = "B") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 14),
    legend.position = "none",
    aspect.ratio = 0.5,
    plot.tag = element_text(face = 'bold', size = 18)
  )

## Results
cols <- paletteer::paletteer_d("ButterflyColors::anteos_menippe", 3)
cols[4] <- "#5D3E99"

## Full model culture
results_df <- fread("./Results/Culture/FullModel_results.csv")
results_df$Predictors <- gsub("LVN", "NEOL", results_df$Predictors)
results_df$Predictors <- factor(results_df$Predictors,
                                levels = rev(c("Intercept", "EHG", "CHG", "NEOL", "Mobility", 
                                               levels(DT$Culture))))
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = c("Intercept", "Ancestry", "Mobility", "Culture"))

# Plot
f4a <- ggplot(data = results_df) +
  geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
                width = 0.2, linewidth = 1) +
  geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
  geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  facet_grid(cols = vars(Variable), rows = vars(Effect),
             scales = "free") +
  theme_bw() +
  labs(tag = "C") +
  # ggtitle("Full model") +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        # aspect.ratio = 0.9,
        # plot.margin = margin(t = 10, b = 0, l = 5, r = 5),
        plot.tag = element_text(face = 'bold', size = 18))

# Model fit
model_fit <- fread("./Results/Culture/Models_fit.csv")

## Plot model fit
model_fit[, deltaWAIC := waic - min(waic), by = "Variable"]
model_fit[, colorWAIC := ifelse(deltaWAIC == 0, pals::ocean.thermal(6)[4], "black")]
model_fit[, deltaDIC := dic - min(dic), by = "Variable"]
model_fit[, colorDIC := ifelse(deltaDIC == 0, pals::ocean.thermal(6)[4], "black")]
model_fit$Variable <- factor(model_fit$Variable, levels = unique(model_fit$Variable))
model_fit$Model <- factor(model_fit$Model, levels = unique(model_fit$Model))

f4b <- ggplot(data = model_fit) +
  geom_line(aes(x = Model, y = deltaWAIC, group = 1)) +
  geom_hline(
    data = model_fit[Model == "baseline"], aes(yintercept = deltaWAIC),
    linetype = 2, color = "black"
  ) +
  geom_point(aes(x = Model, y = deltaWAIC, fill = colorWAIC), shape = 21, size = 2.5) +
  facet_grid(rows = vars(Variable), scales = "free_y") +
  scale_fill_identity() +
  ylab("Delta WAIC") +
  labs(tag = "D") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    # aspect.ratio = 0.3,
    plot.tag = element_text(face = 'bold', size = 18)
  )

p4a <- egg::ggarrange(plots = list(p4.1, p4.2), ncol = 1, nrow = 2, 
                      labels = c("A", "B"), 
                      label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
p4b <- egg::ggarrange(plots = list(f4a, f4b), ncol = 1, nrow = 2, heights = c(1,0.8),
                      labels = c("C", "D"), 
                      label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))

png("./Figures/Supplementary/Figure_S4.png", res = 330, units = 'in', height = 10, width = 18)
plot_grid(p4a, p4b, rel_widths = c(0.8, 1))
dev.off()

## ---------------------------
## FIGURE S5

## Data
bp <- grep("BodyPositioning_", colnames(graves))
bs <- grep("BurialSide_", colnames(graves))
mat <- as.data.frame(cbind(graves[, ..bp], graves[, ..bs]))
mat <- as.matrix(mat)
dim(mat)

## Multiple Correspondence Analysis
mca <- MCA(graves[, .(BurialSide, BodyPositioning)], graph = F)
summary(mca)

## Extract results for individuals
ind <- get_mca_ind(mca)

## Extract results for variables
var <- get_mca_var(mca)

## Eigenvalues
eig.val <- get_eigenvalue(mca)

# Individuals
ind_dims <- as.data.frame(ind$coord)
ind_dims$Individual <- 1:nrow(ind_dims)
ind_dims$Country <- graves$Country
ind_dims$Period <- graves$Period
ind_dims$Culture <- graves$Culture
ind_dims$Sex <- graves$Sex
ind_dims$EHG <- graves$ANCE2_v2
ind_dims$WHG <- graves$ANCE4_v2
ind_dims$LVN <- graves$ANCE6_v2
ind_dims$CHG <- graves$ANCE8_v2
ind_dims$mobility <- graves$mobility_MDS_250y_retrospective
head(ind_dims)
ind_dims <- setDT(ind_dims)

# Variables
var_dims <- as.data.frame(var$coord)
var_dims$Variable <- rownames(var_dims)
var_dims <- var_dims[-which(var_dims$Variable %in% c("BurialSide.NA", "BodyPositioning.NA")),]

# Plot
pS5 <- fviz_screeplot(mca)

png("./Figures/Supplementary/Figure_S5.png", res = 330, units = 'in', height = 8, width = 10)
pS5
dev.off()

## ---------------------------
## FIGURE S6
upset_df <- graves[, .(IndividualID, BodyPositioning_crouched, BodyPositioning_extended, `BodyPositioning_extended, rhomboid`, BodyPositioning_heap, `BodyPositioning_inside a vessel`, BodyPositioning_irregular, BodyPositioning_prone, BodyPositioning_scattered, BodyPositioning_sitting, BurialSide_back, `BurialSide_back/left`, `BurialSide_back/right`, BurialSide_bottom, BurialSide_front, `BurialSide_front/left`, `BurialSide_front/right`, BurialSide_left, `BurialSide_left/back`, `BurialSide_left/front`, BurialSide_right, `BurialSide_right/back`, `BurialSide_right/front`)]
set_size <- data.frame(Set = colnames(upset_df)[-1],
                       Size = colSums(upset_df[,-1], na.rm = T))
sets <- c("BodyPositioning_crouched", 
          "BodyPositioning_extended", 
          "BodyPositioning_heap", 
          "BodyPositioning_inside a vessel", 
          "BodyPositioning_prone", 
          "BodyPositioning_scattered", 
          "BodyPositioning_sitting",
          "BurialSide_bottom", 
          "BurialSide_front", 
          "BurialSide_left", 
          "BurialSide_right", 
          "BurialSide_back")
set_size <- set_size[which(set_size$Set %in% sets),]
set_size$Type <- sapply(strsplit(set_size$Set, "_"), "[", 1)
set_size$Type <- factor(set_size$Type, labels = c("Body Positioning", "Burial Side"))
set_size$Set <- sapply(strsplit(set_size$Set, "_"), "[", 2)

# Plot
pS6 <- ggplot(data = set_size) +
  geom_bar(aes(x = reorder(Set, Size, decreasing = T), y = Size, fill = Type), 
           stat = "identity") +
  scale_fill_manual(values = paletteer_d("beyonce::X6", 6)[3:4]) +
  theme_classic() +
  ylab("Set Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        plot.margin = margin(t = 5, r = 10, l = 20, b = 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

png("./Figures/Supplementary/Figure_S6.png", res = 330, units = 'in', height = 8, width = 14)
pS6
dev.off()

## ---------------------------
## FIGURE S7

## Biplot by PERIOD
dim_period <- ind_dims[, .N, by = c("Period", "Dim 1", "Dim 2")]
dim_period <- dim_period[!is.na(Period)]
dim_period$Period <- factor(dim_period$Period, levels = c("Paleolithic",
                                                          "Mesolithic",
                                                          "Mesolithic/Neolithic",
                                                          "Neolithic",
                                                          "Neolithic/Bronze Age",
                                                          "Eneolithic",
                                                          "Eneolithic/Bronze Age",
                                                          "Bronze Age",
                                                          "Bronze Age/Iron Age",
                                                          "Iron Age",
                                                          "Medieval"))

pS7.a <- ggplot() +
  geom_point(data = dim_period,
             aes(x = `Dim 1`, y = `Dim 2`,
                 fill = Period,
                 size = N),
             position = position_jitter(0.2, 0.2),
             shape = 21,
             alpha = 0.8) +
  scale_fill_paletteer_d("MetBrewer::Hiroshige",
                         guide = guide_legend(override.aes = list(size = 5))) +
  scale_size_continuous("Count", breaks = c(1,50,100,250,500),
                        range = c(3, 20)) +
  geom_point(data = var_dims,
             aes(x = `Dim 1`, y = `Dim 2`),
             shape = 23,
             fill = 'white',
             alpha = 0.5,
             size = 3,
             stroke = 1) +
  geom_label_repel(data = var_dims,
                   aes(x = `Dim 1`, y = `Dim 2`, label = Variable),
                   segment.color = 'black',
                   segment.alpha = 0.5,
                   fill = alpha("white", 0.8),
                   max.overlaps = 50,
                   nudge_x = 1,
                   nudge_y = 0.5,
                   direction = "both") +
  xlab(paste0("Dim 1 (", round(eig.val[1,2], 2), "%)")) +
  ylab(paste0("Dim 2 (", round(eig.val[2,2], 2), "%)")) +
  theme_bw() +
  theme(legend.position = "right",
        aspect.ratio = 1)

## Biplot by CULTURE
dim_culture <- ind_dims[!is.na(Culture)][, .N, by = c("Culture", "Dim 1", "Dim 2")]
culture_keep <- readRDS("/Users/msb290/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Projects/Graves/Data/culture_keep.RDS")
dim_culture <- dim_culture[Culture %in% culture_keep]
dim_culture$Culture <- factor(dim_culture$Culture, 
                              levels = c("Linearbandkeramik",
                                         "Sopot",
                                         "Lažňany",
                                         "Baltic Mesolithic-Neolithic",
                                         "Dnieper-Donets",
                                         "Trichterbecher",
                                         "Baden",
                                         "Corded Ware/Single Grave/Battle Axe/Fatyanovo",
                                         "Bell Beaker complex",
                                         "Únětice"),
                              labels = c("Linearbandkeramik",
                                         "Sopot",
                                         "Lažňany",
                                         "Baltic Meso-Neolithic",
                                         "Dnieper-Donets",
                                         "Trichterbecher",
                                         "Baden",
                                         "Corded Ware",
                                         "Bell Beaker",
                                         "Únětice"))

pS7.b <- ggplot() +
  geom_point(data = dim_culture,
             aes(x = `Dim 1`, y = `Dim 2`,
                 fill = Culture,
                 size = N),
             position = position_jitter(0.2, 0.2),
             shape = 21,
             alpha = 0.8) +
  scale_fill_manual(values = rev(pals::ocean.thermal(12))[2:11],
                    guide = guide_legend(override.aes = list(size = 5))) +
  scale_size_continuous("Count", breaks = c(1,25,50,100,150, 200),
                        range = c(3, 20)) +
  geom_point(data = var_dims,
             aes(x = `Dim 1`, y = `Dim 2`),
             shape = 23,
             fill = 'white',
             alpha = 0.5,
             size = 3,
             stroke = 1) +
  geom_label_repel(data = var_dims,
                   aes(x = `Dim 1`, y = `Dim 2`, label = Variable),
                   segment.color = 'black',
                   segment.alpha = 0.5,
                   fill = alpha("white", 0.8),
                   max.overlaps = 50,
                   nudge_x = 1,
                   nudge_y = 0.5,
                   direction = "both") +
  xlab(paste0("Dim 1 (", round(eig.val[1,2], 2), "%)")) +
  ylab(paste0("Dim 2 (", round(eig.val[2,2], 2), "%)")) +
  theme_bw() +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        aspect.ratio = 1)

## Biplot by sex
dim_sex <- ind_dims[, .N, by = c("Sex", "Dim 1", "Dim 2")]

pS7.c <- ggplot() +
  geom_point(data = dim_sex,
             aes(x = `Dim 1`, y = `Dim 2`,
                 fill = Sex,
                 size = N),
             position = position_jitter(0.2, 0.2),
             shape = 21,
             alpha = 0.8) +
  scale_fill_gradientn("Sex", colors = paletteer_c("grDevices::Tropic", 100)) +
  scale_size_continuous("Count", breaks = c(1,50,100,300,600, 900),
                        range = c(3, 20)) +
  geom_point(data = var_dims,
             aes(x = `Dim 1`, y = `Dim 2`),
             shape = 23,
             fill = 'white',
             alpha = 0.5,
             size = 3,
             stroke = 1) +
  geom_label_repel(data = var_dims,
                   aes(x = `Dim 1`, y = `Dim 2`, label = Variable),
                   segment.color = 'black',
                   segment.alpha = 0.5,
                   fill = alpha("white", 0.8),
                   max.overlaps = 50,
                   nudge_x = 1,
                   nudge_y = 0.5,
                   direction = "both") +
  xlab(paste0("Dim 1 (", round(eig.val[1,2], 2), "%)")) +
  ylab(paste0("Dim 2 (", round(eig.val[2,2], 2), "%)")) +
  theme_bw() +
  theme(legend.position = "right", 
        aspect.ratio = 1)

png("./Figures/Supplementary/Figure_S7.png", width = 16, height = 12, res = 330, units = 'in')
egg::ggarrange(plots = list(pS7.a, pS7.b, pS7.c), ncol = 2, labels = c("A", "B", "C"),
               label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
dev.off()

## ---------------------------
## FIGURE S8

## Biplot by WHG
dim_whg <- ind_dims[, mean(WHG, na.rm = T), by = c("Dim 1", "Dim 2")]
dim_whg <- dim_whg[!is.nan(V1),]
names(dim_whg)[3] <- "WHG"

pS8.a <- ggplot() +
  geom_point(data = dim_whg,
             aes(x = `Dim 1`, y = `Dim 2`,
                 fill = WHG),
             position = position_jitter(0.2, 0.2),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("WHG", colors = pals::parula(100), limits = c(0,1)) +
  geom_point(data = var_dims,
             aes(x = `Dim 1`, y = `Dim 2`),
             shape = 23,
             fill = 'white',
             alpha = 0.5,
             size = 3,
             stroke = 1) +
  geom_label_repel(data = var_dims,
                   aes(x = `Dim 1`, y = `Dim 2`, label = Variable),
                   segment.color = 'black',
                   segment.alpha = 0.5,
                   fill = alpha("white", 0.8),
                   max.overlaps = 50,
                   nudge_x = 1,
                   nudge_y = 0.5,
                   direction = "both") +
  xlab(paste0("Dim 1 (", round(eig.val[1,2], 2), "%)")) +
  ylab(paste0("Dim 2 (", round(eig.val[2,2], 2), "%)")) +
  ggtitle("Western hunter-gatherers ancestry") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

## Biplot by EHG
dim_ehg <- ind_dims[, mean(EHG, na.rm = T), by = c("Dim 1", "Dim 2")]
dim_ehg <- dim_ehg[!is.nan(V1),]
names(dim_ehg)[3] <- "EHG"

pS8.b <- ggplot() +
  geom_point(data = dim_ehg,
             aes(x = `Dim 1`, y = `Dim 2`,
                 fill = EHG),
             position = position_jitter(0.2, 0.2),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("EHG", colors = pals::parula(100), limits = c(0,1)) +
  geom_point(data = var_dims,
             aes(x = `Dim 1`, y = `Dim 2`),
             shape = 23,
             fill = 'white',
             alpha = 0.5,
             size = 3,
             stroke = 1) +
  geom_label_repel(data = var_dims,
                   aes(x = `Dim 1`, y = `Dim 2`, label = Variable),
                   segment.color = 'black',
                   segment.alpha = 0.5,
                   fill = alpha("white", 0.8),
                   max.overlaps = 50,
                   nudge_x = 1,
                   nudge_y = 0.5,
                   direction = "both") +
  xlab(paste0("Dim 1 (", round(eig.val[1,2], 2), "%)")) +
  ylab(paste0("Dim 2 (", round(eig.val[2,2], 2), "%)")) +
  ggtitle("Eastern hunter-gatherers ancestry") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

## Biplot by CHG
dim_chg <- ind_dims[, mean(CHG, na.rm = T), by = c("Dim 1", "Dim 2")]
dim_chg <- dim_chg[!is.nan(V1)]
names(dim_chg)[3] <- "CHG"

pS8.c <- ggplot() +
  geom_point(data = dim_chg,
             aes(x = `Dim 1`, y = `Dim 2`,
                 fill = CHG),
             position = position_jitter(0.2, 0.2),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("CHG", colors = pals::parula(100), limits = c(0,1)) +
  geom_point(data = var_dims,
             aes(x = `Dim 1`, y = `Dim 2`),
             shape = 23,
             fill = 'white',
             alpha = 0.5,
             size = 3,
             stroke = 1) +
  geom_label_repel(data = var_dims,
                   aes(x = `Dim 1`, y = `Dim 2`, label = Variable),
                   segment.color = 'black',
                   segment.alpha = 0.5,
                   fill = alpha("white", 0.8),
                   max.overlaps = 50,
                   nudge_x = 1,
                   nudge_y = 0.5,
                   direction = "both") +
  xlab(paste0("Dim 1 (", round(eig.val[1,2], 2), "%)")) +
  ylab(paste0("Dim 2 (", round(eig.val[2,2], 2), "%)")) +
  ggtitle("Caucasus hunter-gatherers ancestry") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

## Biplot by LVN
dim_lvn <- ind_dims[, mean(LVN, na.rm = T), by = c("Dim 1", "Dim 2")]
dim_lvn <- dim_lvn[!is.nan(V1)]
names(dim_lvn)[3] <- "NEOL"

pS8.d <- ggplot() +
  geom_point(data = dim_lvn,
             aes(x = `Dim 1`, y = `Dim 2`,
                 fill = NEOL),
             position = position_jitter(0.2, 0.2),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("NEOL", colors = pals::parula(100), limits = c(0,1)) +
  geom_point(data = var_dims,
             aes(x = `Dim 1`, y = `Dim 2`),
             shape = 23,
             fill = 'white',
             alpha = 0.5,
             size = 3,
             stroke = 1) +
  geom_label_repel(data = var_dims,
                   aes(x = `Dim 1`, y = `Dim 2`, label = Variable),
                   segment.color = 'black',
                   segment.alpha = 0.5,
                   fill = alpha("white", 0.8),
                   max.overlaps = 50,
                   nudge_x = 1,
                   nudge_y = 0.5,
                   direction = "both") +
  xlab(paste0("Dim 1 (", round(eig.val[1,2], 2), "%)")) +
  ylab(paste0("Dim 2 (", round(eig.val[2,2], 2), "%)")) +
  ggtitle("Neolithic farmers ancestry") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

## Biplot by mobility
dim_mobility <- ind_dims[, mean(mobility, na.rm = T), by = c("Dim 1", "Dim 2")]
dim_mobility <- dim_mobility[!is.nan(V1)]
names(dim_mobility)[3] <- "mobility"

pS8.e <- ggplot() +
  geom_point(data = dim_mobility,
             aes(x = `Dim 1`, y = `Dim 2`,
                 fill = mobility),
             position = position_jitter(0.2, 0.2),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("Mobility", colors = pals::parula(100)) +
  geom_point(data = var_dims,
             aes(x = `Dim 1`, y = `Dim 2`),
             shape = 23,
             fill = 'white',
             alpha = 0.5,
             size = 3,
             stroke = 1) +
  geom_label_repel(data = var_dims,
                   aes(x = `Dim 1`, y = `Dim 2`, label = Variable),
                   segment.color = 'black',
                   segment.alpha = 0.5,
                   fill = alpha("white", 0.8),
                   max.overlaps = 50,
                   nudge_x = 1,
                   nudge_y = 0.5,
                   direction = "both") +
  xlab(paste0("Dim 1 (", round(eig.val[1,2], 2), "%)")) +
  ylab(paste0("Dim 2 (", round(eig.val[2,2], 2), "%)")) +
  ggtitle("Genetically-inferred mobility") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

png("./Figures/Supplementary/Figure_S8.png", width = 14, height = 18, res = 330, units = 'in')
egg::ggarrange(plots = list(pS8.a, pS8.b, pS8.c, pS8.d, pS8.e), ncol = 2, labels = c("A", "B", "C", "D", "E"),
               label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
dev.off()

## ---------------------------
## FIGURE S9

# Data
results_best <- fread("./Results/BestModels_results.csv")
results_best_b5k <- fread("./Results/BestModelsBefore5k_results.csv")
results_best_steppe <- fread("./Results/BestModelsAfter5k_results.csv")

# Full
results_best$Predictors <- factor(results_best$Predictors,
                                levels = rev(c("Intercept", "EHG", "CHG", "LVN", "Mobility")),
                                labels = rev(c("Intercept", "EHG", "CHG", "NEOL", "Mobility")))
results_best$Variable <- factor(results_best$Variable, levels = unique(results_best$Variable))
results_best$Type <- factor(results_best$Type, levels = c("Intercept", "Ancestry", "Mobility"))

# Before 5k
results_best_b5k$Predictors <- factor(results_best_b5k$Predictors,
                                  levels = rev(c("Intercept", "EHG", "CHG", "LVN", "Mobility")),
                                  labels = rev(c("Intercept", "EHG", "CHG", "NEOL", "Mobility")))
results_best_b5k$Variable <- factor(results_best_b5k$Variable, levels = unique(results_best_b5k$Variable))
results_best_b5k$Type <- factor(results_best_b5k$Type, levels = c("Intercept", "Ancestry", "Mobility"))

# After 5k
results_best_steppe$Predictors <- factor(results_best_steppe$Predictors,
                                  levels = rev(c("Intercept", "EHG", "CHG", "LVN", "Steppe ancestry", "Mobility")),
                                  labels = rev(c("Intercept", "EHG", "CHG", "NEOL", "Steppe", "Mobility")))
results_best_steppe$Variable <- factor(results_best_steppe$Variable, levels = unique(results_best_steppe$Variable))
results_best_steppe$Type <- factor(results_best_steppe$Type, levels = c("Intercept", "Ancestry", "Mobility"))

# Plot
cols <- paletteer::paletteer_d("ButterflyColors::anteos_menippe", 3)

pS9.a <- ggplot(data = results_best) +
  geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
                width = 0.2, linewidth = 1) +
  geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
  geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  facet_grid(cols = vars(Variable), scales = "fixed") +
  theme_bw() +
  ggtitle("Full period") +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        # aspect.ratio = 0.9,
        # plot.margin = margin(t = 10, b = 0, l = 5, r = 5),
        plot.tag = element_text(face = 'bold', size = 18))

pS9.b <- ggplot(data = results_best_b5k) +
  geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
                width = 0.2, linewidth = 1) +
  geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
  geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  facet_grid(cols = vars(Variable), scales = "fixed") +
  theme_bw() +
  ggtitle("Before 5k calBP") +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        # aspect.ratio = 0.9,
        # plot.margin = margin(t = 10, b = 0, l = 5, r = 5),
        plot.tag = element_text(face = 'bold', size = 18))

pS9.c <- ggplot(data = results_best_steppe) +
  geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
                width = 0.2, linewidth = 1) +
  geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
  geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  facet_grid(cols = vars(Variable), scales = "free_y") +
  theme_bw() +
  ggtitle("After 5k calBP") +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        # aspect.ratio = 0.9,
        # plot.margin = margin(t = 10, b = 0, l = 5, r = 5),
        plot.tag = element_text(face = 'bold', size = 18))

png("./Figures/Supplementary/Figure_S9.png", width = 12, height = 12, res = 330, units = 'in')
egg::ggarrange(plots = list(pS9.a, pS9.b, pS9.c), nrow = 3, labels = c("A", "B", "C"),
               label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
dev.off()

## ---------------------------
## FIGURE S10

# Data
results_full_b5k <- fread("./Results/FullModelBefore5k_results.csv")
model_fit_b5k <- fread("./Results/ModelsFit_BeforeSteppe.csv")
results_full_b5k$Predictors <- factor(results_full_b5k$Predictors,
                                      levels = rev(c("Intercept",
                                                 "EHG",
                                                 "CHG",
                                                 "LVN",
                                                 "Mobility")),
                                      labels = rev(c("Intercept",
                                                 "EHG",
                                                 "CHG",
                                                 "NEOL",
                                                 "Mobility")))
results_full_b5k$Variable <- factor(results_full_b5k$Variable, levels = c("Left", "Right", "Back"))
results_full_b5k$Type <- factor(results_full_b5k$Type, levels = c("Intercept", "Ancestry", "Mobility"))

model_fit_b5k$Model <- factor(model_fit_b5k$Model,
                              levels = c("baseline", "baseline + no WHG ancestry", "baseline + no WHG ancestry + mobility"),
                              labels = c("baseline", "baseline + ancestry", "baseline + ancestry + mobility"))

model_fit_b5k[, deltaWAIC := waic - min(waic), by = "Variable"]
model_fit_b5k[, colorWAIC := ifelse(deltaWAIC == 0, "red", "black")]
model_fit_b5k[, deltaDIC := dic - min(dic), by = "Variable"]
model_fit_b5k[, colorDIC := ifelse(deltaDIC == 0, "red", "black")]

model_fit_b5k$Variable <- factor(model_fit_b5k$Variable, levels = c("Left", "Right", "Back"))

# Plot
pS10.a <- ggplot(data = results_full_b5k) +
  geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
                width = 0.2, linewidth = 1) +
  geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
  geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  facet_grid(cols = vars(Variable), scales = "free_x") +
  theme_bw() +
  ggtitle("Before 5k calBP") +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        plot.margin = margin(b = 35))

pS10.b <- ggplot(data = model_fit_b5k) +
  geom_line(aes(x = Model, y = deltaWAIC, group = 1)) +
  geom_hline(
    data = model_fit_b5k[Model == "baseline"], aes(yintercept = deltaWAIC),
    linetype = 2, color = "black"
  ) +
  geom_point(aes(x = Model, y = deltaWAIC, fill = colorWAIC), shape = 21, size = 2.5) +
  facet_grid(rows = vars(Variable), scales = "free_y") +
  scale_fill_identity() +
  ylab("Delta WAIC") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(b = 35))


pS10.c <- ggplot(data = model_fit_b5k) +
  geom_line(aes(x = Model, y = deltaDIC, group = 1)) +
  geom_hline(
    data = model_fit_b5k[Model == "baseline"], aes(yintercept = deltaDIC),
    linetype = 2, color = "black"
  ) +
  geom_point(aes(x = Model, y = deltaDIC, fill = colorDIC), shape = 21, size = 2.5) +
  facet_grid(rows = vars(Variable), scales = "free_y") +
  scale_fill_identity() +
  ylab("Delta DIC") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12))

png("./Figures/Supplementary/Figure_S10.png", width = 12, height = 14, res = 330, units = 'in')
egg::ggarrange(plots = list(pS10.a, pS10.b, pS10.c), nrow = 3, labels = c("A", "B", "C"),
               label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
dev.off()

## ---------------------------
## FIGURE S11

# Data
results_full_steppe <- fread("./Results/FullModelSteppe_results.csv")
model_fit_steppe <- fread("./Results/ModelsFit_Steppe.csv")
results_full_steppe$Predictors <- factor(results_full_steppe$Predictors,
                                      levels = rev(c("Intercept",
                                                     "LVN",
                                                     "Steppe",
                                                     "Mobility")),
                                      labels = rev(c("Intercept",
                                                     "NEOL",
                                                     "Steppe",
                                                     "Mobility")))
results_full_steppe$Variable <- factor(results_full_steppe$Variable, levels = c("Left", "Right", "Back"))
results_full_steppe$Type <- factor(results_full_steppe$Type, levels = c("Intercept", "Ancestry", "Mobility"))

model_fit_steppe$Model <- factor(model_fit_steppe$Model,
                              levels = c("baseline", "baseline + mobility + steppe ancestry + LVN", "baseline + steppe ancestry + LVN", "baseline + steppe ancestry"),
                              labels = c("baseline", "baseline + mobility + steppe ancestry + NEOL", "baseline + steppe ancestry + NEOL", "baseline + steppe ancestry"))

model_fit_steppe[, deltaWAIC := waic - min(waic), by = "Variable"]
model_fit_steppe[, colorWAIC := ifelse(deltaWAIC == 0, "red", "black")]
model_fit_steppe[, deltaDIC := dic - min(dic), by = "Variable"]
model_fit_steppe[, colorDIC := ifelse(deltaDIC == 0, "red", "black")]

model_fit_steppe$Variable <- factor(model_fit_steppe$Variable, levels = c("Left", "Right", "Back"))

# Plot
pS11.a <- ggplot(data = results_full_steppe) +
  geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
                width = 0.2, linewidth = 1) +
  geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
  geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  facet_grid(cols = vars(Variable), scales = "free_x") +
  theme_bw() +
  ggtitle("After 5k calBP") +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        plot.margin = margin(b = 35))

pS11.b <- ggplot(data = model_fit_steppe) +
  geom_line(aes(x = Model, y = deltaWAIC, group = 1)) +
  geom_hline(
    data = model_fit_steppe[Model == "baseline"], aes(yintercept = deltaWAIC),
    linetype = 2, color = "black"
  ) +
  geom_point(aes(x = Model, y = deltaWAIC, fill = colorWAIC), shape = 21, size = 2.5) +
  facet_grid(rows = vars(Variable), scales = "free_y") +
  scale_fill_identity() +
  ylab("Delta WAIC") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(b = 35))


pS11.c <- ggplot(data = model_fit_steppe) +
  geom_line(aes(x = Model, y = deltaDIC, group = 1)) +
  geom_hline(
    data = model_fit_steppe[Model == "baseline"], aes(yintercept = deltaDIC),
    linetype = 2, color = "black"
  ) +
  geom_point(aes(x = Model, y = deltaDIC, fill = colorDIC), shape = 21, size = 2.5) +
  facet_grid(rows = vars(Variable), scales = "free_y") +
  scale_fill_identity() +
  ylab("Delta DIC") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12))

png("./Figures/Supplementary/Figure_S11.png", width = 12, height = 14, res = 330, units = 'in')
egg::ggarrange(plots = list(pS11.a, pS11.b, pS11.c), nrow = 3, labels = c("A", "B", "C"),
               label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
dev.off()

## ---------------------------
## FIGURE S12

# Data
model_fit_all <- fread("./Results/Models_fit.csv")
model_fit_all[, deltaWAIC := waic - min(waic), by = "Variable"]
model_fit_all[, colorWAIC := ifelse(deltaWAIC == 0, "red", "black")]
model_fit_all[, deltaDIC := dic - min(dic), by = "Variable"]
model_fit_all[, colorDIC := ifelse(deltaDIC == 0, "red", "black")]
model_fit_all$Variable <- factor(model_fit_all$Variable, levels = unique(model_fit_all$Variable))
model_fit_all$Model <- factor(model_fit_all$Model, levels = unique(model_fit_all$Model))

# Plot
pS12.a <- ggplot(data = model_fit_all) +
  geom_line(aes(x = Model, y = deltaWAIC, group = 1)) +
  geom_hline(
    data = model_fit_all[Model == "baseline"], aes(yintercept = deltaWAIC),
    linetype = 2, color = "black"
  ) +
  geom_point(aes(x = Model, y = deltaWAIC, fill = colorWAIC), shape = 21, size = 2.5) +
  facet_grid(rows = vars(Variable), scales = "free_y") +
  scale_fill_identity() +
  ylab("Delta WAIC") +
  labs(tag = "D") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    # plot.margin = margin(b = 35),
    # aspect.ratio = 0.3,
    plot.tag = element_text(face = 'bold', size = 18)
  )

pS12.b <- ggplot(data = model_fit_all) +
  geom_line(aes(x = Model, y = deltaDIC, group = 1)) +
  geom_hline(
    data = model_fit_all[Model == "baseline"], aes(yintercept = deltaDIC),
    linetype = 2, color = "black"
  ) +
  geom_point(aes(x = Model, y = deltaDIC, fill = colorDIC), shape = 21, size = 2.5) +
  facet_grid(rows = vars(Variable), scales = "free_y") +
  scale_fill_identity() +
  ylab("Delta WAIC") +
  labs(tag = "D") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    # aspect.ratio = 0.3,
    plot.tag = element_text(face = 'bold', size = 18)
  )

png("./Figures/Supplementary/Figure_S12.png", width = 8, height = 12, units = "in", res = 330)
egg::ggarrange(plots = list(pS12.a, pS12.b), nrow = 2, labels = c("A", "B"),
               label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
dev.off()

## ---------------------------
## FIGURE S13

## Create data table
orientation <- graves[, .(Num, BurialOrientationMin, 
                          BurialOrientationMax, 
                          BurialSide_back, 
                          BurialSide_left, 
                          BurialSide_right)]

# Remove NAs in both min and max orientation columns
orientation <- orientation[!(is.na(BurialOrientationMin) & is.na(BurialOrientationMax))]

# Create a table for records that have a single angle value (in BurialOrientationMin)
single_angles <- orientation[is.na(BurialOrientationMax)]
single_angles[, BurialOrientationMean := BurialOrientationMin]

# Remove from overall table
orientation <- orientation[!(is.na(BurialOrientationMin) | is.na(BurialOrientationMax))]

# Average across min and max 
means <- round(sapply(1:nrow(orientation), function(x) mean.circular(circular(c(orientation$BurialOrientationMin[x], orientation$BurialOrientationMax[x]), type = "angles", units = "degrees", rotation = "clock"))),2)
means <- sapply(1:length(means), function(x) ifelse(means[x] < 0, 360 + means[x], means[x]))

orientation[, BurialOrientationMean := means]
# Bind tables together
orientation <- rbind(orientation, single_angles)
orientation <- orientation[order(Num)]

# Calculate cosine and sine from radians
orientation[, BurialOrientationMean_cos := cos(BurialOrientationMean*pi/180)]
orientation[, BurialOrientationMin_cos := cos(BurialOrientationMin*pi/180)]
orientation[, BurialOrientationMax_cos := cos(BurialOrientationMax*pi/180)]
orientation[, BurialOrientationMean_sin := sin(BurialOrientationMean*pi/180)]
orientation[, BurialOrientationMin_sin := sin(BurialOrientationMin*pi/180)]
orientation[, BurialOrientationMax_sin := sin(BurialOrientationMax*pi/180)]

# Add Country, Culture, Period and Ancestry info
orientation$Longitude <- graves[Num %in% orientation$Num]$Longitude
orientation$Latitude <- graves[Num %in% orientation$Num]$Latitude
orientation$Country <- graves[Num %in% orientation$Num]$Country
orientation$Period <- graves[Num %in% orientation$Num]$Period
orientation$Culture <- graves[Num %in% orientation$Num]$Culture
orientation$Sex <- graves[Num %in% orientation$Num]$Sex
orientation$Age <- graves[Num %in% orientation$Num]$gaussianModelMu
orientation$EHG <- graves[Num %in% orientation$Num]$ANCE2_v2
orientation$WHG <- graves[Num %in% orientation$Num]$ANCE4_v2
orientation$LVN <- graves[Num %in% orientation$Num]$ANCE6_v2
orientation$CHG <- graves[Num %in% orientation$Num]$ANCE8_v2
orientation$mobility <- graves[Num %in% orientation$Num]$mobility_MDS_250y_retrospective

#### --------------------------------------------------------------------------------------------------
##### Function to count decimal places #####
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

##### Function for binning when plotting circular histogram #####
circular.bins <- function(v, binwidth = NULL){
  decimals <- max(sapply(v, decimalplaces))
  bounds.zero <- c(360 - binwidth/2, 0 + binwidth/2)
  # bounds.zero <- round(bounds.zero, decimals)
  center.zero <- seq(round(bounds.zero[1], decimals), 360, by = 1/10^decimals)
  center.zero <- center.zero[-length(center.zero)] # remove 360
  center.zero <- append(center.zero, seq(0, bounds.zero[2], by = 1/10^decimals))
  bins <- seq(bounds.zero[2], bounds.zero[1], binwidth)
  cat <- NA
  for(b in 1:(length(bins) - 1)){
    ix <- which(between(v, bins[b], bins[b+1]))
    cat[ix] <- b
  }
  ix2 <- which(v %in% center.zero)
  cat[ix2] <- 0
  return(cat)
}
#### --------------------------------------------------------------------------------------------------
w <- 22.5
orientation[, Bins := circular.bins(BurialOrientationMean, binwidth = w)]
culture_keep <- orientation[!is.na(Culture)][, .N, by = "Culture"][order(N, decreasing = T)][, Culture][1:10]
orientation_culture <- orientation[Culture %in% culture_keep]
orientation_culture$Culture <- factor(orientation_culture$Culture,
                                      labels = c("Linearbandkeramik",
                                                 "Sopot",
                                                 "Lažňany",
                                                 "Baltic Meso-Neolithic",
                                                 "Dnieper-Donets",
                                                 "Trichterbecher",
                                                 "Baden",
                                                 "Corded Ware",
                                                 "Bell Beaker",
                                                 "Únětice"),
                                      levels = c("Linearbandkeramik",
                                                 "Sopot",
                                                 "Lažňany",
                                                 "Baltic Mesolithic-Neolithic",
                                                 "Dnieper-Donets",
                                                 "Trichterbecher",
                                                 "Baden",
                                                 "Corded Ware/Single Grave/Battle Axe/Fatyanovo",
                                                 "Bell Beaker complex",
                                                 "Únětice"))

cols <- rev(pals::ocean.thermal(12))[2:11]

pS13.a <- ggplot(data = orientation, aes(x = Bins)) +
  geom_bar(fill = "grey80", color = "black") +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 2),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")) +
  coord_polar(start = -.2) +
  theme_light() +
  labs(tag = "A") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12),
        plot.tag = element_text(face = "bold", size = 22))

pS13.b <- ggplot(data = orientation_culture,
              aes(x = Bins, fill = Culture)) +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 2),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")) +
  scale_y_continuous() +
  scale_fill_manual("Culture", values = cols) +
  coord_polar(start = -.2) +
  facet_wrap(~Culture, scales = "free_y") +
  theme_light() +
  ylab("Number of individuals") +
  # labs(tag = "B") +
  theme(axis.title = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 22),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 11),
        aspect.ratio = 1,
        panel.spacing.x = unit(1.5, "lines"),
        strip.text = element_text(size = 15))

png("./Figures/Supplementary/Figure_S13.png", res = 330, units = 'in', height = 10, width = 14)
pS13.b
dev.off()

## ---------------------------
## FIGURE S14
left <- orientation[, .(Num, BurialOrientationMean, 
                        BurialOrientationMin, 
                        BurialOrientationMax, 
                        BurialOrientationMean_cos, 
                        BurialOrientationMin_cos, 
                        BurialOrientationMax_cos, 
                        BurialOrientationMean_sin, 
                        BurialOrientationMin_sin, 
                        BurialOrientationMax_sin,
                        Bins), by = "BurialSide_left"][BurialSide_left == 1]
left[, BurialSide := "Left"]
right <- orientation[, .(Num, BurialOrientationMean, 
                         BurialOrientationMin, 
                         BurialOrientationMax, 
                         BurialOrientationMean_cos, 
                         BurialOrientationMin_cos, 
                         BurialOrientationMax_cos, 
                         BurialOrientationMean_sin, 
                         BurialOrientationMin_sin, 
                         BurialOrientationMax_sin,
                         Bins), by = "BurialSide_right"][BurialSide_right == 1]
right[, BurialSide := "Right"]
back <- orientation[, .(Num, BurialOrientationMean, 
                        BurialOrientationMin, 
                        BurialOrientationMax, 
                        BurialOrientationMean_cos, 
                        BurialOrientationMin_cos, 
                        BurialOrientationMax_cos, 
                        BurialOrientationMean_sin, 
                        BurialOrientationMin_sin, 
                        BurialOrientationMax_sin,
                        Bins), by = "BurialSide_back"][BurialSide_back == 1]
back[, BurialSide := "Back"]

# Bind together in a new table
sides <- rbind(left[, -1], right[, -1], back[, -1])

# Plot
pS14 <- ggplot(data = sides, aes(x = Bins, 
                                 fill = factor(BurialSide, levels = c("Left", "Right", "Back")))) +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 2),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")) +
  scale_y_continuous() +
  scale_fill_manual("Burial side", values = pals::brewer.paired(3)) +
  facet_wrap(~factor(BurialSide, levels = c("Left", "Right", "Back")),
             scales = "free_y") +
  ylab("Number of individuals") +
  coord_polar(start = -.2) +
  theme_light() +
  theme(axis.title = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 22),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 11),
        aspect.ratio = 1,
        panel.spacing.x = unit(1.5, "lines"),
        strip.text = element_text(size = 15))

png("./Figures/Supplementary/Figure_S14.png", res = 330, units = 'in', height = 10, width = 14)
pS14
dev.off()

## ---------------------------
## FIGURE S15

# Check number of individuals by period
orientation[, .N, by = "Period"][order(N)]

# Plot
pS15 <- ggplot(data = orientation[Period != "Iron Age"], # only 1 individual
       aes(x = Bins, fill = factor(Period,
                                   levels = c("Mesolithic",
                                              "Mesolithic/Neolithic",
                                              "Neolithic",
                                              "Eneolithic",
                                              "Eneolithic/Bronze Age",
                                              "Bronze Age",
                                              "Bronze Age/Iron Age",
                                              "Iron Age")))) +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 2),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")) +
  scale_y_continuous() +
  scale_fill_paletteer_d("MetBrewer::Hiroshige",
                         guide = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(~factor(Period,
                     levels = c("Mesolithic",
                                "Mesolithic/Neolithic",
                                "Neolithic",
                                "Eneolithic",
                                "Eneolithic/Bronze Age",
                                "Bronze Age",
                                "Bronze Age/Iron Age",
                                "Iron Age")), scales = "free_y") +
  ylab("Number of individuals") +
  coord_polar(start = -.2) +
  theme_light() +
  theme(axis.title = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 22),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 11),
        aspect.ratio = 1,
        panel.spacing.x = unit(1.5, "lines"),
        strip.text = element_text(size = 15))

png("./Figures/Supplementary/Figure_S15.png", res = 330, units = 'in', height = 12, width = 12)
pS15
dev.off()

## ---------------------------
## FIGURE S16

# Group sex
orientation[, Sex := ifelse(Sex < 0.5, 0, ifelse(Sex > 0.5, 1, Sex))]

# Plot
pS16 <- ggplot(data = orientation,
       aes(x = Bins, fill = as.factor(Sex))) +
  geom_bar(color = "black") +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 2),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")) +
  scale_y_continuous() +
  scale_fill_manual("Sex", values = paletteer_c("grDevices::Tropic", 3)) +
  coord_polar(start = -.2) +
  facet_wrap(~as.factor(Sex), scales = "free_y") +
  ylab("Number of individuals") +
  theme_light() +
  theme(axis.title = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 22),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 11),
        aspect.ratio = 1,
        panel.spacing.x = unit(1.5, "lines"),
        strip.text = element_text(size = 15))

png("./Figures/Supplementary/Figure_S16.png", res = 330, units = 'in', height = 5, width = 12)
pS16
dev.off()

## ---------------------------
## FIGURE S17

# Data
bo_all <- fread("./Results/BurialOrientationMeans_AllIndividuals.csv")
bo_all$Culture <- factor(bo_all$Culture, 
                         levels = c("Linearbandkeramik",
                                    "Sopot",
                                    "Lažňany",
                                    "Baltic Meso-Neolithic",
                                    "Dnieper-Donets",
                                    "Trichterbecher",
                                    "Baden",
                                    "Corded Ware",
                                    "Bell Beaker",
                                    "Únětice"))

bo_all <- bo_all[order(Culture)]

## Visualization with circlize package
culture.means <- deg(bo_all$mean)
culture.means <- ifelse(culture.means < 0, 360 + culture.means, culture.means)
culture.ranges.lb <- deg(bo_all$LB)
culture.ranges.lb <- ifelse(culture.ranges.lb < 0, 360 + culture.ranges.lb, culture.ranges.lb)
culture.ranges.ub <- deg(bo_all$UB)
culture.ranges.ub <- ifelse(culture.ranges.ub < 0, 360 + culture.ranges.ub, culture.ranges.ub)

# params
labs <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW", NA)
cols <- rev(pals::ocean.thermal(12))[2:11]
ys <- seq(0, 1, length.out = 10)
texts <- results.dt$Culture

# plot
png("./Figures/Supplementary/Figure_S17.png", width = 13, height = 13, units = 'in', res = 330)
circos.par(gap.degree = 0, cell.padding = c(0, 0, 0, 0), start.degree = 90)
circos.initialize("a", xlim = c(0, 360), ring = T)
circos.track(ylim = c(0.75, 1.05), bg.border = NA)
circos.axis(major.at = seq(0, 360, 45), 
            direction = "outside",
            labels = labs,
            labels.cex = 2)
circos.segments(x0 = culture.ranges.lb[1], 
                y0 = ys[1], 
                x1 = culture.ranges.ub[1], 
                y1 = ys[1], col = cols[1], 
                lwd = 2)
circos.text(x = culture.means[1], 
            y = ys[1] - 0.05,
            labels = texts[1], 
            col = cols[1],
            facing = "bending.outside",
            cex = 1.8)
circos.segments(x0 = culture.ranges.lb[2], 
                y0 = ys[2], 
                x1 = culture.ranges.ub[2], 
                y1 = ys[2], col = cols[2], 
                lwd = 2)
circos.text(x = culture.means[2], 
            y = ys[2] - 0.05,
            labels = texts[2], 
            col = cols[2],
            facing = "bending.inside",
            cex = 1.8)
circos.segments(x0 = culture.ranges.lb[3], 
                y0 = ys[3], 
                x1 = culture.ranges.ub[3], 
                y1 = ys[3], col = cols[3], 
                lwd = 2)
circos.text(x = culture.means[3], 
            y = ys[3] - 0.05,
            labels = texts[3], 
            col = cols[3],
            facing = "bending.outside",
            cex = 1.8)
circos.segments(x0 = culture.ranges.lb[4], 
                y0 = ys[4], 
                x1 = culture.ranges.ub[4], 
                y1 = ys[4], col = cols[4], 
                lwd = 2)
circos.text(x = culture.means[4], 
            y = ys[4] - 0.05,
            labels = texts[4], 
            col = cols[4],
            facing = "bending.outside",
            cex = 1.8)
circos.segments(x0 = culture.ranges.lb[5], 
                y0 = ys[5], 
                x1 = culture.ranges.ub[5], 
                y1 = ys[5], col = cols[5], 
                lwd = 2)
circos.text(x = culture.means[5], 
            y = ys[5] - 0.05,
            labels = texts[5], 
            col = cols[5],
            facing = "bending.outside",
            cex = 1.8)
circos.segments(x0 = culture.ranges.lb[6], 
                y0 = ys[6], 
                x1 = culture.ranges.ub[6], 
                y1 = ys[6], col = cols[6], 
                lwd = 2)
circos.text(x = culture.means[6], 
            y = ys[6] - 0.05,
            labels = texts[6], 
            col = cols[6],
            facing = "bending.outside",
            cex = 1.8)
circos.segments(x0 = culture.ranges.lb[7], 
                y0 = ys[7], 
                x1 = culture.ranges.ub[7], 
                y1 = ys[7], col = cols[7], 
                lwd = 2)
circos.text(x = culture.means[7], 
            y = ys[7] - 0.05,
            labels = texts[7], 
            col = cols[7],
            facing = "bending.outside",
            cex = 1.8)
circos.segments(x0 = culture.ranges.lb[8], 
                y0 = ys[8], 
                x1 = culture.ranges.ub[8], 
                y1 = ys[8], col = cols[8], 
                lwd = 2)
circos.text(x = culture.means[8], 
            y = ys[8] - 0.05,
            labels = texts[8], 
            col = cols[8],
            facing = "bending.outside",
            cex = 1.8)
circos.segments(x0 = culture.ranges.lb[9], 
                y0 = ys[9], 
                x1 = culture.ranges.ub[9], 
                y1 = ys[9], col = cols[9], 
                lwd = 2)
circos.text(x = culture.means[9], 
            y = ys[9] - 0.05,
            labels = texts[9], 
            col = cols[9],
            facing = "bending.outside",
            cex = 1.8)
circos.segments(x0 = culture.ranges.lb[10], 
                y0 = ys[10], 
                x1 = culture.ranges.ub[10], 
                y1 = ys[10], col = cols[10], 
                lwd = 2)
circos.text(x = culture.means[10], 
            y = ys[10] - 0.05,
            labels = texts[10], 
            col = cols[10],
            facing = "bending.outside",
            cex = 1.8)
circos.points(x = culture.means, y = ys, bg = cols, pch = 21, cex = 2)
circos.clear()
dev.off()

## --------------------------------------------------------------------------------------------
## INLA MODELS WITH CULTURE - BEST MODELS

## Results
# cols <- rev(pals::ocean.phase(6)[2:5])
# cols[1] <- "grey50"
cols <- paletteer::paletteer_d("nationalparkcolors::Arches", 4, direction = 1)

## Full model no WHG
results_df <- fread("./Results/Culture/BestModel_results.csv")
results_df$Predictors <- factor(results_df$Predictors,
                                levels = rev(c("Intercept", "EHG", "CHG", "LVN", "Mobility", levels(DT$Culture))))
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = c("Intercept", "Ancestry", "Mobility", "Culture"))

# Plot
p <- ggplot(data = results_df) +
  geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
                width = 0.2, linewidth = 1) +
  geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
  geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  facet_grid(cols = vars(Variable), rows = vars(Effect), scales = "free") +
  theme_bw() +
  # labs(tag = "C") +
  # ggtitle("Full model") +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        # aspect.ratio = 0.9,
        # plot.margin = margin(t = 10, b = 0, l = 5, r = 5),
        plot.tag = element_text(face = 'bold', size = 18))

png("./Figures/Supplementary/Figure_S21.png", res = 330, units = 'in', height = 8, width = 14)
p
dev.off()
