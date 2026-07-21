#### ------- ####
#### FIGURES ####
#### ------- ####

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)
library(factoextra)
library(FactoMineR)
library(UpSetR)
# library(BIADconnect)
library(cowplot)
library(patchwork)
library(circular)
library(CircStats)
library(circlize)
library(bpnreg)
library(sf)
library(terra)

# setwd("~/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Projects/Graves/")
# Sys.setenv(MARIADB_TLS_DISABLE_PEER_VERIFICATION = "1")
# conn <- init.conn()
proj <- "+proj=laea +lon_0=0 +lat_0=49.06 +datum=WGS84 +units=km +no_defs"

## ------------------------------------------------------------------------------------------------------
## Coastlines and continents poly
conts <- rnaturalearth::ne_download(scale = 10, type = "coastline", "physical", returnclass = "sf")
land <- rnaturalearth::ne_download(scale = 10, type = "land", "physical", returnclass = "sf")
countries <- rnaturalearth::ne_download(scale = 10, type = "countries", "cultural", returnclass = "sf")
grat <- rnaturalearth::ne_download(scale = 10, type = "graticules_10", "physical", returnclass = "sf")
conts_prj <- st_transform(conts, crs = proj)
land_prj <- st_transform(land, proj)
countries_prj <- st_transform(countries, crs = proj)
grat_prj <- st_transform(grat, proj)

## Data
graves <- fread("./Data/Burial.csv", na.strings = "")
graves$DepositionType <- as.factor(graves$DepositionType)
graves$BodyPositioning <- as.factor(graves$BodyPositioning)
graves$BurialSide <- as.factor(graves$BurialSide)

# # Sites table (BIAD)
# sites <- setDT(query.database("SELECT * FROM `Sites`", conn))
# 
# # Culture and Period
# culture <- fread("./Data/burial rites_culture_period.csv", na.strings = "\\N")
# culture <- culture[!duplicated(IndividualID)]
# 
# # Add info to graves table
# graves <- merge(graves, sites[, .(SiteID, Country)], by = "SiteID", all.x = T)
# rm(sites)
# graves <- merge(graves, culture, by = "IndividualID", all.x = T)

# Remove outlier
graves <- graves[!is.na(YearBP)]
graves <- graves[YearBP < 15000]
summary(graves$YearBP)

## ------------------------------------------------------------------------------------------------------
## FIGURE 1

# Ancestry and mobility data
ancestry <- fread("./Data/combined.tsv")
ov <- fread("./Data/mobility_estimates_250y_retrospecive_distance.csv")
info <- fread("./Data/Dataset_S1.csv")
mobility <- merge(ov, info, by = "Sample_ID", all.x = T)
rm(ov, info)

# Transform to spatial points
gr_pnts <- project(
  vect(graves[, .(Longitude, Latitude, aDNAID, Num)], geom = c("Longitude", "Latitude"), 
       crs = crs(rast())),
  proj)
ance_pnts <- project(
  vect(ancestry[, .(LON, LAT, IND)], 
       geom = c("LON", "LAT"),
       crs = crs(rast())),
  proj)
mob_pnts <- project(
  vect(mobility[, .(Longitude, Latitude, Sample_ID)], 
       geom = c("Longitude", "Latitude"),
       crs = crs(rast())),
  proj)

# Get common extent
e <- c(min(c(ext(mob_pnts)[1], ext(gr_pnts)[1]), ext(ance_pnts)[1]),
       max(c(ext(mob_pnts)[2], ext(gr_pnts)[2]), ext(ance_pnts)[2]),
       min(c(ext(mob_pnts)[3], ext(gr_pnts)[3]), ext(ance_pnts)[3]),
       max(c(ext(mob_pnts)[4], ext(gr_pnts)[4]), ext(ance_pnts)[4]))

# Add type to points
gr_pnts$Type = "Graves"
mob_pnts$Type = "Mobility"
ance_pnts$Type = "Ancestry"

# Merge points
all_pnts <- rbind(ance_pnts, mob_pnts, gr_pnts)
all_pnts$Type <- factor(all_pnts$Type, 
                        levels = c("Graves", "Mobility", "Ancestry"),
                        labels = c("Graves (BIAD)",
                                   "Genomic sequences (Schmid & Schiffels, 2023)",
                                   "Genomes (Allentoft et al., 2024)"))

# Plot
f1 <- ggplot() +
  geom_sf(data = land_prj, inherit.aes = T, fill = "grey98", color = "black", linewidth = 0.3) +
  geom_sf(data = st_as_sf(all_pnts), aes(fill = Type, 
                                         shape = Type,
                                         size = Type,
                                         alpha = Type)) +
  scale_fill_manual(values = paletteer::paletteer_d("nationalparkcolors::Acadia",
                                                    6)[c(1,5,3)],
                    guide = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  scale_shape_manual(values = c(24, 22, 21)) +
  scale_size_manual(values = c(2.5, 2, 3)) +
  scale_alpha_manual(values = c(0.8, 0.5, 1)) +
  coord_sf(xlim = e[1:2], ylim = e[3:4]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "aliceblue"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12))

png("./Figures/Main/Fig1.png", res = 330, units = 'in', height = 8, width = 10)
f1
dev.off()


## ------------------------------------------------------------------------------------------------------
## FIGURE 3
bp <- grep("BodyPositioning_", colnames(graves))
bs <- grep("BurialSide_", colnames(graves))
mat <- as.data.frame(cbind(graves[, ..bp], graves[, ..bs]))
mat <- as.matrix(mat)
dim(mat)

## Multiple Correspondence Analysis
mca <- MCA(graves[, .(BurialSide, BodyPositioning)], graph = F)
summary(mca)

# Eigenvalues
eig.val <- get_eigenvalue(mca)
# Extract results for individuals
ind <- get_mca_ind(mca)
# Extract results for variables
var <- get_mca_var(mca)

## Create dataframe
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

## FIG 2A
dim_biplot <- ind_dims[, .N, by = c("Dim 1", "Dim 2")]

f2a <- ggplot() +
  geom_point(data = dim_biplot,
             aes(x = `Dim 1`, y = `Dim 2`,
                 size = N),
             position = position_jitter(0.2, 0.2),
             shape = 21,
             fill = "#78CEA3",
             alpha = 0.8) +
  scale_size_continuous("Count", breaks = c(1, 50, 100, 250, 500, 750, 1000),
                        range = c(3,20)) +
  geom_point(data = var_dims,
             aes(x = `Dim 1`, y = `Dim 2`),
             shape = 23,
             fill = 'white',
             alpha = 0.5,
             size = 3,
             stroke = 1) +
  geom_label_repel(data = var_dims,
                   aes(x = `Dim 1`, y = `Dim 2`, label = Variable),
                   segment.color = 'grey30',
                   segment.alpha = 0.5,
                   segment.curvature = 0.05,
                   fill = alpha("white", 0.8),
                   color = "grey20",
                   max.overlaps = 50,
                   nudge_x = 1,
                   nudge_y = 0.5,
                   direction = "both") +
  xlab(paste0("Dim 1 (", round(eig.val[1,2], 2), "%)")) +
  ylab(paste0("Dim 2 (", round(eig.val[2,2], 2), "%)")) +
  theme_bw() +
  theme(legend.position = "right") +
  coord_fixed(1/1.8)

## UpSet Plot
upset_df <- graves[, .(IndividualID, 
                       BodyPositioning_crouched, 
                       BodyPositioning_extended, 
                       BodyPositioning_heap, 
                       `BodyPositioning_inside a vessel`, 
                       BodyPositioning_prone, 
                       BodyPositioning_scattered, 
                       BodyPositioning_sitting, 
                       BurialSide_back, 
                       BurialSide_bottom, 
                       BurialSide_front, 
                       BurialSide_left, 
                       BurialSide_right)]

# Change column names
varnames <- colnames(upset_df)[-1]
vartypes <- sapply(strsplit(varnames, "_"), "[", 1)
vartypes <- gsub("BodyPositioning", "Body positioning", vartypes)
vartypes <- gsub("BurialSide", "Burial side", vartypes)
varnames <- sapply(strsplit(varnames, "_"), "[", 2)
varnames <- sapply(varnames, gtools::capwords)
colnames(upset_df) <- c("IndividualID", varnames)

# Define sets
sets <- varnames

# Metadata table
metadata <- data.frame(sets = varnames, Type = vartypes)

# Generate plot
f2b <- upset(upset_df, sets = sets, 
             order.by = 'freq', 
             nintersects = 100,
             point.size = 3.5,
             mb.ratio = c(0.6, 0.4),
             matrix.dot.alpha = 0.5,
             text.scale = 2,
             set.metadata = list(data = metadata,
                                 plots = list(list(type = "text",
                                                   column = "Type",
                                                   assign = 10,
                                                   colors = c("Body positioning" = "#363873", #"#1283D4"
                                                              "Burial side" = "#4379B6")),
                                              list(type = "matrix_rows",
                                                   column = "Type",
                                                   colors = c("Body positioning" = "#3F4079",
                                                              "Burial side" = "#417498"),
                                                   alpha = 0.8))))



Main_bar_plot = f2b$Main_bar
Matrix_plot = f2b$Matrix
Size_plot = f2b$Sizes
Set_data = f2b$New_data
start_col = f2b$first.col
hratios = f2b$mb.ratio
end_col <- ((start_col + as.integer(length(labels))) - 1)
Set_data <- Set_data[which(rowSums(Set_data[ ,start_col:end_col]) != 0), ]
Main_bar_plot$widths <- Matrix_plot$widths
Matrix_plot$heights <- Size_plot$heights
size_plot_height <- (((hratios[1])+0.01)*100)
matrix_and_mainbar_right <- 100
matrix_and_mainbar_left <- 21
size_bar_right <- 20
size_bar_left <- 1
f2b <- arrangeGrob(Main_bar_plot, Matrix_plot, heights = hratios)

# Plot together
png("./Figures/Main/Fig3a.png", res = 330, units = 'in', height = 8, width = 10)
f2a
dev.off()

png("./Figures/Main/Fig3b.png", res = 330, units = 'in', height = 8, width = 10)
plot_grid(f2b)
dev.off()

png("/Users/msb290/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Projects/Graves/Figures/Main/Fig3.png", res = 330, units = 'in', height = 8, width = 16)
# (f2a + labs(tag = c("(a)")) + theme(plot.tag = element_text(face = "italic", size = 11))) | patchwork::wrap_elements(f2b) + labs(tag = c("(b)"))  + theme(plot.tag = element_text(face = "italic", size = 11))
egg::ggarrange(plots = list(f2a, ggdraw(f2b)), ncol = 2, labels = c("(a)", "(b)"), label.args = list(gp = grid::gpar(font = 3, cex = 1.2)))
dev.off()


## ------------------------------------------------------------------------------------------------------
## FIGURE 4
pnts <- project(
  vect(
    graves[, .(Num, Longitude, Latitude, YearBP,
               Left = BurialSide_left,
               Right = BurialSide_right,
               Back = BurialSide_back,
               mobility = mobility_MDS_250y_retrospective,
               EHG = ANCE2_v2,
               WHG = ANCE4_v2,
               LVN = ANCE6_v2,
               CHG = ANCE8_v2 
    )],
    geom = c("Longitude", "Latitude"),
    crs = crs(rast())
  ),
  proj
)

poly <- vect("./Data/Spatial/IntersectionVector.shp")
pnts <- intersect(pnts, poly)

# Remove data without age and too old sample
pnts <- st_as_sf(pnts)
pnts <- pnts[!is.na(pnts$YearBP),]
pnts <- pnts[which(pnts$YearBP < 15000),]

pnts.dt <- setDT(as.data.frame(pnts))
pnts.dt <- melt(pnts.dt, measure.vars = c("Left", "Right", "Back"))
pnts.dt <- pnts.dt[!is.na(value)]

p3.1 <- ggplot(data = pnts.dt) +
  geom_sf(data = countries_prj, inherit.aes = T, fill = "grey90") +
  geom_sf(
    data = pnts.dt$geometry,
    inherit.aes = T,
    shape = 21,
    fill = "#78CEA3",
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

p3.2 <- ggplot(data = pnts.dt) +
  geom_density(
    data = pnts.dt[value == 1], aes(x = YearBP),
    fill = "grey80",
    color = "transparent"
  ) +
  geom_point(aes(
    x = YearBP, y = 0,
    fill = as.factor(value)
  ), shape = 21, size = 1.5, alpha = 0.8, color = "grey30") +
  scale_fill_manual(values = paletteer::paletteer_d("Redmonder::sPBIYlGn", n = 9, direction = -1)[c(1,9)]) +
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
# cols <- rev(pals::ocean.phase(5)[2:4])
# cols[1] <- "grey50"
cols <- paletteer::paletteer_d("ButterflyColors::anteos_menippe", 3)

## Full model no WHG
results_df <- fread("./Results/FullModelNoWHG_results.csv")
results_df$Predictors <- factor(results_df$Predictors,
                                levels = rev(c("Intercept", "EHG", "CHG", "LVN", "Mobility")),
                                labels = rev(c("Intercept", "EHG", "CHG", "NEOL", "Mobility")))
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = c("Intercept", "Ancestry", "Mobility"))

# Plot
f3a <- ggplot(data = results_df) +
  geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
                width = 0.2, linewidth = 1) +
  geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
  geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  facet_grid(cols = vars(Variable), scales = "free_x") +
  theme_bw() +
  labs(tag = "C") +
  # ggtitle("Full model") +
  theme(plot.title = element_text(size = 11, face = 'bold'),
        axis.title = element_blank(),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size = 11),
        # aspect.ratio = 0.9,
        # plot.margin = margin(t = 10, b = 0, l = 5, r = 5),
        plot.tag = element_text(face = 'italic', size = 11))

# Model fit
model_fit <- fread("./Results/Models_fit.csv")
model_5k_fit <- fread("./Results/ModelsFit_BeforeSteppe.csv")
model_steppe_fit <- fread("./Results/ModelsFit_Steppe.csv")

## Plot model fit
# Whole period
model_fit <- model_fit[!(Model %in% c("baseline + mobility + ancestry",
                                      "baseline + ancestry"))]
model_fit[, deltaWAIC := waic - min(waic), by = "Variable"]
model_fit[, colorWAIC := ifelse(deltaWAIC == 0, "#78CEA3", "black")]
model_fit[, deltaDIC := dic - min(dic), by = "Variable"]
model_fit[, colorDIC := ifelse(deltaDIC == 0, "#78CEA3", "black")]
model_fit$Variable <- factor(model_fit$Variable, levels = unique(model_fit$Variable))
model_fit$Model <- gsub("no WHG ancestry", "ancestry", model_fit$Model)
model_fit$Model <- gsub("LVN", "NEOL", model_fit$Model)
model_fit$Model <- factor(model_fit$Model, levels = unique(model_fit$Model))

f3b <- ggplot(data = model_fit) +
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
    strip.text = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    # aspect.ratio = 0.3,
    plot.tag = element_text(face = 'italic', size = 11)
  )

p3a <- egg::ggarrange(plots = list(p3.1, p3.2), ncol = 1, nrow = 2,
                      labels = c("(a)", "(b)"), 
                      label.args = list(gp = grid::gpar(font = 3, cex = 1.2))
                      )
p3b <- egg::ggarrange(plots = list(f3a, f3b), ncol = 1, nrow = 2, heights = c(.6,1),
                      labels = c("(c)", "(d)"), 
                      label.args = list(gp = grid::gpar(font = 3, cex = 1.2))
                      )

png("/Users/msb290/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Projects/Graves/Figures/Main/Fig4.png", res = 330, units = 'in', height = 8, width = 14)
plot_grid(p3a, p3b, rel_widths = c(0.8, 1))
dev.off()

## ------------------------------------------------------------------------------------------------------
## FIGURE 5

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
# cols <- rev(pals::ocean.phase(6)[2:5])
# cols[1] <- "grey50"
# cols <- paletteer::paletteer_d("nationalparkcolors::Arches", 4, direction = 1)
# cols[4] <- "#5D3E99"

cols <- pals::ocean.thermal(12)[2:11]

## Full model culture
# results_df <- fread("./Results/Culture/FullModel_results.csv")
results_df <- fread("./Results/Culture/Results_culture_test2.csv")
# results_df$Predictors <- factor(results_df$Predictors,
#                                 levels = rev(c("Intercept", "EHG", "CHG", "LVN", "Mobility", levels(DT$Culture))))
results_df$Predictors <- factor(results_df$Predictors,
                                levels = rev(levels(DT$Culture)))
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
# results_df$Type <- factor(results_df$Type, levels = c("Intercept", "Ancestry", "Mobility", "Culture"))

# Plot
f4a <- ggplot(data = results_df) +
  geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Predictors),
                width = 0.2, linewidth = 1) +
  geom_point(aes(x = mean, y = Predictors, fill = Predictors), shape = 21, size = 3) +
  geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  facet_grid(cols = vars(Variable), scales = "free") +
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
# model_fit <- fread("./Results/Culture/Models_fit.csv")
model_fit <- fread("./Results/Culture/Models_fit_test2.csv")

## Plot model fit
model_fit[, deltaWAIC := waic - min(waic), by = "Variable"]
model_fit[, colorWAIC := ifelse(deltaWAIC == 0, pals::ocean.thermal(6)[4], "black")]
model_fit[, deltaDIC := dic - min(dic), by = "Variable"]
model_fit[, colorDIC := ifelse(deltaDIC == 0, pals::ocean.thermal(6)[4], "black")]
model_fit$Variable <- factor(model_fit$Variable, levels = unique(model_fit$Variable))
model_fit$Model <- gsub("no WHG ancestry", "ancestry", model_fit$Model)
model_fit$Model <- gsub("LVN", "NEOL", model_fit$Model)
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
                      labels = c("(a)", "(b)"), 
                      label.args = list(gp = grid::gpar(font = 3, cex = 1.2)))
p4b <- egg::ggarrange(plots = list(f4a, f4b), ncol = 1, nrow = 2, heights = c(1,0.8),
                      labels = c("(c)", "(d)"), 
                      label.args = list(gp = grid::gpar(font = 3, cex = 1.2)))

png("/Users/msb290/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Projects/Graves/Figures/Main/Fig5.png", res = 330, units = 'in', height = 10, width = 18)
plot_grid(p4a, p4b, rel_widths = c(0.8, 1))
dev.off()


## ------------------------------------------------------------------------------------------------------
## FIGURE 6

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

f4b <- ggplot(data = orientation_culture[!is.na(EHG) & !is.na(Age)],
       aes(x = Bins, fill = Culture)) +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 2),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")) +
  scale_fill_manual("Culture", values = cols) +
  coord_polar(start = -.2) +
  facet_wrap(~Culture, scales = "free_y") +
  theme_light() +
  ylab("Number of individuals") +
  # labs(tag = "B") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 15),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 22),
        aspect.ratio = 1,
        panel.spacing.x = unit(1.5, "lines"))

orientation$Type <- "All individuals"
f4a <- ggplot(data = orientation, aes(x = Bins)) +
  geom_bar(fill = "grey80", color = "black") +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 2),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")) +
  coord_polar(start = -.2) +
  facet_wrap(~Type, scales = "free_y") +
  theme_light() +
  # labs(tag = "A") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 15),
        plot.tag = element_text(face = "bold", size = 22),
        aspect.ratio = 1)

png("./Figures/Main/Fig6.png", res = 330, units = 'in', height = 8, width = 12)
# f4b + inset_element(f4a, left = 0.4786, bottom = 0.067, top = 0.3674, right = 0.7889, align_to = "plot")
f4b + inset_element(f4a, left = 0.305, bottom = -0.01, top = 0.3458, right = 0.947, align_to = "plot")
dev.off()

## ------------------------------------------------------------------------------------------------------
## FIGURE 7

## Spatial points
orientation_pnts <- project(
  vect(
    orientation_culture,
    geom = c("Longitude", "Latitude"),
    crs = crs(rast())
  ),
  proj
)

orientation_pnts <- st_as_sf(orientation_pnts)

orientation_pnts.dt <- setDT(as.data.frame(orientation_pnts))

p7.1 <- ggplot(data = orientation_pnts.dt) +
  geom_sf(data = countries_prj, inherit.aes = T, fill = "grey90") +
  geom_sf(
    data = orientation_pnts.dt$geometry,
    inherit.aes = T,
    shape = 21,
    fill = pals::ocean.thermal(6)[2],
    size = 2.5,
    alpha = 0.8
  ) +
  coord_sf(
    xlim = st_bbox(orientation_pnts.dt$geometry)[c(1,3)],
    ylim = st_bbox(orientation_pnts.dt$geometry)[c(2,4)]
  ) +
  # labs(tag = "A") +
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

p7.2 <- ggplot(data = orientation_pnts.dt) +
  geom_density(
    aes(x = Age),
    fill = "grey80",
    color = "transparent"
  ) +
  geom_point(aes(x = Age, y = 0),
             shape = 21, 
             size = 1.5, 
             color = pals::ocean.thermal(6)[1],
             fill = alpha(pals::ocean.thermal(6)[2], 0.5)) +
  scale_x_reverse() +
  xlab("Years (cal BP)") +
  # labs(tag = "B") +
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

## Models list
models.list <- readRDS("./Results/BurialOrientationModels.RDS")
fit.cult <- models.list[[4]]

## Models fit
fit.dt <- fread("./Results/BurialOrientationModelsFit.csv")
fit.dt[, DeltaDIC := DIC - min(DIC)]
fit.dt[, Color := ifelse(DeltaDIC == 0, pals::ocean.thermal(10)[9], "black")]
fit.dt[, DeltaWAIC := WAIC - min(WAIC)]
fit.dt[, Color := ifelse(DeltaWAIC == 0, pals::ocean.thermal(10)[9], "black")]
fit.dt$Model <- gsub(" only", "", fit.dt$Model)


# Means
results.dt <- setDT(data.frame(coef_circ(fit.cult, "categorical", "radians")$Means))[1:10,]
results.dt[, Culture := gsub("Culture", "", rownames(coef_circ(fit.cult, "categorical", "radians")$Means[1:10,]))]
results.dt$Culture <- gsub("\\(Intercept\\)", "Linearbandkeramik", results.dt$Culture)
results.dt$Culture <- factor(results.dt$Culture, levels = levels(orientation_culture$Culture))
results.dt <- results.dt[order(Culture)]
results.dt[, .(Culture, deg(.SD)), .SDcols = c("mean", "mode", "sd", "LB", "UB")]

## Visualization with circlize package
culture.means <- deg(results.dt$mean)
culture.means <- ifelse(culture.means < 0, 360 + culture.means, culture.means)
culture.ranges.lb <- deg(results.dt$LB)
culture.ranges.lb <- ifelse(culture.ranges.lb < 0, 360 + culture.ranges.lb, culture.ranges.lb)
culture.ranges.ub <- deg(results.dt$UB)
culture.ranges.ub <- ifelse(culture.ranges.ub < 0, 360 + culture.ranges.ub, culture.ranges.ub)

# params
labs <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW", NA)
# labs[17] <- "0"
cols <- rev(pals::ocean.thermal(12))[2:11]
ys <- seq(0, 1, length.out = 10)
texts <- results.dt$Culture

# plot
png("./Figures/Main/Fig7a.png", width = 13, height = 13, units = 'in', res = 330)
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

f7d <- ggplot(data = fit.dt) +
  geom_line(aes(x = Model, y = DeltaWAIC, group = 1)) +
  geom_point(aes(x = Model, y = DeltaWAIC, fill = Color), size = 3, shape = 21) +
  scale_fill_identity() +
  ylab("Delta WAIC") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 11),
        # aspect.ratio = 1,
        plot.margin = margin(t = 10))

f7c <- grid::rasterGrob(magick::image_read("./Figures/Main/Fig7a.png"), interpolate = T)

png("/Users/msb290/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Projects/Graves/Figures/Main/Fig7.png", res = 330, units = 'in', height = 9, width = 15)
plot_grid(p7.1, f7c, p7.2, f7d, rel_heights = c(1, 0.8), labels = c("(a)", "(c)", "(b)", "(d)"), label_fontface = "italic", label_size = 11)
dev.off()

## DIFFERENCES TABLE
tbl <- fread("~/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Projects/Graves/Results/BurialOrientationDifferences_culture_ance.csv")

ggplot(data = tbl) +
  geom_errorbar(aes(xmin = LB, xmax = UB, y = Comparison, color = Overlapping),
                width = 0.2, linewidth = 1) +
  # geom_blank(aes(x = -mean, xmin = -LB, xmax = -UB)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 0.9,
        plot.title = element_text(size = 20, face = 'bold'),
        axis.title = element_blank(),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size = 16),
        plot.tag = element_text(face = 'bold', size = 22))

# Calculate overlap
pairs <- t(combn(nrow(results.dt), 2))
pairs <- data.frame(i = pairs[,1], j = pairs[,2])
pairs$overlap <- with(results.dt,
                         pmax(LB[pairs$i], LB[pairs$j]) <=
                           pmin(UB[pairs$i], UB[pairs$j])
)

pairs$Culture1 <- results.dt$Culture[pairs$i]
pairs$Culture2 <- results.dt$Culture[pairs$j]
