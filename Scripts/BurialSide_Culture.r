library(data.table)
library(ggplot2)
library(terra)
library(sf)
library(spatstat)
library(inlabru)
library(INLA)
library(fmesher)
library(pbapply)
library(future)
library(future.apply)
library(furrr)
library(patchwork)
setwd("/maps/projects/racimolab/people/msb290/Graves/")
proj <- "+proj=laea +lon_0=0 +lat_0=49.06 +datum=WGS84 +units=km +no_defs"
proj2 <- "+proj=aea +lon_0=44.296875 +lat_1=43.7864128 +lat_2=69.9512657 +lat_0=56.8688392 +datum=WGS84 +units=m +no_defs"

## ---- ##
## DATA ##
## ---- ##

## Coastlines and continents poly
conts <- rnaturalearth::ne_download(scale = 10, type = "coastline", "physical", returnclass = "sf")
land <- rnaturalearth::ne_download(scale = 10, type = "land", "physical", returnclass = "sf")
countries <- rnaturalearth::ne_download(scale = 50, type = "countries", "cultural", returnclass = "sf")
grat <- rnaturalearth::ne_download(scale = 10, type = "graticules_10", "physical", returnclass = "sf")
conts_prj <- st_transform(conts, crs = proj)
land_prj <- st_transform(land, proj)
countries_prj <- st_transform(countries, crs = proj)
grat_prj <- st_transform(grat, proj)

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

## Spatial points
pnts <- project(
  vect(
    DT,
    geom = c("Longitude", "Latitude"),
    crs = crs(rast())
  ),
  proj
)

poly <- vect("./Spatial/IntersectionVector.shp")
pnts <- intersect(pnts, poly)
pnts <- st_as_sf(pnts)

# Subset the pnts to remove duplicated sites
coords <- st_coordinates(pnts)[, c("X", "Y")]
coords <- coords[!duplicated(coords),]
pnts_new <- vect(coords, crs = proj)
n <- dim(coords)[1]

## -------------- ##
## MODEL SETTINGS ##
## -------------- ##

## Spatial mesh
# Define spatial range and parameters - diagonal (euclidean distance) of the bounding box around the points, in the unit of the coordinates.
bbox_coords <- st_bbox(pnts_new)
summary(dist(coords))
sp.range <- sqrt(sum(c(diff(bbox_coords[c(1,3)])^2, diff(bbox_coords[c(2,4)])^2)))

# Boundary
bound <- fmesher::fm_nonconvex_hull_inla(coords, convex = -0.2, crs = proj)

# Spatial mesh
mesh <-  fm_mesh_2d(loc.domain = coords,
                    boundary = bound,
                    max.edge = c(0.03, 0.1)*sp.range,
                    offset = 0.3*sp.range,
                    cutoff = 0.01*sp.range,
                    crs = proj)

# Plot
ggplot() +
    gg(mesh) +
    geom_sf(data = conts_prj, col = "green4") +
    geom_sf(data = st_as_sf(pnts_new), shape = 21, fill = "blue", size = 1) +
    theme_bw() +
    coord_sf(xlim = range(mesh$loc[, 1]), ylim = range(mesh$loc[, 2]))

# Create a sampler polygon for the model
cover <- st_read("./Spatial/cover.shp")
st_crs(cover) <- proj2
cover <- st_transform(cover, proj)

## Temporal mesh
tknots <- seq(0, 12000, 1000)
mesh.t <- fm_mesh_1d(tknots, boundary = "free")
k <- length(tknots)

## Spatial point pattern
pp <- ppp(coords[,1], coords[,2], bbox_coords[c(1,3)], bbox_coords[c(2,4)], unitname = c("kilometers"))

# Cross Validated Bandwidth Selection
bw <- round(bw.diggle(pp))

## Spatial priors
range.s <- c(bw, 0.01)
sigma.s <- c(1, 0.01)

## Matern
matern <- inla.spde2.pcmatern(mesh,
  prior.sigma = sigma.s,
  prior.range = range.s
)

## Prior for temporal auto-regressive parameter
h.spec <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))

## List of model components
models <- c(
    mc1 = ~ Intercept(1) +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    mc2 = ~ Intercept(1) +
        mobility_scaled +
        EHG_scaled +
        LVN_scaled +
        CHG_scaled +
        culture(Culture,
            model = "factor_contrast"
        ) +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    mc3 = ~ Intercept(1) +
        mobility_scaled +
        EHG_scaled +
        LVN_scaled +
        CHG_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    mc4 = ~ Intercept(1) +
        EHG_scaled +
        LVN_scaled +
        CHG_scaled +
        culture(Culture,
            model = "factor_contrast"
        ) +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    mc5 = ~ Intercept(1) +
        mobility_scaled +
        culture(Culture,
            model = "factor_contrast"
        ) +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    mc6 = ~ Intercept(1) +
        EHG_scaled +
        LVN_scaled +
        CHG_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    mc7 = ~ Intercept(1) +
        mobility_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    mc8 = ~ Intercept(1) +
        culture(Culture,
            model = "factor_contrast"
        ) +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        )
)

## ----------------- ##
## BURIAL SIDE: LEFT ##
## ----------------- ##

# Data
left <- pnts[, c("Num", "YearBP", "Left", "mobility", "EHG", "WHG", "LVN", "CHG", "Culture", "geometry")]
left <- left[!is.na(left$Left),]
nrow(left)

# Jitter the duplicates
duplicates <- left[duplicated(left$geometry),]
duplicates <- st_jitter(duplicates, 0.05)

plot(conts_prj$geometry, xlim = ext(left)[1:2], ylim = ext(left)[3:4]) 
points(left, pch = 1, cex = 1); points(duplicates, pch = 1, col = 'red', cex = 1)

# Add duplicates to the left data
left <- left[!duplicated(left$geometry),]
left <- rbind(left, duplicates)
nrow(left) # same as before

# Divide into time bins
left$YearBin <- plyr:::round_any(left$YearBP, 250)

## Define covariates and scale
namescov <- c("mobility", "EHG", "WHG", "LVN", "CHG")
newnames <- paste0(namescov, "_scaled")
covs <- data.frame(left)[, namescov]
covs <- scale(covs)
colnames(covs) <- newnames
left <- cbind(left, covs)

## FIT MODELS
plan(multisession, workers = 8)
left_models <- future_map(models, function(x) {
    m <- bru(as.formula(paste("Left ~", x)[2]),
        data = left,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T, control.compute = list(cpo = T))
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)

plan(sequential)
save(left_models, file = "./Models/Culture/Left.RData")

## ------------------ ##
## BURIAL SIDE: RIGHT ##
## ------------------ ##

# Data
right <- pnts[, c("YearBP", "Right", "mobility", "EHG", "WHG", "LVN", "CHG", "Culture", "geometry")]
right <- right[!is.na(right$Right),]
nrow(right)

# Jitter the duplicates
duplicates <- right[duplicated(right$geometry),]
duplicates <- st_jitter(duplicates, 0.05)

plot(conts_prj$geometry, xlim = ext(right)[1:2], ylim = ext(right)[3:4]) 
points(right, pch = 1, cex = 1); points(duplicates, pch = 1, col = 'red', cex = 1)

# Add duplicates to the right data
right <- right[!duplicated(right$geometry),]
right <- rbind(right, duplicates)
nrow(right) # same as before

# Divide into time bins
right$YearBin <- plyr:::round_any(right$YearBP, 250)

## Define covariates
namescov <- c("mobility", "EHG", "WHG", "LVN", "CHG")
newnames <- paste0(namescov, "_scaled")
covs <- data.frame(right)[, namescov]
covs <- scale(covs)
colnames(covs) <- newnames
right <- cbind(right, covs)

## FIT MODELS
plan(multisession, workers = 8)

system.time(
right_models <- future_map(models, function(x) {
    m <- bru(as.formula(paste("Right ~", x)[2]),
        data = right,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T, control.compute = list(cpo = T))
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
)

plan(sequential)
save(right_models, file = "./Models/Culture/Right.RData")

## ----------------- ##
## BURIAL SIDE: BACK ##
## ----------------- ##

# Data
back <- pnts[, c("YearBP", "Back", "mobility", "EHG", "WHG", "LVN", "CHG", "Culture", "geometry")]
back <- back[!is.na(back$Back),]
nrow(back)

# Jitter the duplicates
duplicates <- back[duplicated(back$geometry),]
duplicates <- st_jitter(duplicates, 0.05)

plot(conts_prj$geometry, xlim = ext(back)[1:2], ylim = ext(back)[3:4]) 
points(back, pch = 1, cex = 1); points(duplicates, pch = 1, col = 'red', cex = 1)

# Add duplicates to the back data
back <- back[!duplicated(back$geometry),]
back <- rbind(back, duplicates)
nrow(back) # same as before

# Divide into time bins
back$YearBin <- plyr:::round_any(back$YearBP, 250)

## Define covariates
namescov <- c("mobility", "EHG", "WHG", "LVN", "CHG")
newnames <- paste0(namescov, "_scaled")
covs <- data.frame(back)[, namescov]
covs <- scale(covs)
colnames(covs) <- newnames
back <- cbind(back, covs)

## FIT MODELS
plan(multisession, workers = 8)

system.time(
back_models <- future_map(models, function(x) {
    m <- bru(as.formula(paste("Back ~", x)[2]),
        data = back,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T, control.compute = list(cpo = T))
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
)

plan(sequential)
save(back_models, file = "./Models/Culture/Back.RData")

## --------------- ##
## EXTRACT RESULTS ##
## --------------- ##
model_files <- list.files("./Models/Culture/", full.names = T)
for(f in model_files){
    load(f)
}

model_info <- data.frame(
    model = paste0("mc", 1:8),
    name = c(
        "baseline",
        "baseline + mobility + ancestry + culture",
        "baseline + mobility + ancestry",
        "baseline + ancestry + culture",
        "baseline + mobility + culture",
        "baseline + ancestry",
        "baseline + mobility",
        "baseline + culture"
    )
)

## MODEL FIT
left_fit <- list(
    waic = sapply(left_models, function(x) x$waic$waic),
    dic = sapply(left_models, function(x) x$dic$dic),
    mlik = sapply(left_models, function(x) x$mlik[2, 1])
)

right_fit <- list(
    waic = sapply(right_models, function(x) x$waic$waic),
    dic = sapply(right_models, function(x) x$dic$dic),
    mlik = sapply(right_models, function(x) x$mlik[2, 1])
)

back_fit <- list(
    waic = sapply(back_models, function(x) x$waic$waic),
    dic = sapply(back_models, function(x) x$dic$dic),
    mlik = sapply(back_models, function(x) x$mlik[2, 1])
)

fit <- list(left_fit, right_fit, back_fit)
vars <- c("Left", "Right", "Back")
model_fit <- rbindlist(lapply(1:length(fit), function(x){
    d <- as.data.table(do.call(cbind, fit[[x]]))
    d[, Model := model_info$name]
    d[, Variable := vars[x]]
    return(d)
}))

fwrite(model_fit, "./Results/Culture/Models_fit.csv")

## Plot model fit
model_fit$Model <- factor(model_fit$Model, levels = unique(model_fit$Model))
model_fit[, deltaWAIC := waic - min(waic), by = "Variable"]
model_fit[, colorWAIC := ifelse(deltaWAIC == 0, "red", "black")]
model_fit[, deltaDIC := dic - min(dic), by = "Variable"]
model_fit[, colorDIC := ifelse(deltaDIC == 0, "red", "black")]
model_fit$Variable <- factor(model_fit$Variable, levels = unique(model_fit$Variable))

ggplot(data = model_fit) +
    geom_line(aes(x = Model, y = deltaWAIC, group = 1)) +
    geom_hline(
        data = model_fit[Model == "baseline"], aes(yintercept = deltaWAIC),
        linetype = 2, color = "black"
    ) +
    geom_point(aes(x = Model, y = deltaWAIC, fill = colorWAIC), shape = 21, size = 4) +
    facet_grid(rows = vars(Variable), scales = "free_y") +
    scale_fill_identity() +
    ylab("Delta WAIC") +
    ggtitle("Whole period") +
    labs(tag = "A") +
    theme_bw() +
    theme(
        strip.text = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        aspect.ratio = 0.2,
        plot.tag = element_text(face = 'bold', size = 22)
    )

## Extract posterior densities of covariates
full_models <- list(left_models[[2]], right_models[[2]], back_models[[2]])

vf <- lapply(full_models, function(x) x$summary.fixed)
names(vf) <- c("Left", "Right", "Back")
vf <- lapply(1:length(vf), function(x){
    vf[[x]]$Variable <- names(vf)[x]
    vf[[x]]$Predictors <- rownames(vf[[x]])
    return(vf[[x]])
})

vr <- lapply(full_models, function(x) x$summary.random$culture)
names(vr) <- c("Left", "Right", "Back")
vr <- lapply(1:length(vr), function(x){
    vr[[x]]$Variable <- names(vr)[x]
    # vr[[x]]$Predictors <- vr[[x]]$ID
    return(vr[[x]])
})
vr <- rbindlist(vr)
names(vr)[1] <- "Predictors"

results_df <- rbind(rbindlist(vf), vr)
results_df$Type <- sapply(1:nrow(results_df), function(x) {
    ifelse(results_df$Predictors[x] %in% c("EHG_scaled", "CHG_scaled", "LVN_scaled"), "Ancestry",
        ifelse(results_df$Predictors[x] %in% unique(DT$Culture), "Culture",
            ifelse(results_df$Predictors[x] == "Intercept", "Intercept", "Mobility")
        )
    )
})
results_df$Predictors <- factor(results_df$Predictors,
    levels = rev(c("Intercept", "EHG_scaled", "CHG_scaled", "LVN_scaled", "mobility_scaled", levels(DT$Culture))),
    labels = rev(c("Intercept", "EHG", "CHG", "LVN", "Mobility", levels(DT$Culture)))
)
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = c("Intercept", "Ancestry", "Mobility", "Culture"))
results_df$Effect <- ifelse(results_df$Type == "Culture", "Random", "Fixed")
results_df$Effect <- as.factor(results_df$Effect)

# Write to file
fwrite(results_df, "./Results/Culture/FullModel_results.csv")

## Plot
cols <- rev(pals::ocean.phase(6)[2:5])
cols[1] <- "grey50"
ggplot(data = results_df[Effect == "Fixed"]) +
    geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
    width = 0.2, linewidth = 1) +
    geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
    geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
    geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
    scale_color_manual("", values = cols) +
    scale_fill_manual("", values = cols) +
    facet_grid(cols = vars(Variable), scales = "free_x") +
    theme_bw() +
    ggtitle("Full model") +
    theme(aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = 'bold'),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16))

ggplot(data = results_df[Effect == "Random"]) +
    geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
    width = 0.2, linewidth = 1) +
    geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
    geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
    geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
    scale_color_manual("", values = cols) +
    scale_fill_manual("", values = cols) +
    facet_grid(cols = vars(Variable), scales = "free_x") +
    theme_bw() +
    ggtitle("Full model") +
    theme(aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = 'bold'),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16))

## Best models
models <- list(left_models, right_models, back_models)
idx_model <- sapply(model_fit[deltaWAIC == 0]$Model, function(x) model_info[which(model_info$name == x), "model"])

best_models <- lapply(1:length(models), function(i) models[[i]][[idx_model[i]]])

vf <- lapply(best_models, function(x) x$summary.fixed)
names(vf) <- c("Left", "Right", "Back")
vf <- lapply(1:length(vf), function(x){
    vf[[x]]$Variable <- names(vf)[x]
    vf[[x]]$Predictors <- rownames(vf[[x]])
    return(vf[[x]])
})

vr <- lapply(best_models, function(x) x$summary.random$culture)
names(vr) <- c("Left", "Right", "Back")
vr <- lapply(1:length(vr), function(x){
    vr[[x]]$Variable <- names(vr)[x]
    # vr[[x]]$Predictors <- vr[[x]]$ID
    return(vr[[x]])
})
vr <- rbindlist(vr, fill = T)
names(vr)[1] <- "Predictors"

results_df <- rbind(rbindlist(vf), vr)
results_df <- results_df[!is.na(Predictors)]
results_df$Type <- sapply(1:nrow(results_df), function(x) {
    ifelse(results_df$Predictors[x] %in% c("EHG_scaled", "CHG_scaled", "LVN_scaled"), "Ancestry",
        ifelse(results_df$Predictors[x] %in% unique(DT$Culture), "Culture",
            ifelse(results_df$Predictors[x] == "Intercept", "Intercept", "Mobility")
        )
    )
})
results_df$Predictors <- factor(results_df$Predictors,
    levels = rev(c("Intercept", "EHG_scaled", "CHG_scaled", "LVN_scaled", "mobility_scaled", levels(DT$Culture))),
    labels = rev(c("Intercept", "EHG", "CHG", "LVN", "Mobility", levels(DT$Culture)))
)
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = c("Intercept", "Ancestry", "Mobility", "Culture"))
results_df$Effect <- ifelse(results_df$Type == "Culture", "Random", "Fixed")
results_df$Effect <- as.factor(results_df$Effect)

# Write to file
fwrite(results_df, "./Results/Culture/BestModel_results.csv")

ggplot(data = results_df) +
    geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
    width = 0.2, linewidth = 1) +
    geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
    geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
    geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
    scale_color_manual("", values = cols) +
    scale_fill_manual("", values = cols) +
    facet_grid(cols = vars(Variable), scales = "free_x") +
    theme_bw() +
    ggtitle("Full model") +
    theme(aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = 'bold'),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16))

## ---------------------------------------------------------------------------------------------------------------------------

ggplot() +
    geom_sf(data = land_prj, inherit.aes = T, fill = "grey90") +
    geom_sf(
        data = st_as_sf(pnts_new),
        inherit.aes = T,
        shape = 21,
        fill = "grey50",
        size = 4,
        alpha = 0.8
    ) +
    coord_sf(
        xlim = ext(left)[1:2],
        ylim = ext(left)[3:4]
    ) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
    axis.title = element_blank(),
    axis.text = element_blank())

pnts.dt <- melt(DT, measure.vars = c("Left", "Right", "Back"))
pnts.dt <- pnts.dt[!is.na(value)]
ggplot(data = pnts.dt) +
    geom_density(
        data = pnts.dt[value == 1], aes(x = YearBP),
        fill = "grey80",
        color = "transparent"
    ) +
    geom_point(aes(
        x = YearBP, y = 0,
        color = as.factor(value)
    ), shape = 3, size = 3, alpha = 0.8) +
    scale_color_manual(values = rev(paletteer::paletteer_d("nbapalettes::pacers_classic"))) +
    scale_x_reverse() +
    facet_grid(rows = vars(variable)) +
    xlab("Years (cal BP)") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "none",
        aspect.ratio = 0.2
    )
