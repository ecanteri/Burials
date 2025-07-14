## ----- ##
## SETUP ##
## ----- ##
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
setwd("/maps/projects/racimolab/people/msb290/Graves/Scripts/")
proj <- "+proj=laea +lon_0=0 +lat_0=49.06 +datum=WGS84 +units=km +no_defs"
proj2 <- "+proj=aea +lon_0=44.296875 +lat_1=43.7864128 +lat_2=69.9512657 +lat_0=56.8688392 +datum=WGS84 +units=m +no_defs"

## ---- ##
## DATA ##
## ---- ##
graves <- fread("../Data/Burial.csv", na.strings = "")
graves$DepositionType <- as.factor(graves$DepositionType)
graves$BodyPositioning <- as.factor(graves$BodyPositioning)
graves$BurialSide <- as.factor(graves$BurialSide)

## Coastlines and continents poly
conts <- rnaturalearth::ne_download(scale = 10, type = "coastline", "physical", returnclass = "sf")
land <- rnaturalearth::ne_download(scale = 10, type = "land", "physical", returnclass = "sf")
countries <- rnaturalearth::ne_download(scale = 50, type = "countries", "cultural", returnclass = "sf")
grat <- rnaturalearth::ne_download(scale = 10, type = "graticules_10", "physical", returnclass = "sf")
conts_prj <- st_transform(conts, crs = proj)
land_prj <- st_transform(land, proj)
countries_prj <- st_transform(countries, crs = proj)
grat_prj <- st_transform(grat, proj)

# Create spatial points from data
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

poly <- vect("../Spatial/IntersectionVector.shp")
pnts <- intersect(pnts, poly)

# Remove data without age and too old sample
pnts <- st_as_sf(pnts)
pnts <- pnts[!is.na(pnts$YearBP),]
pnts <- pnts[which(pnts$YearBP < 15000),]

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
cover <- st_read("../Spatial/cover.shp")
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
    m1 = ~ Intercept(1) +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    m2 = ~ Intercept(1) +
        mobility_scaled +
        EHG_scaled +
        WHG_scaled +
        LVN_scaled +
        CHG_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    m3 = ~ Intercept(1) +
        EHG_scaled +
        WHG_scaled +
        LVN_scaled +
        CHG_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    m4 = ~ Intercept(1) +
        EHG_scaled +
        LVN_scaled +
        CHG_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    m5 = ~ Intercept(1) +
        EHG_scaled +
        LVN_scaled +
        CHG_scaled +
        mobility_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    m6 = ~ Intercept(1) +
        EHG_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    m7 = ~ Intercept(1) +
        WHG_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    m8 = ~ Intercept(1) +
        CHG_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    m9 = ~ Intercept(1) +
        LVN_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    m10 = ~ Intercept(1) +
        mobility_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        )
)

models_5k <- c(
    mb1 = models[1],
    mb2 = models[4],
    mb3 = models[5]
)
names(models_5k) <- paste0("mb", 1:3)

models_steppe <- c(
    ms1 = ~ Intercept(1) +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    ms2 = ~ Intercept(1) +
        mobility_scaled +
        Steppe_scaled +
        LVN_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),    
    ms3 = ~ Intercept(1) +
        Steppe_scaled +
        LVN_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    ms4 = ~ Intercept(1) +
        Steppe_scaled +
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
left <- pnts[, c("YearBP", "Left", "mobility", "EHG", "WHG", "LVN", "CHG", "geometry")]
left <- left[!is.na(left$Left),]
nrow(left)

# Jitter the duplicates
duplicates <- left[duplicated(left),]
duplicates <- st_jitter(duplicates, 0.05)

plot(conts_prj$geometry, xlim = ext(left)[1:2], ylim = ext(left)[3:4]) 
points(left, pch = 1, cex = 1); points(duplicates, pch = 1, col = 'red', cex = 1)

# Add duplicates to the left data
left <- left[!duplicated(left),]
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
plan(multisession, workers = 5)
left_models <- future_map(models, function(x) {
    m <- bru(as.formula(paste("Left ~", x)[2]),
        data = left,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T, control.compute = list(cpo = T))
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)

plan(sequential)
save(left_models, file = "../Models/Left.RData")

## BEFORE STEPPE
plan(multisession, workers = 3)
left_bs <- future_map(models_5k, function(x) {
    m <- bru(as.formula(paste("Left ~", x)[2]),
        data = left[which(left$YearBP > 5000),], # data older than 5ka
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T, control.compute = list(cpo = T))
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
plan(sequential)
save(left_bs, file = "../Models/Left_BeforeSteppe.RData")

## AFTER STEPPE
plan(multisession, workers = 4)
left_steppe_dt <- left[which(left$YearBP <= 5000),] # data younger than 5ka
left_steppe_dt$Steppe_scaled <- left_steppe_dt$EHG_scaled + left_steppe_dt$CHG_scaled # create steppe ancestry

left_steppe <- future_map(models_steppe, function(x) {
    m <- bru(as.formula(paste("Left ~", x)[2]),
        data = left_steppe_dt,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T, control.compute = list(cpo = T))
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
plan(sequential)
save(left_steppe, file = "../Models/Left_Steppe.RData")

## ------------------ ##
## BURIAL SIDE: RIGHT ##
## ------------------ ##

# Data
right <- pnts[, c("YearBP", "Right", "mobility", "EHG", "WHG", "LVN", "CHG", "geometry")]
right <- right[!is.na(right$Right),]
nrow(right)

# Jitter the duplicates
duplicates <- right[duplicated(right),]
duplicates <- st_jitter(duplicates, 0.05)

plot(conts_prj$geometry, xlim = ext(right)[1:2], ylim = ext(right)[3:4]) 
points(right, pch = 1, cex = 1); points(duplicates, pch = 1, col = 'red', cex = 1)

# Add duplicates to the right data
right <- right[!duplicated(right),]
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
plan(multisession, workers = 5)

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
save(right_models, file = "../Models/Right.RData")

## BEFORE STEPPE
plan(multisession, workers = 3)
right_bs <- future_map(models_5k, function(x) {
    m <- bru(as.formula(paste("Right ~", x)[2]),
        data = right[which(right$YearBP > 5000),],
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T, control.compute = list(cpo = T))
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
plan(sequential)
save(right_bs, file = "../Models/Right_BeforeSteppe.RData")

## AFTER STEPPE
plan(multisession, workers = 4)
right_steppe_dt <- right[which(right$YearBP <= 5000),]
right_steppe_dt$Steppe_scaled <- right_steppe_dt$EHG_scaled + right_steppe_dt$CHG_scaled

right_steppe <- future_map(models_steppe, function(x) {
    m <- bru(as.formula(paste("Right ~", x)[2]),
        data = right_steppe_dt,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T, control.compute = list(cpo = T))
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
plan(sequential)
save(right_steppe, file = "../Models/Right_Steppe.RData")

## ----------------- ##
## BURIAL SIDE: BACK ##
## ----------------- ##

# Data
back <- pnts[, c("YearBP", "Back", "mobility", "EHG", "WHG", "LVN", "CHG", "geometry")]
back <- back[!is.na(back$Back),]
nrow(back)

# Jitter the duplicates
duplicates <- back[duplicated(back),]
duplicates <- st_jitter(duplicates, 0.05)

plot(conts_prj$geometry, xlim = ext(back)[1:2], ylim = ext(back)[3:4]) 
points(back, pch = 1, cex = 1); points(duplicates, pch = 1, col = 'red', cex = 1)

# Add duplicates to the back data
back <- back[!duplicated(back),]
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
plan(multisession, workers = 5)

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
save(back_models, file = "../Models/Back.RData")

## BEFORE STEPPE
plan(multisession, workers = 3)
back_bs <- future_map(models_5k, function(x) {
    m <- bru(as.formula(paste("Back ~", x)[2]),
        data = back[which(back$YearBP > 5000),],
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T, control.compute = list(cpo = T))
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
plan(sequential)
save(back_bs, file = "../Models/Back_BeforeSteppe.RData")

## AFTER STEPPE
plan(multisession, workers = 4)
back_steppe_dt <- back[which(back$YearBP <= 5000),]
back_steppe_dt$Steppe_scaled <- back_steppe_dt$EHG_scaled + back_steppe_dt$CHG_scaled

back_steppe <- future_map(models_steppe, function(x) {
    m <- bru(as.formula(paste("Back ~", x)[2]),
        data = back_steppe_dt,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T, control.compute = list(cpo = T))
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
plan(sequential)
save(back_steppe, file = "../Models/Back_Steppe.RData")

## --------------- ##
## EXTRACT RESULTS ##
## --------------- ##
model_files <- list.files("../Models/", full.names = T)
for(f in model_files){
    load(f)
}

model_info <- data.frame(
    model = paste0("m", 1:10),
    name = c(
        "baseline",
        "baseline + mobility + ancestry",
        "baseline + ancestry",
        "baseline + no WHG ancestry",
        "baseline + no WHG ancestry + mobility",
        "baseline + EHG",
        "baseline + WHG",
        "baseline + CHG",
        "baseline + LVN",
        "baseline + mobility"
    )
)

model_5k_info <- data.frame(
    model = paste0("mb", 1:3),
    name = c(
        "baseline",
        "baseline + no WHG ancestry",
        "baseline + no WHG ancestry + mobility"
    )
)

model_steppe_info <- data.frame(
    model = paste0("ms", 1:4),
    name = c(
        "baseline",
        "baseline + mobility + steppe ancestry + LVN",
        "baseline + steppe ancestry + LVN",
        "baseline + steppe ancestry"
    )
)

## --------------------------------------------------------
## MODEL FIT
# Full period
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

fwrite(model_fit, "../Results/Models_fit.csv")

# Before steppe
left_5k_fit <- list(
    waic = sapply(left_bs, function(x) x$waic$waic),
    dic = sapply(left_bs, function(x) x$dic$dic),
    mlik = sapply(left_bs, function(x) x$mlik[2, 1])
)

right_5k_fit <- list(
    waic = sapply(right_bs, function(x) x$waic$waic),
    dic = sapply(right_bs, function(x) x$dic$dic),
    mlik = sapply(right_bs, function(x) x$mlik[2, 1])
)

back_5k_fit <- list(
    waic = sapply(back_bs, function(x) x$waic$waic),
    dic = sapply(back_bs, function(x) x$dic$dic),
    mlik = sapply(back_bs, function(x) x$mlik[2, 1])
)

fit <- list(left_5k_fit, right_5k_fit, back_5k_fit)
vars <- c("Left", "Right", "Back")
model_5k_fit <- rbindlist(lapply(1:length(fit), function(x){
    d <- as.data.table(do.call(cbind, fit[[x]]))
    d[, Model := model_5k_info$name]
    d[, Variable := vars[x]]
    return(d)
}))

fwrite(model_5k_fit, "../Results/ModelsFit_BeforeSteppe.csv")

# Steppe ancestry
left_steppe_fit <- list(
    waic = sapply(left_steppe, function(x) x$waic$waic),
    dic = sapply(left_steppe, function(x) x$dic$dic),
    mlik = sapply(left_steppe, function(x) x$mlik[2, 1])
)

right_steppe_fit <- list(
    waic = sapply(right_steppe, function(x) x$waic$waic),
    dic = sapply(right_steppe, function(x) x$dic$dic),
    mlik = sapply(right_steppe, function(x) x$mlik[2, 1])
)

back_steppe_fit <- list(
    waic = sapply(back_steppe, function(x) x$waic$waic),
    dic = sapply(back_steppe, function(x) x$dic$dic),
    mlik = sapply(back_steppe, function(x) x$mlik[2, 1])
)

fit <- list(left_steppe_fit, right_steppe_fit, back_steppe_fit)
vars <- c("Left", "Right", "Back")
model_steppe_fit <- rbindlist(lapply(1:length(fit), function(x){
    d <- as.data.table(do.call(cbind, fit[[x]]))
    d[, Model := model_steppe_info$name]
    d[, Variable := vars[x]]
    return(d)
}))

fwrite(model_steppe_fit, "../Results/ModelsFit_Steppe.csv")

## -------------------------------------------------------------------------------------------
## Plot model fit
# Whole period
model_fit$Model <- factor(model_fit$Model, levels = unique(model_fit$Model))
model_fit[, deltaWAIC := waic - min(waic), by = "Variable"]
model_fit[, colorWAIC := ifelse(deltaWAIC == 0, "red", "black")]
model_fit[, deltaDIC := dic - min(dic), by = "Variable"]
model_fit[, colorDIC := ifelse(deltaDIC == 0, "red", "black")]
model_fit$Variable <- factor(model_fit$Variable, levels = unique(model_fit$Variable))

p1 <- ggplot(data = model_fit) +
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
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        strip.text = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 16),
        aspect.ratio = 0.2,
        plot.tag = element_text(face = 'bold', size = 22)
    )

p2 <- ggplot(data = model_fit) +
    geom_line(aes(x = Model, y = deltaDIC, group = 1)) +
    geom_hline(
        data = model_fit[Model == "baseline"], aes(yintercept = deltaDIC),
        linetype = 2, color = "black"
    ) +
    geom_point(aes(x = Model, y = deltaDIC, fill = colorDIC), shape = 21, size = 4) +
    facet_grid(rows = vars(Variable), scales = "free_y") +
    scale_fill_identity() +
    ylab("Delta DIC") +
    labs(tag = "B") +
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

png("../Figures/ModelsFit.png", width = 10, height = 15, res = 330, units = 'in')
p1/p2
dev.off()

# Before steppe
model_5k_fit$Model <- factor(model_5k_fit$Model, levels = unique(model_5k_fit$Model))
model_5k_fit[, deltaWAIC := waic - min(waic), by = "Variable"]
model_5k_fit[, colorWAIC := ifelse(deltaWAIC == 0, "red", "black")]
model_5k_fit[, deltaDIC := dic - min(dic), by = "Variable"]
model_5k_fit[, colorDIC := ifelse(deltaDIC == 0, "red", "black")]
model_5k_fit$Variable <- factor(model_5k_fit$Variable, levels = unique(model_5k_fit$Variable))

p3 <- ggplot(data = model_5k_fit) +
    geom_line(aes(x = Model, y = deltaWAIC, group = 1)) +
    geom_hline(
        data = model_5k_fit[Model == "baseline"], aes(yintercept = deltaWAIC),
        linetype = 2, color = "black"
    ) +
    geom_point(aes(x = Model, y = deltaWAIC, fill = colorWAIC), shape = 21, size = 4) +
    facet_grid(rows = vars(Variable), scales = "free_y") +
    scale_fill_identity() +
    ylab("Delta WAIC") +
    ggtitle("Before 5k BP") +
    labs(tag = "A") +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        strip.text = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 16),
        aspect.ratio = 0.2,
        plot.tag = element_text(face = 'bold', size = 22)
    )

p4 <- ggplot(data = model_5k_fit) +
    geom_line(aes(x = Model, y = deltaDIC, group = 1)) +
    geom_hline(
        data = model_5k_fit[Model == "baseline"], aes(yintercept = deltaDIC),
        linetype = 2, color = "black"
    ) +
    geom_point(aes(x = Model, y = deltaDIC, fill = colorDIC), shape = 21, size = 4) +
    facet_grid(rows = vars(Variable), scales = "free_y") +
    scale_fill_identity() +
    ylab("Delta DIC") +
    labs(tag = "B") +
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

png("../Figures/ModelsFit_BeforeSteppe.png", width = 10, height = 15, res = 330, units = 'in')
p3/p4
dev.off()

# After steppe
model_steppe_fit$Model <- factor(model_steppe_fit$Model, levels = unique(model_steppe_fit$Model))
model_steppe_fit[, deltaWAIC := waic - min(waic), by = "Variable"]
model_steppe_fit[, colorWAIC := ifelse(deltaWAIC == 0, "red", "black")]
model_steppe_fit[, deltaDIC := dic - min(dic), by = "Variable"]
model_steppe_fit[, colorDIC := ifelse(deltaDIC == 0, "red", "black")]
model_steppe_fit$Variable <- factor(model_steppe_fit$Variable, levels = unique(model_steppe_fit$Variable))

p5 <- ggplot(data = model_steppe_fit) +
    geom_line(aes(x = Model, y = deltaWAIC, group = 1)) +
    geom_hline(
        data = model_steppe_fit[Model == "baseline"], aes(yintercept = deltaWAIC),
        linetype = 2, color = "black"
    ) +
    geom_point(aes(x = Model, y = deltaWAIC, fill = colorWAIC), shape = 21, size = 4) +
    facet_grid(rows = vars(Variable), scales = "free_y") +
    scale_fill_identity() +
    ylab("Delta WAIC") +
    ggtitle("Steppe ancestry") +
    labs(tag = "A") +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        strip.text = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 16),
        aspect.ratio = 0.2,
        plot.tag = element_text(face = 'bold', size = 22)
    )

p6 <- ggplot(data = model_steppe_fit) +
    geom_line(aes(x = Model, y = deltaDIC, group = 1)) +
    geom_hline(
        data = model_steppe_fit[Model == "baseline"], aes(yintercept = deltaDIC),
        linetype = 2, color = "black"
    ) +
    geom_point(aes(x = Model, y = deltaDIC, fill = colorDIC), shape = 21, size = 4) +
    facet_grid(rows = vars(Variable), scales = "free_y") +
    scale_fill_identity() +
    ylab("Delta DIC") +
    labs(tag = "B") +
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

png("../Figures/ModelsFit_AfterSteppe.png", width = 10, height = 15, res = 330, units = 'in')
p5/p6
dev.off()

## ------------------------------------------------------------------------------------------------
## Best models
cols <- rev(pals::ocean.phase(5)[2:4])
cols[1] <- "grey50"

models <- list(left_models, right_models, back_models)
idx_model <- sapply(model_fit[deltaWAIC == 0]$Model, function(x) model_info[which(model_info$name == x), "model"])

best_models <- lapply(1:length(models), function(i) models[[i]][[idx_model[i]]])

r <- lapply(best_models, function(x) x$summary.fixed)
names(r) <- c("Left", "Right", "Back")
r <- lapply(1:length(r), function(x){
    r[[x]]$Variable <- names(r)[x]
    r[[x]]$Predictors <- rownames(r[[x]])
    return(r[[x]])
})

results_df <- rbindlist(r)
results_df$Type <- sapply(1:nrow(results_df), function(x){
    ifelse(results_df$Predictors[x] %in% c("WHG_scaled", "EHG_scaled", "CHG_scaled", "LVN_scaled"), "Ancestry",
    ifelse(results_df$Predictors[x] == "Intercept", "Intercept", "Mobility"))
})
results_df$Predictors <- factor(results_df$Predictors,
    levels = rev(c("Intercept", "WHG_scaled", "EHG_scaled", "CHG_scaled", "LVN_scaled")),
    labels = rev(c("Intercept", "WHG", "EHG", "CHG", "LVN"))
)
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = unique(results_df$Type))

pr1 <- ggplot(data = results_df) +
    geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
    width = 0.2, linewidth = 1) +
    geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
    geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
    geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
    scale_color_manual("", values = cols) +
    scale_fill_manual("", values = cols) +
    facet_grid(cols = vars(Variable), scales = "free_x") +
    theme_bw() +
    ggtitle("Whole period") +
    theme(aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = 'bold'),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16))

## Before steppe
models_5k <- list(left_bs, right_bs, back_bs)
idx_model <- sapply(model_5k_fit[deltaWAIC == 0]$Model, function(x) model_5k_info[which(model_5k_info$name == x), "model"])

best_models_5k <- lapply(1:length(models_5k), function(i) models_5k[[i]][[idx_model[i]]])

r <- lapply(best_models_5k, function(x) x$summary.fixed)
names(r) <- c("Left", "Right", "Back")
r <- lapply(1:length(r), function(x){
    r[[x]]$Variable <- names(r)[x]
    r[[x]]$Predictors <- rownames(r[[x]])
    return(r[[x]])
})

results_df <- rbindlist(r)
results_df$Type <- sapply(1:nrow(results_df), function(x){
    ifelse(results_df$Predictors[x] %in% c("WHG_scaled", "EHG_scaled", "CHG_scaled", "LVN_scaled"), "Ancestry",
    ifelse(results_df$Predictors[x] == "Intercept", "Intercept", "Mobility"))
})
results_df$Predictors <- factor(results_df$Predictors,
    levels = rev(c("Intercept", "EHG_scaled", "CHG_scaled", "LVN_scaled", "mobility_scaled")),
    labels = rev(c("Intercept", "EHG", "CHG", "LVN", "Mobility"))
)
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = unique(results_df$Type))

pr2 <- ggplot(data = results_df) +
    geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
    width = 0.2, linewidth = 1) +
    geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
    geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
    geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
    scale_color_manual("", values = cols) +
    scale_fill_manual("", values = cols) +
    facet_grid(cols = vars(Variable), scales = "free_x") +
    theme_bw() +
    ggtitle("Before 5k BP") +
    theme(aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = 'bold'),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16))

## After steppe
models_steppe <- list(left_steppe, right_steppe, back_steppe)
idx_model <- sapply(model_steppe_fit[deltaWAIC == 0]$Model, function(x) model_steppe_info[which(model_steppe_info$name == x), "model"])

best_models_steppe <- lapply(1:length(models_steppe), function(i) models_steppe[[i]][[idx_model[i]]])

r <- lapply(best_models_steppe, function(x) x$summary.fixed)
names(r) <- c("Left", "Right", "Back")
r <- lapply(1:length(r), function(x){
    r[[x]]$Variable <- names(r)[x]
    r[[x]]$Predictors <- rownames(r[[x]])
    return(r[[x]])
})

results_df <- rbindlist(r)
results_df$Type <- sapply(1:nrow(results_df), function(x){
    ifelse(results_df$Predictors[x] %in% c("WHG_scaled", "EHG_scaled", "CHG_scaled", "LVN_scaled", "Steppe_scaled"), "Ancestry",
    ifelse(results_df$Predictors[x] == "Intercept", "Intercept", "Mobility"))
})
results_df$Predictors <- factor(results_df$Predictors,
    levels = rev(c("Intercept", "LVN_scaled", "Steppe_scaled")),
    labels = rev(c("Intercept", "LVN", "Steppe ancestry"))
)
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = unique(results_df$Type))

pr3 <- ggplot(data = results_df) +
    geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
    width = 0.2, linewidth = 1) +
    geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
    geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
    geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
    scale_color_manual("", values = cols) +
    scale_fill_manual("", values = cols) +
    facet_grid(cols = vars(Variable), scales = "free_x") +
    theme_bw() +
    ggtitle("After 5k BP") +
    theme(aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = 'bold'),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16))

png("../Figures/Models_Results.png", width = 10, height = 11, res = 330, units = 'in')
pr1/pr2/pr3
dev.off()

## ------------------------------------------------------------------------------------------------
## Full model
cols <- rev(pals::ocean.phase(5)[2:4])
cols[1] <- "grey50"
full_models <- list(left_models[[2]], right_models[[2]], back_models[[2]])

r <- lapply(full_models, function(x) x$summary.fixed)
names(r) <- c("Left", "Right", "Back")
r <- lapply(1:length(r), function(x){
    r[[x]]$Variable <- names(r)[x]
    r[[x]]$Predictors <- rownames(r[[x]])
    return(r[[x]])
})

results_df <- rbindlist(r)
results_df$Type <- sapply(1:nrow(results_df), function(x){
    ifelse(results_df$Predictors[x] %in% c("WHG_scaled", "EHG_scaled", "CHG_scaled", "LVN_scaled"), "Ancestry",
    ifelse(results_df$Predictors[x] == "Intercept", "Intercept", "Mobility"))
})
results_df$Predictors <- factor(results_df$Predictors,
    levels = rev(c("Intercept", "WHG_scaled", "EHG_scaled", "CHG_scaled", "LVN_scaled", "mobility_scaled")),
    labels = rev(c("Intercept", "WHG", "EHG", "CHG", "LVN", "Mobility"))
)
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = c("Intercept", "Ancestry", "Mobility"))

p_full <- ggplot(data = results_df) +
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

png("../Figures/FullModels_Results.png", width = 10, height = 6, res = 330, units = 'in')
p_full
dev.off()


## Full model no WHG
noWHG_models <- list(left_models[[5]], right_models[[5]], back_models[[5]])

r <- lapply(noWHG_models, function(x) x$summary.fixed)
names(r) <- c("Left", "Right", "Back")
r <- lapply(1:length(r), function(x){
    r[[x]]$Variable <- names(r)[x]
    r[[x]]$Predictors <- rownames(r[[x]])
    return(r[[x]])
})

results_df <- rbindlist(r)
results_df$Type <- sapply(1:nrow(results_df), function(x){
    ifelse(results_df$Predictors[x] %in% c("WHG_scaled", "EHG_scaled", "CHG_scaled", "LVN_scaled"), "Ancestry",
    ifelse(results_df$Predictors[x] == "Intercept", "Intercept", "Mobility"))
})
results_df$Predictors <- factor(results_df$Predictors,
    levels = rev(c("Intercept", "WHG_scaled", "EHG_scaled", "CHG_scaled", "LVN_scaled", "mobility_scaled")),
    labels = rev(c("Intercept", "WHG", "EHG", "CHG", "LVN", "Mobility"))
)
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = c("Intercept", "Ancestry", "Mobility"))

p_noWHG <- ggplot(data = results_df) +
    geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
    width = 0.2, linewidth = 1) +
    geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
    geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
    geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
    scale_color_manual("", values = cols) +
    scale_fill_manual("", values = cols) +
    facet_grid(cols = vars(Variable), scales = "free_x") +
    theme_bw() +
    ggtitle("Full model without WHG") +
    theme(aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = 'bold'),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16))

png("../Figures/noWHG_Results.png", width = 10, height = 6, res = 330, units = 'in')
p_noWHG
dev.off()

## Full model before 5k BP
models_before5k <- list(left_bs[[3]], right_bs[[3]], back_bs[[3]])

r <- lapply(models_before5k, function(x) x$summary.fixed)
names(r) <- c("Left", "Right", "Back")
r <- lapply(1:length(r), function(x){
    r[[x]]$Variable <- names(r)[x]
    r[[x]]$Predictors <- rownames(r[[x]])
    return(r[[x]])
})

results_df <- rbindlist(r)
results_df$Type <- sapply(1:nrow(results_df), function(x){
    ifelse(results_df$Predictors[x] %in% c("EHG_scaled", "CHG_scaled", "LVN_scaled"), "Ancestry",
    ifelse(results_df$Predictors[x] == "Intercept", "Intercept", "Mobility"))
})
results_df$Predictors <- factor(results_df$Predictors,
    levels = rev(c("Intercept", "EHG_scaled", "CHG_scaled", "LVN_scaled", "mobility_scaled")),
    labels = rev(c("Intercept", "EHG", "CHG", "LVN", "Mobility"))
)
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = c("Intercept", "Ancestry", "Mobility"))

p_5k <- ggplot(data = results_df) +
    geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
    width = 0.2, linewidth = 1) +
    geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
    geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
    geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
    scale_color_manual("", values = cols) +
    scale_fill_manual("", values = cols) +
    facet_grid(cols = vars(Variable), scales = "free_x") +
    theme_bw() +
    ggtitle("Full model before 5k BP") +
    theme(aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = 'bold'),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16))

png("../Figures/FullModel_before5k_Results.png", width = 10, height = 6, res = 330, units = 'in')
p_5k
dev.off()

## Full model after steppe
models_st <- list(left_steppe[[2]], right_steppe[[2]], back_steppe[[2]])

r <- lapply(models_st, function(x) x$summary.fixed)
names(r) <- c("Left", "Right", "Back")
r <- lapply(1:length(r), function(x){
    r[[x]]$Variable <- names(r)[x]
    r[[x]]$Predictors <- rownames(r[[x]])
    return(r[[x]])
})

results_df <- rbindlist(r)
results_df$Type <- sapply(1:nrow(results_df), function(x){
    ifelse(results_df$Predictors[x] %in% c("Steppe_scaled", "LVN_scaled"), "Ancestry",
    ifelse(results_df$Predictors[x] == "Intercept", "Intercept", "Mobility"))
})
results_df$Predictors <- factor(results_df$Predictors,
    levels = rev(c("Intercept", "Steppe_scaled", "LVN_scaled", "mobility_scaled")),
    labels = rev(c("Intercept", "Steppe", "LVN", "Mobility"))
)
results_df$Variable <- factor(results_df$Variable, levels = unique(results_df$Variable))
results_df$Type <- factor(results_df$Type, levels = c("Intercept", "Ancestry", "Mobility"))

p_steppe <- ggplot(data = results_df) +
    geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`, x = mean, y = Predictors, color = Type),
    width = 0.2, linewidth = 1) +
    geom_point(aes(x = mean, y = Predictors, fill = Type), shape = 21, size = 3) +
    geom_blank(aes(x = -mean, xmin = -`0.025quant`, xmax = -`0.975quant`)) +
    geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.5) +
    scale_color_manual("", values = cols) +
    scale_fill_manual("", values = cols) +
    facet_grid(cols = vars(Variable), scales = "free_x") +
    theme_bw() +
    ggtitle("Full model after 5k BP") +
    theme(aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = 'bold'),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16))

png("../Figures/FullModel_Steppe_Results.png", width = 10, height = 6, res = 330, units = 'in')
p_steppe
dev.off()

## ------------------------------------------------------------------------------------------------
## Predictions
# Predicting dataframe
pix <- fm_pixels(mesh, mask = cover)
pred.t <- seq(4000, 8000, 500)
df <- fm_cprod(pix, data.frame(YearBin = pred.t))

# Extract ancestry values
ance <- readRDS("/projects/racimolab/people/msb290/MesoNeo/SelectedAncestries_100ystep.RDS")
time <- raster::getZ(ance[[1]])
layers <- which(raster::getZ(ance[[1]]) %in% sort(pred.t*-1))
ance <- pblapply(ance, rast)
ance <- pblapply(ance, project, proj)
ance_pred <- lapply(ance, function(a){
    pblapply(layers, function(r){
        extract(a[[r]], st_coordinates(pix))
    })
})
ance_pred <- lapply(1:length(ance_pred), function(a){
    pblapply(1:length(ance_pred[[a]]), function(d){
        df <- cbind(st_coordinates(pix), ance_pred[[a]][[d]], data.frame(YearBin = sort(pred.t, decreasing = T)[d]))
        colnames(df)[3] <- names(ance_pred[a])
        return(as.data.table(df))
    })
})
ance_pred <- pblapply(ance_pred, rbindlist)

data_pred <- cbind(ance_pred[[1]][, c(1,2,4)], ance_pred[[1]][, 3], ance_pred[[2]][, 3], ance_pred[[3]][,3], ance_pred[[4]][,3])
data_pred <- st_as_sf(vect(data_pred, geom=c("X", "Y"), crs = proj))
colnames(data_pred)[2:5] <- newnames[-1]

## Best models
fixed <- lapply(best_models, function(x) x$names.fixed)
formulas <- lapply(fixed, function(f) as.formula(paste("~ 1/(1 + exp(-(", paste(c(f, "spt"), collapse = "+"), ")))")))
names(formulas) <- c("Left", "Right", "Back")

# Predict
plan(multisession, workers = 3)
pred_best <- future_map(1:length(best_models), function(i) {
    pr <- predict(
        best_models[[i]],
        data_pred,
        formulas[[i]]
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
plan(sequential)
names(pred_best) <- names(formulas)
save(pred_best, file = "../Results/BestModels_predictions.RData")

pred_df <- lapply(pred_best, function(x) cbind(as.data.frame(x), st_coordinates(x)))
plot_pred <- lapply(1:length(pred_df), function(d) {
    pl <- ggplot() +
        geom_tile(data = pred_df[[d]], aes(X, Y, fill = mean)) +
        geom_sf(data = cover, inherit.aes = F, fill = "transparent", alpha = 0.1) +
        scale_fill_gradientn("Mean",
            colors = rev(pals::brewer.rdylbu(100)),
            breaks = seq(0, 1, 0.2),
            limits = c(0, 1)
        ) +
        facet_wrap(~ -YearBin) +
        ggtitle(names(pred_df)[d]) +
        theme_bw() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_blank()
        )
    return(pl)
})
plot_pred[[1]]
plot_pred[[2]]
plot_pred[[3]]

## Full models
fixed <- lapply(full_models, function(x) x$names.fixed)
formulas <- lapply(fixed, function(f) as.formula(paste("~ 1/(1 + exp(-(", paste(c(f, "spt"), collapse = "+"), ")))")))
names(formulas) <- c("Left", "Right", "Back")

plan(multisession, workers = 3)
pred_full <- future_map(1:length(full_models), function(i) {
    pr <- predict(
        full_models[[i]],
        data_pred,
        formulas[[i]]
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
plan(sequential)
names(pred_full) <- names(formulas)
save(pred_full, file = "../Results/FullModels_predictions.RData")

pred_df <- lapply(pred_full, function(x) cbind(as.data.frame(x), st_coordinates(x)))
plot_pred <- lapply(1:length(pred_df), function(d) {
    pl <- ggplot() +
        geom_tile(data = pred_df[[d]], aes(X, Y, fill = mean)) +
        geom_sf(data = cover, inherit.aes = F, fill = "transparent", alpha = 0.1) +
        scale_fill_gradientn("Mean",
            colors = rev(pals::brewer.rdylbu(100)),
            breaks = seq(0, 1, 0.2),
            limits = c(0, 1)
        ) +
        facet_wrap(~ -YearBin) +
        ggtitle(names(pred_df)[d]) +
        theme_bw() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_blank()
        )
    return(pl)
})
plot_pred[[1]]
plot_pred[[2]]
plot_pred[[3]]

## No WHG models
fixed <- lapply(noWHG_models, function(x) x$names.fixed)
formulas <- lapply(fixed, function(f) as.formula(paste("~ 1/(1 + exp(-(", paste(c(f, "spt"), collapse = "+"), ")))")))
names(formulas) <- c("Left", "Right", "Back")

plan(multisession, workers = 3)
pred_noWHG <- future_map(1:length(noWHG_models), function(i) {
    pr <- predict(
        noWHG_models[[i]],
        data_pred,
        formulas[[i]]
    )
}, .options = furrr_options(seed = 123), .progress = TRUE)
plan(sequential)
names(pred_noWHG) <- names(formulas)
save(pred_noWHG, file = "../Results/noWHGModels_predictions.RData")

pred_df <- lapply(pred_noWHG, function(x) cbind(as.data.frame(x), st_coordinates(x)))
plot_pred <- lapply(1:length(pred_df), function(d) {
    pl <- ggplot() +
        geom_tile(data = pred_df[[d]], aes(X, Y, fill = mean)) +
        geom_sf(data = cover, inherit.aes = F, fill = "transparent", alpha = 0.1) +
        scale_fill_gradientn("Mean",
            colors = rev(pals::brewer.rdylbu(100)),
            breaks = seq(0, 1, 0.2),
            limits = c(0, 1)
        ) +
        facet_wrap(~ -YearBin) +
        ggtitle(names(pred_df)[d]) +
        theme_bw() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_blank()
        )
    return(pl)
})
plot_pred[[1]]
plot_pred[[2]]
plot_pred[[3]]

## ---------------------------------------------------------------------------------------------
## Check model performance - Permutation.R file