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
components <- c(
    c1 = ~ Intercept(1) +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        ),
    c2 = ~ Intercept(1) +
        EHG_scaled +
        LVN_scaled +
        CHG_scaled +
        mobility_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        )
)

## Model info
models_info <- data.frame(
    model = paste0("c", 1:2),
    name = c(
        "baseline",
        "baseline + no WHG ancestry + mobility"
    )
)

## Prediction formulas
formulas <- list(
    f1 = as.formula("~ 1 / (1 + exp(-(Intercept + spt)))"),
    f2 = as.formula("~ 1 / (1 + exp(-(Intercept + EHG_scaled + LVN_scaled + CHG_scaled + mobility_scaled + spt)))")
)

## ------------ ##
## BURIAL SIDES ##
## ------------ ##

## LEFT
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

## RIGHT
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

## BACK
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

## Number of sample points for training and validation
all.equal(nrow(left), nrow(right))
all.equal(nrow(back), nrow(right))
nr <- nrow(left)
n <- nr/2

#### ------------------------------------ ####
#### BASELINE AND FULL MODEL PERMUTATIONS ####
#### ------------------------------------ ####

plan(multisession, workers = 6)

left_results <- list()
right_results <- list()
back_results <- list()

pb <- txtProgressBar(min = 1, max = 100, style = 3)
for (i in 1:100) {
    ## Data
    smpl <- sample(1:nr, n)
    # Left
    train_left <- left[smpl, ]
    valid_left <- left[-smpl, ]
    pred_df_left <- valid_left[, which(names(valid_left) != "Left")]
    # Right
    train_right <- right[smpl, ]
    valid_right <- right[-smpl, ]
    pred_df_right <- valid_right[, which(names(valid_right) != "Right")]
    # Back
    train_back <- back[smpl, ]
    valid_back <- back[-smpl, ]
    pred_df_back <- valid_back[, which(names(valid_back) != "Back")]
    ## Models
    # Left
    left_m <- future_map(components, function(c) {
        m <- bru(as.formula(paste("Left ~", c)[2]),
            data = train_left,
            domain = list(geometry = mesh, YearBin = mesh.t),
            family = "binomial",
            options = list(safe = T)
        )
    }, .options = furrr_options(seed = 123))
    # Right
    right_m <- future_map(components, function(c) {
        m <- bru(as.formula(paste("Right ~", c)[2]),
            data = train_right,
            domain = list(geometry = mesh, YearBin = mesh.t),
            family = "binomial",
            options = list(safe = T)
        )
    }, .options = furrr_options(seed = 123))
    # Back
    back_m <- future_map(components, function(c) {
        m <- bru(as.formula(paste("Back ~", c)[2]),
            data = train_back,
            domain = list(geometry = mesh, YearBin = mesh.t),
            family = "binomial",
            options = list(safe = T)
        )
    }, .options = furrr_options(seed = 123))
    ## Predictions
    # Left
    left_pred <- future_map(1:length(formulas), function(f) {
        p <- predict(
            left_m[[f]],
            pred_df_left,
            formulas[[f]]
        )
    }, .options = furrr_options(seed = 123))
    # Right
    right_pred <- future_map(1:length(formulas), function(f) {
        p <- predict(
            right_m[[f]],
            pred_df_right,
            formulas[[f]]
        )
    }, .options = furrr_options(seed = 123))
    # Back
    back_pred <- future_map(1:length(formulas), function(f) {
        p <- predict(
            back_m[[f]],
            pred_df_back,
            formulas[[f]]
        )
    }, .options = furrr_options(seed = 123))
    ## Results
    # Left
    left_results[[i]] <- lapply(1:length(left_pred), function(x) {
        d <- data.table(
            Run = i,
            Model = models_info$name[x],
            Observed = valid_left$Left,
            Predicted = ifelse(left_pred[[x]]$mean >= 0.5, 1, 0)
        )
        d[, Delta := Observed - Predicted]
        return(d)
    })
    left_results[[i]] <- rbindlist(left_results[[i]])
    # Right
    right_results[[i]] <- lapply(1:length(right_pred), function(x) {
        d <- data.table(
            Run = i,
            Model = models_info$name[x],
            Observed = valid_right$Right,
            Predicted = ifelse(right_pred[[x]]$mean >= 0.5, 1, 0)
        )
        d[, Delta := Observed - Predicted]
        return(d)
    })
    right_results[[i]] <- rbindlist(right_results[[i]])
    # Back
    back_results[[i]] <- lapply(1:length(back_pred), function(x) {
        d <- data.table(
            Run = i,
            Model = models_info$name[x],
            Observed = valid_back$Back,
            Predicted = ifelse(back_pred[[x]]$mean >= 0.5, 1, 0)
        )
        d[, Delta := Observed - Predicted]
        return(d)
    })
    back_results[[i]] <- rbindlist(back_results[[i]])
    setTxtProgressBar(pb, i)
}
close(pb)

plan(sequential)

## Save data table
# Left
left_perm <- rbindlist(left_results)
fwrite(left_perm, "../Results/Left_Permutations.csv")

# Right
right_perm <- rbindlist(right_results)
fwrite(right_perm, "../Results/Right_Permutations.csv")

# Back
back_perm <- rbindlist(back_results)
fwrite(back_perm, "../Results/Back_Permutations.csv")

## THE BEST MODEL FOR THE BACK IS THE BASELINE SO WE ALREADY HAVE THE RESULTS

## --------------------------------------------------------------------------------
## BEST MODELS PERMUTATIONS

## Models fit
model_fit <- fread("../Results/Models_fit.csv")
best_models_names <- model_fit[, Model[which.min(waic)], by = "Variable"]$V1

plan(multisession, workers = 10)
bm_perm <- future_map(1:100, function(i) {
    ## Data
    smpl <- sample(1:nr, n)
    # Left
    train_left <- left[smpl, ]
    valid_left <- left[-smpl, ]
    pred_df_left <- valid_left[, which(names(valid_left) != "Left")]
    # Right
    train_right <- right[smpl, ]
    valid_right <- right[-smpl, ]
    pred_df_right <- valid_right[, which(names(valid_right) != "Right")]
    ## Components
    components_left <- ~ Intercept(1) +
        LVN_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        )
    components_right <- ~ Intercept(1) +
        WHG_scaled +
        EHG_scaled +
        LVN_scaled +
        CHG_scaled +
        mobility_scaled +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        )
    ## Model
    left_best_m <- bru(as.formula(paste("Left ~", components_left)[2]),
        data = train_left,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T)
    )
    right_best_m <- bru(as.formula(paste("Right ~", components_right)[2]),
        data = train_right,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = T)
    )
    ## Predictions
    # Left
    left_bestmodel_pred <- predict(
        left_best_m,
        pred_df_left,
        ~ 1 / (1 + exp(-(Intercept + LVN_scaled + spt)))
    )
    # Right
    right_bestmodel_pred <- predict(
        right_best_m,
        pred_df_right,
        ~ 1 / (1 + exp(-(Intercept + WHG_scaled + EHG_scaled + LVN_scaled + CHG_scaled +
            mobility_scaled + spt)))
    )
    ## Results
    # Left
    left_bestmodel_results <- data.table(
        Run = i,
        Model = best_models_names[1],
        Observed = valid_left$Left,
        Predicted = ifelse(left_bestmodel_pred$mean >= 0.5, 1, 0)
    )
    left_bestmodel_results[, Delta := Observed - Predicted]

    # Right
    right_bestmodel_results <- data.table(
        Run = i,
        Model = best_models_names[2],
        Observed = valid_right$Right,
        Predicted = ifelse(right_bestmodel_pred$mean >= 0.5, 1, 0)
    )
    right_bestmodel_results[, Delta := Observed - Predicted]

    ## Return
    return(list(left_bestmodel_results, right_bestmodel_results))
}, .options = furrr_options(seed = 123), .progress = T)

plan(sequential)
saveRDS(bm_perm, "../Results/BestModels_Permutations.RDS")

# Divide between Left and Right
bm_perm_left <- sapply(bm_perm, "[", 1)
bm_perm_left <- rbindlist(bm_perm_left)
fwrite(bm_perm_left, "../Results/BestModelPerm_Left.csv")

bm_perm_right <- sapply(bm_perm, "[", 2)
bm_perm_right <- rbindlist(bm_perm_right)
fwrite(bm_perm_right, "../Results/BestModelPerm_Right.csv")

# Extract for back
back_perm <- fread("../Results/Back_Permutations.csv")
bm_perm_back <- back_perm[Model == "baseline"]
fwrite(bm_perm_back, "../Results/BestModelPerm_Back.csv")

## --------------------------------------------------------------------------------------------------------  
## PERFORMANCE
left_perm <- fread("../Results/Left_Permutations.csv")
right_perm <- fread("../Results/Right_Permutations.csv")

## --------------- ##
## BASELINE MODELS ##
## --------------- ##

# Left
left_null <- left_perm[Model == "baseline"][, sum(Delta == 0)/.N, by = Run]
mean(left_null$V1)
sd(left_null$V1)

# Right
right_null <- right_perm[Model == "baseline"][, sum(Delta == 0)/.N, by = Run]
mean(right_null$V1)
sd(right_null$V1)

# Back
back_null <- back_perm[Model == "baseline"][, sum(Delta == 0)/.N, by = Run]
mean(back_null$V1)
sd(back_null$V1)

## ----------- ##
## FULL MODELS ##
## ----------- ##

# Left
left_full <- left_perm[Model != "baseline"][, sum(Delta == 0)/.N, by = Run]
mean(left_full$V1)
sd(left_full$V1)

# Right
right_full <- right_perm[Model != "baseline"][, sum(Delta == 0)/.N, by = Run]
mean(right_full$V1)
sd(right_full$V1)

# Back
back_full <- back_perm[Model != "baseline"][, sum(Delta == 0)/.N, by = Run]
mean(back_full$V1)
sd(back_full$V1)

## ----------- ##
## BEST MODELS ##
## ----------- ##

# Left
left_best <- bm_perm_left[, sum(Delta == 0)/.N, by = Run]
mean(left_best$V1)
sd(left_best$V1)

# Right
right_best <- bm_perm_right[, sum(Delta == 0)/.N, by = Run]
mean(right_best$V1)
sd(right_best$V1)

# Back
back_best <- bm_perm_back[, sum(Delta == 0)/.N, by = Run]
mean(back_best$V1)
sd(back_best$V1)

## PUT TOGETHER
performance.dt <- list(left_null, right_null, back_null, left_full, right_full, back_full, left_best, right_best, back_best)
performance.dt <- rbindlist(performance.dt)
names(performance.dt)[2] <- "Performance"
performance.dt[, Model := c(
    rep("Baseline", 300),
    rep("Full model", 300),
    rep("Best model", 300)
)]
performance.dt[, Variable := rep(
    c(
        rep("Left", 100),
        rep("Right", 100),
        rep("Back", 100)
    ),
    3
)]
performance.dt$Model <- factor(performance.dt$Model, levels = unique(performance.dt$Model))
performance.dt$Variable <- factor(performance.dt$Variable, levels = unique(performance.dt$Variable))
fwrite(performance.dt, "../Results/ModelPerformance.csv")

summary.performance.dt <- data.table()
summary.performance.dt[, Mean := c(
    mean(left_null$V1),
    mean(right_null$V1),
    mean(back_null$V1),
    mean(left_full$V1),
    mean(right_full$V1),
    mean(back_full$V1),
    mean(left_best$V1),
    mean(right_best$V1),
    mean(back_best$V1)
)]

summary.performance.dt[, SD := c(
    mean(left_null$V1),
    sd(right_null$V1),
    sd(back_null$V1),
    sd(left_full$V1),
    sd(right_full$V1),
    sd(back_full$V1),
    sd(left_best$V1),
    sd(right_best$V1),
    sd(back_best$V1)
)]

summary.performance.dt[, Variable := rep(c("Left", "Right", "Back"), 3)]

summary.performance.dt[, Model := c(
    rep("Baseline", 3),
    rep("Full model", 3),
    rep("Best model", 3)
)]

summary.performance.dt$Variable <- factor(summary.performance.dt$Variable,
    levels = unique(summary.performance.dt$Variable)
)
summary.performance.dt$Model <- factor(summary.performance.dt$Model,
    levels = unique(performance.dt$Model)
)
fwrite(summary.performance.dt, "../Results/ModelPerformanceSummary.csv")

## PLOT
# cols <- c(
#     rev(pals::brewer.blues(4))[c(1,4)],
#     rev(pals::brewer.rdpu(4))[c(1,4)],
#     rev(pals::brewer.ylgn(4))[c(1,4)]
# )
cols <- c("#125eaa", "#daecfc", "#33870e", "#d7f6cd", "#68228b", "#edd9ff")

ggplot(data = summary.performance.dt) +
    geom_errorbar(aes(y = Model, xmin = Mean - SD, xmax = Mean + SD, color = Model),
    width = 0.25) +
    geom_point(aes(y = Model, x = Mean, fill = Model, shape = Model)) +
    scale_shape_manual(values = c(21, 23, 24)) +
    facet_wrap(~Variable, nrow = 3, scales = "free_x") +
    theme_light() +
    theme(axis.title.y = element_blank(),
    legend.position = 'none')

performance.dt[, Color := as.character(NA)]
performance.dt[Variable == "Left" & Model == "Best model"]$Color <- cols[1]
performance.dt[Variable == "Left" & Model != "Best model"]$Color <- cols[2]
performance.dt[Variable == "Right" & Model == "Baseline"]$Color <- cols[3]
performance.dt[Variable == "Right" & Model != "Baseline"]$Color <- cols[4]
performance.dt[Variable == "Back" & Model == "Baseline"]$Color <- cols[5]
performance.dt[Variable == "Back" & Model == "Best model"]$Color <- cols[5]
performance.dt[is.na(Color)]$Color <- cols[6]

png("../Plots/Permutations.png", width = 10, height = 6, res = 330, units = 'in')
ggplot(data = performance.dt) +
    geom_violin(aes(y = Performance, x = Model, fill = Color),
    width = 0.8) +
    geom_boxplot(aes(y = Performance, x = Model, fill = Color),
        width = 0.25,
        linewidth = 1,
        color = "grey60"
    ) +
    scale_fill_identity() +
    facet_wrap(~Variable, ncol = 3, scales = "fixed") +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        legend.position = "none",
    )
dev.off()


cols <- pals::brewer.paired(3)

png("./Plots/Permutations.png", width = 10, height = 8, res = 330, units = 'in')
ggplot(data = performance.dt) +
  geom_violin(aes(y = Performance, x = Model, fill = Variable),
              width = 0.6) +
  geom_boxplot(aes(y = Performance, x = Model, fill = Variable),
               width = 0.15,
               linewidth = 0.8,
               color = "grey30"
  ) +
  scale_fill_manual(values = cols) +
  facet_wrap(~Variable, ncol = 3, scales = "fixed") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
  )
dev.off()