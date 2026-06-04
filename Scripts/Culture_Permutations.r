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
components <- c(
    cb = ~ Intercept(1) +
        culture(Culture,
            model = "factor_full"
        ) -1 +
        spt(geometry,
            model = matern,
            group = YearBin,
            group_mapper = bru_mapper(mesh.t, indexed = TRUE),
            control.group = list(model = "ar1", hyper = h.spec)
        )
)

## Model info
models_info <- data.frame(
    model = "cb",
    name = c(
        "baseline + culture"
    )
)

## Prediction formulas
formulas <- list(
    f1 = as.formula("~ 1 / (1 + exp(-(Intercept + spt + culture)))")
)

## ------------ ##
## BURIAL SIDES ##
## ------------ ##

## LEFT
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

## RIGHT
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

## BACK
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

## Number of sample points for training and validation
all.equal(nrow(left), nrow(right))
all.equal(nrow(back), nrow(right))
nr <- nrow(left)
n <- floor(nr/2)

#### ------------ ####
#### PERMUTATIONS ####
#### ------------ ####

library(furrr)
library(future)
library(data.table)
library(inlabru)

# Set up your multi-session backend
plan(multisession, workers = 10)

culture_perm <- future_map(1:100, function(i) {
    ## Data Sampling
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
    left_culture_m <- bru(as.formula(paste("Left ~", components$cb)[2]),
        newdata = train_left,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = TRUE)
    )
    
    right_culture_m <- bru(as.formula(paste("Right ~", components$cb)[2]),
        data = train_right,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = TRUE)
    )
    
    back_culture_m <- bru(as.formula(paste("Back ~", components$cb)[2]),
        data = train_back,
        domain = list(geometry = mesh, YearBin = mesh.t),
        family = "binomial",
        options = list(safe = TRUE)
    )
    
    ## Predictions
    left_culture_pred  <- predict(left_culture_m,  pred_df_left,  formulas$f1)
    right_culture_pred <- predict(right_culture_m, pred_df_right, formulas$f1)
    back_culture_pred  <- predict(back_culture_m,  pred_df_back,  formulas$f1)
    
    ## Results Collection
    # Left
    left_culture_results <- data.table(
        Run = i,
        Model = models_info$name,
        Observed = valid_left$Left,
        Predicted = ifelse(left_culture_pred$mean >= 0.5, 1, 0)
    )
    left_culture_results[, Delta := Observed - Predicted]

    # Right
    right_culture_results <- data.table(
        Run = i,
        Model = models_info$name,
        Observed = valid_right$Right,
        Predicted = ifelse(right_culture_pred$mean >= 0.5, 1, 0)
    )
    right_culture_results[, Delta := Observed - Predicted]

    # Back (Fixed the assignment variable typo here)
    back_culture_results <- data.table(
        Run = i,
        Model = models_info$name,
        Observed = valid_back$Back,
        Predicted = ifelse(back_culture_pred$mean >= 0.5, 1, 0)
    )
    back_culture_results[, Delta := Observed - Predicted]

    ## Return
    return(list(left_culture_results, right_culture_results, back_culture_results))
    
}, .options = furrr_options(
    seed = 123,
    # CRITICAL: Manually export the variables hidden from automated static text analysis
    globals = c("matern", "mesh", "mesh.t", "components", "formulas", "models_info", "nr", "n", "left", "right", "back", "h.spec"),
    packages = c("inlabru", "INLA", "data.table")
), .progress = TRUE)

plan(sequential)
saveRDS(culture_perm, "./Results/Culture_Permutations.RDS")

## ------- ##
## RESULTS ##
## ------- ##

# Left
culture_perm_left <- sapply(culture_perm, "[", 1)
culture_perm_left <- rbindlist(culture_perm_left)
fwrite(culture_perm_left, "./Results/CulturePerm_Left.csv")

left_perms <- culture_perm_left[, sum(Delta == 0)/.N, by = Run]
mean(left_perms$V1)
sd(left_perms$V1)

# Right
culture_perm_right <- sapply(culture_perm, "[", 2)
culture_perm_right <- rbindlist(culture_perm_right)
fwrite(culture_perm_right, "./Results/CulturePerm_Right.csv")

right_perms <- culture_perm_right[, sum(Delta == 0)/.N, by = Run]
mean(right_perms$V1)
sd(right_perms$V1)

# Back
culture_perm_back <- sapply(culture_perm, "[", 3)
culture_perm_back <- rbindlist(culture_perm_back)
fwrite(culture_perm_back, "./Results/CulturePerm_Back.csv")

back_perms <- culture_perm_back[, sum(Delta == 0)/.N, by = Run]
mean(back_perms$V1)
sd(back_perms$V1)

## PUT TOGETHER
culture_performance.dt <- list(left_perms, right_perms, back_perms)
culture_performance.dt <- rbindlist(culture_performance.dt)
names(culture_performance.dt)[2] <- "Performance"
culture_performance.dt[, Model := "baseline + culture"]
culture_performance.dt[, Variable := 
    c(
        rep("Left", 100),
        rep("Right", 100),
        rep("Back", 100)
    )
]
culture_performance.dt$Variable <- factor(culture_performance.dt$Variable, levels = unique(culture_performance.dt$Variable))
fwrite(culture_performance.dt, "./Results/CulturePerformance.csv")

## --------------------------------------------------------------------------------------
## PERCENTAGE OF EXPLAINED VARIANCE
load("/projects/racimolab/people/msb290/Graves/Models/Culture/Left.RData")
load("/projects/racimolab/people/msb290/Graves/Models/Culture/Right.RData")
load("/projects/racimolab/people/msb290/Graves/Models/Culture/Back.RData")

left_culture <- left_models[[8]]
right_culture <- right_models[[8]]
back_culture <- back_models[[8]]

## BUILD A COMPREHENSIVE TABLE OF NAKAGAWA, EFRON, and GELMAN R-SQUARED METRICS

# --- 1. Set Up the Parallel Infrastructure ---
# We use 3 workers since we have exactly 3 burial sides to process
plan(multisession, workers = 3)

# Define our list structures
model_list <- list(
  Left  = left_culture,
  Right = right_culture,
  Back  = back_culture
)

data_list <- list(
  Left  = left,
  Right = right,
  Back  = back
)

# Shared logistic variance constant (pi^2 / 3)
var_binomial <- (pi^2) / 3

# --- 2. Parallel Processing across Sides ---
parallel_r2_results <- future_map(names(model_list), function(side) {
    # Extract the correct model and data objects for this background worker
    current_model <- model_list[[side]]
    current_data <- data_list[[side]]
    observed <- current_data[[side]]
    n_obs <- length(observed)

    # Baseline total sum of squares for Efron calculation
    mean_obs <- mean(observed)
    ss_total <- sum((observed - mean_obs)^2)

    # Draw 1000 posterior realizations inside the worker environment
    set.seed(123)
    post_samples <- inlabru::generate(
        object = current_model,
        newdata = current_data,
        formula = ~ list(
            eta_full = Intercept + culture + spt,
            culture_link = culture,
            spatial_link = spt,
            p_full = INLA::inla.link.invlogit(Intercept + culture + spt),
            p_culture = INLA::inla.link.invlogit(Intercept + culture),
            p_spatial = INLA::inla.link.invlogit(Intercept + spt)
        ),
        n.samples = 1000
    )

    # Initialize tracking arrays for all frameworks
    nak_culture <- numeric(1000)
    nak_spatial <- numeric(1000)
    nak_full <- numeric(1000)
    efr_culture <- numeric(1000)
    efr_spatial <- numeric(1000)
    efr_full <- numeric(1000)
    gel_culture <- numeric(1000)
    gel_spatial <- numeric(1000)
    gel_full <- numeric(1000)

    for (s in 1:1000) {
        # Extract Link Scale Vectors
        v_nak_full <- as.numeric(post_samples[[s]]$eta_full)
        v_nak_culture <- as.numeric(post_samples[[s]]$culture_link)
        v_nak_spatial <- as.numeric(post_samples[[s]]$spatial_link)

        # --- A. NAKAGAWA FRAMEWORK ---
        var_c_link <- var(v_nak_culture)
        var_s_link <- var(v_nak_spatial)
        total_model_var <- var_c_link + var_s_link
        total_denom <- total_model_var + var_binomial

        nak_culture[s] <- max(0, min(1, var_c_link / total_denom))
        nak_spatial[s] <- max(0, min(1, var_s_link / total_denom))
        nak_full[s] <- max(0, min(1, total_model_var / total_denom))

        # Extract Probability Scale Vectors
        pi_full <- as.numeric(post_samples[[s]]$p_full)
        pi_culture <- as.numeric(post_samples[[s]]$p_culture)
        pi_spatial <- as.numeric(post_samples[[s]]$p_spatial)

        # --- B. EFRON FRAMEWORK ---
        ss_res_full <- sum((observed - pi_full)^2)
        ss_res_culture <- sum((observed - pi_culture)^2)
        ss_res_spatial <- sum((observed - pi_spatial)^2)

        efr_culture[s] <- max(0, min(1, 1 - (ss_res_culture / ss_total)))
        efr_spatial[s] <- max(0, min(1, 1 - (ss_res_spatial / ss_total)))
        efr_full[s] <- max(0, min(1, 1 - (ss_res_full / ss_total)))

        # --- C. GELMAN FRAMEWORK ---
        # Full Model
        var_fit_full <- var(pi_full)
        var_res_full <- mean(pi_full * (1 - pi_full))
        gel_full[s] <- var_fit_full / (var_fit_full + var_res_full)

        # Culture Only
        var_fit_cult <- var(pi_culture)
        var_res_cult <- mean(pi_culture * (1 - pi_culture))
        gel_culture[s] <- var_fit_cult / (var_fit_cult + var_res_cult)

        # Spatial Only
        var_fit_spat <- var(pi_spatial)
        var_res_spat <- mean(pi_spatial * (1 - pi_spatial))
        gel_spatial[s] <- var_fit_spat / (var_fit_spat + var_res_spat)
    }

    # Package all 3 frameworks neatly for this side
    side_table <- data.table::data.table(
        Side = side,
        Framework = c(
            rep("Nakagawa (Latent Log-Odds)", 3),
            rep("Efron (Observed Probability)", 3),
            rep("Gelman (Model-Based Probability)", 3)
        ),
        Effect = rep(c("Culture Effect Alone", "Spatiotemporal Process Alone", "Full Model Fit"), 3),
        Mean = c(
            mean(nak_culture), mean(nak_spatial), mean(nak_full),
            mean(efr_culture), mean(efr_spatial), mean(efr_full),
            mean(gel_culture), mean(gel_spatial), mean(gel_full)
        ) * 100,
        Lower = c(
            quantile(nak_culture, 0.025), quantile(nak_spatial, 0.025), quantile(nak_full, 0.025),
            quantile(efr_culture, 0.025), quantile(efr_spatial, 0.025), quantile(efr_full, 0.025),
            quantile(gel_culture, 0.025), quantile(gel_spatial, 0.025), quantile(gel_full, 0.025)
        ) * 100,
        Upper = c(
            quantile(nak_culture, 0.975), quantile(nak_spatial, 0.975), quantile(nak_full, 0.975),
            quantile(efr_culture, 0.975), quantile(efr_spatial, 0.975), quantile(efr_full, 0.975),
            quantile(gel_culture, 0.975), quantile(gel_spatial, 0.975), quantile(gel_full, 0.975)
        ) * 100
    )

    return(side_table)
},
.options = furrr_options(
    seed = 123,
    packages = c("inlabru", "INLA", "data.table"),
    globals = c("model_list", "data_list", "var_binomial")
),
.progress = TRUE
)

# --- 3. Bind and Format Master Combined Output ---
comprehensive_master_table <- data.table::rbindlist(parallel_r2_results)

comprehensive_master_table[, `:=`(
  Mean  = round(Mean, 2),
  Lower = round(Lower, 2),
  Upper = round(Upper, 2)
)]

# Revert to standard execution and print master compiled metrics
plan(sequential)
print(comprehensive_master_table)

# Save to file
fwrite(comprehensive_master_table, file = "./Results/Culture_ExplainedVariance.csv")