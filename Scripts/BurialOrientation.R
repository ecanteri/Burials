library(data.table)
library(ggplot2)
library(FactoMineR)
library(circular)
library(CircStats)
library(BIADconnect)
conn <- init.conn()

#### DATA ####
graves <- fread("./Data/Burial.csv", na.strings = "")
graves$DepositionType <- as.factor(graves$DepositionType)
graves$BodyPositioning <- as.factor(graves$BodyPositioning)
graves$BurialSide <- as.factor(graves$BurialSide)
graves[order(IndividualID)]
str(graves)

# Sites table (BIAD)
sites <- setDT(query.database("SELECT * FROM Sites", conn))

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
hist(graves$YearBP, breaks = 100)

## Create data table
orientation <- graves[, .(Num, BurialOrientationMin, BurialOrientationMax, BurialSide_back, BurialSide_left, BurialSide_right)]

## We can't just average between the min and max angles because for some angles the average will be 0, so we will not get burials facing south.
## We create a function to transform the angles that works around this issue.
#----------------------------------------------------------------------
find.mid.deg <- function(min.deg, max.deg){
  
  difference <- max.deg - min.deg
  i.accute <- difference<180
  i.obtuse <- difference>180
  i.bad <- difference==180
  
  mid <- numeric()
  mid[i.accute] <- max.deg[i.accute] - difference[i.accute]/2
  mid[i.obtuse] <- max.deg[i.obtuse] + (360 - difference[i.obtuse])/2
  mid[i.bad] <- NA
  
  i <- mid>=360
  mid[i] <- mid[i]-360
  
  i <- mid<0
  mid[i] <- mid[i]+360           
  return(mid)}
#----------------------------------------------------------------------
min.deg <- c(5,15,170,270)
max.deg <- c(356,300,185,360)
find.mid.deg(min.deg, max.deg)
#----------------------------------------------------------------------

# Remove NAs in both min and max orientation columns
orientation <- orientation[!(is.na(BurialOrientationMin) & is.na(BurialOrientationMax))]

# Create a table for records that have a single angle value (in BurialOrientationMin)
single_angles <- orientation[is.na(BurialOrientationMax)]
single_angles[, BurialOrientationMean := BurialOrientationMin]

# Remove from overall table
orientation <- orientation[!(is.na(BurialOrientationMin) | is.na(BurialOrientationMax))]

# Average across min and max 
orientation[, BurialOrientationMean := find.mid.deg(BurialOrientationMin, BurialOrientationMax)]

## Using circular package
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

#### HISTOGRAMS ####
png("./Figures/BurialOrientation/BurialOrientation_Mean.png", width = 10, height = 5, units = 'in', res = 330)
ggplot(data = sides) +
  geom_histogram(aes(x = BurialOrientationMean), bins = 30, fill = "white", col = "black") +
  xlab("Burial orientation (mean)") +
  facet_wrap(~BurialSide) +
  theme_bw()
dev.off()

ggplot(data = sides) +
  geom_histogram(aes(x = BurialOrientationMin), bins = 30, fill = "white", col = "black") +
  xlab("Burial orientation (min)") +
  facet_wrap(~BurialSide) +
  theme_bw()

ggplot(data = sides) +
  geom_histogram(aes(x = BurialOrientationMax), bins = 30, fill = "white", col = "black") +
  xlab("Burial orientation (max)") +
  facet_wrap(~BurialSide) +
  theme_bw()

png("./Figures/BurialOrientation/BurialOrientation_cosMean.png", width = 10, height = 5, units = 'in', res = 330)
ggplot(data = sides) +
  geom_histogram(aes(x = BurialOrientationMean_cos), bins = 30, fill = "white", col = "black") +
  xlab("Burial orientation cos(mean)") +
  facet_wrap(~BurialSide) +
  theme_bw()
dev.off()

ggplot(data = sides) +
  geom_histogram(aes(x = BurialOrientationMin_cos), bins = 30, fill = "white", col = "black") +
  xlab("Burial orientation cos(min)") +
  facet_wrap(~BurialSide) +
  theme_bw()

ggplot(data = sides) +
  geom_histogram(aes(x = BurialOrientationMax_cos), bins = 30, fill = "white", col = "black") +
  xlab("Burial orientation cos(max)") +
  facet_wrap(~BurialSide) +
  theme_bw()

png("./Figures/BurialOrientation/BurialOrientation_sinMean.png", width = 10, height = 5, units = 'in', res = 330)
ggplot(data = sides) +
  geom_histogram(aes(x = BurialOrientationMean_sin), bins = 30, fill = "white", col = "black") +
  xlab("Burial orientation sin(mean)") +
  facet_wrap(~BurialSide) +
  theme_bw()
dev.off()

ggplot(data = sides) +
  geom_histogram(aes(x = BurialOrientationMin_sin), bins = 30, fill = "white", col = "black") +
  xlab("Burial orientation sin(min)") +
  facet_wrap(~BurialSide) +
  theme_bw()

ggplot(data = sides) +
  geom_histogram(aes(x = BurialOrientationMax_sin), bins = 30, fill = "white", col = "black") +
  xlab("Burial orientation sin(max)") +
  facet_wrap(~BurialSide) +
  theme_bw()

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
  bounds.zero <- round(bounds.zero, decimals)
  center.zero <- seq(bounds.zero[1], 360, by = 1/10^decimals)
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

#### PLOTS ####
##### All individuals ######
dim_angles <- orientation[, .N, by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]

ggplot(data = dim_angles) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, size = N), 
             shape = 21) +
  scale_size_continuous("Count", breaks = c(1, 10, 50, 100, 150, 200, 250),
                        range = c(3,10)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

# Circular histogram
w <- 22.5
orientation[, Bins := circular.bins(BurialOrientationMean, binwidth = w)]
ggplot(data = orientation, aes(x = Bins)) +
  geom_bar(fill = "grey80", color = "black") +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  coord_polar(start = -.19) +
  theme_light() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12))

# w <- 22.5
# ggplot(data = orientation, aes(x = BurialOrientationMean)) +
#   stat_bin(binwidth = w, fill = "grey80", color = "black") +
#   # geom_histogram(fill = "grey80", color = "black") +
#   scale_x_continuous(expand = expansion(0, 0),
#                      limits = c(0, 360),
#                      breaks = seq(0, 360, by = w)) +
#   scale_y_continuous(expand = expansion(c(0, 0.05))) +
#   coord_polar() +
#   # coord_radial() +
#   theme_light() +
#   theme(axis.title = element_blank(),
#         axis.text = element_text(size = 12))


##### Burial side ######
## Separate data based on burial side
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
# Group by count within same burial side and same angles
# Count number of entries per burial side
dim_side <- sides[, .N, by = c("BurialSide", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]

png("./Figures/BurialOrientation/BurialOrientation_BurialSide.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_side) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = BurialSide, size = N),
             position = position_jitter(0.05, 0.05),
             shape = 21) +
  scale_fill_manual("Burial side", values = pals::brewer.paired(3)) +
  scale_size_continuous("Count", breaks = c(1, 10, 50, 100, 150),
                        range = c(3,10)) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Circular histogram
cols <- pals::brewer.paired(3)
w <- 22.5
ggplot(data = sides, aes(x = Bins,
                         fill = factor(BurialSide, levels = c("Left", "Right", "Back"))),
       alpha = 0.5) +
  geom_bar(color = "black") +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Burial side", values = pals::brewer.paired(3)) +
  coord_polar(start = -.19) +
  ylab("Count") +
  theme_light() +
  theme(axis.text = element_text(size = 12))

# ggplot(data = sides, aes(x = BurialOrientationMean, fill = factor(BurialSide, levels = c("Left", "Right", "Back"))), alpha = 0.5) +
#   stat_bin(bins = 30, boundary = 0, color = 'black') +
#   scale_x_continuous(expand = expansion(0, 0), 
#                      limits = c(0, 360),
#                      breaks = seq(0, 360, by = 45)) +
#   scale_y_continuous(expand = expansion(c(0, 0.05))) +
#   scale_fill_manual("Burial side", values = pals::brewer.paired(3)) +
#   ylab("Count") +
#   coord_polar() +
#   theme_light() +
#   theme(axis.text = element_text(size = 12))

# Single plots
ggplot(data = sides, aes(x = Bins, 
                         fill = factor(BurialSide, levels = c("Left", "Right", "Back")))) +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Burial side", values = pals::brewer.paired(3)) +
  facet_wrap(~factor(BurialSide, levels = c("Left", "Right", "Back")),
             scales = "free_y") +
  ylab("Count") +
  coord_polar(start = -.19) +
  theme_light() +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none")

##### Country #####
dim_country <- orientation[, .N, by = c("Country", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_country$Country <- factor(dim_country$Country, levels = c("Netherlands",
                                                              "France",
                                                              "Germany",
                                                              "Denmark",
                                                              "Sweden",
                                                              "Latvia",
                                                              "Poland",
                                                              "Czechia / Czech Republic",
                                                              "Slovakia",
                                                              "Austria",
                                                              "Hungary",
                                                              "Croatia",
                                                              "Romania",
                                                              "Serbia",
                                                              "Bulgaria",
                                                              "Ukraine",
                                                              "Russia",
                                                              "Armenia",
                                                              "Kazakhstan"))

png("./Figures/BurialOrientation/BurialOrientation_Country.png", width = 13, height = 10, units = 'in', res = 330)
ggplot(data = dim_country) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = Country, size = N),
             position = position_jitter(0.05, 0.05),
             shape = 21) +
  scale_fill_manual("Country", values = pals::parula(19)) +
  scale_size_continuous("Count", breaks = c(1, 10, 25, 50, 75, 100),
                        range = c(3,10)) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Circular histogram
w <- 22.5
ggplot(data = orientation,
       aes(x = Bins, fill = Country),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Country", values = pals::tol.rainbow(19)) +
  coord_polar(start = -.19) +
  theme_light() +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank())


# Single plots
ggplot(data = orientation, 
       aes(x = Bins, fill = Country)) +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Country", values = pals::tol.rainbow(19)) +
  facet_wrap(~Country, scales = "free_y") +
  ylab("Count") +
  coord_polar(start = -.19) +
  theme_light() +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none")

##### Period ######
dim_period <- orientation[, .N, by = c("Period", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_period$Period <- factor(dim_period$Period, levels = c("Mesolithic",
                                                          "Mesolithic/Neolithic",
                                                          "Neolithic",
                                                          "Eneolithic",
                                                          "Eneolithic/Bronze Age",
                                                          "Bronze Age",
                                                          "Bronze Age/Iron Age",
                                                          "Iron Age"))

png("./Figures/BurialOrientation/BurialOrientation_Period.png", width = 13, height = 10, units = 'in', res = 330)
ggplot(data = dim_period) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = Period, size = N),
             position = position_jitter(0.05, 0.05),
             shape = 21) +
  scale_fill_manual("Period", values = pals::parula(8)) +
  scale_size_continuous("Count", breaks = c(1, 10, 50, 100, 200),
                        range = c(3,10)) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Circular histogram
w <- 22.5
ggplot(data = orientation,
       aes(x = Bins, fill = factor(Period,
                                   levels = c("Mesolithic",
                                              "Mesolithic/Neolithic",
                                              "Neolithic",
                                              "Eneolithic",
                                              "Eneolithic/Bronze Age",
                                              "Bronze Age",
                                              "Bronze Age/Iron Age",
                                              "Iron Age"))),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Period", values = pals::parula(8)) +
  coord_polar(start = -.19) +
  theme_light() +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12))

# Single plots
ggplot(data = orientation, 
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
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Period", values = pals::parula(8)) +
  facet_wrap(~factor(Period,
                     levels = c("Mesolithic",
                                "Mesolithic/Neolithic",
                                "Neolithic",
                                "Eneolithic",
                                "Eneolithic/Bronze Age",
                                "Bronze Age",
                                "Bronze Age/Iron Age",
                                "Iron Age")), scales = "free_y") +
  ylab("Count") +
  coord_polar(start = -.19) +
  theme_light() +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none")


##### Culture #####
dim_culture <- orientation[, .N, by = c("Culture", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_culture <- dim_culture[!is.na(Culture)]
# Plot only the 10 most frequent cultures
culture_keep <- orientation[!is.na(Culture)][, .N, by = "Culture"][order(N, decreasing = T)][, Culture][1:10]

png("./Figures/BurialOrientation/BurialOrientation_Culture.png", width = 15, height = 10, units = 'in', res = 330)
ggplot(data = dim_culture[Culture %in% culture_keep]) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = Culture, size = N),
             position = position_jitter(0.05, 0.05),
             shape = 21) +
  scale_fill_manual("Culture", values = pals::parula(10)) +
  scale_size_continuous("Count", breaks = c(1, 10, 25, 50, 75), range = c(3,10)) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Circular histogram
w <- 22.5
ggplot(data = orientation[Culture %in% culture_keep],
       aes(x = Bins, fill = Culture),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Culture", values = pals::parula(10)) +
  coord_polar(start = -.19) +
  theme_light() +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12))

# Single plots
ggplot(data = orientation[Culture %in% culture_keep],
       aes(x = Bins, fill = Culture),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Culture", values = pals::parula(10)) +
  coord_polar(start = -.19) +
  facet_wrap(~Culture, scales = "free_y") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

##### Sex #####
dim_sex <- orientation[, .N, by = c("Sex", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]

png("./Figures/BurialOrientation/BurialOrientation_Sex.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_sex) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = Sex, size = N),
             position = position_jitter(0.08, 0.08),
             shape = 21) +
  scale_fill_gradientn("Sex", colors = rev(pals::brewer.piyg(100))) +
  scale_size_continuous("Count", breaks = c(1, 30, 60, 90, 120), range = c(3,10)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Circular histogram
w <- 22.5
ggplot(data = orientation,
       aes(x = Bins, fill = as.factor(Sex)),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Sex", values = rev(pals::brewer.piyg(9))) +
  coord_polar(start = -.19) +
  theme_light() +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12))

# Single plots
ggplot(data = orientation,
       aes(x = Bins, fill = as.factor(Sex)),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Sex", values = rev(pals::brewer.piyg(9))) +
  coord_polar(start = -.19) +
  facet_wrap(~as.factor(Sex), scales = "free_y") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

##### Age #####
dim_age <- orientation[, .N, by = c("Age", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]

png("./Figures/BurialOrientation/BurialOrientation_Age.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_age) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = Age, size = N),
             position = position_jitter(0.1, 0.1),
             shape = 21) +
  scale_fill_gradientn("Age", colors = pals::coolwarm(100)) +
  scale_size_continuous("Count", breaks = c(1, 20, 50, 75), range = c(3,10)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Age interval between 8k and 2k BP
ggplot(data = dim_age[between(Age, 2000, 8000)]) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = Age, size = N),
             position = position_jitter(0.1, 0.1),
             shape = 21) +
  scale_fill_gradientn("Age", colors = pals::coolwarm(100), limits = c(2000, 8000),
                       breaks = seq(2000, 8000, 2000)) +
  scale_size_continuous("Count", breaks = c(1, 25, 50, 75, 100), range = c(3,10)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

# Histograms
w <- 22.5
orientation$YearBin <- plyr::round_any(orientation$Age, 500)
ggplot(data = orientation[between(Age, 3000, 8000)],
       aes(x = Bins, fill = factor(YearBin, levels = as.character(rev(seq(3000, 8000, 500))))),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("YearBin", values = rev(pals::brewer.piyg(11))) +
  coord_polar(start = -.19) +
  ylab("Count") +
  facet_wrap(~factor(YearBin, levels = as.character(rev(seq(3000, 8000, 500)))), scales = "free_y") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")


##### Ancestry: WHG ##### 
png("/Users/msb290/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Graves/Figures/BurialOrientation/BurialOrientation_WHG.png", width = 12, height = 10, units = 'in', res = 330)
dim_whg <- orientation[, mean(WHG, na.rm = T), by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_whg <- dim_whg[!is.nan(V1),]
names(dim_whg)[3] <- "WHG"

ggplot(data = dim_whg) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = WHG),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("WHG", colors = pals::parula(100), limits = c(0,1)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Histograms
w <- 22.5
orientation$WHG_bin <- plyr::round_any(orientation$WHG, 0.05)
ggplot(data = orientation[!is.na(WHG)],
       aes(x = Bins, fill = as.factor(WHG_bin)),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("WHG", values = pals::parula(15)) +
  coord_polar(start = -.19) +
  ylab("Count") +
  facet_wrap(~as.factor(WHG_bin), scales = "free_y") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

##### Ancestry: EHG ##### 
dim_ehg <- orientation[, mean(EHG, na.rm = T), by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_ehg <- dim_ehg[!is.nan(V1),]
names(dim_ehg)[3] <- "EHG"

png("/Users/msb290/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Graves/Figures/BurialOrientation/BurialOrientation_EHG.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_ehg) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = EHG),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("EHG", colors = pals::parula(100), limits = c(0,1)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Histograms
w <- 22.5
orientation$EHG_bin <- plyr::round_any(orientation$EHG, 0.05)
ggplot(data = orientation[!is.na(EHG)],
       aes(x = Bins, fill = as.factor(EHG_bin)),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("EHG", values = pals::parula(14)) +
  coord_polar(start = -.19) +
  ylab("Count") +
  facet_wrap(~as.factor(EHG_bin), scales = "free_y") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

##### Ancestry: CHG ##### 
dim_chg <- orientation[, mean(CHG, na.rm = T), by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_chg <- dim_chg[!is.nan(V1),]
names(dim_chg)[3] <- "CHG"

png("/Users/msb290/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Graves/Figures/BurialOrientation/BurialOrientation_CHG.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_chg) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = CHG),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("CHG", colors = pals::parula(100), limits = c(0,1)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Histograms
w <- 22.5
orientation$CHG_bin <- plyr::round_any(orientation$CHG, 0.05)
ggplot(data = orientation[!is.na(CHG)],
       aes(x = Bins, fill = as.factor(CHG_bin)),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("CHG", values = pals::parula(14)) +
  coord_polar(start = -.19) +
  ylab("Count") +
  facet_wrap(~as.factor(CHG_bin), scales = "free_y") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

##### Ancestry: LVN ##### 
dim_lvn <- orientation[, mean(LVN, na.rm = T), by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_lvn <- dim_lvn[!is.nan(V1),]
names(dim_lvn)[3] <- "LVN"

png("/Users/msb290/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Graves/Figures/BurialOrientation/BurialOrientation_LVN.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_lvn) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = LVN),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("LVN", colors = pals::parula(100), limits = c(0,1)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Histograms
w <- 22.5
orientation$LVN_bin <- plyr::round_any(orientation$LVN, 0.05)
ggplot(data = orientation[!is.na(LVN)],
       aes(x = Bins, fill = as.factor(LVN_bin)),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  coord_polar(start = -.19) +
  ylab("Count") +
  facet_wrap(~as.factor(LVN_bin), scales = "free_y") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

##### Mobility #####
dim_mob <- orientation[, mean(mobility, na.rm = T), by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_mob <- dim_mob[!is.nan(V1),]
names(dim_mob)[3] <- "mobility"

png("/Users/msb290/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Graves/Figures/BurialOrientation/BurialOrientation_mobility.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_mob) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = mobility),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("Mobility", colors = pals::parula(100)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

# Histograms
w <- 22.5
orientation$mob_bin <- plyr::round_any(orientation$mobility, 50)
ggplot(data = orientation[!is.na(mobility)],
       aes(x = Bins, fill = as.factor(mob_bin)),
       color = 'black') +
  geom_bar(color = 'black') +
  scale_x_continuous(breaks = seq(0, length(seq(0, 360, w))-1, 1),
                     labels = seq(0, 360, w)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_fill_manual("Mobility", values = pals::parula(28)) +
  coord_polar(start = -.19) +
  ylab("Count") +
  facet_wrap(~as.factor(mob_bin), scales = "free_y") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")
