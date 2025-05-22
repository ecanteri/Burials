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

# Bind tables together
orientation <- rbind(orientation, single_angles)
orientation <- orientation[order(Num)]

## Separate data based on burial side
left <- orientation[, .(Num, BurialOrientationMean, BurialOrientationMin, BurialOrientationMax), by = "BurialSide_left"][BurialSide_left == 1]
left[, BurialSide := "Left"]
right <- orientation[, .(Num, BurialOrientationMean, BurialOrientationMin, BurialOrientationMax), by = "BurialSide_right"][BurialSide_right == 1]
right[, BurialSide := "Right"]
back <- orientation[, .(Num, BurialOrientationMean, BurialOrientationMin, BurialOrientationMax), by = "BurialSide_back"][BurialSide_back == 1]
back[, BurialSide := "Back"]

# Bind together in a new table
sides <- rbind(left[, -1], right[, -1], back[, -1])

# Calculate cosine and sine from radians
sides[, BurialOrientationMean_cos := cos(BurialOrientationMean*pi/180)]
sides[, BurialOrientationMin_cos := cos(BurialOrientationMin*pi/180)]
sides[, BurialOrientationMax_cos := cos(BurialOrientationMax*pi/180)]
sides[, BurialOrientationMean_sin := sin(BurialOrientationMean*pi/180)]
sides[, BurialOrientationMin_sin := sin(BurialOrientationMin*pi/180)]
sides[, BurialOrientationMax_sin := sin(BurialOrientationMax*pi/180)]

# Add Country, Culture, Period and Ancestry info
sides$Country <- graves[Num %in% sides$Num]$Country
sides$Period <- graves[Num %in% sides$Num]$Period
sides$Culture <- graves[Num %in% sides$Num]$Culture
sides$Sex <- graves[Num %in% sides$Num]$Sex
sides$Age <- graves[Num %in% sides$Num]$gaussianModelMu
sides$EHG <- graves[Num %in% sides$Num]$ANCE2_v2
sides$WHG <- graves[Num %in% sides$Num]$ANCE4_v2
sides$LVN <- graves[Num %in% sides$Num]$ANCE6_v2
sides$CHG <- graves[Num %in% sides$Num]$ANCE8_v2
sides$mobility <- graves[Num %in% sides$Num]$mobility_MDS_250y_retrospective

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

#### PLOTS ####
##### All individuals ######
dim_angles <- sides[, .N, by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]

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


##### Burial side ######
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

##### Country ######
dim_country <- sides[, .N, by = c("Country", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]
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

##### Period ######
dim_period <- sides[, .N, by = c("Period", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]
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

##### Culture #####
dim_culture <- sides[, .N, by = c("Culture", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_culture <- dim_culture[!is.na(Culture)]
# Plot only the 10 most frequent cultures
culture_keep <- dim_culture[, .N, by = "Culture"][order(N, decreasing = T)][, Culture][1:10]

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

##### Sex #####
dim_sex <- sides[, .N, by = c("Sex", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]

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

##### Age #####
dim_age <- sides[, .N, by = c("Age", "BurialOrientationMean_sin", "BurialOrientationMean_cos")]

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
  scale_fill_gradientn("Age", colors = pals::coolwarm(100)) +
  scale_size_continuous("Count", breaks = c(1, 25, 50, 75, 100), range = c(3,10)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

##### Ancestry: WHG ##### 
png("./Figures/BurialOrientation/BurialOrientation_WHG.png", width = 12, height = 10, units = 'in', res = 330)
dim_whg <- sides[, mean(WHG, na.rm = T), by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_whg <- dim_whg[!is.nan(V1),]
names(dim_whg)[3] <- "WHG"

ggplot(data = dim_whg) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = WHG),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("WHG", colors = pals::parula(100), limits = c(0,1)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

##### Ancestry: EHG ##### 
dim_ehg <- sides[, mean(EHG, na.rm = T), by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_ehg <- dim_ehg[!is.nan(V1),]
names(dim_ehg)[3] <- "EHG"

png("./Figures/BurialOrientation/BurialOrientation_EHG.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_ehg) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = EHG),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("EHG", colors = pals::parula(100), limits = c(0,1)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

##### Ancestry: CHG ##### 
dim_chg <- sides[, mean(CHG, na.rm = T), by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_chg <- dim_chg[!is.nan(V1),]
names(dim_chg)[3] <- "CHG"

png("./Figures/BurialOrientation/BurialOrientation_CHG.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_chg) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = CHG),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("CHG", colors = pals::parula(100), limits = c(0,1)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

##### Ancestry: LVN ##### 
dim_lvn <- sides[, mean(LVN, na.rm = T), by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_lvn <- dim_lvn[!is.nan(V1),]
names(dim_lvn)[3] <- "LVN"

png("./Figures/BurialOrientation/BurialOrientation_LVN.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_lvn) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = LVN),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("LVN", colors = pals::parula(100), limits = c(0,1)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

##### Mobility #####
dim_mob <- sides[, mean(mobility, na.rm = T), by = c("BurialOrientationMean_sin", "BurialOrientationMean_cos")]
dim_mob <- dim_mob[!is.nan(V1),]
names(dim_mob)[3] <- "mobility"

png("./Figures/BurialOrientation/BurialOrientation_mobility.png", width = 12, height = 10, units = 'in', res = 330)
ggplot(data = dim_mob) +
  geom_point(aes(BurialOrientationMean_sin, BurialOrientationMean_cos, fill = mobility),
             shape = 21,
             size = 5) +
  scale_fill_gradientn("Mobility", colors = pals::parula(100)) +
  xlab(expression(sin(theta))) +
  ylab(expression(cos(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()
