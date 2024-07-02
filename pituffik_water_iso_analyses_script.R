#===============================================================================
#=================Pituffik Surface Waters Script===================================
#===============================================================================
# This script is focused on the water isotopes of lakes and streams in Pituffik, Greenland, based
  # on samples taken in 2018 and 2019.
# Written Pete D Akers, March 2021-August 2023, revised March-July 2024

# Loading required packages
library(raster)
library(ncdf4)
library(Rmisc)
library(gridExtra)
library(cowplot)
library(tidyverse)
library(clock)
library(broom)

#### Additional Functions ####
lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

darken <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col - (col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
#### END Additional Functions ####

# Reading in data
water_iso_full <- as_tibble(read.csv("pituffik_H2O_iso_2018_2019.csv", header=TRUE, sep=","))
water_iso_full$date <- as.Date(water_iso_full$date, format= "%d/%m/%Y")
water_iso_full$yr <- year(water_iso_full$date) # Extracting year
water_iso_full$doy <- yday(water_iso_full$date)
water_iso_full$type <- factor(water_iso_full$type, levels=c("lake", "pool", "stream", "surface flow", # Making so the order is always correct
                                                            "snow or ice", "precipitation rain", "precipitation snow"))
water_iso <- water_iso_full[water_iso_full$qc_flag != 1, ] # Removing data flagged for quality
water_iso$laketype <- factor(water_iso$laketype, levels=c("megapool", "headwater", "downstream", # Making so the order is always correct
                                                          "vale", "proglacial", "altered"))

# Pulling out lakes and pools with observations in all three periods
water_iso_periodcompare <- water_iso %>%
  filter(type_number < 3) %>% # Restricting to lake and pool data
  group_by(site_name) %>%
  filter(timing > 0) %>%
  filter(n() == 3) %>% # Only sites with all three periods
  ungroup()
#write_csv(water_iso_periodcompare, "water_iso_periodcompare.csv") # Remove comment marker to output file

# Pituffik GNIP precip data for LMWL data
pfk_gnip <- as_tibble(read.csv("pfk_gnip.csv", header=TRUE, sep=","))
pfk_gnip <- pfk_gnip[pfk_gnip$qc_flag != 1, ] # Removing data flagged for quality
pfk_gnip$date_start <- as.Date(pfk_gnip$date_start, format="%d/%m/%Y")
pfk_gnip$date_end <- as.Date(pfk_gnip$date_end, format="%d/%m/%Y")
pfk_lmwl_prcp <- lm(pfk_gnip$d2H~pfk_gnip$d18O)

# Adding GNIP data to water iso data for plotting comparisons
water_iso_plus_gnip <- bind_rows(water_iso, pfk_gnip[c("site", "d18O", "d2H", "dxs", "type")]) # Adding GNIP iso data
water_iso_plus_gnip$type <- factor(water_iso_plus_gnip$type, levels=c("lake", "pool", "stream", # Making so the order is always correct
                                                                      "surface flow", "snow or ice",
                                                                      "precipitation rain", "precipitation snow", "GNIP"))

#Finding unique sampled sites in database
unique_sites_full <- water_iso_full %>%
  distinct(latitude, longitude, .keep_all = T)
unique_sites_bygroup_full <- water_iso_full %>%
  group_by(type) %>%
  distinct(latitude, longitude, .keep_all = T)
unique_sites <- water_iso %>%
  distinct(latitude, longitude, .keep_all = T)
unique_sites_bygroup <- water_iso %>%
  group_by(type) %>%
  distinct(latitude, longitude, .keep_all = T)

# Setting plotting and analysis variables
vbl_clm <- c("d18O", "d2H", "dxs") # columns with plotting data
vbl_clm_index <- NA
for (i in 1:length(vbl_clm)) {
  vbl_clm_index[i] <- which(colnames(water_iso) == vbl_clm[i])
}

# Creating color-iso variable table
iso_color_names <- c("turquoise3", "firebrick", "mediumorchid4")
iso_color_index <- data.frame(colnames(water_iso[vbl_clm_index]), iso_color_names[seq(1,length(vbl_clm_index))])
colnames(iso_color_index) <- c("iso", "color")

lake <- water_iso[water_iso$type == "lake", ]
pool <- water_iso[water_iso$type == "pool", ]
surf_flow <- water_iso[water_iso$type == "surface flow", ]
stream <- water_iso[water_iso$type == "stream", ]
precip_rain <- water_iso[water_iso$type == "precipitation rain", ]
precip_snow <- water_iso[water_iso$type == "precipitation snow", ]

####============================####
####===CLIMATE DATA AND STATS===####
####============================####
# Pituffik weather data from Giovanni Muscari and USAF and water vapor isotopes from Akers 2020
pfk_isowx_day <- as_tibble(read.csv("pfk_iso_wx_day.csv", header=TRUE))
pfk_isowx_day$daybreak <- as.Date(pfk_isowx_day$daybreak, format= "%d/%m/%Y")
pfk_isowx_day <- pfk_isowx_day %>% rename(date = daybreak)
pfk_isowx_day$doy <- yday(pfk_isowx_day$date)
pfk_isowx_day$yr <- year(pfk_isowx_day$date)

# Pituffik daily weather obs 2000-2021 from USAF
pfk_wx_usaf_day <- as_tibble(read.csv("pituffik_USAF_wx_daily_QC.csv", header=TRUE)) # T in C, Prcp in mm
pfk_wx_usaf_day$date <- as.Date(pfk_wx_usaf_day$date, format= "%d/%m/%Y")
pfk_wx_usaf_day$doy <- yday(pfk_wx_usaf_day$date)
pfk_wx_usaf_day$wk <- week(pfk_wx_usaf_day$date)

# Loading PET data extracted for Pituffik from Singer et al., 2021
# Data extraction point coordinates (76.5 N 68.5 W), bilinear interpolation 
pfk_pet_1981_2023 <- as_tibble(read.csv("pituffik_pet_1981_2023.csv", header=TRUE)) # PET in mm/day
pfk_pet_1981_2023$date <- as.Date(pfk_pet_1981_2023$date, format="%Y-%m-%d")
pfk_isowx_day <- pfk_isowx_day %>%
  left_join(pfk_pet_1981_2023 %>% select(pet, date), by = "date")

# Monthly climate norms
pfk_wx_yrmo <- pfk_wx_usaf_day %>%
  group_by(yr, mo) %>% # This step needed due to summing prcp by month
  summarize(tavg_mo = mean(tavg, na.rm=FALSE), prcp_mo = mean(prcp, na.rm=FALSE), prcp_sum_mo = sum(prcp, na.rm=FALSE))
pfk_mo_clim <- pfk_wx_yrmo %>%
  group_by(mo) %>%
  summarize(tavg = mean(tavg_mo, na.rm=TRUE), prcp = mean(prcp_mo, na.rm=TRUE), prcp_sum = mean(prcp_sum_mo, na.rm=TRUE))
midmonth_doy <- c(15,46,74,105,135,166,196,227,258,288,319,349) # DOY for middle of month, for plotting alignment
pfk_mo_clim$doy <- midmonth_doy
pfk_pet_mo_clim <- pfk_pet_1981_2023 %>% # Calculating monthly means
  group_by(mo) %>%
  summarize(pet_mo = mean(pet))
pfk_mo_clim <- pfk_mo_clim %>% # Adding to Pituffik monthly climate norms
  left_join(pfk_pet_mo_clim, by="mo") %>%
  rename(pet = pet_mo)

# Day of year climate norms
pfk_wx_yrdoy <- pfk_wx_usaf_day %>%
  group_by(yr, doy) %>% # This step needed due to summing prcp by month
  summarize(tavg_doy = mean(tavg, na.rm=FALSE), prcp_doy = mean(prcp, na.rm=FALSE), prcp_sum_doy = sum(prcp, na.rm=FALSE))
pfk_doy_clim <- pfk_wx_yrdoy %>%
  group_by(doy) %>%
  summarize(tavg = mean(tavg_doy, na.rm=TRUE), prcp = mean(prcp_doy, na.rm=TRUE), prcp_sum = mean(prcp_sum_doy, na.rm=TRUE))
pfk_pet_doy_means <- pfk_pet_1981_2023 %>% # Calculating monthly means
  group_by(doy) %>%
  summarize(pet_doy = mean(pet))
pfk_doy_clim <- pfk_doy_clim %>% # Adding to Pituffik monthly climate norms
  left_join(pfk_pet_doy_means, by="doy") %>%
  rename(pet = pet_doy)

# Calculating sinusoidal climate norms to make infer a curve that smooths stochastic noise in means
# Sinusoidal regression function adapted from StackOverflow Glen_b https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
sinusoid_lm <- function(sinusoid_x, sinusoid_y, plot_graph, harmonic, period) {
  #ssp <- spectrum(sinusoid_y, plot=FALSE)  #Used if do not know period
  #per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
  if(period == "monthly") {
    per <- 12 #Period in months of year
  } else if (period == "weekly") {
    per <- 52.18
  } else if (period == "daily") {
    per <- 365.25
  }
  sinusoid_lm_1harm <- lm(sinusoid_y ~ sin(2*pi/per*sinusoid_x)+cos(2*pi/per*sinusoid_x)) # Basic sinusoidal regression
  sinusoid_range <- diff(range(sinusoid_y))
  sinusoid_lm_2harm <- lm(sinusoid_y ~ sin(2*pi/per*sinusoid_x)+cos(2*pi/per*sinusoid_x) + # Adding 2nd harmonic to improve the fit
                            sin(4*pi/per*sinusoid_x)+cos(4*pi/per*sinusoid_x))
  sinusoid_lm_3harm <- lm(sinusoid_y ~ sin(2*pi/per*sinusoid_x)+cos(2*pi/per*sinusoid_x) + # Adding another harmonic to improve the fit
                            sin(4*pi/per*sinusoid_x)+cos(4*pi/per*sinusoid_x) + 
                            sin(6*pi/per*sinusoid_x)+cos(6*pi/per*sinusoid_x))
  if(plot_graph == TRUE) {
    windows() # Plotting sinusoidal regression
    plot(sinusoid_y~sinusoid_x,ylim=c(min(sinusoid_y)-0.1*sinusoid_range,max(sinusoid_y)+0.1*sinusoid_range), type="l")
    lines(fitted(sinusoid_lm_1harm)~sinusoid_x,col="dodgerblue", lty=2)   # dashed blue line is sin fit
    lines(fitted(sinusoid_lm_2harm)~sinusoid_x,col="orange", lty=4)    # dot dash green line is periodic with second harmonic
    lines(fitted(sinusoid_lm_3harm)~sinusoid_x, col="violet")    # solid violet line is periodic with third harmonic
    legend("topleft", legend=c("Sine", "Second Harmonic", "Third Harmonic"),
           col=c("dodgerblue", "orange", "violet"), lty=c(2,4,1), bty="n")
  }
  
  if(harmonic==1) {
    return(sinusoid_lm_1harm)
  }
  if(harmonic==2) {
    return(sinusoid_lm_2harm)
  }
  if(harmonic==3) {
    return(sinusoid_lm_3harm)
  }
}

pfk_mo_clim_sinusoid <- tibble(mo = pfk_mo_clim$mo, # Monthly sinusoidal 
                               tavg_clim = fitted(sinusoid_lm(pfk_mo_clim$mo, pfk_mo_clim$tavg,
                                                              plot_graph = FALSE, harmonic = 3, period = "monthly")),
                               prcp_clim = exp(fitted(sinusoid_lm(pfk_mo_clim$mo, log(pfk_mo_clim$prcp),
                                                                  plot_graph = FALSE, harmonic = 3, period = "monthly"))),
                               pet_clim = exp(fitted(sinusoid_lm(pfk_mo_clim$mo, log(pfk_mo_clim$pet),
                                                                 plot_graph = FALSE, harmonic = 3, period = "monthly"))))

pfk_doy_clim_sinusoid <- tibble(doy = pfk_doy_clim$doy, # Daily sinusoidal
                                tavg_clim = fitted(sinusoid_lm(pfk_doy_clim$doy, pfk_doy_clim$tavg,
                                                               plot_graph = FALSE, harmonic = 3, period = "daily")),
                                #prcp_clim = exp(fitted(sinusoid_lm(pfk_doy_clim$doy, log(pfk_doy_clim$prcp), # Unsolvable regression due to data spread
                                #                                          plot_graph = TRUE, harmonic = 3, period = "daily))),
                                pet_clim = exp(fitted(sinusoid_lm(pfk_doy_clim$doy, log(pfk_doy_clim$pet),
                                                                  plot_graph = FALSE, harmonic = 3, period = "daily"))))

# Long term climate stats
pfk_wx_usaf_day %>% # 95% confidence intervals for climate data
  group_by(yr) %>%
  summarize(tavg_yr = mean(tavg, na.rm=TRUE), prcp_yr = mean(prcp, na.rm=TRUE),
            prcp_sum_yr = sum(prcp, na.rm=TRUE)) %>%
  ungroup() %>%
  reframe(across(tavg_yr:prcp_sum_yr, CI))

pfk_pet_1981_2023 %>% # PET totals in summer for each year in the database
  filter(sn == "sum") %>%
  group_by(yr) %>%
  summarize(pet_sum_yr = sum(pet, na.rm=TRUE))
  

# Climate comparisons between summers 2018 and 2019
yr_index <- c(2018, 2019)
wxmeans_index <- c("tavg", "pres", "rh_ice", "nao", "ao", "seaice", "Eness", "wind_sp")
wxsums_index <- c("prcp_H2O", "prcp_snow", "pet")

pfk_isowx_day %>% # Mean summer values for selected climate variables
  filter(yr %in% yr_index) %>%
  filter(sn == "sum") %>%
  group_by(yr) %>%
  summarize(across(all_of(wxmeans_index), \(x) mean(x, na.rm = TRUE)))

pfk_isowx_day %>% # Median summer values for selected climate variables
  filter(yr %in% yr_index) %>%
  filter(sn == "sum") %>%
  group_by(yr) %>%
  summarize(across(all_of(wxmeans_index), \(x) median(x, na.rm = TRUE)))

pfk_isowx_day %>% # Summed summer values for precipitation variables
  filter(yr %in% yr_index) %>%
  filter(sn == "sum") %>%
  group_by(yr) %>%
  summarize(across(all_of(wxsums_index), \(x) sum(x, na.rm = TRUE)))

####========END CLIMATE======####


####============================####
####========BASIC STATS=========####
####============================####
nrow(water_iso_full) # Number of total samples, before QC removal
water_iso_full %>% # Number of samples by category, before QC removal
  count(type)
nrow(unique_sites_full) # Number of uniquely sampled sites, before QC removal
unique_sites_bygroup_full %>% # Number of unique sampling sites by category, before QC removal
  count(type)

nrow(water_iso) # Number of total samples, after QC removal
water_iso %>% # Number of samples by category, after QC removal
  count(type)
nrow(unique_sites) # Number of uniquely sampled sites, after QC removal
unique_sites_bygroup %>% # Number of unique sampling sites by category, after QC removal
  count(type)

# Total isotope means and CI
water_iso %>%
  reframe(across(all_of(vbl_clm), CI))

# Total isotope 1*SD
water_iso %>%
  summarize(across(all_of(vbl_clm), sd))

# Total isotope 2*SD
water_iso %>%
  summarize(across(all_of(vbl_clm), sd))*2

# Sorting sample sites by number of samples
print(water_iso %>%
        count(site_name) %>%
        arrange(desc(n)),
      n=20)

# Local water line parameters and stat values
lwl_bytype <- water_iso_plus_gnip %>% # Running regression by type
  group_by(type) %>%
  nest() %>%
  mutate(lwl = map(data, ~lm(d2H~d18O, data=.)))

lwl_slope <- 
  lwl_bytype %>%
  mutate(tidyit = map(lwl, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "d18O") %>%
  dplyr::select(type, estimate, std.error)
colnames(lwl_slope) <- c("type", "slope", "se_slope")
lwl_slope$conf_int_slope <- lwl_slope$se_slope*qnorm(0.975)

lwl_intercept <- 
  lwl_bytype %>%
  mutate(tidyit = map(lwl, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(type, estimate, std.error)
colnames(lwl_intercept) <- c("type", "intercept", "se_intercept")
lwl_intercept$conf_int_intercept <- lwl_intercept$se_intercept*qnorm(0.975)

lwl_r2 <- 
  lwl_bytype %>%
  mutate(glanceit = map(lwl, broom::glance)) %>%
  unnest(glanceit) %>%
  dplyr::select(type, r.squared, p.value, nobs)

lwl <- lwl_slope %>% # Combining all data
  inner_join(lwl_intercept, by="type") %>%
  inner_join(lwl_r2, by="type")
print(lwl)
#write.csv(lwl,"lwl_regression.csv", row.names=FALSE) # Remove comment marker to output file

# GNIP weighted mean isotopic values
gnip_wtmean <- as_tibble(t(c(weighted.mean(pfk_gnip$d18O, pfk_gnip$precip_amt, na.rm=TRUE),
                             weighted.mean(pfk_gnip$d2H, pfk_gnip$precip_amt, na.rm=TRUE))), .name_repair = "unique")
colnames(gnip_wtmean) <- c("d18O", "d2H")
print(gnip_wtmean)

# Finding intersections between LEL and LMWL
type_index <- c("lake", "pool", "stream", "surface flow", "snow or ice")
lwl_intersect <- as_tibble(matrix(nrow=length(type_index),ncol=2), .name_repair = "unique")
lwl_intersect <- as_tibble(cbind(type_index, lwl_intersect))
colnames(lwl_intersect) <- c("type", "d18O", "d2H")

for (i in 1:length(type_index)) {
  coef_iter <- rbind(coef(lwl_bytype$lwl[lwl_bytype$type=="GNIP"][[1]]), # Coefficient matrix
                     coef(lwl_bytype$lwl[lwl_bytype$type==type_index[i]][[1]]))
  lwl_intersect[i,c(2,3)] <- t(c(-solve(cbind(coef_iter[,2],-1)) %*% coef_iter[,1]))
}
print(lwl_intersect)

####========END BASIC STATS======####

####============================####
####=====Spatial Analyses ======####
####============================####
#==Lake type analysis
# Making subsets of lake data grouped by sample period and main lake zone
lakeset <- water_iso %>%
  filter(type == "lake") %>% # Getting only lakes
  filter(!is.na(timing)) # Keeping only samples from 3 set periods

lakeset_all_byperiod <- water_iso %>%
  filter(type == "lake") %>% # Getting only lakes
  filter(!is.na(timing)) %>% # Keeping only samples from 3 set periods
  group_by(timing) %>% # Grouping by sampling period
  nest()

lakeset_main_byperiod <- water_iso %>%
  filter(type == "lake") %>% # Getting only lakes
  filter(!is.na(timing)) %>% # Keeping only samples from 3 set periods
  filter(main_lakes == 1) %>% # Keeping only main lakes
  group_by(timing) %>% # Grouping by sampling period
  nest()

lakeset_byperiod <- rbind(lakeset_all_byperiod, lakeset_main_byperiod)

# Lake LEL by laketype for lakes with >=3 samples
multilake <- water_iso %>% # All lakes with >=3 samples taken
  filter(type=="lake") %>%
  group_by(site_name) %>%
  filter(n() >= 3) %>%
  filter(is.na(laketype) == FALSE)

lel_lake_bylake <- multilake %>% # Running regression by type
  group_by(site_name) %>%
  nest() %>%
  mutate(lel = map(data, ~lm(d2H~d18O, data=.)))

lel_slope <- 
  lel_lake_bylake %>%
  mutate(tidyit = map(lel, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "d18O") %>%
  dplyr::select(site_name, estimate, std.error)
colnames(lel_slope) <- c("site_name", "slope", "se_slope")
lel_slope$conf_int_slope <- lel_slope$se_slope*qnorm(0.975)

lel_intercept <- 
  lel_lake_bylake %>%
  mutate(tidyit = map(lel, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(site_name, estimate, std.error)
colnames(lel_intercept) <- c("site_name", "intercept", "se_intercept")
lel_intercept$conf_int_intercept <- lel_intercept$se_intercept*qnorm(0.975)

lel_r2 <- 
  lel_lake_bylake %>%
  mutate(glanceit = map(lel, broom::glance)) %>%
  unnest(glanceit) %>%
  dplyr::select(site_name, r.squared, p.value, nobs)

lel_lake_bylake_params <- lel_slope %>% # Combining all data
  inner_join(lel_intercept, by="site_name") %>%
  inner_join(lel_r2, by="site_name")

lel_lake_bylake_params$laketype <- distinct(multilake, site_name, .keep_all = TRUE) %>%
  ungroup() %>%
  dplyr::pull(laketype)
print(lel_lake_bylake_params) # Lake LEL parameters by each lake

lel_lake_bylake_bylaketype_slope <- lel_lake_bylake_params %>%
  group_by(laketype) %>%
  summarize(count=n(),mean_slope=mean(slope, na.rm=TRUE), conf_int_slope=sd(slope)/sqrt(n())*qnorm(0.975))
print(lel_lake_bylake_bylaketype_slope) # Lake LEL slopes averaged by laketype

# Comparing overall mean LEL to means of individual lakes
mean_LEL_bylake <- lel_lake_bylake_params %>% # Mean LEL slope of individual lake LELs
  ungroup() %>%
  summarize(mean_slope=mean(slope, na.rm=TRUE),
            conf_int_slope=sd(slope, na.rm=TRUE)/sqrt(sum(!is.na(slope)))*qnorm(0.975))

multilake_naomit <- multilake %>%
  filter (site_name != "Interramp Proglacial Lake") # Removing to match bylake slope dataset
single_LEL_lm <- summary(lm(multilake_naomit$d2H~multilake_naomit$d18O)) # Single LEL for all lake data combined
single_LEL <- tibble(mean_slope = single_LEL_lm$coefficients[2],
                     conf_int_slope = single_LEL_lm$coefficients[4]*qnorm(0.975))
                        

LEL_method_compare <- rbind(mean_LEL_bylake, single_LEL)
LEL_method_compare$method <- c("Mean of individual lake LELs", "Single LEL of all lake data combined")
print(LEL_method_compare) # Comparison of two LEL methods

# Finding intersections between lake LEL and LMWL
lake_index <- lel_lake_bylake$site_name
lake_lel_intersect <- as_tibble(matrix(nrow=length(lake_index),ncol=2), .name_repair = "unique")
lake_lel_intersect <- as_tibble(cbind(lake_index, lake_lel_intersect))
colnames(lake_lel_intersect) <- c("site_name", "d18O", "d2H")

for (i in 1:length(lake_index)) {
  coef_iter <- rbind(coef(lwl_bytype$lwl[lwl_bytype$type=="GNIP"][[1]]), # Coefficient matrix
                     coef(lel_lake_bylake$lel[lel_lake_bylake$site_name==lake_index[i]][[1]]))
  lake_lel_intersect[i,c(2,3)] <- t(c(-solve(cbind(coef_iter[,2],-1)) %*% coef_iter[,1]))
}
lake_lel_intersect$laketype <- lel_lake_bylake_params$laketype

lake_lel_intersect_bylaketype <- lake_lel_intersect %>%
  group_by(laketype) %>%
  summarize(count=n(), mean_d18O=mean(d18O, na.rm=TRUE)*1000, conf_int_d18O=sd(d18O)/sqrt(n())*1000*qnorm(0.975),
            mean_d2H=mean(d2H, na.rm=TRUE)*1000, conf_int_d2H=sd(d2H)/sqrt(n())*1000*qnorm(0.975))
print(lake_lel_intersect_bylaketype)

#==Lake E/I ratios and source mixing modeling
library(isoWater)
library(beepr)
# NOTE: This code is included to show how the E/I and source data was created.
#   Because it is Bayesian modeling, slight differences in the output may occur from
#   run to run. Thus we load one file that was the output of one run to ensure
#   consistent results and stats in later sections.

# # Setting a meteoric water line for Pituffik
# pfk_gnip_mwl_data <- na.omit(data.frame(pfk_gnip$d2H, pfk_gnip$d18O)) # Based off of GNIP data
# MWL_pfk <- mwl(pfk_gnip_mwl_data)
# 
# # Pulling lake LEL data to estimate source water isotopic composition
# lel_input <- lel_lake_bylake_params %>%
#   mutate(sd_slope = se_slope*sqrt(nobs), sd_intercept = se_intercept*sqrt(nobs)) %>%
#   select(site_name, slope, sd_slope, intercept, sd_intercept)
# lel_lake_iso <- na.omit(inner_join(water_iso %>% select(site_name, date, d18O, d2H),lel_input, by="site_name"))
# 
# # Looping through individual lake samples to estimate source water isotopic composition (Bowen et al 2018)
# mwl_output <- as_tibble(matrix(nrow=nrow(lel_lake_iso), ncol=5))
# mwl_output[1] <- lel_lake_iso$site_name
# for(i in 1:nrow(lel_lake_iso)) {
#   mwl_sample_iter <- mwlSource(iso(lel_lake_iso$d2H[i], lel_lake_iso$d18O[i], 0.6, 0.1), MWL_pfk, slope=c(lel_lake_iso$slope[i], lel_lake_iso$sd_slope[i]), ngens = 1e6)
#   mwl_output[i,2] <- mwl_sample_iter$summary[c(4),c(1)]
#   mwl_output[i,3] <- mwl_sample_iter$summary[c(4),c(2)]
#   mwl_output[i,4] <- mwl_sample_iter$summary[c(5),c(1)]
#   mwl_output[i,5] <- mwl_sample_iter$summary[c(5),c(2)]
# }
# colnames(mwl_output) <- c("site_name", "d18O_inflow", "d18O_inflow_sd", "d2H_inflow", "d2H_inflow_sd")
# 
# iso_inflow_bayesian <- mwl_output %>% # Organizing the Bowen model output to get isotopic estimates
#   group_by(site_name) %>%
#   summarize(d18O_inflow = mean(d18O_inflow),
#             d18O_inflow_sd = sd(mapply(rnorm, n=100000, mean=d18O_inflow, sd=d18O_inflow_sd)), # Generating distributions for each lake sample based on mean and SD, then finding SD of combined distributions
#             d2H_inflow = mean(d2H_inflow),
#             d2H_inflow_sd = sd(mapply(rnorm, n=100000, mean=d2H_inflow, sd=d2H_inflow_sd)),
#             dxs_inflow = mean(d2H_inflow-8*d18O_inflow),
#             dxs_inflow_sd = sd(mapply(rnorm, n=100000, mean=d2H_inflow, sd=d2H_inflow_sd)-
#                                  8*mapply(rnorm, n=100000, mean=d18O_inflow, sd=d18O_inflow_sd))) # Note that dxs sd is large
# beep(4)
#  
# #==== Estimating E/I for select lakes with LELs (Gonfiantini 1986 model, following Gibson et al 2016)
# # For each sampling date, calculate the PET-weighted mean tavg, rh, d18Oatm, d2Hatm, dxsatm from Akers et al 2020 data
# lake_sample_isowx <- tibble(lake_iso_inflow %>% select(site_name, date)) %>% # Pulling out sampling dates for each lake
#   add_column(pet_tavg = NA, pet_rh_ice = NA, pet_d18Oatm = NA, pet_d2Hatm = NA, pet_dxsatm = NA)
# for(i in 1:nrow(lake_sample_isowx)) { # Looping for each lake isotope sample
#   lake_sample_isowx[i,3:7] <- pfk_isowx_day %>%
#     filter(sn=="sum") %>% 
#     filter(yr==year(lake_sample_isowx$date[i]) & date <= lake_sample_isowx$date[i]) %>% # Pulling only summer dates of same year up through sampling date
#     mutate(dailyfrac_pet = (pet)/sum(pet)) %>% # Calculating the fraction of total PET per day
#     mutate(pet_tavg = dailyfrac_pet*tavg, pet_rh_ice = dailyfrac_pet*rh_ice, # Caculating daily fractional contributions of select variable to the overall weighted mean
#            pet_d18O = dailyfrac_pet*d18O, pet_d2H = dailyfrac_pet*d2H,
#            pet_dxs = dailyfrac_pet*dxs) %>%
#     summarize(pet_tavg = sum(pet_tavg), pet_rh_ice = sum(pet_rh_ice), pet_d18Oatm = sum(pet_d18O), # Determining PET-weighted means for the summer up to the date of sampling
#               pet_d2Hatm = sum(pet_d2H), pet_dxsatm = sum(pet_dxs))
# }
# # Creating a working tibble for model data input and results for iso mass balance model (imbm)
# lake_iso_imbm <- bind_cols(lake_iso_inflow, lake_sample_isowx %>% select(-site_name,-date)) # Joining PET weighted means to lake inflow data
# 
# # Organizing parameters for iso mass balance model
# T_summer <- lake_iso_imbm$pet_tavg + 273.15 # Mean air temperature weighted by daily PET, from 01 Jun through sample date
# rel_hum <- lake_iso_imbm$pet_rh_ice/100 #  Mean relative humidity weighted by daily PET, from 01 Jun through sample date
# iso_atm_d18O <- lake_iso_imbm$pet_d18Oatm #  Mean atmospheric vapor d18O weighted by daily PET, from 01 Jun through sample date
# iso_atm_d2H <- lake_iso_imbm$pet_d2Hatm #  Mean atmospheric vapor d2H weighted by daily PET, from 01 Jun through sample date
# 
# # Horita Wesolowski 1994
# alpha_plus_d18O <- exp(-7.685/10^3 + 6.7123/T_summer - 1666.4/T_summer^2 + 350410/T_summer^3)
# alpha_plus_d2H <- exp(1158.8 * (T_summer^3/(10^12)) - 1620.1*((T_summer^2)/(10^(9))) + 794.84*(T_summer/10^6) - 161.04/10^3 + 2999200/T_summer^3)
# 
# eps_plus_d18O <- (alpha_plus_d18O-1)*1000
# eps_plus_d2H <- (alpha_plus_d2H-1)*1000
# 
# eps_k_d18O <- (1-rel_hum) * 14.2   # Gat 1995, Gonfiantini 1986
# eps_k_d2H <- (1-rel_hum)^2 * 124.66 - (1-rel_hum) * 9.5737 + 0.6326 # Estimated from Fig 5 in Zuber 1983, following Biggs 2015 and Cui 2016
# 
# m_d18O <- (rel_hum-(10^-3)*(eps_k_d18O+eps_plus_d18O/alpha_plus_d18O)) / (1-rel_hum + (10^-3)*eps_k_d18O)
# d18O_star <- (rel_hum*iso_atm_d18O*1000+eps_k_d18O+eps_plus_d18O/alpha_plus_d18O) / (rel_hum-(10^-3)*(eps_k_d18O+eps_plus_d18O/alpha_plus_d18O))
# iso_lake_d18O <- lake_iso_imbm %>% pull(d18O)
# iso_inflow_d18O <- lake_iso_imbm %>% pull(d18O_inflow)
# iso_inflow_d18O_sd <- lake_iso_imbm %>% pull(d18O_inflow_sd)
# evap_d18O <- -1 * m_d18O*(d18O_star-iso_lake_d18O*1000)+iso_lake_d18O*1000 # Additional calculation for dxs E/I
# E_I_d18O_bayes <- (iso_lake_d18O*1000 - iso_inflow_d18O*1000) / (m_d18O*(d18O_star-iso_lake_d18O*1000))
# E_I_d18O_1sd_bayes <- abs(((iso_lake_d18O*1000 - (iso_inflow_d18O+iso_inflow_d18O_sd)*1000) / (m_d18O*(d18O_star-iso_lake_d18O*1000)))-E_I_d18O_bayes)
# 
# m_d2H <- (rel_hum-(10^-3)*(eps_k_d2H+eps_plus_d2H/alpha_plus_d2H)) / (1-rel_hum + (10^-3)*eps_k_d2H)
# d2H_star <- (rel_hum*iso_atm_d2H*1000+eps_k_d2H+eps_plus_d2H/alpha_plus_d2H) / (rel_hum-(10^-3)*(eps_k_d2H+eps_plus_d2H/alpha_plus_d2H))
# iso_lake_d2H <- lake_iso_imbm %>% pull(d2H)
# iso_inflow_d2H <- lake_iso_imbm %>% pull(d2H_inflow)
# iso_inflow_d2H_sd <- lake_iso_imbm %>% pull(d2H_inflow_sd)
# evap_d2H <- -1 * m_d2H*(d2H_star-iso_lake_d2H*1000)+iso_lake_d2H*1000 # Additional calculation for dxs E/I
# E_I_d2H_bayes <- (iso_lake_d2H*1000 - iso_inflow_d2H*1000) / (m_d2H*(d2H_star-iso_lake_d2H*1000))
# E_I_d2H_1sd_bayes <- abs(((iso_lake_d2H*1000 - (iso_inflow_d2H+iso_inflow_d2H_sd)*1000) / (m_d2H*(d2H_star-iso_lake_d2H*1000)))-E_I_d2H_bayes)
# 
# dxs_star <- d2H_star - d18O_star*8  
# iso_lake_dxs <- lake_iso_imbm %>% pull(dxs)
# iso_inflow_dxs <- lake_iso_imbm %>% pull(dxs_inflow)
# iso_inflow_dxs_sd <- lake_iso_imbm %>% pull(dxs_inflow_sd)
# evap_dxs <- evap_d2H - 8 * evap_d18O
# E_I_dxs_bayes <- (iso_inflow_dxs*1000 - iso_lake_dxs*1000) / ((evap_dxs-iso_lake_dxs*1000))
# E_I_dxs_1sd_bayes <- abs((((iso_inflow_dxs+iso_inflow_dxs_sd)*1000 - iso_lake_dxs*1000) / (evap_dxs-iso_lake_dxs*1000))-E_I_dxs_bayes)
# 
# lake_e_i <- tibble(lake_iso_inflow %>% select(sample_id, site_name, date),
#                    E_I_d18O_bayes, E_I_d18O_1sd_bayes, E_I_d2H_bayes,
#                    E_I_d2H_1sd_bayes, E_I_dxs_bayes, E_I_dxs_1sd_bayes)
# 
# lake_e_i_round <- tibble(lake_e_i%>%select(sample_id, site_name, date), round(lake_e_i%>%select(-sample_id, -site_name, -date),3)) #Rounding to 3 digits for output csv
# print(lake_e_i_round, n=nrow(lake_iso_inflow))
# 
# lake_e_i_fulldata <- inner_join(water_iso, lake_e_i_round, by=c("sample_id", "date", "site_name"))
# #write_csv(lake_e_i_fulldata, "pituffik_lake_iso_e_i.csv")
# 
# # Calculating seasonal precipitation inputs from inflow isotopic composition
# library(simmr)
# 
# season_index <- tibble(c(1:12), # Tibble of months with seasonal and winter-summer assignments
#                        c(rep("win",2),rep("spr",3), rep("sum",3), rep("aut",3), rep("win",1)),
#                        c(rep("freeze",5),rep("thaw",3), rep("freeze",4))) # Freeze = months predominantly below freezing air temperature (09 Sep-04 Jun)
# colnames(season_index) <- c("mo", "sn", "freezethaw")
# 
# pfk_gnip_nm <- pfk_gnip %>% # Creating a tibble of GNIP data with seasonal identifiers
#   mutate(mo = month(date_start)) %>%
#   drop_na(d18O, d2H, dxs) %>%
#   inner_join(season_index, by="mo")
# 
# pfk_gnip_sn_iso_mean <- pfk_gnip_nm %>% # Seasonal GNIP mean isotopic values (Not used in main study but here for alternative examination)
#   group_by(sn) %>%
#   summarize(across(c(d18O, d2H, dxs), mean))
# 
# pfk_gnip_sn_iso_sd <- pfk_gnip_nm %>% # Seasonal GNIP standard deviations of isotopic values (Not used in main study but here for alternative examination)
#   group_by(sn) %>%
#   summarize(across(c(d18O, d2H, dxs), sd))
# 
# pfk_gnip_freezethaw_iso_mean <- pfk_gnip_nm %>% # Freeze-thaw GNIP mean isotopic values
#   group_by(freezethaw) %>%
#   summarize(across(c(d18O, d2H, dxs), mean))
# 
# pfk_gnip_freezethaw_iso_sd <- pfk_gnip_nm %>% # Freeze-thaw GNIP standard deviations of isotopic values
#   group_by(freezethaw) %>%
#   summarize(across(c(d18O, d2H, dxs), sd))
# 
# # Executing source attribution mixing model for lake inflow isotopic compositions in freeze-thaw seasons
# pfk_freezethaw_simmr_in <- simmr_load( # Loading in input data
#   mixtures = iso_inflow_bayesian %>% select(d18O_inflow, d2H_inflow, dxs_inflow),
#   source_names = pfk_gnip_freezethaw_iso_mean %>% pull(freezethaw),
#   source_means = pfk_gnip_freezethaw_iso_mean %>% select(d18O, d2H, dxs),
#   source_sds = pfk_gnip_freezethaw_iso_sd %>% select(d18O, d2H, dxs),
#   group = iso_inflow_bayesian %>% pull(site_name)
# )
# pfk_freezethaw_simmr_out <- simmr_mcmc(pfk_freezethaw_simmr_in) # Running mixing model
# 
# # Organizing output data
# pfk_iso_inflow_freezethaw_mixmodel_full <- list()
# pfk_iso_inflow_freezethaw_mixmodel <- as_tibble(matrix(ncol=5, nrow=length(pfk_freezethaw_simmr_out$output)))
# for(i in 1:length(pfk_freezethaw_simmr_out$output)){
#   pfk_iso_inflow_freezethaw_mixmodel_full[[i]] <- pfk_freezethaw_simmr_out$output[[i]]$BUGSoutput$summary
#   pfk_iso_inflow_freezethaw_mixmodel[i,2] <- pfk_iso_inflow_freezethaw_mixmodel_full[[i]][2,1]
#   pfk_iso_inflow_freezethaw_mixmodel[i,3] <- pfk_iso_inflow_freezethaw_mixmodel_full[[i]][2,2]
#   pfk_iso_inflow_freezethaw_mixmodel[i,4] <- pfk_iso_inflow_freezethaw_mixmodel_full[[i]][3,1]
#   pfk_iso_inflow_freezethaw_mixmodel[i,5] <- pfk_iso_inflow_freezethaw_mixmodel_full[[i]][3,2]
# }
# names(pfk_iso_inflow_freezethaw_mixmodel_full) <- iso_inflow_bayesian %>% pull(site_name)
# pfk_iso_inflow_freezethaw_mixmodel[,1] <- iso_inflow_bayesian %>% select(site_name)
# colnames(pfk_iso_inflow_freezethaw_mixmodel) <- c("site_name", "freeze_frac", "freeze_frac_sd", "thaw_frac", "thaw_frac_sd")
# 
# pfk_iso_inflow_freezethaw_mixmodel_round <- tibble(pfk_iso_inflow_freezethaw_mixmodel%>%select(site_name), #Rounding to 3 digits for output csv
#                                                    round(pfk_iso_inflow_freezethaw_mixmodel%>%select(-site_name),3))
# print(pfk_iso_inflow_freezethaw_mixmodel_round, n=nrow(lake_iso_inflow))
# iso_inflow_bayesian_round <- tibble(iso_inflow_bayesian%>%select(site_name), #Rounding to 4 digits for output csv
#                                     round(iso_inflow_bayesian%>%select(-site_name),4))
# pfk_iso_inflow <- inner_join(iso_inflow_bayesian_round, pfk_iso_inflow_freezethaw_mixmodel_round, by= "site_name")
# #write_csv(pfk_iso_inflow, "pituffik_lake_iso_inflow_bylake.csv", col_names=TRUE)

# Loading in tibbles of data output by code commented out above
pfk_iso_inflow <- read_csv("pituffik_lake_iso_inflow_bylake.csv", show_col_types = FALSE)
pfk_iso_inflow <- pfk_iso_inflow %>% # QC removing two lakes without LELs and poor inflow results
  filter(site_name != "Lake Tuto" & site_name != "Ice Ramp Base Pond")
pfk_iso_inflow_fulldata <- inner_join(water_iso, pfk_iso_inflow, by="site_name")
lake_e_i_fulldata <- read_csv("pituffik_lake_iso_e_i.csv", show_col_types = FALSE)
lake_e_i_fulldata$laketype <- factor(lake_e_i_fulldata$laketype, levels=c("megapool", "headwater", "downstream", # Making so the order is always correct
                                                          "vale", "proglacial", "altered"))
lake_e_i_fulldata <- lake_e_i_fulldata %>% # QC removing two lakes without LELs and poor inflow results
  filter(site_name != "Lake Tuto" & site_name != "Ice Ramp Base Pond")
lake_e_i_fulldata <- lake_e_i_fulldata %>%
  left_join(pfk_iso_inflow, by="site_name")

# Mean and CI summaries of inflow and E/I data
inflow_stat_clm <- c("d18O_inflow", "d2H_inflow", "dxs_inflow", "freeze_frac", "thaw_frac")
pfk_iso_inflow %>%
           reframe(across(all_of(inflow_stat_clm), CI))
e_i_stat_clm <- c("E_I_d18O_bayes", "E_I_d2H_bayes", "E_I_dxs_bayes")
lake_e_i_fulldata %>%
  reframe(across(all_of(e_i_stat_clm), median))

# Linear regressions between estimated source waters and sampled lake water isotopes.
inflow_clm <- c("d18O_inflow", "d2H_inflow", "dxs_inflow")
inflow_obs_lake_param <- list()
for (i in 1:length(vbl_clm)){
  inflow_obs_lake_iter <- pfk_iso_inflow_fulldata %>% # Running regression by type
    filter(!is.na(timing)) %>% # Removing Lake Potato and Power Lake samples in between main sampling periods
    group_by(timing) %>%
    nest() %>%
    mutate(inflow_obs = map(data, ~lm(.[[vbl_clm[i]]]~.[[inflow_clm[i]]])))
  
  inflow_obs_slope_iter <- 
    inflow_obs_lake_iter %>%
    mutate(tidyit = map(inflow_obs, broom::tidy)) %>%
    unnest(tidyit) %>%
    filter(term == ".[[inflow_clm[i]]]") %>%
    dplyr::select(timing, estimate, std.error)
  colnames(inflow_obs_slope_iter) <- c("timing", "slope", "se_slope")
  inflow_obs_slope_iter$conf_int_slope <- inflow_obs_slope_iter$se_slope*qnorm(0.975)
  
  inflow_obs_intercept_iter <- 
    inflow_obs_lake_iter %>%
    mutate(tidyit = map(inflow_obs, broom::tidy)) %>%
    unnest(tidyit) %>%
    filter(term == "(Intercept)") %>%
    dplyr::select(timing, estimate, std.error)
  colnames(inflow_obs_intercept_iter) <- c("timing", "intercept", "se_intercept")
  inflow_obs_intercept_iter$conf_int_intercept_iter <- inflow_obs_intercept_iter$se_intercept*qnorm(0.975)
  
  inflow_obs_r2_iter <- 
    inflow_obs_lake_iter %>%
    mutate(glanceit = map(inflow_obs, broom::glance)) %>%
    unnest(glanceit) %>%
    dplyr::select(timing, r.squared, p.value, nobs)
  
  inflow_obs_lake_params_iter <- inflow_obs_slope_iter %>% # Combining all data
    inner_join(inflow_obs_intercept_iter, by="timing") %>%
    inner_join(inflow_obs_r2_iter, by="timing")
  
  inflow_obs_lake_param[[i]] <- inflow_obs_lake_params_iter
}
names(inflow_obs_lake_param) <- vbl_clm
print(inflow_obs_lake_param) # Regression results

# Individual and multiple regression of source water isotopes vs. environmental factors
lake_inflow_regr_data <- pfk_iso_inflow_fulldata
lake_inflow_regr_data$log_area <- log(lake_inflow_regr_data$lake_surface_area) # Log transformation of lake surface area
lake_inflow_regr_data$log_lakeshed_area <- log(lake_inflow_regr_data$lakeshed_area) # Log transformation of lake watershed area
envparam_names <- c("elevation", "log_area", "log_lakeshed_area", "dist_gris", "dist_ocean", "doy") # Selected environmental parameters
inflow_clm <- c("d18O_inflow", "d2H_inflow", "dxs_inflow")

inflow_data <- lake_inflow_regr_data %>% filter(timing==3)
corr_inflow_r_iter <- as_tibble(matrix(ncol=length(inflow_clm), nrow=length(envparam_names), dimnames=list(NULL,inflow_clm)))
corr_inflow_pval_iter <- as_tibble(matrix(ncol=length(inflow_clm), nrow=length(envparam_names), dimnames=list(NULL,inflow_clm)))
for(j in 1:length(inflow_clm)){
  for(k in 1:length(envparam_names)){
    corr_inflow_r_iter[k,j] <- cor.test(inflow_data[[inflow_clm[j]]], inflow_data[[envparam_names[k]]])$estimate
    corr_inflow_pval_iter[k,j] <- cor.test(inflow_data[[inflow_clm[j]]], inflow_data[[envparam_names[k]]])$p.value
  }
}

# Organizing correlation output
iso_inflow_corr <- tibble(envparam = envparam_names, n_lakes = nrow(inflow_data),
                          r_d18O_inflow = corr_inflow_r_iter$d18O_inflow, pval_d18O_inflow = corr_inflow_pval_iter$d18O_inflow,
                          r_d2H_inflow = corr_inflow_r_iter$d2H_inflow, pval_d2H_inflow = corr_inflow_pval_iter$d2H_inflow,
                          r_dxs_inflow = corr_inflow_r_iter$dxs_inflow, pval_dxs_inflow = corr_inflow_pval_iter$dxs_inflow)
print(iso_inflow_corr)

# Multiple regression results for correlated parameters
multlm_final_envparam <- list(d18O_inflow=c("dist_ocean", "log_area", "log_lakeshed_area"), # Important parameters based on correlations
                              d2H_inflow=c("dist_ocean", "log_area", "log_lakeshed_area"),
                              dxs_inflow=c("dist_ocean", "log_area", "log_lakeshed_area"))
multlm_final_inflow_coef <- list()
multlm_final_inflow_strength <- as_tibble(matrix(nrow=length(inflow_clm), ncol=3, dimnames=list(NULL,c("fstat", "fstat_pval", "adj_r2"))))
for (i in 1:length(inflow_clm)) {
  iso_iter <- inflow_data[[inflow_clm[i]]]
  multlm_formula_iter <- paste("iso_iter~",paste(multlm_final_envparam[[inflow_clm[i]]],collapse="+"), sep="") # formula for multiple regressions
  multlm_final_inflow_iter <- summary(lm(as.formula(multlm_formula_iter), data=inflow_data))
  multlm_final_inflow_coef[[i]] <- multlm_final_inflow_iter$coefficients
  multlm_final_inflow_coef[[i]] <- bind_cols(signif(multlm_final_inflow_coef[[i]][,1:3],2), signif(multlm_final_inflow_coef[[i]][,4],1))
  names(multlm_final_inflow_coef[[i]]) <- c("coef_value", "st_err", "tval", "pval")
  multlm_final_inflow_coef[[i]] <- multlm_final_inflow_coef[[i]] %>%
    add_column(parameter = c("intercept", multlm_final_envparam[[i]]), .before = 1)
  multlm_final_inflow_strength[i,] <- as.list(c(multlm_final_inflow_iter$fstatistic[[1]],
                                                pf(multlm_final_inflow_iter$fstatistic[[1]], 
                                                   multlm_final_inflow_iter$fstatistic[[2]],
                                                   multlm_final_inflow_iter$fstatistic[[3]],
                                                   lower.tail = FALSE),
                                                multlm_final_inflow_iter$adj.r.squared))
}
names(multlm_final_inflow_coef) <- inflow_clm
multlm_final_inflow_strength <- multlm_final_inflow_strength %>%
  add_column(iso = inflow_clm, .before = 1)

print(multlm_final_inflow_coef)
print(multlm_final_inflow_strength)


#==Dendrogram analysis of lake types
library(ggdendro)
library(dendextend)
library(ape)

# Hierarchical clustering of period 3 lakes based on d18O and dxs
lake_per3 <- lakeset %>% # Selecting only lakes from period 3
  filter(timing == 3) 
lake_isoz_per3 <-lake_per3 %>%
  dplyr::select(d18O, dxs)
lake_isoz_per3 <- as.data.frame(lake_isoz_per3) # Necessary because tibbles can't have rownames
rownames(lake_isoz_per3) <- lake_per3$site_name
lake_isoz_per3$d18O <- scale(lake_isoz_per3$d18O, center=mean(lake_isoz_per3$d18O), scale=sd(lake_isoz_per3$d18O))
lake_isoz_per3$dxs <- scale(lake_isoz_per3$dxs, center=mean(lake_isoz_per3$dxs), scale=sd(lake_isoz_per3$dxs))
lake_isoz_per3_hclust <- hclust(dist(lake_isoz_per3))
lake_isoz_per3_dendro <- as.dendrogram(lake_isoz_per3_hclust)

# Lake parameter means and confidence interval per lake type
lake_per3_means_bylaketype <- lake_per3 %>%
  select(elevation, d18O, d2H, dxs, laketype, lake_surface_area, lakeshed_area, dist_gris, dist_ocean, doy) %>%
  group_by(laketype) %>%
  summarize(across(everything(), list(mean=mean, conf_int=~sd(.)/sqrt(length(.))*qnorm(0.975))))
print(lake_per3_means_bylaketype)

#==Lake isotope environmental parameter regression
# This section uses LASSO regression first to see which environmental parameters
#   are strongest predictors of isotopic variability in Pituffik lakes, then also
#   looks at linear regressions for more insight.
library(glmnet)
library(gridGraphics)
library(reshape2)

# Creating list of different headwater-downstream lake sample sets, varying based
#   on whether from sample period 3 only, from any of the 3 sample periods, or
#   from all samples, and also whether from the main lake set or full lake set.
base_headdown <- water_iso %>%
  filter(type == "lake" | type == "pond") %>%
  filter(laketype == "headwater" | laketype == "downstream")
base_headdown$log_area <- log(base_headdown$lake_surface_area) # Log transformation of lake surface area
base_headdown$log_lakeshed_area <- log(base_headdown$lakeshed_area) # Log transformation of lake watershed area
envparam_names <- c("elevation", "log_area", "log_lakeshed_area", "dist_gris", "dist_ocean", "doy") # Selected environmental parameters

headdown_list <- list()
headdown_list_names <- c("Main Period 3", "Main All Periods", "Main Alltime",
                         "Full Period 3", "Full All Periods", "Full Alltime")
headdown_list[[1]] <- base_headdown %>%
  filter(timing == 3 & main_lakes == 1)
headdown_list[[2]] <- base_headdown %>%
  filter(!is.na(timing) & main_lakes == 1)
headdown_list[[3]] <- base_headdown %>%
  filter(main_lakes == 1)
headdown_list[[4]] <- base_headdown %>%
  filter(timing == 3)
headdown_list[[5]] <- base_headdown %>%
  filter(!is.na(timing))
headdown_list[[6]] <- base_headdown
names(headdown_list) <- headdown_list_names

# Performing analyses where output gives results of LASSO regression, linear regressions
#   by individual parameters, and multiple regression.
lasso_headdown_fit <- list()
lasso_headdown_cvfit <- list()
lasso_headdown_cvfit_coef <- list()
lm_headdown_r2 <- list()
lm_headdown_pval <- list()
multlm_formula <- paste("iso_iter~",paste(envparam_names,collapse="+"), sep="") # Formula for multiple regressions
multlm_headdown <- list()
multlm_headdown_beta <- list()
multlm_headdown_pval <- list()
for(i in 1:length(headdown_list)) {
  headdown_envparam_iter <- headdown_list[[i]] %>% # Selecting environmental parameters for regression
    dplyr::select(all_of(envparam_names))
  lasso_headdown_fit_iter <- list()
  lasso_headdown_cvfit_iter <- list()
  lasso_headdown_cvfit_coef_iter <- list()
  lm_headdown_r2_iter <- as_tibble(matrix(ncol=length(vbl_clm), nrow=length(envparam_names), dimnames=list(NULL,vbl_clm)))
  lm_headdown_pval_iter <- as_tibble(matrix(ncol=length(vbl_clm), nrow=length(envparam_names), dimnames=list(NULL,vbl_clm)))
  multlm_headdown_iter <- list()
  multlm_headdown_beta_iter <- as_tibble(matrix(ncol=length(vbl_clm), nrow=length(envparam_names), dimnames=list(NULL,vbl_clm)))
  multlm_headdown_pval_iter <- as_tibble(matrix(ncol=length(vbl_clm), nrow=length(envparam_names), dimnames=list(NULL,vbl_clm)))
  for(j in 1:length(vbl_clm)){
    iso_iter <- headdown_list[[i]][[vbl_clm[j]]] # Selecting isotopic variable for regression
    # LASSO
    lasso_headdown_fit_iter[[j]] <- glmnet(headdown_envparam_iter, iso_iter)
    lasso_headdown_cvfit_iter[[j]] <- cv.glmnet(makeX(headdown_envparam_iter), iso_iter)
    lasso_headdown_cvfit_coef_iter[[j]] <- coef(lasso_headdown_cvfit_iter[[j]], s="lambda.1se")
    
    # Linear regression (single)
    for(k in 1:length(envparam_names)) {
      lm_headdown_r2_iter[k,j] <- summary(lm(iso_iter~headdown_envparam_iter[[k]]))$r.squared
      lm_headdown_pval_iter[k,j] <- summary(lm(iso_iter~headdown_envparam_iter[[k]]))$coefficients[8]
    }
    
    # Multiple regression
    multlm_headdown_iter[[j]] <- lm(as.formula(multlm_formula), data=headdown_envparam_iter)
    multlm_headdown_beta_iter[j] <- summary(multlm_headdown_iter[[j]])$coefficients[2:(2+length(envparam_names)-1)]
    multlm_headdown_pval_iter[j] <- summary(multlm_headdown_iter[[j]])$coefficients[23:(23+length(envparam_names)-1)]
  }
  names(lasso_headdown_fit_iter) <- vbl_clm
  names(lasso_headdown_cvfit_iter) <- vbl_clm
  names(lasso_headdown_cvfit_coef_iter) <- vbl_clm
  lm_headdown_r2_iter <- signif(lm_headdown_r2_iter,2)
  lm_headdown_pval_iter <- signif(lm_headdown_pval_iter,1)
  lm_headdown_r2_iter <- lm_headdown_r2_iter %>%
    add_column(envparam = envparam_names, .before = 1)
  lm_headdown_pval_iter <- lm_headdown_pval_iter %>%
    add_column(envparam = envparam_names, .before = 1)
  names(multlm_headdown_iter) <- vbl_clm
  multlm_headdown_beta_iter <- signif(multlm_headdown_beta_iter,2)
  multlm_headdown_pval_iter <- signif(multlm_headdown_pval_iter,1)
  multlm_headdown_beta_iter <- multlm_headdown_beta_iter %>%
    add_column(envparam = envparam_names, .before = 1)
  multlm_headdown_pval_iter <- multlm_headdown_pval_iter %>%
    add_column(envparam = envparam_names, .before = 1)
  
  # Putting iso variables together in list group by lakeset
  lasso_headdown_fit[[i]] <- lasso_headdown_fit_iter
  lasso_headdown_cvfit[[i]] <- lasso_headdown_cvfit_iter
  lasso_headdown_cvfit_coef[[i]] <- lasso_headdown_cvfit_coef_iter
  lm_headdown_r2[[i]] <- lm_headdown_r2_iter
  lm_headdown_pval[[i]] <- lm_headdown_pval_iter
  multlm_headdown[[i]] <- multlm_headdown_iter
  multlm_headdown_beta[[i]] <- multlm_headdown_beta_iter
  multlm_headdown_pval[[i]] <- multlm_headdown_pval_iter
}
names(lasso_headdown_fit) <- headdown_list_names
names(lasso_headdown_cvfit) <- headdown_list_names
names(lasso_headdown_cvfit_coef) <- headdown_list_names
names(lm_headdown_r2) <- headdown_list_names
names(lm_headdown_pval) <- headdown_list_names
names(multlm_headdown) <- headdown_list_names
names(multlm_headdown_beta) <- headdown_list_names
names(multlm_headdown_pval) <- headdown_list_names

# Option to print all six lakeset outputs
#print(lasso_headdown_cvfit_coef) # LASSO cross-validated coefficients
#print(lm_headdown_r2) # LM r2
#print(lm_headdown_pval) # LM p-values
#print(multlm_headdown_beta) # MultLM beta coefficients
#print(multlm_headdown_pval) # MultLM p-values

# Pulling out LASSO and linear regression output from only period 3 main lake samples 
print(lasso_headdown_cvfit_coef[["Main Period 3"]]) # LASSO coefficients
print(tibble(envparam = lm_headdown_r2[["Main Period 3"]][["envparam"]], # Linear regression r2 and fstat p-values
             d18O_r2=lm_headdown_r2[["Main Period 3"]][["d18O"]], d18O_pval=lm_headdown_pval[["Main Period 3"]][["d18O"]],
             d2H_r2=lm_headdown_r2[["Main Period 3"]][["d2H"]], d2H_pval=lm_headdown_pval[["Main Period 3"]][["d2H"]],
             dxs_r2=lm_headdown_r2[["Main Period 3"]][["dxs"]], dxs_pval=lm_headdown_pval[["Main Period 3"]][["dxs"]]))
print(tibble(envparam = multlm_headdown_beta[["Main Period 3"]][["envparam"]], # Multiple regression beta coefficients and fstat p-values
             d18O_beta=multlm_headdown_beta[["Main Period 3"]][["d18O"]], d18O_pval=multlm_headdown_pval[["Main Period 3"]][["d18O"]],
             d2H_beta=multlm_headdown_beta[["Main Period 3"]][["d2H"]], d2H_pval=multlm_headdown_pval[["Main Period 3"]][["d2H"]],
             dxs_beta=multlm_headdown_beta[["Main Period 3"]][["dxs"]], dxs_pval=multlm_headdown_pval[["Main Period 3"]][["dxs"]]))

# Extracting multiple regression information from final environmental parameter set
multlm_final_envparam <- list(d18O=c("elevation", "log_area", "log_lakeshed_area"), # Important parameters based on LASSO & multiple regression
                              d2H=c("elevation", "log_lakeshed_area"),
                              dxs=c("elevation", "log_area", "log_lakeshed_area"))
multlm_final_headdown_coef <- list()
multlm_final_headdown_strength <- as_tibble(matrix(nrow=length(vbl_clm), ncol=3, dimnames=list(NULL,c("fstat", "fstat_pval", "adj_r2"))))
for (i in 1:length(vbl_clm)) {
  iso_iter <- headdown_list[["Main Period 3"]][[vbl_clm[i]]]
  multlm_formula_iter <- paste("iso_iter~",paste(multlm_final_envparam[[vbl_clm[i]]],collapse="+"), sep="") # formula for multiple regressions
  multlm_final_headdown_iter <- summary(lm(as.formula(multlm_formula_iter), data=headdown_list[["Main Period 3"]]))
  multlm_final_headdown_coef[[i]] <- multlm_final_headdown_iter$coefficients
  multlm_final_headdown_coef[[i]] <- bind_cols(signif(multlm_final_headdown_coef[[i]][,1:3],2), signif(multlm_final_headdown_coef[[i]][,4],1))
  names(multlm_final_headdown_coef[[i]]) <- c("coef_value", "st_err", "tval", "pval")
  multlm_final_headdown_coef[[i]] <- multlm_final_headdown_coef[[i]] %>%
    add_column(parameter = c("intercept", multlm_final_envparam[[i]]), .before = 1)
  multlm_final_headdown_strength[i,] <- as.list(c(multlm_final_headdown_iter$fstatistic[[1]],
                                                  pf(multlm_final_headdown_iter$fstatistic[[1]], 
                                                     multlm_final_headdown_iter$fstatistic[[2]],
                                                     multlm_final_headdown_iter$fstatistic[[3]],
                                                     lower.tail = FALSE),
                                                  multlm_final_headdown_iter$adj.r.squared))
}
names(multlm_final_headdown_coef) <- vbl_clm
multlm_final_headdown_strength <- multlm_final_headdown_strength %>%
  add_column(iso = vbl_clm, .before = 1)

print(multlm_final_headdown_coef)
print(multlm_final_headdown_strength)

#==Stream basin analysis
# Stream LEL by basin for main sampling sites
stream_spatial_index <- c("North River Mouth", "North River Shelter 2 South Fork", "North River Shelter 5.8",
                          "South River Mouth", "Fox Canyon Bridge North Fork", "Fox Canyon Bridge South Fork",
                          "Amitsuarsuk")
stream_subset_northsouth <- water_iso %>% # Extracting out stream data
  filter(type == "stream") %>%
  filter(site_name %in% stream_spatial_index)
amitsuarsuk_subset<- water_iso %>% # Extracting out stream data
  filter(type == "stream") %>%
  filter(basin_name == "Amitsuarsuk")
amitsuarsuk_subset$site_name <- "Amitsuarsuk" # Renaming all sites along stream to be simply Amitsuarsuk
stream_subset <- rbind(stream_subset_northsouth, amitsuarsuk_subset) # Putting all desired stream data in one tibble
stream_subset$site_name <- factor(stream_subset$site_name, levels=stream_spatial_index) # Making so the order is always correct

mean_iso_pituffik_sioraq <- stream_subset %>% # Calculating the mean isotopic values of all Pituffik and Sioraq River samples
  filter (basin_name != "Amitsuarsuk") %>%
  summarize(count=n(),mean_d18O=mean(d18O, na.rm=TRUE)*1000, conf_int_d18O=sd(d18O, na.rm=TRUE)/sqrt(sum(!is.na(d18O)))*qnorm(0.975)*1000,
            mean_d2H=mean(d2H, na.rm=TRUE)*1000, conf_int_d2H=sd(d2H, na.rm=TRUE)/sqrt(sum(!is.na(d2H)))*qnorm(0.975)*1000,
            mean_dxs=mean(dxs, na.rm=TRUE)*1000, conf_int_dxs=sd(dxs, na.rm=TRUE)/sqrt(sum(!is.na(dxs)))*qnorm(0.975)*1000)
print(mean_iso_pituffik_sioraq)

lel_stream_bystream <- stream_subset %>% # Running regression by type
  group_by(site_name) %>%
  nest() %>%
  mutate(lel = map(data, ~lm(d2H~d18O, data=.)))

lel_slope <- 
  lel_stream_bystream %>%
  mutate(tidyit = map(lel, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "d18O") %>%
  dplyr::select(site_name, estimate, std.error)
colnames(lel_slope) <- c("site_name", "slope", "se_slope")
lel_slope$conf_int_slope <- lel_slope$se_slope*qnorm(0.975)

lel_intercept <- 
  lel_stream_bystream %>%
  mutate(tidyit = map(lel, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(site_name, estimate, std.error)
colnames(lel_intercept) <- c("site_name", "intercept", "se_intercept")
lel_intercept$conf_int_intercept <- lel_intercept$se_intercept*qnorm(0.975)

lel_r2 <- 
  lel_stream_bystream %>%
  mutate(glanceit = map(lel, broom::glance)) %>%
  unnest(glanceit) %>%
  dplyr::select(site_name, r.squared, p.value, nobs)

lel_stream_bystream_params <- lel_slope %>% # Combining all data
  inner_join(lel_intercept, by="site_name") %>%
  inner_join(lel_r2, by="site_name")

lel_stream_bystream_params$basin_name <- distinct(stream_subset, site_name, .keep_all = TRUE) %>%
  ungroup() %>%
  dplyr::pull(basin_name)
print(lel_stream_bystream_params)

lel_stream_bystream_bybasin_slope <- lel_stream_bystream_params %>%
  group_by(basin_name) %>%
  summarize(count=n(),mean_slope=mean(slope, na.rm=TRUE), conf_int_slope=sd(slope, na.rm=TRUE)/sqrt(sum(!is.na(slope)))*qnorm(0.975))
print(lel_stream_bystream_bybasin_slope)

####========END SPATIAL ANALYSIS======####


####==================================####
####=======Temporal Analyses==========####
####==================================####
#==Lake temporal analyses
# Local water line parameters and stat values for 2018 Lake Potato/Power Lake
lwl_lake18_bytype <- water_iso_plus_gnip %>% # Running regression by type
  filter(yr == 2018) %>%
  filter(site_name == "Lake Potato" | site_name == "Power Lake") %>%
  group_by(site_name) %>%
  nest() %>%
  mutate(lwl = map(data, ~lm(d2H~d18O, data=.)))

lwl_lake18_slope <- 
  lwl_lake18_bytype %>%
  mutate(tidyit = map(lwl, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "d18O") %>%
  dplyr::select(site_name, estimate, std.error)
colnames(lwl_lake18_slope) <- c("site_name", "slope", "se_slope")
lwl_lake18_slope$conf_int_slope <- lwl_lake18_slope$se_slope*qnorm(0.975)

lwl_lake18_intercept <- 
  lwl_lake18_bytype %>%
  mutate(tidyit = map(lwl, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(site_name, estimate, std.error)
colnames(lwl_lake18_intercept) <- c("site_name", "intercept", "se_intercept")
lwl_lake18_intercept$conf_int_intercept <- lwl_lake18_intercept$se_intercept*qnorm(0.975)

lwl_lake18_r2 <- 
  lwl_lake18_bytype %>%
  mutate(glanceit = map(lwl, broom::glance)) %>%
  unnest(glanceit) %>%
  dplyr::select(site_name, r.squared, p.value, nobs)

lwl_lake18 <- lwl_lake18_slope %>% # Combining all data
  inner_join(lwl_lake18_intercept, by="site_name") %>%
  inner_join(lwl_lake18_r2, by="site_name")
print(lwl_lake18)
#write.csv(lwl,"lwl_lake18_regression.csv", row.names=FALSE) # Remove comment marker to output file

# Finding intersections between LEL and LMWL for 2018 Lake Potato/Power Lake
lake18_index <- c("Lake Potato", "Power Lake")
lwl_lake18_intersect <- as_tibble(matrix(nrow=length(lake18_index),ncol=2), .name_repair = "unique")
lwl_lake18_intersect <- as_tibble(cbind(lake18_index, lwl_lake18_intersect))
colnames(lwl_lake18_intersect) <- c("lake", "d18O", "d2H")

for (i in 1:length(lake18_index)) {
  coef_iter <- rbind(coef(lwl_bytype$lwl[lwl_bytype$type=="GNIP"][[1]]), # Coefficient matrix
                     coef(lwl_lake18_bytype$lwl[lwl_lake18_bytype$site_name==lake18_index[i]][[1]]))
  lwl_lake18_intersect[i,c(2,3)] <- t(c(-solve(cbind(coef_iter[,2],-1)) %*% coef_iter[,1]))
}
print(lwl_lake18_intersect)

# Comparing how the isotopic values change between different time periods
lm_iso_doy <- list()
lm_iso_doy_means <- list()
isodiff <- list()
water_iso_periodcompare_bytiming <- water_iso_periodcompare %>% # Splitting into a list by timing period
  group_split(timing)
period_index <- c("x2018a_2018b", "x2018a_2019", "x2018b_2019")
timing_avoid_index <- c(3, 2, 1)
diff_index_a <- c(2,3,3)
diff_index_b <- c(1,1,2)
# Linear regressions to be used in loop
iso_doy_func <- list()
iso_doy_func[[1]] <- function(data) {
  lm(d18O~doy, data=data)
}
iso_doy_func[[2]] <- function(data) {
  lm(d2H~doy, data=data)
}
iso_doy_func[[3]] <- function(data) {
  lm(dxs~doy, data=data)
}
# Loop to calculate slopes between day of year and isotopic values between different periods
for (j in 1:length(period_index)) { # Looping between time periods comparison pairs
  lm_iso_doy_iter <- list()
  lm_iso_doy_means_iter <- as_tibble(data.frame(matrix(nrow=length(vbl_clm), ncol=3)))
  isodiff_iter <- as_tibble(data.frame(matrix(nrow=length(vbl_clm), ncol=3)))
  for (i in 1:length(vbl_clm)){ # Looping through isotopic variables
    lm_iso_doy_iter[[i]] <- water_iso_periodcompare %>%
      filter(timing != timing_avoid_index[j]) %>%
      group_by(site_name) %>%
      nest() %>%
      mutate(iso_lm = lapply(data, iso_doy_func[[i]])) %>%
      mutate(tidyit = map(iso_lm, broom::tidy)) %>%
      unnest(tidyit) 
    lm_iso_doy_means_iter[i,1] <- vbl_clm[i] # Getting the mean slope and intercept per period compare
    lm_iso_doy_means_iter[i,2] <- lm_iso_doy_iter[[i]] %>%
      ungroup() %>%
      filter(term == "doy") %>%
      summarize(mean=mean(estimate))
    lm_iso_doy_means_iter[i,3] <- lm_iso_doy_iter[[i]] %>%
      ungroup() %>%
      filter(term == "(Intercept)") %>%
      summarize(mean=mean(estimate))
    isodiff_iter[i,1] <- vbl_clm[i] # Getting the difference in mean isotopic value pairwise per period compare
    isodiff_iter[i,2] <- mean(water_iso_periodcompare_bytiming[[diff_index_a[j]]][[vbl_clm_index[i]]] - 
                                water_iso_periodcompare_bytiming[[diff_index_b[j]]][[vbl_clm_index[i]]])
    isodiff_iter[i,3] <- sd(water_iso_periodcompare_bytiming[[diff_index_a[j]]][[vbl_clm_index[i]]] - 
                              water_iso_periodcompare_bytiming[[diff_index_b[j]]][[vbl_clm_index[i]]])
  }
  names(lm_iso_doy_iter) <- vbl_clm
  colnames(lm_iso_doy_means_iter) <- c("iso", "slope", "yint")
  colnames(isodiff_iter) <- c("iso", "mean", "sd")
  lm_iso_doy[[j]] <- lm_iso_doy_iter
  lm_iso_doy_means[[j]] <- lm_iso_doy_means_iter
  isodiff[[j]] <- isodiff_iter
}
names(lm_iso_doy) <- period_index
names(lm_iso_doy_means) <- period_index
names(isodiff) <- period_index

# Predicting what the isotopic values in 2018 would be based on the day of year of 2019 sampling
predict_iso_2018 <- as_tibble(matrix(nrow=nrow(water_iso_periodcompare_bytiming[[3]]), ncol=2*length(vbl_clm)+1))
predict_iso_2018[1] <- water_iso_periodcompare_bytiming[[3]]$site_name
for (i in 1:length(vbl_clm)){ # Looping through isotopic variables
  predict_iso_2018[i+1] <- water_iso_periodcompare_bytiming[[3]]$doy *
    lm_iso_doy[[1]][[i]]$estimate[lm_iso_doy[[1]][[i]]$term == "doy"] +
    lm_iso_doy[[1]][[i]]$estimate[lm_iso_doy[[1]][[i]]$term == "(Intercept)"]
  predict_iso_2018[i+4] <- water_iso_periodcompare_bytiming[[3]][[vbl_clm[i]]]
}
colnames(predict_iso_2018) <- c("site_name", "d18O_predicted", "d2H_predicted", "dxs_predicted",
                                "d18O_2019", "d2H_2019", "dxs_2019")

t.test(predict_iso_2018$d18O_predicted, predict_iso_2018$d18O_2019)
t.test(predict_iso_2018$d2H_predicted, predict_iso_2018$d2H_2019)
t.test(predict_iso_2018$dxs_predicted, predict_iso_2018$dxs_2019)

#### END TEMPORAL ANALYSES ####



#=====================================================================#
####==================Plotting Figures=============================####
#=====================================================================#

####==================================####
####====Climatology plots====####
####==================================####
# This script sets out the Pituffik climatology figures
wx2018_mo <- pfk_isowx_day %>% # Calculating monthly means/sums of 2018 wx data
  filter(yr == 2018) %>%
  group_by(mo) %>%
  summarize(tavg = mean(tavg, na.rm=TRUE), prcp = mean(prcp_H2O, na.rm=TRUE),
            d18O = mean(d18O, na.rm=TRUE), d2H = mean(d2H, na.rm=TRUE), 
            dxs = mean(dxs, na.rm=TRUE), pet = mean(pet, na.rm=TRUE)) %>%
  left_join(pfk_mo_clim%>%select(mo, doy), by="mo")

wx2019_mo <- pfk_isowx_day %>% # Calculating monthly means/sums of 2019 wx data
  filter(yr == 2019) %>%
  group_by(mo) %>%
  summarize(tavg = mean(tavg, na.rm=TRUE), prcp = mean(prcp_H2O, na.rm=TRUE),
            d18O = mean(d18O, na.rm=TRUE), d2H = mean(d2H, na.rm=TRUE), 
            dxs = mean(dxs, na.rm=TRUE), pet = mean(pet, na.rm=TRUE)) %>%
  left_join(pfk_mo_clim%>%select(mo, doy), by="mo")

pfk_gnip_mo <- pfk_gnip %>%  # Calculating monthly means of GNIP data
  mutate(mo = month(date_start)) %>%
  group_by(mo) %>%
  summarize(d18O = mean(d18O, na.rm=TRUE), d2H = mean(d2H, na.rm=TRUE), 
            dxs = mean(dxs, na.rm=TRUE)) %>%
  left_join(pfk_mo_clim %>% select(mo, doy), by="mo")

wx2018_wk <- pfk_isowx_day %>% # Calculating weekly means/sums of 2018 wx data
  filter(yr == 2018) %>%
  mutate(wk = week(date)) %>%
  group_by(wk) %>%
  summarize(tavg = mean(tavg, na.rm=TRUE), prcp = mean(prcp_H2O, na.rm=TRUE),
            d18O = mean(d18O, na.rm=TRUE), d2H = mean(d2H, na.rm=TRUE), 
            dxs = mean(dxs, na.rm=TRUE), pet = mean(pet, na.rm=TRUE)) %>%
  mutate(doy = wk*7-4)

wx2019_wk <- pfk_isowx_day %>% # Calculating weekly means/sums of 2019 wx data
  filter(yr == 2019) %>%
  mutate(wk = week(date)) %>%
  group_by(wk) %>%
  summarize(tavg = mean(tavg, na.rm=TRUE), prcp = mean(prcp_H2O, na.rm=TRUE),
            d18O = mean(d18O, na.rm=TRUE), d2H = mean(d2H, na.rm=TRUE), 
            dxs = mean(dxs, na.rm=TRUE), pet = mean(pet, na.rm=TRUE)) %>%
  mutate(doy = wk*7-4)

# Weekly + Monthly temperature plot
tavg_clim_plot <- ggplot() +
  theme_classic() +
  geom_hline(yintercept=0, color="gray80") +
  geom_line(aes(x=wx2018_wk$doy, # 2018 SMtn wx obs
                y=wx2018_wk$tavg), color="dodgerblue", linetype="solid", linewidth=0.5) +
  geom_line(aes(x=wx2019_wk$doy, # 2019 SMtn wx obs
                y=wx2019_wk$tavg), color="deeppink2", linetype="solid", linewidth=0.5) +
  geom_line(aes(x=pfk_doy_clim_sinusoid$doy, # Climatology data from USAF data
                  y=pfk_doy_clim_sinusoid$tavg_clim), color="gray30", linewidth=1) +
  scale_y_continuous(name="Air Temperature (°C)", position = "left", limits = c(-30,15)) +
  scale_x_continuous(name=NULL, position = "bottom", breaks=pfk_mo_clim$doy,
                     labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))


# Monthly precipitation plot (Weekly data is hard to view and interpret)
prcp_clim_plot <- ggplot() +
  theme_classic() +
  geom_hline(yintercept=0, color="gray80") +
  geom_col(aes(x=wx2018_mo$mo, # 2018 USAF THU obs
               y=wx2018_mo$prcp), fill="dodgerblue") +
  geom_col(aes(x=wx2019_mo$mo, # 2019 USAF THU 
               y=wx2019_mo$prcp), fill="deeppink2") + # MISSING AUG PRECIP IN WEEKLY
  geom_step(aes(x=pfk_mo_clim$mo,  # USAF THU 2000-2021 means
                y=pfk_mo_clim$prcp), color="gray30", linewidth=1) +
  scale_y_continuous(name="Precipitation (mm d-1)", position = "left", limits = c(0,1)) +
  scale_x_continuous(name=NULL, position = "bottom", breaks=seq(1,12),
                     labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

# Weekly + Monthly PET plot
pet_clim_plot <- ggplot() +
  theme_classic() +
  geom_hline(yintercept=0, color="gray80") +
  geom_line(aes(x=wx2018_wk$doy, # Weekly mean rate in 2018
                y=wx2018_wk$pet), color="dodgerblue", linetype="solid", linewidth=0.5) +
  geom_line(aes(x=wx2019_wk$doy, # Weekly mean rate in 2019
                y=wx2019_wk$pet), color="deeppink2", linetype="solid", linewidth=0.5) +
  geom_line(aes(x=pfk_doy_clim_sinusoid$doy, # 1981-2023 daily mean
                  y=pfk_doy_clim_sinusoid$pet_clim), color="gray30", linewidth=1) +
  scale_y_continuous(name="PET (mm d-1)", position = "left", limits = c(0,3)) +
  scale_x_continuous(name=NULL, position = "bottom", breaks=pfk_mo_clim$doy,
                     labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

# Monthly precipitation and water vapor d18O plot
d18O_clim_plot <- ggplot() +
  theme_classic() +
  geom_step(aes(x=wx2018_mo$mo, # Water vapor isotopes 2018
                y=wx2018_mo$d18O*1000), color="dodgerblue", linetype="dashed", linewidth=1) +
  geom_step(aes(x=wx2019_mo$mo, # Water vapor isotopes 2019
                y=wx2019_mo$d18O*1000), color="deeppink2", linetype="dashed", linewidth=1) +
  geom_step(aes(x=pfk_gnip_mo$mo, # GNIP precipitation isotopes
                y=pfk_gnip_mo$d18O*1000), color="gray30", linewidth=1) +
  scale_y_continuous(name="d18O (‰)", position = "left") +
  scale_x_continuous(name=NULL, position = "bottom", breaks=seq(1,12),
                     labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

# Monthly precipitation and water vapor d2H plot (not used in composite)
d2H_clim_plot <- ggplot() +
  theme_classic() +
  geom_step(aes(x=wx2018_mo$mo, # Water vapor isotopes 2018
                y=wx2018_mo$d2H*1000), color="dodgerblue", linetype="dashed", linewidth=1) +
  geom_step(aes(x=wx2019_mo$mo, # Water vapor isotopes 2019
                y=wx2019_mo$d2H*1000), color="deeppink2", linetype="dashed", linewidth=1) +
  geom_step(aes(x=pfk_gnip_mo$mo, # GNIP precipitation
                y=pfk_gnip_mo$d2H*1000), color="gray30", linewidth=1) +
  scale_y_continuous(name="d2H (‰)", position = "left") +
  scale_x_continuous(name=NULL, position = "bottom", breaks=seq(1,12),
                     labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

# Monthly precipitation and water vapor dxs plot
dxs_clim_plot <- ggplot() +
  theme_classic() +
  geom_step(aes(x=wx2018_mo$mo, # Water vapor isotopes 2018
                y=wx2018_mo$dxs*1000), color="dodgerblue", linetype="dashed", linewidth=1) +
  geom_step(aes(x=wx2019_mo$mo, # Water vapor isotopes 2019
                y=wx2019_mo$dxs*1000), color="deeppink2",linetype="dashed", linewidth=1) +
  geom_step(aes(x=pfk_gnip_mo$mo, # GNIP precipitation isotopes
                y=pfk_gnip_mo$dxs*1000), color="gray30", linewidth=1) +
  scale_y_continuous(name="dxs (‰)", position = "left") +
  scale_x_continuous(name=NULL, position = "bottom", breaks=seq(1,12),
                     labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

# Compiling all plots into composite plot
ptfk_clim_plot <- plot_grid(plotlist=list(tavg_clim_plot, prcp_clim_plot, pet_clim_plot,
                                          d18O_clim_plot, d2H_clim_plot, dxs_clim_plot),
                            ncol=3, align = "hv")
windows(height=14, width=21)
#pdf("Figures/Pituffik_climatology.pdf", height=14, width=21)
ptfk_clim_plot
#ggsave("Figures/Pituffik_climatology.png", height=14, width=21, dpi=600)
#dev.off()

#### END CLIMATOLOGY PLOTS ####


####==================================####
####====Overall data results plots====####
####==================================####
# This makes a plot of the general structure of data sets by the different sample origins.
overall_mean_iso <- water_iso %>%
  summarize(across(vbl_clm, mean))
overall_sd_iso <- water_iso %>%
  summarize(across(vbl_clm, sd))
iso_violin_plot <- list()
colorset <- c("deepskyblue3", "darkgoldenrod2", "violetred4", "springgreen4",
              "cyan3", "slateblue4", "lightpink2", "gray70")
type_labels <- c("lake", "pool", "stream", "surface flow", "snow/ice",
                 "prcp rain", "prcp snow", "GNIP 1966-1971")
for (i in 1:length(vbl_clm)) 
  local({ #Have to do local or else plot later only does each plot in most recent i
    i <- i
    iso_violin_plot[[i]] <<- ggplot(data=water_iso_plus_gnip) +
      theme_classic() +
      geom_violin(aes(x=type, y=.data[[vbl_clm[i]]]*1000, color=type, fill=type), scale="width", linetype='dotted', draw_quantiles = c(0.25, 0.75)) +
      geom_violin(aes(x=type, y=.data[[vbl_clm[i]]]*1000, color=type), scale="width", fill='transparent', draw_quantiles = c(0.5)) +
      geom_hline(yintercept=overall_mean_iso[[i]]*1000, linetype="dashed", color="gray60") +
      geom_hline(yintercept=overall_mean_iso[[i]]*1000+overall_sd_iso[[i]]*1000, linetype="dotted", color="gray60") +
      geom_hline(yintercept=overall_mean_iso[[i]]*1000-overall_sd_iso[[i]]*1000, linetype="dotted", color="gray60") +
      scale_color_manual(values=colorset) +
      scale_fill_manual(values=alpha(colorset, 0.4)) +
      scale_y_continuous(name=paste(vbl_clm[i],"(‰)"), position = "left") +
      scale_x_discrete(labels=toupper(type_labels), position = "bottom", name=NULL) +
      theme(legend.position = "none",
            axis.title.x = element_text(size=16, color='gray40'),
            axis.text.x =  element_text(size=10, color=colorset, angle=30, hjust=1),
            axis.ticks.x = element_line(color='gray40'),
            axis.line.x = element_line(color='gray40'),
            axis.title.y = element_text(size=16, color='gray40'),
            axis.text.y = element_text(size=14, color='gray40'),
            axis.ticks.y = element_line(color='gray40'),
            axis.line.y = element_line(color='gray40'))
  })

windows(height=8, width=4)
#pdf("Figures/Iso_Violin_Plot_ByGroup.pdf", height=8, width=4)
plot_grid(plotlist = iso_violin_plot, ncol=1, align = "hv")
#ggsave("Figures/Iso_Violin_Plot_ByGroup.png", height=8, width=4, dpi=600)
#dev.off()

# This makes a plot comparing the LWL of samples by groups
water_iso_lwl <- subset(water_iso, type != "precipitation rain" & type != "precipitation snow")
water_iso_lwl$type <- factor(water_iso_lwl$type, levels=c("lake", "pool", "stream", # Making so the order is always correct
                                                            "surface flow", "snow or ice"))
colorset_inv <- rev(colorset[seq(1:length(unique(water_iso_lwl$type)))])

# Full LMWL plot
iso_lwl_plot <- ggplot() +
  theme_classic() +
  geom_abline(aes(slope=8, intercept=10), linetype="solid", lwd=1, color="gray40") + # GMWL
  geom_abline(aes(slope=7.46, intercept=-3.30), linetype="dashed", lwd=1, color="gray40") + # GNIP LMWL
  geom_abline(aes(slope=6.95, intercept=-18.24), linetype="dotted", lwd=1, color="gray40") + # Water Vapor Akers 2020 LMWL
  geom_point(aes(x=d18O*1000, y=d2H*1000, color=fct_rev(type)), data=water_iso_lwl, alpha=0.3) +
  geom_smooth(aes(x=d18O*1000, y=d2H*1000, color=fct_rev(type), fill=fct_rev(type)), data=water_iso_lwl, method="lm") +
  geom_rect(aes(xmin=-23, xmax=-15, ymin=-170, ymax=-125), color="gray60", fill=alpha("gray", 0))+
  scale_color_manual(values=colorset_inv, guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=alpha(colorset_inv, 0.4), guide=guide_legend(reverse=TRUE)) +
  scale_y_continuous(name="d2H (‰)") +
  scale_x_continuous(name="d18O (‰)") +
  coord_cartesian(xlim=c(-31,-4), ylim=c(-230,-80)) +
  theme(legend.position = c(0.95,0.05),
        legend.justification = c("right", "bottom"),
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40'),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

# Zoomed in LMWL plot
iso_lwl_plot_zoom <- ggplot() +
  theme_classic() +
  geom_abline(aes(slope=8, intercept=10), linetype="solid", lwd=1, color="gray40") + # GMWL
  geom_abline(aes(slope=7.46, intercept=-3.30), linetype="dashed", lwd=1, color="gray40") + # GNIP LMWL
  geom_abline(aes(slope=6.95, intercept=-18.24), linetype="dotted", lwd=1, color="gray40") + # Water Vapor Akers 2020 LMWL
  geom_point(aes(x=d18O*1000, y=d2H*1000, color=fct_rev(type)), data=water_iso_lwl, alpha=0.3) +
  geom_smooth(aes(x=d18O*1000, y=d2H*1000, color=fct_rev(type), fill=fct_rev(type)), data=water_iso_lwl, method="lm") +
  scale_color_manual(values=colorset_inv, guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=alpha(colorset_inv, 0.4), guide=guide_legend(reverse=TRUE)) +
  scale_y_continuous(name="d2H (‰)") +
  scale_x_continuous(name="d18O (‰)") +
  coord_cartesian(xlim=c(-23,-15), ylim=c(-170,-125)) +
  theme(legend.position = c(0.95,0.05),
        legend.justification = c("right", "bottom"),
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40'),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

windows(height=6, width=12)
#pdf("Figures/pfk_LMWL_plots.pdf", height=6, width=12)
plot_grid(plotlist = list(iso_lwl_plot, iso_lwl_plot_zoom), ncol=2, align = "hv")
#ggsave("Figures/pfk_LMWL_plots.png", height=6, width=12, dpi = 600)
#dev.off()

#### END OVERALL PLOTS ####

####==================================####
####======Spatial analyses plots======####
####==================================####
#====Violin and LEL plots of lakes by the different lake types====#
laketype_labels <- setNames(c("megapool", "headwater", "downstream", "vale", "proglacial", "altered"), levels(water_iso$laketype))
colorset <- setNames(c("violetred", "turquoise4", "yellow3", "orange2", "mediumpurple1", "gray70"), levels(water_iso$laketype))
period_label <- c("E. Summer 2018: All", "L. Summer 2018: All", "M. Summer 2019: All",
                  "E. Summer 2018: Main", "L. Summer 2018: Main", "M. Summer 2019: Main") #Seems to be a flaw/bug in getting labeling to work with violin plots

# Calculating overall y limits per variable for group plotting
ylims <- list()
for (i in 1:length(vbl_clm)) {
  ylims[[i]] <- range(lakeset[vbl_clm[i]])
}
names(ylims) <- vbl_clm

# Building list of violin plots of lake iso by laketype (Period 1 vs 2 vs 3)
lake_iso_violin_plot <- list()
lake_iso_violin_plot_iter <- list()
for (j in 1:nrow(lakeset_byperiod)) {
  data_iter <- lakeset_byperiod[[2]][[j]]
  for (i in 1:length(vbl_clm)) 
    local({ #Have to do local or else plot later only does each plot in most recent i
      i <- i
      lake_iso_violin_plot_iter[[i]] <<- ggplot(data=data_iter) +
        theme_classic() +
        geom_violin(aes(x=laketype, y=.data[[vbl_clm[i]]]*1000, color=laketype, fill=laketype), scale="width", linetype='dotted', draw_quantiles = c(0.25, 0.75)) +
        geom_violin(aes(x=laketype, y=.data[[vbl_clm[i]]]*1000, color=laketype), scale="width", fill='transparent', draw_quantiles = c(0.5)) +
        scale_color_manual(values=colorset) +
        scale_fill_manual(values=alpha(colorset, 0.4)) +
        scale_y_continuous(name=paste(vbl_clm[i],"(‰)"), position = "left", limits=ylims[[vbl_clm[i]]]*1000) +
        scale_x_discrete(labels=toupper(laketype_labels), position = "bottom", name=NULL) +
        theme(legend.position = "none",
              axis.title.x = element_text(size=16, color='gray40'),
              axis.text.x =  element_text(size=10, color=colorset[unique(data_iter$laketype)[order(unique(data_iter$laketype))]], angle=30, hjust=1),
              axis.ticks.x = element_line(color='gray40'),
              axis.line.x = element_line(color='gray40'),
              axis.title.y = element_text(size=16, color='gray40'),
              axis.text.y = element_text(size=14, color='gray40'),
              axis.ticks.y = element_line(color='gray40'),
              axis.line.y = element_line(color='gray40'))
    })
  lake_iso_violin_plot[[j]] <- plot_grid(plotlist=lake_iso_violin_plot_iter, ncol=1,
                                         align = "hv")
}

# Violin plots of isotopes by laketype, all 3 periods, all + main lakes, used for broader context
windows(height=10, width=16)
#pdf("Figures/Lake_Iso_Violin_Plot_ByType.pdf", height=10, width=16)
plot_grid(plotlist = lake_iso_violin_plot, ncol=nrow(lakeset_byperiod), align = "hv",
          labels=c(period_label,period_label), label_fontface="plain", label_size = 12, label_y = 1.03) +
  theme(plot.margin=unit(c(1.2,0,0,0), "cm"))
#ggsave("Figures/Lake_Iso_Violin_Plot_ByType.png", height=10, width=16, dpi=600)
#dev.off()

# Violin plots of isotopes by laketype, only period 3 (midsummer 2019), all lakes
lake_per3_iso_violin_plot <- lake_iso_violin_plot[[3]]

# A lot of Period 3 lakes all in one violin, used to add to violin plot in final version made in Illustrator
lake_total_iso_violin_plot_iter <- list()
for (i in 1:length(vbl_clm)) 
  local({ #Have to do local or else plot later only does each plot in most recent i
    i <- i
    lake_total_iso_violin_plot_iter[[i]] <<- ggplot(data=lakeset %>% filter(timing == 3)) +
      theme_classic() +
      geom_violin(aes(x=1,y=.data[[vbl_clm[i]]]*1000), scale="width", color="gray30", fill="gray50", linetype='dotted', draw_quantiles = c(0.25, 0.75)) +
      geom_violin(aes(x=1,y=.data[[vbl_clm[i]]]*1000), scale="width", color="gray30", fill='transparent', draw_quantiles = c(0.5)) +
      scale_y_continuous(name=paste(vbl_clm[i],"(‰)"), position = "left") +
      scale_x_discrete(labels=toupper(laketype_labels), position = "bottom", name=NULL) +
      theme(legend.position = "none",
            axis.title.x = element_text(size=16, color='gray40'),
            axis.text.x =  element_text(size=10, color='gray40', angle=30, hjust=1),
            axis.ticks.x = element_line(color='gray40'),
            axis.line.x = element_line(color='gray40'),
            axis.title.y = element_text(size=16, color='gray40'),
            axis.text.y = element_text(size=14, color='gray40'),
            axis.ticks.y = element_line(color='gray40'),
            axis.line.y = element_line(color='gray40'))
  })

lake_total_iso_violin_plot <- plot_grid(plotlist = lake_total_iso_violin_plot_iter, ncol=1, align = "hv")


# Plotting E/I values per lake type for 3+ sampled lakes
lake_e_i_fulldata_plot <- lake_e_i_fulldata %>%  # Adjusting values out of bounds to fit bounds
  mutate(E_I_d18O_bayes = replace(E_I_d18O_bayes, E_I_d18O_bayes > 1, 1), # Values > 1 represent unstable lakes with high evap
         E_I_d18O_bayes = replace(E_I_d18O_bayes, E_I_d18O_bayes < 0, 0)) # Values < 1 come from a lake with direct rapid glacial melt supply and no evap

data_iter <- lake_e_i_fulldata_plot %>% filter(timing == 3) # Summer 2019 E/I values (note, 1 altered lake later excluded due to poor LEL input data)
lake_e_i_2019_plot <- ggplot(data=data_iter) +
  theme_classic() +
  geom_jitter(aes(x=laketype, y=E_I_d18O_bayes, color=laketype), position=position_jitter(width=0), size=6, alpha=0.4, shape=19) +
  scale_color_manual(values=colorset) +
  scale_y_continuous(name="E/I", position = "left") +
  scale_x_discrete(labels=toupper(laketype_labels), position = "bottom", name=NULL) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=10, color=colorset[unique(data_iter$laketype)[order(unique(data_iter$laketype))]], angle=30, hjust=1),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))


# Plotting LELs for lakes with more than 3 samples
colorset <- setNames(c("violetred", "turquoise4", "yellow3", "orange2", "mediumpurple1", "gray70"), levels(water_iso$laketype))
lake_lel_bylake_plot <- ggplot() +
  theme_classic() +
  geom_abline(aes(slope=8, intercept=10), linetype="solid", lwd=1, color="gray40") + # GMWL
  geom_abline(aes(slope=7.46, intercept=-3.30), linetype="dashed", lwd=1, color="gray40") + # GNIP LMWL
  geom_abline(aes(slope=6.95, intercept=-18.24), linetype="dotted", lwd=1, color="gray40") + # Water Vapor Akers 2020 LMWL
  geom_point(aes(x=d18O*1000, y=d2H*1000, group=site_name, color=laketype), data=multilake, alpha=0.3) +
  geom_smooth(aes(x=d18O*1000, y=d2H*1000, group=site_name, color=laketype), fill=NA, data=multilake, method="lm") +
  scale_color_manual(values=colorset) +
  scale_y_continuous(name="d2H (‰)") +
  scale_x_continuous(name="d18O (‰)") +
  coord_cartesian(xlim=c(-23,-4), ylim=c(-175,-80)) +
  theme(legend.position = c(0.95,0.05),
        legend.justification = c("right", "bottom"),
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40'),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

windows(height=8, width=16)
#pdf("Figures/lake_violin_lel_bylake_plots.pdf", height=8, width=16)
plot_grid(plotlist = list(lake_total_iso_violin_plot, lake_per3_iso_violin_plot,
                          lake_e_i_2019_plot, lake_lel_bylake_plot),
          ncol=4, align = "hv",
          rel_widths = c(0.3,0.4,0.3,1))   
#ggsave("Figures/lake_violin_lel_bylake_plots.png", height=8, width=16)
#dev.off()

#=========Dendrogram plots=========#
# Plotting clustering results in dendrogram
colorbytype <- colorset[lake_per3$laketype_number]
labels_colors(lake_isoz_per3_dendro) <- colorbytype[order.dendrogram(lake_isoz_per3_dendro)]
windows(height=8, width=12)
#pdf("Figures/ptfk_lake_dendro1.pdf", height=8, width=12)
#png("Figures/ptfk_lake_dendro1.png", height=800, width=1200)
par(cex=0.7, mar=c(2,2,2,10))
plot(lake_isoz_per3_dendro, horiz=TRUE)
#dev.off()

# Plotting same results as unrooted tree
lake_isoz_per3_hclust_unroot <- lake_isoz_per3_hclust
labels(lake_isoz_per3_hclust_unroot) <- rep("•", nrow(lake_isoz_per3))
labels(lake_isoz_per3_hclust_unroot) <- rep("o", nrow(lake_isoz_per3))
labels(lake_isoz_per3_hclust_unroot) <- rep("–", nrow(lake_isoz_per3))
windows(height=10, width=10)
#pdf("Figures/ptfk_lake_dendro2.pdf", height=8, width=8)
#png("Figures/ptfk_lake_dendro2.png", height=800, width=800)
plot(as.phylo(lake_isoz_per3_hclust_unroot), type = "unrooted", cex = 1.5,
     no.margin = TRUE, show.tip.label = FALSE, rotate.tree=78)
tiplabels(pch=20, col=colorbytype)
#dev.off()

### Regression plots for Main Lakes, Sample Period 3 ###
lm_headdown_plots <- list()
multlm_headdown_addvar_iter <- list()
multlm_headdown_addvar <- list()

# Function for making added variable plots
#   From stackoverflow user Ben (3 Dec 2019), with some updating by PDA
avPlots.invis <- function(MODEL, ...) {
  
  ff <- tempfile()
  png(filename = ff)
  OUT <- car::avPlots(MODEL, ...)
  dev.off()
  unlink(ff)
  OUT }

gg_avplots  <- function(model_input, ylabel = NULL) {
  
  #Extract the information for AV plots
  avplots_model <- avPlots.invis(model_input)
  avplots_length <- length(avplots_model)
  avplots_yrange <- range(bind_rows(lapply(avplots_model, FUN=data.frame))[[2]])
  #Create the added variable plots using ggplot
  avplots_ggplot <- vector('list', avplots_length)
  for (i in 1:avplots_length) {
    DATA <- data.frame(avplots_model[[i]])
    avplots_ggplot[[i]] <- ggplot2::ggplot(aes_string(x = colnames(!!DATA)[1],
                                                      y = colnames(!!DATA)[2]),
                                           data = DATA) +
      scale_y_continuous(limits=avplots_yrange, labels=function(y) y*1000) +
      theme_classic() +
      geom_point(colour = 'mediumorchid4', alpha=0.5) + 
      geom_smooth(method = 'lm', color = 'seagreen4', fill='seagreen4', alpha=0.3,
                  formula = y ~ x, linetype = 'dashed') +
      xlab(paste0('Predictor Residual \n (', 
                  names(DATA)[1], ' | others)')) +
      ylab(paste0('Response Residual \n (',
                  ifelse(is.null(ylabel), 
                         paste0(names(DATA)[2], ' | others'), ylabel), ')')) }
  #Return output object
  avplots_ggplot }
# END ADDED VARIABLE PLOT FUNCTION

# Individual and partial regression plot creation
for (i in 1:length(vbl_clm)) 
  local({ #Have to do local or else plot later only does each plot in most recent i
    i <- i
    # Individual regressions of environmental parameters (no multiple regression weighting)
    lm_headdown_plots[[i]] <<- multlm_headdown[["Main Period 3"]][[vbl_clm[i]]] %>%
      augment() %>%
      melt(measure.vars = envparam_names, variable.name = c("envparam")) %>%
      ggplot(., aes(value, iso_iter)) +
      theme_classic() +
      geom_point(color="dodgerblue4", alpha=0.5) +
      geom_smooth(method = "lm", color="firebrick", fill="firebrick", alpha=0.3) +
      scale_y_continuous(name=vbl_clm[i]) +
      theme(plot.title = element_text(hjust=0.5)) +
      facet_wrap(~envparam, scales = "free_x", ncol=length(envparam_names))
    
    # Partial regression (added variable) plots
    multlm_headdown_addvar_iter <- gg_avplots(multlm_headdown[["Main Period 3"]][[vbl_clm[i]]], ylabel=vbl_clm[i])
    multlm_headdown_addvar[[i]] <<- plot_grid(plotlist = multlm_headdown_addvar_iter, ncol=length(envparam_names), align = "hv")
  })

title_lm <- ggdraw() + 
  draw_label("Individual parameter linear regression plots", x = 0.5, hjust = 0.5, size=18) +
  theme(plot.margin = margin(0, 0, 0, 4))
title_multlm <- ggdraw() + 
  draw_label("Added variable plots of multiple regression", x = 0.5, hjust = 0.5, size=18) +
  theme(plot.margin = margin(0, 0, 0, 4))

# Individual regression plots for all isotopes
windows(height=10, width=16)
#pdf("Figures/lm_headdown_individual_plots.pdf", height=10, width=16)
plot_grid(title_lm, plotlist = lm_headdown_plots, ncol=1, align = "hv",
          rel_heights = c(0.1,seq(1,1,length=length(lm_headdown_plots))))   
#ggsave("Figures/lm_headdown_individual_plots.png", height=10, width=16)
#dev.off()

# Partial regressions for all isotopes
windows(height=10, width=16)
#pdf("Figures/multlm_headdown_partial_plots.pdf", height=10, width=16)
plot_grid(title_multlm, plotlist = multlm_headdown_addvar, ncol=1, align = "hv",
          rel_heights = c(0.1,seq(1,1,length=length(lm_headdown_plots))))   
#ggsave("Figures/multlm_headdown_partial_plots.png", height=10, width=16)
#dev.off()

#=====Stream spatial plots=====#
# This makes a violin plot of the general structure of stream isotopes by the different streams.
overall_mean_iso <- stream_subset %>%
  summarize(across(vbl_clm, mean))
overall_sd_iso <- stream_subset %>%
  summarize(across(vbl_clm, sd))
stream_violin_plot_iter <- list()
colorset <- c(rep("deepskyblue3",3), rep("darkgoldenrod2",3), "violetred4")
stream_spatial_labels <- c("Pituffik: Mouth", "Pituffik: Snoutwash", "Pituffik: Ice Wall",
                           "Sioraq: Mouth", "Sioraq: Tuto", "Sioraq: Pingorsuit",
                           "Amitsuarsuk")

for (i in 1:length(vbl_clm)) 
  local({ #Have to do local or else plot later only does each plot in most recent i
    i <- i
    stream_violin_plot_iter[[i]] <<- ggplot(data=stream_subset) +
      theme_classic() +
      geom_violin(aes(x=site_name, y=.data[[vbl_clm[i]]]*1000, color=site_name, fill=site_name), scale="width", linetype='dotted', draw_quantiles = c(0.25, 0.75)) +
      geom_violin(aes(x=site_name, y=.data[[vbl_clm[i]]]*1000, color=site_name), scale="width", fill='transparent', draw_quantiles = c(0.5)) +
      geom_hline(yintercept=overall_mean_iso[[i]]*1000, linetype="dashed", color="gray60") +
      geom_hline(yintercept=overall_mean_iso[[i]]*1000+overall_sd_iso[[i]]*1000, linetype="dotted", color="gray60") +
      geom_hline(yintercept=overall_mean_iso[[i]]*1000-overall_sd_iso[[i]]*1000, linetype="dotted", color="gray60") +
      scale_color_manual(values=colorset) +
      scale_fill_manual(values=alpha(colorset, 0.4)) +
      scale_y_continuous(name=paste(vbl_clm[i],"(‰)"), position = "left") +
      scale_x_discrete(labels=toupper(stream_spatial_labels), position = "bottom", name=NULL) +
      theme(legend.position = "none",
            axis.title.x = element_text(size=16, color='gray40'),
            axis.text.x =  element_text(size=10, color=colorset, angle=30, hjust=1),
            axis.ticks.x = element_line(color='gray40'),
            axis.line.x = element_line(color='gray40'),
            axis.title.y = element_text(size=16, color='gray40'),
            axis.text.y = element_text(size=14, color='gray40'),
            axis.ticks.y = element_line(color='gray40'),
            axis.line.y = element_line(color='gray40'))
  })

stream_violin_plot <- plot_grid(plotlist = stream_violin_plot_iter, ncol=1, align = "hv")

# Plotting LELs for stream sites with more than 3 samples
stream_subset$basin_name <- factor(stream_subset$basin_name, levels = c("Pituffik", "Sioraq", "Amitsuarsuk"))
colorset <- setNames(c("deepskyblue3", "darkgoldenrod2", "violetred4"), levels(stream_subset$basin_name))
stream_lel_bystream_plot <- ggplot() +
  theme_classic() +
  geom_abline(aes(slope=8, intercept=10), linetype="solid", lwd=1, color="gray40") + # GMWL
  geom_abline(aes(slope=7.46, intercept=-3.30), linetype="dashed", lwd=1, color="gray40") + # GNIP LMWL
  geom_abline(aes(slope=6.95, intercept=-18.24), linetype="dotted", lwd=1, color="gray40") + # Water Vapor Akers 2020 LMWL
  geom_point(aes(x=d18O*1000, y=d2H*1000, group=site_name, color=basin_name), data=stream_subset, alpha=0.3) +
  geom_smooth(aes(x=d18O*1000, y=d2H*1000, group=site_name, color=basin_name), fill=NA, data=stream_subset, method="lm") +
  scale_color_manual(values=colorset) +
  scale_y_continuous(name="d2H (‰)") +
  scale_x_continuous(name="d18O (‰)") +
  coord_cartesian(xlim=c(-25,-16), ylim=c(-182,-138)) +
  theme(legend.position = c(0.95,0.05),
        legend.justification = c("right", "bottom"),
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40'),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

windows(height=8, width=12)
#pdf("Figures/Stream_Spatial_Violin_Plot.pdf", height=8, width=12)
plot_grid(plotlist = list(stream_violin_plot, stream_lel_bystream_plot),
          ncol=2, align = "hv", rel_widths = c(0.5,1))
#ggsave("Figures/Stream_Spatial_Violin_Plot.png", height=8, width=12, dpi=600)
#dev.off()

#### END SPATIAL ANALYSES PLOTS ####

####===================================================####
####=====Times Series of Lake and Stream Isotopes======####
####===================================================####
#====Lake temporal plots====#
  # Selected set of lakes with early and late samples
  iso_ts_plot_lakeset <- list()
  ylims <- data.frame(c(-23,-180,-25),c(-8,-100,10))
  colnames(ylims) <- c("ymin", "ymax")
  for (i in 1:length(vbl_clm_index)) 
    local({ #Have to do local or else plot later only does each plot in most recent i
      i <- i
      color_select <- as.character(iso_color_index[iso_color_index$iso == colnames(water_iso_periodcompare[vbl_clm_index[i]]),2])
      iso_ts_plot_lakeset[[i]] <<- ggplot() +
        theme_classic() +
        geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                      y=.data[[vbl_clm[i]]]*1000, group = factor(site_name)),
                  color=color_select, data=water_iso_periodcompare %>% filter(timing != 3)) +
        geom_point(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                       y=.data[[vbl_clm[i]]]*1000, group = factor(site_name), size=lake_surface_area), shape=21,
                   color=color_select, fill=alpha(color_select, 0.3), data=water_iso_periodcompare %>% filter(timing != 3)) +
        scale_size(range = c(1,5)) +
        scale_y_continuous(name=paste(colnames(water_iso_periodcompare[vbl_clm_index[i]]),"(‰)"), position = "left",
                           limits = c(ylims[i,1],ylims[i,2])) +
        scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                     date_breaks="2 weeks", date_labels = "%d %b") +
        theme(legend.position = "none",
              axis.title.x = element_text(size=16, color='gray40'),
              axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
              axis.ticks.x = element_line(color='gray40'),
              axis.line.x = element_line(color='gray40'),
              axis.title.y = element_text(size=16, color=color_select),
              axis.text.y = element_text(size=14, color=color_select),
              axis.ticks.y = element_line(color=color_select),
              axis.line.y = element_line(color=color_select))
    })
  
  # Focused time series of Lake Potato and Power Lake
  lake_ts_index <- c("Lake Potato", "Power Lake")
  lake_ts_plot <- list()
  title_iter <- list()
  #ylims <- data.frame(c(-23,-180,-10),c(-14,-120,10))
  ylims <- data.frame(c(-23,-180,-25),c(-8,-100,10))
  colnames(ylims) <- c("ymin", "ymax")
  
  data_iter <- water_iso %>%
    filter(site_name %in% lake_ts_index) %>%
    filter(yr == 2018)
  iso_ts_plot_lake <- list()
  #title_iter <- ggdraw() + 
  #  draw_label("Lakes", x = 0, hjust = 0, size=14) +
  #  theme(plot.margin = margin(0, 0, 0, 4))
  for (i in 1:length(vbl_clm_index)) 
    local({ #Have to do local or else plot later only does each plot in most recent i
      i <- i
      color_select <- as.character(iso_color_index[iso_color_index$iso == colnames(data_iter[vbl_clm_index[i]]),2])
      iso_ts_plot_lake[[i]] <<- ggplot(data=data_iter) +
        theme_classic() +
        geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                      y=.data[[vbl_clm[i]]]*1000, group = factor(site_name), linetype=factor(site_name)), color=color_select) +
        geom_point(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                       y=.data[[vbl_clm[i]]]*1000, group = factor(site_name), shape=factor(site_name)),
                   color=color_select, fill=alpha(color_select, 0.3), size=2) + # Scaling to match other plot
        scale_shape_manual(values=c(21, 23)) +
        scale_y_continuous(name=paste(colnames(data_iter[vbl_clm_index[i]]),"(‰)"), position = "left",
                           limits = c(ylims[i,1],ylims[i,2])) +
        scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                     date_breaks="2 weeks", date_labels = "%d %b") +
        theme(legend.position = "none",
              axis.title.x = element_text(size=16, color='gray40'),
              axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
              axis.ticks.x = element_line(color='gray40'),
              axis.line.x = element_line(color='gray40'),
              axis.title.y = element_text(size=16, color=color_select),
              axis.text.y = element_text(size=14, color=color_select),
              axis.ticks.y = element_line(color=color_select),
              axis.line.y = element_line(color=color_select))
    })
  
  # E/I Lake Time Series Plots
  e_i_periodcompare <- lake_e_i_fulldata_plot %>%
    filter(!is.na(timing) & laketype != "altered") %>%
    filter(site_name %in% unique(water_iso_periodcompare$site_name))
  
  e_i_ts_plot_lakeset <- ggplot() +
    theme_classic() +
    geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                  y=E_I_d18O_bayes, group = factor(site_name)),
              color="brown4", data=e_i_periodcompare %>% filter(timing != 3)) +
    geom_point(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                   y=E_I_d18O_bayes, group = factor(site_name), size=lake_surface_area), shape=21,
               color="brown4", fill=alpha("brown4", 0.3), data=e_i_periodcompare %>% filter(timing != 3)) +
    scale_size(range = c(1,5)) +
    scale_y_continuous(name="E/I", position = "left",
                       limits = c(0,1)) +
    scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                 date_breaks="2 weeks", date_labels = "%d %b") +
    theme(legend.position = "none",
          axis.title.x = element_text(size=16, color='gray40'),
          axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
          axis.ticks.x = element_line(color='gray40'),
          axis.line.x = element_line(color='gray40'),
          axis.title.y = element_text(size=16, color="brown4"),
          axis.text.y = element_text(size=14, color="brown4"),
          axis.ticks.y = element_line(color="brown4"),
          axis.line.y = element_line(color="brown4"))
  
  data_iter <- lake_e_i_fulldata_plot %>%
    filter(site_name %in% lake_ts_index) %>%
    filter(yr == 2018)
  e_i_ts_plot_lake <- ggplot(data=data_iter) +
    theme_classic() +
    geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                  y=E_I_d18O_bayes, group = factor(site_name), linetype=factor(site_name)), color="brown4") +
    geom_point(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                   y=E_I_d18O_bayes, group = factor(site_name), shape=factor(site_name)),
               color="brown4", fill=alpha("brown4", 0.3), size=2) + # Scaling to match other plot
    scale_shape_manual(values=c(21, 23)) +
    scale_y_continuous(name="E/I", position = "left",
                       limits = c(0,1)) +
    scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                 date_breaks="2 weeks", date_labels = "%d %b") +
    theme(legend.position = "none",
          axis.title.x = element_text(size=16, color='gray40'),
          axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
          axis.ticks.x = element_line(color='gray40'),
          axis.line.x = element_line(color='gray40'),
          axis.title.y = element_text(size=16, color="brown4"),
          axis.text.y = element_text(size=14, color="brown4"),
          axis.ticks.y = element_line(color="brown4"),
          axis.line.y = element_line(color="brown4"))
  
  # Temperature record 2018
  color_select <- "seagreen4"
  wxdata_select <- pfk_isowx_day %>%
    filter(yr == 2018)
  tavg_ts_plot_2018 <- ggplot(data=wxdata_select) +
    theme_classic() +
    geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                  y=tavg), color=color_select) +
    scale_y_continuous(name="Air Temperature (°C)", position = "left", limits = c(-2,15)) +
    scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                 date_breaks="2 weeks", date_labels = "%d %b") +
    theme(legend.position = "none",
          axis.title.x = element_text(size=16, color='gray40'),
          axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
          axis.ticks.x = element_line(color='gray40'),
          axis.line.x = element_line(color='gray40'),
          axis.title.y = element_text(size=16, color=color_select),
          axis.text.y = element_text(size=14, color=color_select),
          axis.ticks.y = element_line(color=color_select),
          axis.line.y = element_line(color=color_select))
  
  # Precipitation record 2018
  color_select <- "navyblue"
  wxdata_select <- pfk_isowx_day %>%
    filter(yr == 2018)
  prcp_ts_plot_2018 <- ggplot(data=wxdata_select) +
    theme_classic() +
    geom_col(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                 y=prcp_H2O), fill=color_select) +
    scale_y_continuous(name="Precipitation (mm)", position = "left", limits = c(0,10.5)) +
    scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                 date_breaks="2 weeks", date_labels = "%d %b") +
    theme(legend.position = "none",
          axis.title.x = element_text(size=16, color='gray40'),
          axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
          axis.ticks.x = element_line(color='gray40'),
          axis.line.x = element_line(color='gray40'),
          axis.title.y = element_text(size=16, color=color_select),
          axis.text.y = element_text(size=14, color=color_select),
          axis.ticks.y = element_line(color=color_select),
          axis.line.y = element_line(color=color_select))
  
  # PET record 2018 (Singer et al. 2021 data)
  color_select <- "darkorange"
  wxdata_select <- pfk_isowx_day %>%
    filter(yr == 2018)
  pet_ts_plot_2018 <- ggplot(data=wxdata_select) +
    theme_classic() +
    geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                  y=pet), color=color_select) +
    scale_y_continuous(name="PET (mm d-1)", position = "left", limits = c(0,3)) +
    scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                 date_breaks="2 weeks", date_labels = "%d %b") +
    theme(legend.position = "none",
          axis.title.x = element_text(size=16, color='gray40'),
          axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
          axis.ticks.x = element_line(color='gray40'),
          axis.line.x = element_line(color='gray40'),
          axis.title.y = element_text(size=16, color=color_select),
          axis.text.y = element_text(size=14, color=color_select),
          axis.ticks.y = element_line(color=color_select),
          axis.line.y = element_line(color=color_select))
  
  iso_ts_plot_lakeset_wx <- append(iso_ts_plot_lakeset, list(e_i_ts_plot_lakeset, tavg_ts_plot_2018, prcp_ts_plot_2018, pet_ts_plot_2018)) # Joining time series to wx data
  iso_ts_plot_lake_wx <- append(iso_ts_plot_lake, list(e_i_ts_plot_lake, tavg_ts_plot_2018, prcp_ts_plot_2018, pet_ts_plot_2018)) # Joining time series to wx data
  windows(height=12, width=6)
  #pdf("Figures/TS_lakes_2018.pdf", height=12, width=6)
  plot_grid(plotlist = append(iso_ts_plot_lakeset_wx, iso_ts_plot_lake_wx),
            ncol=2, align = "v", byrow = FALSE,
            rel_heights = c(seq(1,1,length=length(iso_ts_plot_lake)), 1, 0.5, 0.5, 0.5))
  #ggsave("Figures/TS_lakes_2018.png", height=12, width=6, dpi=600)
  #dev.off()
      
#======Stream temporal plots======#
river_ts_north_index <- c("North River Mouth", "North River Shelter 2 South Fork", "North River Shelter 5.8")
river_ts_south_index <- c("South River Mouth", "Fox Canyon Bridge North Fork", "Fox Canyon Bridge South Fork")
river_ts_index <- list(river_ts_north_index, river_ts_south_index)
river_ts_label_index <- c("Pituffik River", "Sioraq River")
yr_index <- c(2018, 2019)
river_ts_plot <- list()
river_ts_plot_byyear <- list()
title_iter <- list()
ylims <- data.frame(c(-25,-185,6),c(-19,-140,14))
colnames(ylims) <- c("ymin", "ymax")

for (k in 1:length(yr_index)) {
  for (j in 1:length(river_ts_index)) {
    data_iter <- water_iso %>%
      filter(yr == yr_index[k]) %>%
      filter(site_name %in% river_ts_index[[j]])
    data_iter$site_name <- factor(data_iter$site_name, levels = river_ts_index[[j]])
    iso_ts_plot_stream <- list()
    title_iter <- ggdraw() + 
      draw_label(river_ts_label_index[j], x = 0, hjust = 0, size=14) +
      theme(plot.margin = margin(0, 0, 0, 4))
    for (i in 1:length(vbl_clm_index)) 
      local({ #Have to do local or else plot later only does each plot in most recent i
        i <- i
        color_select <- as.character(iso_color_index[iso_color_index$iso == colnames(data_iter[vbl_clm_index[i]]),2])
        iso_ts_plot_stream[[i]] <<- ggplot(data=data_iter) +
          theme_classic() +
          geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                        y=.data[[vbl_clm[i]]]*1000, group = factor(site_name), linetype=factor(site_name)), color=color_select) +
          geom_point(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                         y=.data[[vbl_clm[i]]]*1000, group = factor(site_name), shape=factor(site_name)), color=color_select) +
          scale_y_continuous(name=paste(colnames(data_iter[vbl_clm_index[i]]),"(‰)"), position = "left",
                             limits = c(ylims[i,1],ylims[i,2])) +
          scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                       date_breaks="2 weeks", date_labels = "%d %b") +
          theme(legend.position = "none",
                axis.title.x = element_text(size=16, color='gray40'),
                axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
                axis.ticks.x = element_line(color='gray40'),
                axis.line.x = element_line(color='gray40'),
                axis.title.y = element_text(size=16, color=color_select),
                axis.text.y = element_text(size=14, color=color_select),
                axis.ticks.y = element_line(color=color_select),
                axis.line.y = element_line(color=color_select))
      })
    
    # Temperature record
    color_select <- "seagreen4"
    wxdata_select <- pfk_isowx_day %>%
      filter(yr == yr_index[k])
    tavg_ts_plot <- ggplot(data=wxdata_select) +
      theme_classic() +
      geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                    y=tavg), color=color_select) +
      scale_y_continuous(name="Air Temperature (°C)", position = "left", limits = c(-2,15)) +
      scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                   date_breaks="2 weeks", date_labels = "%d %b") +
      theme(legend.position = "none",
            axis.title.x = element_text(size=16, color='gray40'),
            axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
            axis.ticks.x = element_line(color='gray40'),
            axis.line.x = element_line(color='gray40'),
            axis.title.y = element_text(size=16, color=color_select),
            axis.text.y = element_text(size=14, color=color_select),
            axis.ticks.y = element_line(color=color_select),
            axis.line.y = element_line(color=color_select))
    
    # Precipitation record
    color_select <- "navyblue"
    wxdata_select <- pfk_isowx_day %>%
      filter(yr == yr_index[k])
    prcp_ts_plot <- ggplot(data=wxdata_select) +
      theme_classic() +
      geom_col(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                   y=prcp_H2O), fill=color_select) +
      scale_y_continuous(name="Precipitation (mm)", position = "left", limits = c(0,10.5)) +
      scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                   date_breaks="2 weeks", date_labels = "%d %b") +
      theme(legend.position = "none",
            axis.title.x = element_text(size=16, color='gray40'),
            axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
            axis.ticks.x = element_line(color='gray40'),
            axis.line.x = element_line(color='gray40'),
            axis.title.y = element_text(size=16, color=color_select),
            axis.text.y = element_text(size=14, color=color_select),
            axis.ticks.y = element_line(color=color_select),
            axis.line.y = element_line(color=color_select))
    
    # PET record 2018 (Singer et al. 2021 data)
    color_select <- "darkorange"
    wxdata_select <- pfk_isowx_day %>%
      filter(yr == yr_index[k])
    pet_ts_plot <- ggplot(data=wxdata_select) +
      theme_classic() +
      geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes
                    y=pet), color=color_select) +
      scale_y_continuous(name="PET (mm d-1)", position = "left", limits = c(0,4)) +
      scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
                   date_breaks="2 weeks", date_labels = "%d %b") +
      theme(legend.position = "none",
            axis.title.x = element_text(size=16, color='gray40'),
            axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
            axis.ticks.x = element_line(color='gray40'),
            axis.line.x = element_line(color='gray40'),
            axis.title.y = element_text(size=16, color=color_select),
            axis.text.y = element_text(size=14, color=color_select),
            axis.ticks.y = element_line(color=color_select),
            axis.line.y = element_line(color=color_select))
    
    river_ts_plot[[(k-1)*length(river_ts_index)+j]] <- plot_grid(title_iter,
                                                                 plotlist=append(iso_ts_plot_stream, list(tavg_ts_plot, prcp_ts_plot, pet_ts_plot)),
                                                                 ncol=1, align = "v",
                                                                 rel_heights = c(0.2, seq(1,1,length=length(iso_ts_plot_stream)), 0.5, 0.5, 0.5))
    
  }
}
windows(height=10, width=14)
#pdf("Figures/TS_rivers.pdf", height=10, width=14)
plot_grid(plotlist = river_ts_plot, ncol=4, align = "hv")
#ggsave("Figures/TS_rivers.png", height=10, width=14, dpi=600)
#dev.off()


#======Lake Interannual Variability plots======#
# Subplot comparing sample values from all three sampling periods in 2018 and 2019
iso_lake_3period_plot <- list()
ylims <- data.frame(c(-20,-155,-43),c(-4,-80,5))
colnames(ylims) <- c("ymin", "ymax")
for (i in 1:length(vbl_clm_index)) 
  local({ #Have to do local or else plot later only does each plot in most recent i
    i <- i
    color_select <- as.character(iso_color_index[iso_color_index$iso == colnames(water_iso_periodcompare[vbl_clm_index[i]]),2])
    iso_lake_3period_plot[[i]] <<- ggplot() +
      theme_classic() +
      geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                    y=.data[[vbl_clm[i]]]*1000, group = factor(site_name)),
                color=lighten(color_select, factor=0.7), data=water_iso_periodcompare %>% filter(timing != 3)) +
      geom_point(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                     y=.data[[vbl_clm[i]]]*1000, group = factor(site_name), shape=as.factor(timing), size=lake_surface_area),
                 color=lighten(color_select, factor=0.7), fill=alpha(lighten(color_select, factor=0.7), 0.3),
                 data=water_iso_periodcompare %>% filter(timing != 3)) +
      scale_shape_manual(values=c(1,21)) +
      scale_size(range = c(1,5)) +
      geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                    y=.data[[vbl_clm[i]]]*1000, group = factor(site_name)),
                linetype = "dashed", color=darken(color_select, factor=0.3), data=water_iso_periodcompare %>% filter(timing != 2)) +
      geom_point(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                     y=.data[[vbl_clm[i]]]*1000, size=lake_surface_area),
                 shape = 22, color=darken(color_select, factor=0.3), fill=alpha(darken(color_select, factor=0.3), 0.3),
                 data=water_iso_periodcompare %>% filter(timing == 3)) +
      scale_y_continuous(name=paste(colnames(water_iso_periodcompare[vbl_clm_index[i]]),"(‰)"), position = "left",
                         limits = c(ylims[i,1],ylims[i,2])) +
      scale_x_date(name=NULL, position = "bottom", date_labels = "%d %b") +
      theme(legend.position = "none",
            axis.title.x = element_text(size=16, color='gray40'),
            axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
            axis.ticks.x = element_line(color='gray40'),
            axis.line.x = element_line(color='gray40'),
            axis.title.y = element_text(size=16, color=color_select),
            axis.text.y = element_text(size=14, color=color_select),
            axis.ticks.y = element_line(color=color_select),
            axis.line.y = element_line(color=color_select))
  })

# Subplot comparing sample values from the end of summer 2018 to midsummer 2019
iso_lake_endsum_violinplot <- list()
for (i in 1:length(vbl_clm_index)) 
  local({ #Have to do local or else plot later only does each plot in most recent i
    i <- i
    color_select <- as.character(iso_color_index[iso_color_index$iso == colnames(water_iso_periodcompare[vbl_clm_index[i]]),2])
    iso_lake_endsum_violinplot[[i]] <<- ggplot(data = water_iso_periodcompare %>% filter(timing != 1)) +
      theme_classic() +
      geom_violin(aes(x=timing, y=.data[[vbl_clm[i]]]*1000, group=timing, color=as.factor(timing), fill=as.factor(timing)), 
                  scale="width", draw_quantiles = c(0.25, 0.5, 0.75)) +
      geom_line(aes(timing, .data[[vbl_clm[i]]]*1000, group=site_name), color=color_select, linetype = "dashed") + # Lines connecting 2018 to 2019
      geom_point(aes(x=timing, y=.data[[vbl_clm[i]]]*1000, shape=as.factor(timing),
                     color=as.factor(timing), fill=as.factor(timing), size=lake_surface_area)) + # Individual lake points
      scale_shape_manual(values=c(21,22)) +
      scale_size(range = c(1,5)) +
      scale_color_manual(values=c(lighten(color_select, factor=0.7), darken(color_select, factor=0.3))) +
      scale_fill_manual(values=c(alpha(lighten(color_select, factor=0.7), 0.3), alpha(darken(color_select, factor=0.3), 0.3))) +
      scale_x_continuous(name=NULL, labels=c("2018", "2019"), position = "bottom", breaks=c(2,3)) +
      scale_y_continuous(name=paste(colnames(water_iso_periodcompare[vbl_clm_index[i]]),"(‰)"), position = "left",
                         limits = c(ylims[i,1],ylims[i,2])) +
      theme(legend.position = "none",
            axis.title.x = element_text(size=16, color='gray40'),
            axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
            axis.ticks.x = element_line(color='gray40'),
            axis.line.x = element_line(color='gray40'),
            axis.title.y = element_text(size=16, color=color_select),
            axis.text.y = element_text(size=14, color=color_select),
            axis.ticks.y = element_line(color=color_select),
            axis.line.y = element_line(color=color_select))
  })

# E I plots of same nature as above
e_i_lake_3period_plot <- ggplot() +
  theme_classic() +
  geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                y=E_I_d18O_bayes, group = factor(site_name)),
            color=lighten("brown4", factor=0.7), data=e_i_periodcompare %>% filter(timing != 3)) +
  geom_point(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                 y=E_I_d18O_bayes, group = factor(site_name), shape=as.factor(timing), size=lake_surface_area),
             color=lighten("brown4", factor=0.7), fill=alpha(lighten("brown4", factor=0.7), 0.3),
             data=e_i_periodcompare %>% filter(timing != 3)) +
  scale_shape_manual(values=c(1,21)) +
  scale_size(range = c(1,5)) +
  geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                y=E_I_d18O_bayes, group = factor(site_name)),
            linetype = "dashed", color=darken("brown4", factor=0.3), data=e_i_periodcompare %>% filter(timing != 2)) +
  geom_point(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                 y=E_I_d18O_bayes, size=lake_surface_area),
             shape = 22, color=darken("brown4", factor=0.3), fill=alpha(darken("brown4", factor=0.3), 0.3),
             data=e_i_periodcompare %>% filter(timing == 3)) +
  scale_y_continuous(name="E/I", position = "left", limits = c(0,1)) +
  scale_x_date(name=NULL, position = "bottom", date_labels = "%d %b") +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color="brown4"),
        axis.text.y = element_text(size=14, color="brown4"),
        axis.ticks.y = element_line(color="brown4"),
        axis.line.y = element_line(color="brown4"))

e_i_lake_endsum_violinplot <- ggplot(data = e_i_periodcompare %>% filter(timing != 1)) +
  theme_classic() +
  geom_violin(aes(x=timing, y=E_I_d18O_bayes, group=timing, color=as.factor(timing), fill=as.factor(timing)), 
              scale="width", draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_line(aes(x=timing, y=E_I_d18O_bayes, group=site_name), color="brown4", linetype = "dashed") + # Lines connecting 2018 to 2019
  geom_point(aes(x=timing, y=E_I_d18O_bayes, shape=as.factor(timing),
                 color=as.factor(timing), fill=as.factor(timing), size=lake_surface_area)) + # Individual lake points
  scale_shape_manual(values=c(21,22)) +
  scale_size(range = c(1,5)) +
  scale_color_manual(values=c(lighten("brown4", factor=0.7), darken("brown4", factor=0.3))) +
  scale_fill_manual(values=c(alpha(lighten("brown4", factor=0.7), 0.3), alpha(darken("brown4", factor=0.3), 0.3))) +
  scale_x_continuous(name=NULL, labels=c("2018", "2019"), position = "bottom", breaks=c(2,3)) +
  scale_y_continuous(name="E/I", position = "left",limits = c(0,1)) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40', angle=45, hjust=1),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color="brown4"),
        axis.text.y = element_text(size=14, color="brown4"),
        axis.ticks.y = element_line(color="brown4"),
        axis.line.y = element_line(color="brown4"))

windows(height=14, width=8)
#pdf("Figures/Lake_compare_2018-2019.pdf", height=14, width=7)
lake_3period_plot <- append(iso_lake_3period_plot, list(e_i_lake_3period_plot))
lake_ensum_violinplot <- append(iso_lake_endsum_violinplot, list(e_i_lake_endsum_violinplot))
plot_grid(plotlist = append(lake_3period_plot, lake_ensum_violinplot), ncol=2, align = "hv", byrow = FALSE)
#ggsave("Figures/Lake_compare_2018-2019.png", height=14, width=7, dpi = 600)
#dev.off()

####==========END PLOTTING============####



#=====================================================================#
####==================Supplemental=================================####
#=====================================================================#

#======LWL of samples by groups plot======#
water_iso_lwl <- subset(water_iso, type == "precipitation rain" | type == "precipitation snow")
colorset_inv <- rev(c("slateblue4", "lightpink2"))

# Full LMWL plot
iso_lwl_plot <- ggplot() +
  theme_classic() +
  geom_abline(aes(slope=8, intercept=10), linetype="solid", lwd=1, color="gray40") + # GMWL
  geom_abline(aes(slope=7.46, intercept=-3.30), linetype="dashed", lwd=1, color="gray40") + # GNIP LMWL
  geom_abline(aes(slope=6.95, intercept=-18.24), linetype="dotted", lwd=1, color="gray40") + # Water Vapor Akers 2020 LMWL
  geom_point(aes(x=d18O*1000, y=d2H*1000, color=fct_rev(type)), data=water_iso_lwl, alpha=0.3) +
  geom_smooth(aes(x=d18O*1000, y=d2H*1000, color=fct_rev(type), fill=fct_rev(type)), data=water_iso_lwl, method="lm") +
  geom_rect(aes(xmin=-23, xmax=-15, ymin=-170, ymax=-125), color="gray60", fill=alpha("gray", 0))+
  scale_color_manual(values=colorset_inv, guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=alpha(colorset_inv, 0.4), guide=guide_legend(reverse=TRUE)) +
  scale_y_continuous(name="d2H (‰)") +
  scale_x_continuous(name="d18O (‰)") +
  coord_cartesian(xlim=c(-40,-4), ylim=c(-300,-80)) +
  theme(legend.position = c(0.95,0.05),
        legend.justification = c("right", "bottom"),
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40'),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

# Zoomed in LMWL plot
iso_lwl_plot_zoom <- ggplot() +
  theme_classic() +
  geom_abline(aes(slope=8, intercept=10), linetype="solid", lwd=1, color="gray40") + # GMWL
  geom_abline(aes(slope=7.46, intercept=-3.30), linetype="dashed", lwd=1, color="gray40") + # GNIP LMWL
  geom_abline(aes(slope=6.95, intercept=-18.24), linetype="dotted", lwd=1, color="gray40") + # Water Vapor Akers 2020 LMWL
  geom_point(aes(x=d18O*1000, y=d2H*1000, color=fct_rev(type)), data=water_iso_lwl, alpha=0.3) +
  geom_smooth(aes(x=d18O*1000, y=d2H*1000, color=fct_rev(type), fill=fct_rev(type)), data=water_iso_lwl, method="lm") +
  scale_color_manual(values=colorset_inv, guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=alpha(colorset_inv, 0.4), guide=guide_legend(reverse=TRUE)) +
  scale_y_continuous(name="d2H (‰)") +
  scale_x_continuous(name="d18O (‰)") +
  coord_cartesian(xlim=c(-23,-15), ylim=c(-170,-125)) +
  theme(legend.position = c(0.95,0.05),
        legend.justification = c("right", "bottom"),
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40'),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

windows(height=6, width=12)
#pdf("Figures/pfk_LMWL_prcp_plots.pdf", height=6, width=12)
plot_grid(plotlist = list(iso_lwl_plot, iso_lwl_plot_zoom), ncol=2, align = "hv")
#ggsave("Figures/pfk_LMWL_prcp_plots.png", height=6, width=12, dpi = 600)
#dev.off()

#======Lake LEL analyses and plots by sample period timing=======#
lel_lake_bytiming <- water_iso_periodcompare %>% # Running regression by type
  group_by(timing) %>%
  nest() %>%
  mutate(lel = map(data, ~lm(d2H~d18O, data=.)))

lel_slope <- 
  lel_lake_bytiming %>%
  mutate(tidyit = map(lel, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "d18O") %>%
  dplyr::select(timing, estimate, std.error)
colnames(lel_slope) <- c("timing", "slope", "se_slope")
lel_slope$conf_int_slope <- lel_slope$se_slope*qnorm(0.975)

lel_intercept <- 
  lel_lake_bytiming %>%
  mutate(tidyit = map(lel, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(timing, estimate, std.error)
colnames(lel_intercept) <- c("timing", "intercept", "se_intercept")
lel_intercept$conf_int_intercept <- lel_intercept$se_intercept*qnorm(0.975)

lel_r2 <- 
  lel_lake_bytiming %>%
  mutate(glanceit = map(lel, broom::glance)) %>%
  unnest(glanceit) %>%
  dplyr::select(timing, r.squared, p.value, nobs)

lel_lake_bytiming_params <- lel_slope %>% # Combining all data
  inner_join(lel_intercept, by="timing") %>%
  inner_join(lel_r2, by="timing")
print(lel_lake_bytiming_params)

# This makes a plot of the LEL of lakes by the three sampling periods
colorset <- c("deepskyblue3", "darkgoldenrod2", "violetred4")
water_iso_periodcompare$timing <- factor(water_iso_periodcompare$timing)
colorset_inv <- rev(colorset[seq(1:length(unique(water_iso_periodcompare$timing)))])

lake_lel_bytiming_plot <- ggplot() +
  theme_classic() +
  geom_abline(aes(slope=8, intercept=10), linetype="solid", lwd=1, color="gray40") + # GMWL
  geom_abline(aes(slope=7.46, intercept=-3.30), linetype="dashed", lwd=1, color="gray40") + # GNIP LMWL
  geom_abline(aes(slope=6.95, intercept=-18.24), linetype="dotted", lwd=1, color="gray40") + # Water Vapor Akers 2020 LMWL
  geom_point(aes(x=d18O*1000, y=d2H*1000, color=fct_rev(timing)), data=water_iso_periodcompare, alpha=0.3) +
  geom_smooth(aes(x=d18O*1000, y=d2H*1000, color=fct_rev(timing), fill=fct_rev(timing)), data=water_iso_periodcompare, method="lm") +
  scale_color_manual(values=colorset_inv, guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=alpha(colorset_inv, 0.4), guide=guide_legend(reverse=TRUE)) +
  scale_y_continuous(name="d2H (‰)") +
  scale_x_continuous(name="d18O (‰)") +
  coord_cartesian(xlim=c(-20,-4), ylim=c(-155,-80)) +
  theme(legend.position = c(0.95,0.05),
        legend.justification = c("right", "bottom"),
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40'),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

windows(height=6, width=6)
#pdf("Figures/lake_lel_bytiming_plot.pdf", height=6, width=6)
plot(lake_lel_bytiming_plot)
#ggsave("Figures/lake_lel_bytiming_plot.png", height=6, width=6, dpi = 600)
#dev.off()

# Stream dxs plotted by day of year for three basins
colorset <- c(rep("deepskyblue3",3), rep("darkgoldenrod2",3), "violetred4")
stream_spatial_labels <- c("Pituffik: Mouth", "Pituffik: Snoutwash", "Pituffik: Ice Wall",
                           "Sioraq: Mouth", "Sioraq: Tuto", "Sioraq: Pingorsuit",
                           "Amitsuarsuk")
stream_dxs_doy_plot <- ggplot() +
  theme_classic() +
  geom_point(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                 y=dxs*1000, color=site_name, shape=as.character(yr)), data=stream_subset, alpha=0.6, size=2) +
  geom_line(aes(x=as.Date(doy, origin = "2017-12-31"), # Day of Year all put on 2018 for plotting purposes,
                y=dxs*1000, color=site_name, group=interaction(as.character(yr), site_name)), data=stream_subset) +
  scale_color_manual(values=colorset) +
  scale_y_continuous(name="dxs (‰)") +
  scale_x_date(name=NULL, position = "bottom", limits = as.Date(c("2018-06-10", "2018-08-25")),
               date_breaks="2 weeks", date_labels = "%d %b") +
  theme(legend.position = c(0.05,0.05),
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill='transparent'),
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40'),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

windows(height=8, width=8)
#pdf("Figures/stream_dxs_doy_plot.pdf", height=6, width=6)
stream_dxs_doy_plot
#ggsave("Figures/stream_dxs_doy_plot.png", height=6, width=6, dpi=600)
#dev.off()
