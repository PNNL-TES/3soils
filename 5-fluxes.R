# Produce diagnostics, analysis plots, and stats for 3Soils Picarro data
# 
# Ben Bond-Lamberty February 2018

source("0-functions.R")

SCRIPTNAME  	<- "5-fluxes.R"
PROBLEM       <- FALSE


# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE) # open log

printlog("Welcome to", SCRIPTNAME)

printlog("Reading in summarized data...")

sdata <- read_csv(SUMMARYDATA_CLEAN_FILE, na = c("NA", "#VALUE!", "NO HS"), guess_max = 1e6)

# Fill in missing headspace data (TODO)
# Cores that have missing values get assigned the median of all good cores
sdata %>% 
  group_by(SampleID) %>% 
  summarise(hs = mean(Headspace_height_cm, na.rm = TRUE)) %>% 
  summarise(hs = median(hs, na.rm = TRUE)) %>% 
  pull(hs) ->
  mean_headspace
sdata %>% 
  replace_na(list(Headspace_height_cm = mean_headspace)) ->
  sdata

# -----------------------------------------------------------------------------
# Flux computation

# At this point, `sdata` has slopes (CO2 ppm/s and CH4 ppb/s).
# We want to convert this to mg C/s, using
# A = dC/dt * V/M * Pa/RT (cf. Steduto et al. 2002), where
# 	A is the flux (µmol/g/s)
#	  dC/dt is raw respiration as above (mole fraction/s)
# 	V is total chamber volume (cm3)
#	  M is [dry] soil mass (g)
#	  Pa is atmospheric pressure (kPa)
#	  R is universal gas constant (8.3 x 10^3 cm3 kPa mol-1 K-1)
#	  T is air temperature (K)

# The instrument tubing is 455 cm long by ID 1/16"
V_tubing <- (1/16 * 2.54 / 2 ) ^ 2 * pi * 455
# Headspace is in each is the total volume of the sleeve minus the soil volume
V_headspace <- (7.5 / 2) ^ 2 * pi * sdata$Headspace_height_cm
# Replace missing headspace values with the mean
V_headspace[is.na(V_headspace)] <- mean(V_headspace, na.rm = TRUE)

# Internal volume of Picarro? 
V_picarro <- 9 # Assume same as PP-Systems

sdata$V_cm3 <- V_tubing + V_headspace + V_picarro

Pa 			<- 101						# kPa				(Richland is ~120 m asl)
R 			<- 8.3145e+3			# cm3 kPa K−1 mol−1
Tair    <- 273.1 + 20 # TODO: fluxdata$Temperature     # C -> K

# Calculate mass-corrected respiration, µmol/s
CO2_flux_µmol_g_s <- 
  with(sdata,
       CO2_ppm_s / 1 * # from ppm/s to µmol/s
         V_cm3 * Pa / (R * Tair)) # ideal gas law
CH4_flux_µmol_g_s <- 
  with(sdata,
       CH4_ppb_s / 1000 * # from ppb/s to µmol/s
         V_cm3 * Pa / (R * Tair)) # ideal gas law

# Calculate flux of mg C/hr
sdata$CO2_flux_mgC_hr <- CO2_flux_µmol_g_s /
  1e6 * # to mol 
  12 *  # to g C
  1000 * # to mg C
  60 * 60 # to /hr
sdata$CH4_flux_mgC_hr <- CH4_flux_µmol_g_s /
  1e6 * # to mol 
  16 *  # to g C
  1000 *  # to mg C
  60 * 60 # to /hr


# Cumulative emissions
sdata %>% 
  filter(!is.na(DATETIME)) %>% 
  arrange(DATETIME) %>% 
  group_by(SampleID, PHASE) %>% 
  # Compute incubation time
  mutate(inctime_hours = as.numeric(difftime(DATETIME, min(DATETIME), units = "hours")),
         # interpolate missing fluxes
         CO2_flux_mgC_hr_interp = approx(inctime_hours, CO2_flux_mgC_hr, xout = inctime_hours, rule = 2)[['y']],
         CH4_flux_mgC_hr_interp = approx(inctime_hours, CH4_flux_mgC_hr, xout = inctime_hours, rule = 2)[['y']],
         # ...and compute cumulative emissions
         delta_hrs = (inctime_hours - lag(inctime_hours)),
         CO2_flux_mgC = CO2_flux_mgC_hr_interp * delta_hrs,
         cumCO2_flux_mgC_gSoil = c(0, cumsum(CO2_flux_mgC[-1] / DryMass_SoilOnly_g[-1])),
         CH4_flux_mgC = CH4_flux_mgC_hr_interp * delta_hrs,
         cumCH4_flux_mgC_gSoil = c(0, cumsum(CH4_flux_mgC[-1] / DryMass_SoilOnly_g[-1])),
         label = if_else(inctime_hours == max(inctime_hours), SampleID, "")) %>% 
  ungroup ->
  fluxdata

save_data(fluxdata, fn = FLUXDATA_FILE, scriptfolder = FALSE)

p_collar_co2 <- ggplot(fluxdata, aes(inctime_hours, cumCO2_flux_mgC_gSoil, group = SampleID, color = Site)) + 
  geom_point() + geom_line() + geom_text(aes(x = inctime_hours * 1.1, label = label), size = 3) +
  facet_wrap(~PHASE, scales = "free") +
  ggtitle("Cumulative CO2 emissions by core")
print(p_collar_co2)
save_plot("cumulative_co2_inc")

p_collar_ch4 <- ggplot(fluxdata, aes(inctime_hours, cumCH4_flux_mgC_gSoil, group = SampleID, color = Site)) + 
  geom_point() + geom_line() + geom_text(aes(x = inctime_hours * 1.1, label = label), size = 3) +
  facet_wrap(~PHASE, scales = "free") +
  ggtitle("Cumulative CH4 emissions by core")
print(p_collar_ch4)
save_plot("cumulative_ch4_inc")

# Plot final (end of incubation) totals
fluxdata %>% 
  group_by(Site, PHASE, SampleID) %>% 
  summarise(inctime_hours = last(inctime_hours),
            cumCO2_flux_mgC_gSoil = last(cumCO2_flux_mgC_gSoil),
            cumCH4_flux_mgC_gSoil = last(cumCH4_flux_mgC_gSoil)) ->
  fluxdata_final

p_final_co2 <- ggplot(fluxdata_final, aes(Site, cumCO2_flux_mgC_gSoil, color = Site)) +
  geom_boxplot() + geom_point() +
  facet_wrap(~PHASE, scales = "free")
print(p_final_co2)
save_plot("final_co2_inc")
p_final_ch4 <- ggplot(fluxdata_final, aes(Site, cumCH4_flux_mgC_gSoil, color = Site)) +
  geom_boxplot() + geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~PHASE, scales = "free")
print(p_final_ch4)
save_plot("final_ch4_inc")

# DON'T KNOW IF ANYTHING WORKS BEYOND HERE!!!!!
stop("# DON'T KNOW IF ANYTHING WORKS BEYOND HERE!!!!!")

# Site-level means
fluxdata %>%
  filter(!is.na(cumCO2_flux_mgC_gC)) %>% 
  # We need to first interpolate the fluxes to every hour
  group_by(PHASE, Site, SampleID) %>% 
  filter(n() > 2) %>%   # make sure enough samples to interpolate
  group_by(PHASE, Site, SampleID) %>% 
  do(tibble(label = first(.$label),
            inctime_hours = seq_len(max(.$inctime_hours)),
            cumCH4_flux_mgC_gC = approx(.$inctime_hours, .$cumCH4_flux_mgC_gC, xout = inctime_hours, rule = 2)[['y']],
            cumCO2_flux_mgC_gC = approx(.$inctime_hours, .$cumCO2_flux_mgC_gC, xout = inctime_hours, rule = 2)[['y']])) %>% 
  # Now can average for each hour
  group_by(PHASE, Site, inctime_hours) %>% 
  summarise(cumCO2_flux_mgC_gC_sd = sd(cumCO2_flux_mgC_gC, na.rm = TRUE),
            cumCO2_flux_mgC_gC = mean(cumCO2_flux_mgC_gC, na.rm = TRUE),
            cumCH4_flux_mgC_gC_sd = sd(cumCH4_flux_mgC_gC, na.rm = TRUE),
            cumCH4_flux_mgC_gC = mean(cumCH4_flux_mgC_gC, na.rm = TRUE)) ->
  fluxdata_by_site

fd1 <- filter(fluxdata_by_site,
              PHASE %in% c("DROUGHT_PRE_INCUBATION", "DROUGHT_INCUBATION", "FIELD_MOIST_INCUBATION", "SATURATION_INCUBATION"))
fd2 <- filter(fluxdata_by_site,
              !PHASE %in% c("DROUGHT_PRE_INCUBATION", "DROUGHT_INCUBATION", "FIELD_MOIST_INCUBATION", "SATURATION_INCUBATION"))

p_site_co2 <- ggplot(fd1, aes(inctime_hours, cumCO2_flux_mgC_gC, color = Site)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = I(0.5), color = NA,
              aes(fill = Site,
                  ymin = cumCO2_flux_mgC_gC - cumCO2_flux_mgC_gC_sd,
                  ymax = cumCO2_flux_mgC_gC + cumCO2_flux_mgC_gC_sd)) + 
  facet_wrap(~PHASE, scales = "free") + xlab("Incubation (hrs)") +
  ggtitle("Cumulative CO2 emissions by site")
print(p_site_co2)
save_plot("cumulative_co2_site")

p_site_ch4 <- ggplot(fd1, aes(inctime_hours, cumCH4_flux_mgC_gC, color = Site)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = I(0.5), color = NA,
              aes(fill = Site,
                  ymin = cumCH4_flux_mgC_gC - cumCH4_flux_mgC_gC_sd,
                  ymax = cumCH4_flux_mgC_gC + cumCH4_flux_mgC_gC_sd)) + 
  facet_wrap(~PHASE, scales = "free") + xlab("Incubation (hrs)") +
  ggtitle("Cumulative CH4 emissions by site")
print(p_site_ch4)
save_plot("cumulative_ch4_site")


# -----------------------------------------------------------------------------
# Join with C/N data



# -----------------------------------------------------------------------------
# Done! 

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
