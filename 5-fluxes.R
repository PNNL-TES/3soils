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

sdata <- read_csv(SUMMARYDATA_FILE, na = c("NA", "#VALUE!", "NO HS"), guess_max = 1e6)

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
V_headspace <- (7.5 / 2) ^ 2 * pi * 30 - sdata$Headspace_height_cm
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

# Normalize by soil mass
sdata$CO2_flux_mgC_hr_gsoil <- sdata$CO2_flux_mgC_hr / sdata$Soil_dry_weight_equivalent_g
sdata$CH4_flux_mgC_hr_gsoil <- sdata$CH4_flux_mgC_hr / sdata$Soil_dry_weight_equivalent_g

# Cumulative emissions
sdata %>% 
  arrange(inctime_hours) %>% 
  group_by(PHASE, SampleID) %>%
  mutate(CO2_cum_mgC_gsoil = cumsum(CO2_flux_mgC_hr_gsoil)) ->
  fluxdata

qplot(inctime_hours, CO2_cum_mgC_gsoil, data = fluxdata, group = SampleID, geom = "line", color = Site) + 
  facet_wrap(~PHASE) +  theme(strip.text.x = element_text(size = 5))
save_plot("cumulative_co2")

# Site-level means
fluxdata %>%
  filter(PHASE %in% c("DROUGHT_PRE_INCUBATION", "DROUGHT_INCUBATION", "FIELD_MOIST_INCUBATION", "SATURATION_INCUBATION")) %>% 
  group_by(PHASE, Site, round(inctime_hours/12)) %>% 
  summarise(CO2_cum_mgC_gsoil_sd = sd(CO2_cum_mgC_gsoil),
            CO2_cum_mgC_gsoil = mean(CO2_cum_mgC_gsoil)) ->
  fluxdata_by_site

p <- ggplot(fluxdata_by_site, aes(`round(inctime_hours/12)` * 12, CO2_cum_mgC_gsoil, color = Site)) +
  geom_point() + geom_smooth(method = "loess") + facet_wrap(~PHASE)
p <- ggplot(fluxdata_by_site, aes(`round(inctime_hours/12)` * 12, CO2_cum_mgC_gsoil, color = Site)) +
  geom_point() + geom_line() +
  geom_errorbar(alpha = I(0.5),
                aes(ymin = CO2_cum_mgC_gsoil - CO2_cum_mgC_gsoil_sd,
                    ymax = CO2_cum_mgC_gsoil + CO2_cum_mgC_gsoil_sd)) + 
  facet_wrap(~PHASE) + xlab("Incubation (hrs)")
print(p)
save_plot("cumulative_co2_site")


# -----------------------------------------------------------------------------
# Join with C/N data



# -----------------------------------------------------------------------------
# Done! 

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
