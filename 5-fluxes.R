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

sdata <- read_csv(SUMMARYDATA_CLEAN_FILE,
                  #                  col_types = "cccddddddcdcdddddddTdd",
                  na = c("NA", "#VALUE!", "NO HS"), guess_max = 1e6)

# Fill in missing headspace data (TODO)
# Cores that have missing values get assigned the median of all good cores
sdata %>% 
  group_by(SampleID) %>% 
  summarise(hs = mean(HeadSpace_Ht_cm, na.rm = TRUE)) %>% 
  summarise(hs = mean(hs, na.rm = TRUE)) %>% 
  pull(hs) ->
  mean_headspace
sdata %>% 
  replace_na(list(HeadSpace_Ht_cm = mean_headspace)) ->
  sdata

# -----------------------------------------------------------------------------
# Water content computation
sdata %>% 
  replace_na(list(Additional_Wt_toRemove = 0)) %>% 
  mutate(sm_gravimetric = (Total_core_mass_pot_pie_pans_g - DryMass_SoilOnly_g - DryMass_NONsoil_ALL_g - Additional_Wt_toRemove) / DryMass_SoilOnly_g,
         sm_volumetric = (Total_core_mass_pot_pie_pans_g - DryMass_SoilOnly_g - DryMass_NONsoil_ALL_g - Additional_Wt_toRemove) / VolumeSoil_cm3,
         bd = DryMass_SoilOnly_g / VolumeSoil_cm3,
         sm_gravimetric = round(sm_gravimetric, 3),
         sm_volumetric = round(sm_volumetric, 3),
         bd = round(bd, 3)) -> 
  sdata

# Summary plots and table
qplot(sm_gravimetric, SampleID, data = sdata, color=Treatment)
save_plot("wc_gravimetric_sample", height = 6, width = 4)
qplot(DATETIME, sm_gravimetric, data=sdata, color = Treatment)
save_plot("wc_gravimetric_time", height = 4, width = 6)
qplot(sm_volumetric, SampleID, data = sdata, color=Treatment)
save_plot("wc_volumetric_sample", height = 6, width = 4)
qplot(DATETIME, sm_volumetric, data=sdata, color = Treatment)
save_plot("wc_volumetric_time", height = 4, width = 6)

sdata %>% 
  group_by(TREATMENT_PHASE, Site) %>% 
  summarise(sm_grav = mean(sm_gravimetric, na.rm = TRUE),
            sm_grav_sd = sd(sm_gravimetric, na.rm = TRUE),
            sm_vol = mean(sm_volumetric, na.rm = TRUE),
            sm_vol_sd = sd(sm_volumetric, na.rm = TRUE),
            BD = mean(bd, na.rm = TRUE),
            BD_sd = sd(bd, na.rm = TRUE)) ->
  sd_wc_summary

print(sd_wc_summary, n = 50)

qplot(TREATMENT_PHASE, sm_grav, color = Site, data = sd_wc_summary) + 
  geom_errorbar(aes(ymin = sm_grav - sm_grav_sd, ymax = sm_grav + sm_grav_sd)) + 
  facet_wrap(~Site) + coord_flip()
save_plot("wc_grav_summary", width = 6, height = 4)
qplot(TREATMENT_PHASE, sm_vol, color = Site, data = sd_wc_summary) + 
  geom_errorbar(aes(ymin = sm_vol - sm_vol_sd, ymax = sm_vol + sm_vol_sd)) + 
  facet_wrap(~Site) + coord_flip()
save_plot("wc_vol_summary", width = 6, height = 4)
qplot(TREATMENT_PHASE, BD, color = Site, data = sd_wc_summary) + 
  geom_errorbar(aes(ymin = BD - BD_sd, ymax = BD + BD_sd)) + 
  facet_wrap(~Site) + coord_flip()
save_plot("bd_summary", width = 6, height = 4)

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
V_headspace <- (7.5 / 2) ^ 2 * pi * sdata$HeadSpace_Ht_cm
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
  group_by(SampleID, TREATMENT_PHASE) %>% 
  # Compute incubation time
  mutate(inctime_hours = as.numeric(difftime(DATETIME, min(DATETIME), units = "hours")) %>% round(2),
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

# Peyton's cutoffs from 2018-12-04 email
fluxdata %>% 
  select(Site, TREATMENT_PHASE, SampleID, DATETIME, inctime_hours,
         sm_gravimetric, sm_volumetric, 
         `CO2 flux (mgC/hr)` = CO2_flux_mgC_hr, 
         `CH4 flux (mgC/hr)` = CH4_flux_mgC_hr, 
         `cumCO2 flux (mgC/gSoil)` = cumCO2_flux_mgC_gSoil, 
         `cumCH4 flux (mgC/gSoil)` = cumCH4_flux_mgC_gSoil, label) %>% 
  arrange(Site, TREATMENT_PHASE, SampleID, inctime_hours) %>% 
  mutate(cutoff = if_else(TREATMENT_PHASE %in% c("SATURATION_IMMEDIATE", 
                                                 "SATURATION_SATURATION", 
                                                 "FIELD_MOIST_SATURATION", 
                                                 "DROUGHT_SATURATION"), 7.5, 999),
         cutoff = if_else(TREATMENT_PHASE %in% c("DROUGHT_INCUBATION", 
                                                 "FIELD_MOIST_INCUBATION", 
                                                 "SATURATION_INCUBATION"), 600, cutoff)) %>% 
  filter(inctime_hours <= cutoff) ->
  fluxdata_cutoff
save_data(fluxdata_cutoff)

p_collar_co2 <- ggplot(fluxdata_cutoff, aes(inctime_hours, `cumCO2 flux (mgC/gSoil)`, group = SampleID, color = Site)) + 
  geom_point() + geom_line() + geom_text(aes(x = inctime_hours * 1.1, label = label), size = 3) +
  ggtitle("Cumulative CO2 emissions by core")
p_collar_co2_1 <- p_collar_co2 +
  facet_wrap(~TREATMENT_PHASE, scales = "free") + theme(strip.text = element_text(size = 7))
print(p_collar_co2_1)
save_plot("cumulative_co2_by_core", width = 8, height = 6)
p_collar_co2_2 <- p_collar_co2 +
  facet_wrap(~TREATMENT_PHASE) + theme(strip.text = element_text(size = 7))
print(p_collar_co2_2)
save_plot("cumulative_co2_by_core_samescale", width = 8, height = 6)

p_collar_ch4 <- ggplot(fluxdata_cutoff, aes(inctime_hours, `cumCH4 flux (mgC/gSoil)`, group = SampleID, color = Site)) + 
  geom_point() + geom_line() + geom_text(aes(x = inctime_hours * 1.1, label = label), size = 3) +
  ggtitle("Cumulative CH4 emissions by core")
p_collar_ch4_1 <- p_collar_ch4 +
  facet_wrap(~TREATMENT_PHASE, scales = "free") + theme(strip.text = element_text(size = 7))
print(p_collar_ch4_1)
save_plot("cumulative_ch4_by_core", width = 8, height = 6)
p_collar_ch4_2 <- p_collar_ch4 +
  facet_wrap(~TREATMENT_PHASE) + theme(strip.text = element_text(size = 7))
print(p_collar_ch4_2)
save_plot("cumulative_ch4_by_core_samescale", width = 8, height = 6)

# Plot final (end of incubation) totals
fluxdata_cutoff %>% 
  group_by(Site, TREATMENT_PHASE, SampleID) %>% 
  summarise(inctime_hours = last(inctime_hours),
            `cumCO2 flux (mgC/gSoil)` = last(`cumCO2 flux (mgC/gSoil)`),
            `cumCH4 flux (mgC/gSoil)` = last(`cumCH4 flux (mgC/gSoil)`)) ->
  fluxdata_final
save_data(fluxdata_final)

fluxdata_final %>% 
  group_by(Site, TREATMENT_PHASE) %>%  # should be already grouped...
  summarise(`cumCO2 flux (mgC/gSoil) s.d.` = sd(`cumCO2 flux (mgC/gSoil)`),
            `cumCO2 flux (mgC/gSoil)` = mean(`cumCO2 flux (mgC/gSoil)`),
            `cumCH4 flux (mgC/gSoil) s.d.` = sd(`cumCH4 flux (mgC/gSoil)`),
            `cumCH4 flux (mgC/gSoil)` = mean(`cumCH4 flux (mgC/gSoil)`)) ->
  fluxdata_final_by_site
save_data(fluxdata_final_by_site)


p_final_co2 <- ggplot(fluxdata_final, aes(Site, `cumCO2 flux (mgC/gSoil)`, color = Site)) +
  geom_boxplot() + geom_point() +
  facet_wrap(~TREATMENT_PHASE, scales = "free") + theme(strip.text = element_text(size = 7))
print(p_final_co2)
save_plot("final_co2_inc", width = 8, height = 6)
p_final_ch4 <- ggplot(fluxdata_final, aes(Site, `cumCH4 flux (mgC/gSoil)`, color = Site)) +
  geom_boxplot() + geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~TREATMENT_PHASE, scales = "free") + theme(strip.text = element_text(size = 7))
print(p_final_ch4)
save_plot("final_ch4_inc", width = 8, height = 6)


# -----------------------------------------------------------------------------
# Done! 

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
