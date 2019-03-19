# Process Picarro data for 3soils lab experiment
# This workhorse script summarizes individual (raw) Picarro observations to 
# summaries of "samples" (groups of consecutive observations made from a given 
# core at a point in time). It computes gas concentration changes, performs 
# some QC, merges the Picarro data with valve map and other ancillary data,
# and writes SUMMARYDATA_FILE.
# 
# Ben Bond-Lamberty November 2017

source("0-functions.R")

SCRIPTNAME  	<- "3-summarize.R"
PROBLEM       <- FALSE


qc_keydata <- function(keydata) {
  dupes <- which(duplicated(keydata$`Field_#`) |
                   duplicated(keydata$`Core_#`))
  if(length(dupes)) {
    flaglog("There are duplicate core and/or field numbers in the key data:")
    print(keydata[dupes,])
    PROBLEM <<- TRUE
  }
}

qc_valvemap <- function(valvemap) {
  valvemap %>% 
    filter(!is.na(Valve)) %>%
    group_by(Date, Valve) %>% 
    summarise(n = n()) %>%
    summarise(dupeflag = any(n > 1)) %>%
    filter(dupeflag) ->
    dupes
  
  if(nrow(dupes)) {
    flaglog("There are duplicate core number in the valvemap data on these dates:")
    print(dupes)
    PROBLEM <<- TRUE
  }
}


# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE) # open log

printlog("Welcome to", SCRIPTNAME)

printlog("Reading in raw data...")
read_csv(RAWDATA_FILE, col_types = cols(filename = col_character(),
                                        DATE = col_date(format = ""),
                                        TIME = col_time(format = ""),
                                        ALARM_STATUS = col_integer(),
                                        MPVPosition = col_integer(),
                                        CH4_dry = col_double(),
                                        CO2_dry = col_double(),
                                        h2o_reported = col_double()
)) %>%
  # Convert date/time to POSIXct
  mutate(DATETIME = ymd_hms(paste(DATE, TIME))) %>%
  select(-DATE, -TIME) %>%
  arrange(DATETIME) %>%
  print_dims("rawdata") ->
  rawdata
print(summary(rawdata))
printlog("First timestamp:")
print(min(rawdata$DATETIME))
printlog("Last timestamp:")
print(max(rawdata$DATETIME))

# -----------------------------------------------------------------------------
# Prep work: data cleaning, dates, sample numbers, elapsed time

# Assign a different sample number to each sample group 
# (we know we're on a new sample when MPVPosition changes)
printlog("Assigning sample numbers and computing elapsed time...")
rawdata %>%
  mutate(newsample = MPVPosition != lag(MPVPosition)) %>%
  replace_na(list(newsample = FALSE)) %>% 
  mutate(samplenum = cumsum(newsample)) %>%
  select(-newsample) %>%
  group_by(samplenum) %>%
  mutate(elapsed_seconds = as.double(difftime(DATETIME, min(DATETIME), units = "secs"))) ->
  rawdata_samples

printlog("Removing ambient and very long samples...")
AMBIENT_VALVE <- 16
rawdata_samples %>%
  filter(MPVPosition != AMBIENT_VALVE, 
         elapsed_seconds < MAX_MEASUREMENT_TIME) %>% 
  group_by(samplenum) %>%
  print_dims("rawdata_samples") ->
  rawdata_samples

printlog("Visualizing...")
rawdata_samples %>% 
  group_by(MPVPosition, elapsed_seconds) %>% 
  summarise_at(c("CO2_dry", "CH4_dry"), funs(min, max, mean, sd)) -> 
  rawdata_summary

ggplot(rawdata_summary, aes(elapsed_seconds, CO2_dry_mean)) +
  geom_line() + 
  geom_ribbon(aes(ymin = CO2_dry_mean - CO2_dry_sd,
                                ymax = CO2_dry_mean + CO2_dry_sd), alpha = I(0.5)) +
  facet_wrap(~MPVPosition) + 
  ggtitle("CO2 concentration by valve")
save_plot("raw_co2_by_valve", ptype = ".png")

ggplot(rawdata_summary, aes(elapsed_seconds, CH4_dry_mean)) +
  geom_line() + 
  geom_ribbon(aes(ymin = CH4_dry_mean - CH4_dry_sd,
                  ymax = CH4_dry_mean + CH4_dry_sd), alpha = I(0.5)) +
  facet_wrap(~MPVPosition) + 
  ggtitle("CH4 concentration by valve") +
  coord_cartesian(ylim = c(0, 5))
save_plot("raw_ch4_by_valve", ptype = ".png")

# -----------------------------------------------------------------------------
# Load and QC the key and valvemap data

# The 'valvemap' data maps Picarro valve numbers to sample IDs
printlog(SEPARATOR)
printlog("Reading valve and core mapping data...")
read_csv(VALVEMAP_FILE, na = c("NA", "#VALUE!", "NO VALVE", "NO WEIGHT", "", "DATA?")) %>% 
  mutate(rownum = row_number()) %>% 
  filter(!is.na(SampleID)) %>%
  mutate(Picarro_start = mdy_hm(Start_Date_Time, tz = "America/Los_Angeles"),
         Picarro_stop = mdy_hm(Stop_Date_Time, tz = "America/Los_Angeles"),
         sequence_valve = as.numeric(sequence_valve)) %>% 
  select(rownum, SampleID, TREATMENT_PHASE, sequence_valve,
         Picarro_start,
         Picarro_stop, 
         Total_core_mass_pot_pie_pans_g, Additional_Wt_toRemove) %>% 
  arrange(Picarro_start) ->
  valvemap

# The `gs_key` file maps SampleID to (at the moment) core dry mass and pH
read_csv(KEY_FILE) %>% 
  select(SampleID, Site, Treatment, HeadSpace_Ht_cm, 
         soil_pH_water_airdried, 
         DryMass_SoilOnly_g, DryMass_NONsoil_ALL_g, VolumeSoil_cm3) %>% 
  right_join(valvemap, by = "SampleID") ->
  valvemap

# valvemap diagnostic
valvemap %>% 
  ggplot(aes(SampleID, DryMass_SoilOnly_g, color = Site)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))
save_plot("diag_DryMass_SoilOnly_g", width = 8, height = 4)


# -----------------------------------------------------------------------------
# Compute concentration changes and match the Picarro data with valvemap data

printlog( "Computing summary statistics for each sample..." )
rawdata_samples %>% 
  ungroup %>% 
  # we only look at first 45 seconds, after first few
  filter(elapsed_seconds <= MAX_MAXCONC_TIME,
         elapsed_seconds >= MIN_MEASUREMENT_TIME) %>% 
  # find max CO2 time for each sample
  group_by(samplenum) %>% 
  mutate(N = n(),
         max_co2_time = nth(elapsed_seconds, which.max(CO2_dry))) %>%
  # filter for at least 3 data points and for max CO2 time window
  filter(N >= 3, 
         elapsed_seconds <= max_co2_time) %>%
  # now compute the slope for each group
  # stupid to fit each model twice, and I could use do(), but this is easy...
  group_by(samplenum) %>%
  summarise(CO2_ppm_s = lm(CO2_dry ~ elapsed_seconds)$coefficients["elapsed_seconds"],
            CH4_ppb_s = lm(CH4_dry ~ elapsed_seconds)$coefficients["elapsed_seconds"],
            r2_CO2 = summary(lm(CO2_dry ~ elapsed_seconds))$adj.r.squared,
            r2_CH4 = summary(lm(CH4_dry ~ elapsed_seconds))$adj.r.squared,
            DATETIME = mean(DATETIME),
            N = n(),
            MPVPosition	= mean(MPVPosition),
            h2o_reported = mean(h2o_reported)) %>% 
  print_dims("summarydata") ->
  summarydata

printlog("Time zone conversion...")
summarydata <- mutate(summarydata, DATETIME = with_tz(DATETIME, "America/Los_Angeles"))

printlog("Joining key data and Picarro output...")
newdata <- list()
valvemap$picarro_records <- 0
#summarydata$SampleID <- NA_character_
for(i in seq_len(nrow(valvemap))) {
  d <- filter(summarydata, DATETIME >= valvemap$Picarro_start[i], DATETIME <= valvemap$Picarro_stop[i])
  valvemap$picarro_records[i] <- nrow(d)
  newdata[[i]] <- left_join(valvemap[i,], d, by = c("sequence_valve" = "MPVPosition"))
}
bind_rows(newdata) %>% 
  select(-Picarro_start, -Picarro_stop) ->
  summarydata_clean


summarydata$matched <- summarydata$samplenum %in% newdata$samplenum
summarydata %>%
  filter(!matched) %>% 
  save_data(fn = "outputs/unmatched_picarro_data.csv", scriptfolder = FALSE)

# -----------------------------------------------------------------------------
# Done! 

save_data(summarydata, fn = SUMMARYDATA_FILE, scriptfolder = FALSE)
save_data(summarydata_clean, fn = SUMMARYDATA_CLEAN_FILE, scriptfolder = FALSE)
save_data(rawdata_samples, fn = RAWDATA_SAMPLES_FILE, scriptfolder = FALSE)
save_data(valvemap, fn = "valvemap_combined.csv", scriptfolder = FALSE)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
