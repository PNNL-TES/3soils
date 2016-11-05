# Process Picarro data for 3soils lab experiment
# This workhorse script summarizes individual (raw) Picarro observations to 
# summaries of "samples" (groups of consecutive observations made from a given 
# core at a point in time). It computes gas concentration changes, performs 
# some QC, merges the Picarro data with valve map and other ancillary data,
# and writes SUMMARYDATA_FILE.
# 
# Ben Bond-Lamberty November 2016

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
read_csv(RAWDATA_FILE, col_types = "ccddddiiiddddddddc") %>%
  # immediately discard columns we don't need
  select(DATE, TIME, MPVPosition, CH4_dry, CO2_dry, h2o_reported) %>%
  # Fractional solenoid values mean that the analyzer was shifting
  # between two samples. Discard these.
  filter(MPVPosition == trunc(MPVPosition)) %>%
  # Convert date/time to POSIXct
  mutate(DATETIME = ymd_hms(paste(DATE, TIME))) %>%
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
printlog("Assigning sample numbers...")
oldsampleflag <- with(rawdata, c(FALSE, MPVPosition[-length(MPVPosition)] == MPVPosition[-1]))
rawdata$samplenum <- cumsum(!oldsampleflag)

printlog("Computing elapsed seconds...")
rawdata %>%
  group_by(samplenum) %>%
  mutate(elapsed_seconds = difftime(DATETIME, min(DATETIME), units = "secs")) ->
  rawdata_samples

printlog("Identifying and removing ambient samples...")
rawdata_samples %>%
  group_by(samplenum) %>%
  summarise(max_seconds = max(as.numeric(elapsed_seconds))) %>%
  left_join(rawdata_samples, by = "samplenum") %>%
  filter(max_seconds <= MAX_MEASUREMENT_TIME) %>%
  select(-max_seconds) %>%
  print_dims("rawdata_samples") ->
  rawdata_samples

printlog("Visualizing...")
p <- ggplot(rawdata_samples, aes(x = elapsed_seconds, color = DATE, group = samplenum)) +
  facet_wrap(~MPVPosition) + 
  ggtitle("CO2 concentration by date and valve")
print(p + geom_line(aes(y = CO2_dry)))
save_plot("co2_by_valve")
print(p + geom_line(aes(y = CH4_dry)))
save_plot("ch4_by_valve")

# -----------------------------------------------------------------------------
# Load and QC the key and valvemap data

printlog("Loading key data...")
read_csv(KEY_FILE, col_types = "ciiccc", na = c("NA", "", "na")) %>%
  filter(!is.na(`Core_#`)) %>%
  print_dims("keydata") ->
  keydata

qc_keydata(keydata)

printlog("Loading valve map data...")
read_csv(VALVEMAP_FILE, col_types = "cccicidcc", na = c("NA", "", "na")) %>%
  mutate(Date = mdy(Date)) %>%
  mutate(Time = lubridate::hm(Start_time_PDT)) %>%
  arrange(Date, Time) %>%
  print_dims("valvemap") ->
  valvemap

qc_valvemap(valvemap)


# -----------------------------------------------------------------------------
# Compute concentration changes and match the Picarro data with valvemap data

printlog( "Computing summary statistics for each sample..." )

# We want to apply different criteria here, so three different pipelines
# to compute the min and max gas concentrations
rawdata_samples %>%
  filter(elapsed_seconds <= MAX_MINCONC_TIME) %>%
  group_by(samplenum) %>%
  summarise(min_CO2 = min(CO2_dry),
            min_CO2_time = nth(elapsed_seconds, which.min(CO2_dry)),
            min_CH4 = min(CH4_dry),
            min_CH4_time = nth(elapsed_seconds, which.min(CH4_dry))) ->
  summarydata_min

# Now we want to look for the max concentration AFTER the minimum
rawdata_samples %>%
  left_join(summarydata_min, by = "samplenum") ->
  rawdata_temp

rawdata_temp %>%
  filter(elapsed_seconds > min_CO2_time & elapsed_seconds < MAX_MAXCONC_TIME) %>%
  group_by(samplenum) %>%
  summarise(max_CO2 = max(CO2_dry),
            max_CO2_time = nth(elapsed_seconds, which.max(CO2_dry))) ->
  summarydata_maxCO2
rawdata_temp %>%
  filter(elapsed_seconds > min_CH4_time & elapsed_seconds < MAX_MAXCONC_TIME) %>%
  group_by(samplenum) %>%
  summarise(max_CH4 = max(CH4_dry),
            max_CH4_time = nth(elapsed_seconds, which.max(CH4_dry))) ->
  summarydata_maxCH4

# Final pipeline: misc other data, and match up with valve map entries
rawdata_samples %>%
  group_by(samplenum) %>%
  summarise(
    DATETIME = mean(DATETIME),
    N = n(),
    MPVPosition	= mean(MPVPosition),
    h2o_reported = mean(h2o_reported)) ->
  summarydata_other

# Merge pieces together to form final summary data set
printlog("Removing N=1 and MPVPosition=0 data, and merging...")
summarydata_other %>%
  filter(N > 1) %>% # N=1 observations are...? Picarro quirk
  filter(MPVPosition > 0) %>% # ? Picarro quirk
  left_join(summarydata_min, by = "samplenum") %>%
  left_join(summarydata_maxCO2, by = "samplenum") %>% 
  left_join(summarydata_maxCH4, by = "samplenum") ->
  summarydata

printlog("Merging Picarro and mapping data...")
summarydata %>%
  left_join(valvemap, by = c("MPVPosition" = "Valve")) ->
  summarydata

printlog("Reading and merging key data...")
left_join(summarydata, keydata, by = "Core_#") -> 
  summarydata

# The treatment codes should match. Warn if not, then proceed
# with the treatment defined in the key data
if(nrow(subset(summarydata, tolower(Treatment.x) != tolower(Treatment.y)))) {
  printlog("WARNING - some treatments in the valve data don't match those given in key data")
  printlog("This probably doesn't matter but may indicate a problem")
}
summarydata %>%
  select(-Treatment.x) %>%
  rename(Treatment = Treatment.y) ->
  summarydata

printlog("Computing per-second rates...")
summarydata %>%
  mutate(CO2_ppm_s = (max_CO2 - min_CO2) / (max_CO2_time - min_CO2_time),
         CH4_ppb_s = (max_CH4 - min_CH4) / (max_CH4_time - min_CH4_time),
         inctime_days = 1 + as.numeric(difftime(DATETIME, min(DATETIME), units = "days"))) ->
  summarydata

# -----------------------------------------------------------------------------
# Done! 

save_data(summarydata, fn = SUMMARYDATA_FILE, scriptfolder = FALSE)
save_data(rawdata_samples, fn = RAWDATA_SAMPLES_FILE, scriptfolder = FALSE)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
