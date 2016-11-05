# Process Picarro data for 3soils lab experiment
# This workhorse script summarizes individual (raw) Picarro observations to 
# summaries of "samples" (groups of consecutive observations made from a given 
# core at a point in time). It computes gas concentration changes, performs 
# some QC, merges the Picarro data with valve map and other ancillary data,
# and writes `outputs/summarydata.csv`.
# 
# Ben Bond-Lamberty November 2016

source("0-functions.R")

SCRIPTNAME  	<- "3-summarize.R"
PROBLEM       <- FALSE


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
  arrange(DATETIME) ->
  rawdata
print_dims(rawdata)
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

printlog("Visualizing...")
p <- ggplot(rawdata_samples, aes(x = elapsed_seconds, color = DATE, group = samplenum)) +
  xlim(c(0,150)) + 
  facet_wrap(~MPVPosition) + 
  ggtitle("CO2 concentration by date and valve")
print(p + geom_line(aes(y = CO2_dry)))
save_plot("co2_by_valve")
print(p + geom_line(aes(y = CH4_dry)))
save_plot("ch4_by_valve")

stop()

# -----------------------------------------------------------------------------
# Load and QC the key data

printlog("Loading valve map (key) data...")
read_csv(KEY_FILE, col_types = "ciiccc", na = c("NA", "", "na")) %>%
  filter(!is.na(`Core_#`)) ->
  keydata
# printlog( "Converting date/time info to POSIXct..." )
# keydata$StartDateTime <- mdy_hm(paste(keydata$Date, keydata$Time_set_start_UTC))
# keydata <- arrange(keydata, StartDateTime)
# keydata$valvemaprow <- seq_len(nrow(keydata))
qc_valvemap(keydata)

# -----------------------------------------------------------------------------
# Compute concentration changes and match the Picarro data with valvemap data

# Function to match up Picarro data with mapping file data
# This is done by date and valve number (see plot saved above)
matchfun <- function(DATETIME, MPVPosition) {
  DATETIME <- as.POSIXct(DATETIME, origin = lubridate::origin, tz="UTC")
  rowmatches <- which(DATETIME >= valvemap$StartDateTime & 
                        yday(DATETIME) == yday(valvemap$StartDateTime) &
                        MPVPosition == valvemap$MPVPosition)
  if(length(rowmatches) == 0) rowmatches <- NA
  max(rowmatches)  # return latest time match
}

printlog( "Computing summary statistics for each sample..." )

# We want to apply different criteria here, so three different pipelines
# to compute the min and max gas concentrations
summarydata_min <- rawdata_samples %>%
  filter(elapsed_seconds <= MAX_MINCONC_TIME) %>%
  group_by(samplenum) %>%
  summarise(min_CO2 = min(CO2_dry),
            min_CO2_time = nth(elapsed_seconds, which.min(CO2_dry)),
            min_CH4 = min(CH4_dry),
            min_CH4_time = nth(elapsed_seconds, which.min(CH4_dry)))

# Now we want to look for the max concentration AFTER the minimum
rawdata_temp <- rawdata_samples %>%
  left_join(summarydata_min, by = "samplenum") 

summarydata_maxCO2 <- rawdata_temp %>%
  filter(elapsed_seconds > min_CO2_time & elapsed_seconds < MAX_MAXCONC_TIME) %>%
  summarise(max_CO2 = max(CO2_dry),
            max_CO2_time = nth(elapsed_seconds, which.max(CO2_dry))
  )
summarydata_maxCH4 <- rawdata_temp %>%
  filter(elapsed_seconds > min_CH4_time & elapsed_seconds < MAX_MAXCONC_TIME) %>%
  summarise(max_CH4 = max(CH4_dry),
            max_CH4_time = nth(elapsed_seconds, which.max(CH4_dry)))

# Final pipeline: misc other data, and match up with valve map entries
summarydata_other <- rawdata_samples %>%
  group_by(samplenum) %>%
  summarise(
    DATETIME = mean(DATETIME),
    N = n(),
    MPVPosition	= mean(MPVPosition),
    h2o_reported = mean(h2o_reported),
    valvemaprow = matchfun(DATETIME, MPVPosition)) 

# Merge pieces together to form final summary data set
printlog("Removing N=1 and MPVPosition=0 data, and merging...")
summarydata <- summarydata_other %>%
  filter(N > 1) %>% # N=1 observations are...? Picarro quirk
  filter(MPVPosition > 0) %>% # ? Picarro quirk
  left_join(summarydata_min, by = "samplenum") %>%
  left_join(summarydata_maxCO2, by = "samplenum") %>% 
  left_join(summarydata_maxCH4, by = "samplenum")

printlog("Merging Picarro and mapping data...")
summarydata <- left_join(summarydata, valvemap, by = c("MPVPosition", "valvemaprow"), all.x=TRUE)

printlog("Reading and merging treatment data...")
trtdata <- read_csv(TREATMENTS, skip = 1)
summarydata <- left_join(summarydata, trtdata, by = "Core")

printlog("Computing per-second rates...")
summarydata <- summarydata %>%
  mutate(CO2_ppm_s = (max_CO2 - min_CO2) / (max_CO2_time - min_CO2_time),
         CH4_ppb_s = (max_CH4 - min_CH4) / (max_CH4_time - min_CH4_time),
         inctime_days = 1 + as.numeric(difftime(DATETIME, min(DATETIME), units = "days")))

printlog("Saving a comparison of MPVPosition sequence in Picarro data and valvemap")
checkdata <- select(summarydata, DATETIME, MPVPosition)
checkdata$sequence <- seq_len(nrow(checkdata))
vdata <- data.frame(sequence = seq_len(nrow(valvemap)),
                    DATETIME_valvemap = paste(valvemap$Date, valvemap$Time_set_start_UTC),
                    MPVPosition_valvemap = valvemap$MPVPosition)
MPVPosition_checkdata <- left_join(vdata, checkdata, by = "sequence")
save_data(MPVPosition_checkdata)

# -----------------------------------------------------------------------------
# Done! 

save_data(summarydata, fn = SUMMARYDATA_FILE, scriptfolder = FALSE)
save_data(rawdata_samples, fn = RAWDATA_SAMPLES_FILE, scriptfolder = FALSE)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
