# Process Picarro data for the 3soils lab experiment
# This script reads all available Picarro outputs in `data/picarro/`,
# concatenating and writing to an `outputs/rawdata.csv.gz` file.
# Ben Bond-Lamberty July 2015

source("0-functions.R")

SCRIPTNAME  	<- "1-googlesheets.R"
PROBLEM       <- FALSE

library(googlesheets) # 0.2.1

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)

printlog("Welcome to", SCRIPTNAME)

# Register file - need to be authenticated to Google
# Not storing an OAuth token here
# File is "CPCRW Soil Cores Key, Weights"

# Note that "registration by key is the safest, long-run strategy"
# https://cran.r-project.org/web/packages/googlesheets/vignettes/basic-usage.html
KEY <- "1vYmwhZaymQNjkl-IPZPxPe1-aX8CQWyUQVX7Ub8yZSs"
gap <- gs_key(KEY)

print(gap)

sheet_key <- "Key"
printlog("Downloading", sheet_key, "to", KEY_FILE)

gap %>%
  gs_download(ws = sheet_key, to = KEY_FILE, overwrite = TRUE)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
