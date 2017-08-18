# Process Picarro data for the 3soils lab experiment
# This script reads all available Picarro outputs in `data/picarro/`,
# concatenating and writing to an `outputs/rawdata.csv.gz` file.
# Ben Bond-Lamberty July 2015

source("0-functions.R")

SCRIPTNAME  	<- "1-googlesheets.R"
PROBLEM       <- FALSE

library(googlesheets) # 0.2.2

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)

printlog("Welcome to", SCRIPTNAME)

# Register file - need to be authenticated to Google
# Not storing an OAuth token here
# File is "DWP_Picarro_data_collection"

# Note that "registration by key is the safest, long-run strategy"
# https://cran.r-project.org/web/packages/googlesheets/vignettes/basic-usage.html
KEY <- "15iM3XSFchLNm010D9n9pTpFapw0hX7unVqh35rtqpaI"
gap <- gs_key(KEY)

print(gap)
keydata <- list()
sheets <- c("Drought", "Field_moist", "Saturation_II")
for(sheet_key in sheets) {
  
  printlog("Downloading", sheet_key, "to", KEY_FILE)
  gap %>%
    gs_download(ws = sheet_key, to = KEY_FILE, overwrite = TRUE)
  keydata[[sheet_key]] <- read_csv(KEY_FILE, na = c("", "NA", "n/a"),
                                   col_types = cols(Headspace_height = col_double(),
                                                    Weight_upon_saturation = col_double()))
  
}

keydata <- bind_rows(keydata)
save_data(keydata, fn = KEY_FILE, scriptfolder = FALSE)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
