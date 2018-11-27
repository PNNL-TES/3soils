# Download the valve_map and sampleID_key files from Google Sheets
# Ben Bond-Lamberty November 2017

source("0-functions.R")

SCRIPTNAME  	<- "1-googlesheets.R"
PROBLEM       <- FALSE

library(googlesheets) # 0.2.2
library(googledrive)

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)

printlog("Welcome to", SCRIPTNAME)

# Notes   November 27, 2018
# The do_download below downloads all Picarro files but note it's starting in my GDrive folder
# The code below that is broken--they shifted which files were using for key information
# (per Peyton email 20 November 2018)
# For now I've just downloaded the relevant file (3Soils_CPCRW_SR_DWP_PicarroLog) directly


printlog("Getting list of Picarro files...")
folder_id <- as_id()
filelist <- drive_ls("~/3Soils_Picarro/")

do_download <- function(filelist, depth = 1) {
  for(f in seq_len(nrow(filelist))) {
    cat(depth, filelist$name[f], "\t")
    if(is_folder(filelist[f,])) {
      cat("folder\n")
      do_download(drive_ls(as_id(filelist$id[f])), depth + 1)
    } else {
      localfile <- file.path("data/picarro/", filelist$name[f])
      if(file.exists(localfile)) {
        cat("already exists\n")
      } else {
        cat("downloading", filelist$name[f], "\n")
        drive_download(filelist[f,], path = localfile)
      }
    }
  }
}

suppressWarnings(do_download(filelist))



# Register file - need to be authenticated to Google
# Not storing an OAuth token here
# File is "3Soils_CPCRW_SR_DWP_PicarroLog"

# Note that "registration by key is the safest, long-run strategy"
# https://cran.r-project.org/web/packages/googlesheets/vignettes/basic-usage.html
KEY <- "1wsI3tldbhhMDSDoRejmS2jE2c9U-75sht_Dyum3QUBY"
valvemap <- gs_key(KEY)
print(valvemap)

vmdata <- list()
tf <- tempfile(fileext = ".csv")
sheet <- "valve_map"
printlog("Downloading", sheet)
valvemap %>%
  gs_download(ws = sheet, to = tf, overwrite = TRUE)
vmdata <- read_csv(tf, na = c("", "NA", "n/a", "na"))

sheet <- "sampleID_key"
printlog("Downloading", sheet)
valvemap %>%
  gs_download(ws = sheet, to = tf, overwrite = TRUE)
keydata <- read_csv(tf, na = c("", "NA", "n/a", "na"))

save_data(vmdata, fn = VALVEMAP_FILE, scriptfolder = FALSE)
save_data(keydata, fn = KEY_FILE, scriptfolder = FALSE)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
