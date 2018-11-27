# Process Picarro data for the 3soils lab experiment
# This script reads all available Picarro outputs in `data/picarro/`,
# concatenating and writing to an `outputs/rawdata.csv.gz` file.
# Ben Bond-Lamberty July 2015

source("0-functions.R")

SCRIPTNAME  	<- "2-data.R"
PROBLEM       <- FALSE

PICARRO_DATA_DIR      <- file.path(DATA_DIR, "picarro/")

# -----------------------------------------------------------------------------
# read a single output file, returning data frame
read_outputfile <- function(fqfn) {
  printlog("Reading", fqfn)
  stopifnot(file.exists(fqfn))
  
  f <- fqfn
  if(grepl(".gz$", fqfn)) {
    f <- gzfile(fqfn)
  } else if(grepl(".zip$", fqfn)) {
    f <- unz(fqfn)
  }
  
  d <- read.table(f, header = TRUE)
  
  print_dims(d)
  
  # Add ancillary data
  d$file <- basename(fqfn)
  #  d$dir <- dirname(fqfn)
  
  return(d)
} # read_outputfile

# -----------------------------------------------------------------------------
# scan a directory and process all files in it, returning tempfile names
process_directory <- function(input_path) {
  printlog("process_directory", input_path)
  filelist <- list.files(path = input_path, 
                         pattern = "dat$|dat.gz$|dat.zip$", 
                         recursive = TRUE,
                         full.names = TRUE)
  filedata <- list()
  for(f in filelist) {
    printlog("Reading", f)
    tibble::as_tibble(read.table(f, header = TRUE, stringsAsFactors = FALSE)) %>%
      # select only the columns we need, and discard any fractional valve numbers
      select(DATE, TIME, ALARM_STATUS, MPVPosition, CH4_dry, CO2_dry, h2o_reported) %>%
      filter(MPVPosition == floor(MPVPosition)) ->
      filedata[[basename(f)]]
  }
  filedata %>%
    bind_rows(.id = "filename")
}

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)

printlog("Welcome to", SCRIPTNAME)
printlog("Data directory is", PICARRO_DATA_DIR)

rawdata <- process_directory(PICARRO_DATA_DIR)

printlog("Writing output file...")
save_data(rawdata, fn = RAWDATA_FILE, scriptfolder = FALSE)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
