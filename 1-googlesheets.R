
library(googlesheets)

# Register file - need to be authenticated to Google
# Not storing an OAuth token here
#gap <- gs_title("CPCRW Soil Cores Key, Weights-BEN")

# Note that "registration by key is the safest, long-run strategy"
# https://cran.r-project.org/web/packages/googlesheets/vignettes/basic-usage.html
KEY <- "1vYmwhZaymQNjkl-IPZPxPe1-aX8CQWyUQVX7Ub8yZSs"
gap <- gs_key(KEY)

print(gap)

gap %>%
  gs_read(ws = "Key") -> 
  keydata

print(keydata)
