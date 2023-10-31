# Function that takes in full path to CSV file of FANC IDs and writes to a 
#  specified folder all the change logs for those IDs
# Expects no header for CSV file

csvFancChangeLog <- function(csvFilePath, outPath) {
  # read in CSV file
  csvData <- read.csv(csvFilePath,header=FALSE, colClasses='character')
  
  # set wd
  setwd(outPath)
  
  # for each of these IDs, get change log
  for (i in csvData$V1) {
    thisChangeLog <- fanc_change_log(i)
    # file name is ID
    thisFileName <- paste(i,".csv",sep="")
    # write as csv and save to output folder
    write.csv(thisChangeLog, thisFileName, quote = FALSE, row.names = FALSE)
  }
}