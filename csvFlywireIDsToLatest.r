# Function that takes in full path to CSV file of FlyWire IDs and writes to a CSV
#  file in the same folder the latest FlyWire IDs
# Expects no header for CSV file

csvFlywireIDsToLatest <- function(csvFilePath) {
  # read in CSV file
  csvData <- read.csv(csvFilePath,header=FALSE, colClasses='character')
  
  # for each of these IDs, get latest IDs
  newIDs <- flywire_latestid(csvData$V1)
  
  # write these new IDs to csv file
  newCSVFilePath <- paste(substr(csvFilePath,1,nchar(csvFilePath)-4),"_newIDs.csv",sep="")
  write.csv(newIDs, newCSVFilePath, quote = FALSE, row.names = FALSE)
}