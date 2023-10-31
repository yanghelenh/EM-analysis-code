# csvFancGetNumInputSyn.r
#
# Function that takes in CSV file of FANC IDs and returns the total number of
#  synapses that neuron receives
# Writes output CSV file with FANC IDs and number of synapses (name is same 
#  with _synNum appended)

csvFancGetNumInputSyn <- function(csvFilePath) {
  # read in CSV file
  csvData <- read.csv(csvFilePath,header=FALSE, colClasses='character')
  
  # for each of these IDs, get latest IDs
  newIDs <- fanc_latestid(csvData$V1)
  
  # initialize
  synNum <- NULL
  
  for (i in newIDs) {
    thisInputs <- fanc_partner_summary(i,"inputs")
    synNum[NROW(synNum)+1] = sum(thisInputs$weight)
  }
  
  # convert to data frame (use old FANC IDs, not latest, for matching)
  df <- data.frame(id = csvData$V1, numSyn = synNum)
  
  # write these IDs and synapse counts to CSV
  newCSVFilePath <- paste(substr(csvFilePath,1,nchar(csvFilePath)-4),"_synNum.csv",sep="")
  write.table(df, newCSVFilePath, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
}