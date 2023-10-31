# getNTpredFlywire.r
#
# Function that generates text files of neurotransmitter summary, given csv file
#  of IDs (e.g. processedNeurons.csv from genPresynCSVsFlyWire.r)
# 

getNTpredFlywire <- function(neuronsCSVpath, cleftThresh, outPath) {
  # read in csv file
  csvData <- read.csv(neuronsCSVpath,header=TRUE, colClasses='character')
  
  # set wd
  setwd(outPath)
  
  # newIDs <- fafbseg::flywire_latestid(csvData$V2)
  
  # loop through all neurons
  for (i in csvData$x) {
    # get NT prediction
    thisNTpred <-flywire_ntpred(i,cleft.threshold = cleftThresh)
    # file name is neuron ID
    thisFileName <- paste(i,".txt",sep="")
    
    # write print summary to text file
    sink(thisFileName)
    print(thisNTpred)
    sink()
  }
}