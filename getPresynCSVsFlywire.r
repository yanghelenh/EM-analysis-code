# getPresynCSVsFlywire.r
#
# Function that generates CSV files of input neurons to neurons of interest and
#  the number of synapses they make
# 

getPresynCSVsFlywire <- function(rootid, thresh, cleftThresh, csvFilePath) {
  thisID <- fafbseg::flywire_latestid(rootid)
  
  theseInputs <- fafbseg::flywire_partner_summary(rootid, partners = "inputs", 
                                                  threshold = thresh, 
                                                  cleft.threshold = cleftThresh)
  
  # convert to data frame
  df <- data.frame(pre_id = theseInputs$pre_id, weight = theseInputs$weight)
  
  # write to csv file
  write.table(df, csvFilePath, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
}