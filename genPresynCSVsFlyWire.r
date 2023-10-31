# genPresynCSVsFlyWire.r
# 
# Script that generates CSV files of presynaptic neurons given starting FlyWire ID
#  Will return 2 layers of inputs: direct and their inputs
#  Writes CSV files, named by current FlyWire ID, into specified directory
#  Set threshold of minimum number of synapses

# make sure fafbseg library is loaded
# library(fafbseg)

# Directory to save CSV files
setwd("/Users/hyang/Dropbox (Personal)/Wilson Lab/EM reconstruction/Brain/FlyWirePresynCSVs/230903")

# Threshold for minimum number of synapses
MIN_SYN <- 5

# Cleft threshold
CLEFT_THRESH <- 80

# starting neuron
# DNa02 ID - note, check if ID is up to date
dna02lID <- "720575940604737708"
dna02rID <- "720575940629327659"
# DNg13 ID - note, check if ID is up to date
# dng13lID <- "720575940633270497"
dng13lID <- "720575940606112940"
dng13rID <- "720575940616471052"

#startNeuron <- dna02lID
# startNeuron <- dna02rID
# startNeuron <- dng13lID
startNeuron <- dng13rID

# get postsynaptic targets of start neuron as data.frame
primTarg <- flywire_partner_summary(startNeuron, "inputs", MIN_SYN, TRUE, CLEFT_THRESH)

# save primary targets as CSV file
write.csv(primTarg, paste(startNeuron,".csv",sep=""),quote = FALSE)

# running list of neurons processed
# processedNeurons <- startNeuron
# processedNeurons <- NULL

# running list of neurons that throw errors
errorNeurons <- NULL

# loop over list of postsynaptic targets of start neuron
for (i in primTarg$pre_id) {
  # check if in processedNeurons list
  if (!(i %in% processedNeurons)) {
    
    # try catch b/c getting errors for some IDs
    secTarg <- tryCatch(
      {
        flywire_partner_summary(i, "inputs", MIN_SYN, TRUE, CLEFT_THRESH)
      },
      error=function(cond) {
        # print ID that was invalid
        print(i)
        # add this ID to list of neurons that throw errors
        errorNeurons[NROW(errorNeurons) + 1] = i
        # print error message
        message(cond)
        # return NULL for secTarg value
        return(NULL)
      }
    )
    
    # only do this if successfully generated list of postsynaptic targets
    if (!(is.null(secTarg))) {
      # save these targets into CSV file
      write.csv(secTarg, paste(i,".csv",sep=""),quote = FALSE)
      
      # add this neuron to processedNeurons list
      processedNeurons[NROW(processedNeurons) + 1] = i
      
      # loop over list of postsynaptic targets of secondary neuron
      for (j in secTarg$pre_id) {
        if (!(j %in% processedNeurons)) {
          # try catch b/c some IDs throwing errors
          terTarg <- tryCatch(
            {
              flywire_partner_summary(j, "inputs", MIN_SYN, TRUE, CLEFT_THRESH)
            },
            error=function(cond) {
              # print ID that was invalid
              print(j)
              # add this ID to list of neurons that throw errors
              errorNeurons[NROW(errorNeurons) + 1] = j
              # add this ID to list of processed neurons (so don't try again)
              processedNeurons[NROW(processedNeurons)+1] = j
              # print error message
              message(cond)
              # return null for terTarg value
              return(NULL)
            }
          )

          # only if successfully generated list of tertiary targets
          if (!(is.null(terTarg))) {
            # save these targets into CSV file
            write.csv(terTarg, paste(j,".csv",sep=""),quote = FALSE)

            # add this ID to list of processed neurons
            processedNeurons[NROW(processedNeurons)+1] = j
          }
        }
      }
    }
  }
}

# write processedNeurons to csv file
write.csv(processedNeurons, "processedNeurons.csv", quote = FALSE)

# write errorNeurons to csv file
write.csv(errorNeurons,"errorNeurons.csv", quote = FALSE)
