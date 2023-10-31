# genPostsynCSVsFANC.r
# 
# Script that generates CSV files of postsynaptic neurons given starting FANC ID
#  Will return 3 layers of targets: direct, their targets, and their targets
#  Writes CSV files, named by current FANC ID, into specified directory
#  Set threshold of minimum number of synapses

# make sure fancr library is loaded
# library(fancr)

# Directory to save CSV files
setwd("/Users/hyang/Dropbox (Personal)/Wilson Lab/EM reconstruction/VNC/PostsynCSVs/230903")

# Threshold for minimum number of synapses
MIN_SYN <- 5

# Cleft threshold
CLEFT_THRESH <- 50

# starting neuron
# DNa02 ID - note, check if ID is up to date
dna02ID <- "648518346478550356"
# DNg13 ID - note, check if ID is up to date
dng13ID <- "648518346497743463"
# DNa01 IDs (4 cells)
dna01aID <- "648518346476472438"
# dna01bID <- "648518346494596515"
dna01bID <- "648518346475464576"
dna01cID <- "648518346481256591"
dna01dID <- "648518346494759562"

# startNeuron <- dna02ID
# startNeuron <- dng13ID
# startNeuron <- dna01aID
# startNeuron <- dna01bID
# startNeuron <- dna01cID
startNeuron <- dna01dID

# get postsynaptic targets of start neuron as data.frame
primTarg <- fanc_partner_summary(startNeuron, "outputs", MIN_SYN, TRUE)

# save primary targets as CSV file
write.csv(primTarg, paste(startNeuron,".csv",sep=""),quote = FALSE)

# running list of neurons processed
# processedNeurons <- startNeuron
# processedNeurons <- NULL

# running list of neurons that throw errors
errorNeurons <- NULL

# loop over list of postsynaptic targets of start neuron
for (i in primTarg$post_id) {
  # check if in processedNeurons list
  if (!(i %in% processedNeurons)) {
  
    # try catch b/c getting errors for some IDs
    secTarg <- tryCatch(
      {
        fanc_partner_summary(i, "outputs", MIN_SYN, TRUE)
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
      for (j in secTarg$post_id) {
        if (!(j %in% processedNeurons)) {
          # try catch b/c some IDs throwing errors
          terTarg <- tryCatch(
            {
              fanc_partner_summary(j, "outputs", MIN_SYN, TRUE)
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
