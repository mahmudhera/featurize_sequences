library(DNAshapeR)
library("Biostrings")
library(dplyr)
library(tibble)

fa <- readDNAStringSet("./locs.fasta")
sequence_names <- names(fa)
# prediction
pred_inter <- getShape("./locs.fasta", 
                       shapeType = c('HelT', 'Rise', 'Roll', 'Shift', 'Slide', 'Tilt' # inter bp
                       ))
pred_intra <- getShape("./locs.fasta", 
                       shapeType = c(
                         'Buckle', 'Opening', 'ProT', 'Shear', 'Stagger', 'Stretch', "MGW" #intra bp
                       ))
# pred into list of data frames
pred_inter_list <- 
  lapply(pred_inter, function(X){
    X[is.na(X)] = 0
    rownames(X) = sequence_names
    return(X)
  })
pred_intra_list <- 
  lapply(pred_intra, function(X){
    X[is.na(X)] = 0
    rownames(X) = sequence_names
    return(X)
  })
# add colnames to data frame in the list
pred_inter_list_wi_colnames <- 
  lapply(names(pred_inter_list), function(X){
    out = pred_inter_list[[X]]
    colnames(out) = paste0(X, ".nt", 1:ncol(out))
    return(out)
    })
pred_intra_list_wi_colnames <- 
  lapply(names(pred_intra_list), function(X){
    out = pred_intra_list[[X]]
    colnames(out) = paste0(X, ".nt", 1:ncol(out))
    return(out)
  })
# reduce list to df
pred_inter_df <- Reduce(`cbind`, pred_inter_list_wi_colnames)
pred_intra_df <- Reduce(`cbind`, pred_intra_list_wi_colnames)

pred <- cbind(pred_inter_df, pred_intra_df)

# write out
write.csv(pred, "DNAshapeR.csv", quote = F)

## remove proc files
procfiles <- list.files('.', pattern = "locs.fasta..*")
file.remove(procfiles)
