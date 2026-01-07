setwd("/projectsp/f_ak1833_1/jliu/ASD/02.Features/03.CalcDiff")
library(dplyr)
library(data.table)
library(Biostrings)
library(RecordLinkage)

# import raw --- 
## features -----
features_odd.raw <- 
  fread("~/ASD/02.Features/02.FeatureSummary/summary_odd.csv") %>% as.data.frame
features_even.raw <-
  fread("~/ASD/02.Features/02.FeatureSummary/summary_even.csv") %>% as.data.frame
## log2fc file -----
log2fc <- read.delim("~/ASD/00.raw/RefAlt_resout.tsv", check.names = F, stringsAsFactors = F)
## fasta -----
fasta_odd = readDNAStringSet("~/ASD/02.Features/01.Featurizer/odd/locs.fasta")
fasta_even = readDNAStringSet("~/ASD/02.Features/01.Featurizer/even/locs.fasta")
## fimo -----
fimo <- 
  rbind(
    fread("~/ASD/02.Features/01.Featurizer/even/fimo_hg19/fimo.txt"), 
    fread("~/ASD/02.Features/01.Featurizer/even/fimo_encode/fimo.txt"), 
    fread("~/ASD/02.Features/01.Featurizer/odd/fimo_hg19/fimo.txt"), 
    fread("~/ASD/02.Features/01.Featurizer/odd/fimo_encode/fimo.txt")
  ) 

# data cleaning -----
## fimo -----
fimo <- fimo %>% select(-motif_alt_id) %>% .[complete.cases(.), ]
## log2fc -----
log2fc$sequence_name_odd = sapply(stringr::str_split(log2fc$Names, '\\|'), function(X) X[1])
log2fc$sequence_name_even = sapply(stringr::str_split(log2fc$Names, '\\|'), function(X) X[2])
## features list -----
features.list <- 
  list(odd = features_odd.raw, even = features_even.raw)
features.list <- 
  lapply(features.list, function(X){
    X$regions = stringr::str_extract(X$sequence_name, "chr[0-9]+\\:[0-9]+\\-[0-9]+")
    X$regions[is.na(X$regions)] = stringr::str_extract(X$sequence_name[is.na(X$regions)], "chrX\\:[0-9]+\\-[0-9]+")
    return(X)
  })

## get features -----
features_even <- features_even.raw[match(log2fc$sequence_name_even, features_even.raw$sequence_name), ]
features_odd <- features_odd.raw[match(log2fc$sequence_name_odd, features_odd.raw$sequence_name), ]
## sanity check
#all.equal(features_even$sequence_name, log2fc$sequence_name_even) # should be TRUE
#all.equal(features_odd$sequence_name, log2fc$sequence_name_odd) # should be TRUE
#all.equal(colnames(features_odd), colnames(features_even)) # should be TRUE
## calc diff -----
diff.df <- 
  data.frame(matrix(vector(), 
                    nrow = nrow(features_even), ncol = ncol(features_even)-1,
                    dimnames=list(c(), colnames(features_even)[-1])),
             stringsAsFactors=F)
diff.df = features_odd[, -1] - features_even[, -1]
diff.df$sequence_name_even = features_even$sequence_name
diff.df$sequence_name_odd = features_odd$sequence_name
diff.df = diff.df %>% relocate(sequence_name_even, .before = everything())
diff.df = diff.df %>% relocate(sequence_name_odd, .before = everything())

# calc similarity -----
## fasta odd and even -----
fasta_odd = fasta_odd[log2fc$sequence_name_odd]
fasta_even = fasta_even[log2fc$sequence_name_even]
s_odd <- data.frame(sequence_name_odd = names(fasta_odd), seq_odd = paste(fasta_odd))
s_even <- data.frame(sequence_name_even = names(fasta_even), seq_even = paste(fasta_even))

s_even_odd = cbind(s_odd, s_even)
s_even_odd = 
  s_even_odd %>% mutate(levenshteinSim = levenshteinSim(seq_odd, seq_even))
similarity <- s_even_odd %>% select(sequence_name_odd, sequence_name_even, levenshteinSim)

## update diff.df -----
if(all.equal(diff.df$sequence_name_odd, similarity$sequence_name_odd)){
  diff.df$levenshteinSim = similarity$levenshteinSim
}

# calc n.motif -----
n.motif = fimo %>% select(sequence_name, motif_id)
n.motif = table(n.motif) %>% as.data.frame(., stringsAsFactors = F)
n.motif = dcast(n.motif, sequence_name~motif_id, fill = 0)
n.motif_even = n.motif[match(log2fc$sequence_name_even, n.motif$sequence_name), ]
n.motif_odd = n.motif[match(log2fc$sequence_name_odd, n.motif$sequence_name), ]
# all.equal(n.motif_even$sequence_name, diff.df$sequence_name_even)
# all.equal(n.motif_odd$sequence_name, diff.df$sequence_name_odd)
delta.n.motif = n.motif_odd[, -1] - n.motif_even[, -1]
## update diff.df -----
diff.df = cbind(diff.df, delta.n.motif)

# save diff -----
fwrite(diff.df, "diff.txt", sep = "\t", quote = F, row.names = F)
