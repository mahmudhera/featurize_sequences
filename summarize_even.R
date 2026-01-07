setwd("~/ASD/02.Features/02.FeatureSummary")

library(dplyr)
library(magrittr)
library(data.table)

dir.input <- "~/ASD/02.Features/01.Featurizer/even/"

# import raw -----
features <- list.files(dir.input, pattern = ".csv|.txt|predictions.tsv", recursive=T)
print(features)
#[1] "5mer.csv"             "deepbind.csv"         "dna_shape.csv"        "fimo_encode/fimo.txt"
#[5] "fimo_hg19/fimo.txt"   "fimo_summary.csv"     "polyA_polyT_GC.csv" 
kmer5 <- read.csv(paste0(dir.input, "5mer.csv"))
deepbind <- read.csv(paste0(dir.input, "deepbind.csv"))
deepsea <- fread(paste0(dir.input, "deepsea/261912e4-352b-4821-91ea-81d5c11442be_fa_variants.even_w0_predictions.tsv"), drop = 1) %>% as.data.frame
#dna_shape <- read.csv(paste0(dir.input, "dna_shape.csv"))
dna_shapeR <- read.csv(paste0(dir.input, "DNAshapeR.csv"))
fimo_encode <- read.delim(paste0(dir.input, "fimo_encode/fimo.txt"))
fimo_hg19 <- read.delim(paste0(dir.input, "fimo_hg19/fimo.txt"))
fimo_summary <- read.csv(paste0(dir.input, "fimo_summary.csv"))
polyA_polyT <- read.csv(paste0(dir.input, "polyA_polyT_GC.csv"))
## TFDB -----
TFDB <- read.delim("~/data/AnimalTFDB/Homo_sapiens_TF_Family.main.txt")

# data cleaning -----
## change the 1st column name as "sequence name" -----
colnames(kmer5)[1] = "sequence_name"
colnames(deepbind)[1] = "sequence_name"
colnames(deepsea)[1] = "sequence_name"
colnames(dna_shapeR)[1] = "sequence_name"
colnames(fimo_encode)[3] = "sequence_name"
colnames(fimo_hg19)[3] = "sequence_name"
colnames(fimo_summary)[1] = "sequence_name"
colnames(polyA_polyT)[1] = "sequence_name"
## fimo
fimo_encode <- fimo_encode %>% select(-motif_alt_id) %>% .[complete.cases(.), ]
fimo_hg19 <- fimo_hg19 %>% select(-motif_alt_id) %>% .[complete.cases(.), ]

## deepbind_top -----
deepbind_matrix = deepbind %>% `rownames<-`(deepbind$sequence_name) %>% .[, -1] %>% as.matrix
cutoff.deepbind <- quantile(deepbind_matrix, probs = c(0.9)) 
deepbind_90p_bl <- deepbind_matrix >= cutoff.deepbind
deepbind_top <- rowSums(deepbind_90p_bl)
deepbind_top <- data.frame(sequence_name = names(deepbind_top), n.deepbind_top = deepbind_top)
## deepsea_top -----
deepsea_matrix = deepsea %>% `rownames<-`(deepsea$sequence_name) %>% .[, -1] %>% as.matrix
cutoff.deepsea <- quantile(deepsea_matrix, probs = c(0.9)) 
deepsea_90p_bl <- deepsea_matrix >= cutoff.deepsea
deepsea_top <- rowSums(deepsea_90p_bl)
deepsea_top <- data.frame(sequence_name = names(deepsea_top), n.deepsea_top = deepsea_top)
## n.5mer -----
kmer5_matrix = kmer5 %>% `rownames<-`(kmer5$sequence_name) %>% .[, -1] %>% as.matrix
kmer5_matrix_distinct <- kmer5_matrix !=  0
n.5mer <- rowSums(kmer5_matrix_distinct)
n.5mer <- data.frame(sequence_name = names(n.5mer), n.5mer = n.5mer)
## polyA_polyT -----
polyA_polyT <- polyA_polyT %>% select(-GC) %>% `colnames<-`(c("sequence_name", "n.polyA", "n.polyT"))
## nMotifs, MotifDensity -----
fimo_summary = fimo_summary %>% `colnames<-`(c("sequence_name", "n.Motifs_ENCODE", "MotifDensity_ENCODE", "n.Motifs_hg19", "MotifDensity_hg19"))
## n activators, n repressors, TFDB families
### encode
fimo_encode$regions = stringr::str_extract(fimo_encode$sequence_name, "chr[0-9]+\\:[0-9]+\\-[0-9]+")
fimo_encode$regions[is.na(fimo_encode$regions)] = stringr::str_extract(fimo_encode$sequence_name[is.na(fimo_encode$regions)], "chrX\\:[0-9]+\\-[0-9]+")
fimo_encode$Symbol <- gsub("\\_.*", "", fimo_encode$motif_id)
fimo_encode <- merge(fimo_encode, TFDB %>% select(Symbol, Family, Family.main), 
                     by = "Symbol", 
                     all.x = T)
fimo_encode$Family.main[is.na(fimo_encode$Family.main)] = "unknown"
fimo_encode$Family[is.na(fimo_encode$Family)] = "unknown"
Family.freq_fimo_encode = fimo_encode %>% select(sequence_name, Family.main) %>% table %>% as.data.frame.matrix()
colnames(Family.freq_fimo_encode) = paste0("n.", colnames(Family.freq_fimo_encode), "_encode")
Family.freq_fimo_encode = Family.freq_fimo_encode[fimo_summary$sequence_name, ]
if(all.equal(rownames(Family.freq_fimo_encode), fimo_summary$sequence_name)){
  fimo_summary = cbind(fimo_summary, Family.freq_fimo_encode)
}
### hg19
fimo_hg19$regions = stringr::str_extract(fimo_hg19$sequence_name, "chr[0-9]+\\:[0-9]+\\-[0-9]+")
fimo_hg19$regions[is.na(fimo_hg19$regions)] = stringr::str_extract(fimo_hg19$sequence_name[is.na(fimo_hg19$regions)], "chrX\\:[0-9]+\\-[0-9]+")
fimo_hg19$Symbol <- gsub("\\_.*", "", fimo_hg19$motif_id)
fimo_hg19 <- merge(fimo_hg19, TFDB %>% select(Symbol, Family, Family.main), 
                   by = "Symbol", 
                   all.x = T)
fimo_hg19$Family.main[is.na(fimo_hg19$Family.main)] = "unknown"
fimo_hg19$Family[is.na(fimo_hg19$Family)] = "unknown"
Family.freq_fimo_hg19 = fimo_hg19 %>% select(sequence_name, Family.main) %>% table %>% as.data.frame.matrix()
colnames(Family.freq_fimo_hg19) = paste0("n.", colnames(Family.freq_fimo_hg19), "_hg19")
Family.freq_fimo_hg19 = Family.freq_fimo_hg19[fimo_summary$sequence_name, ]
if(all.equal(rownames(Family.freq_fimo_hg19), fimo_summary$sequence_name)){
  fimo_summary = cbind(fimo_summary, Family.freq_fimo_hg19)
}

## merge -----
summary <- 
  Reduce(function(x, y) merge(x,y, by = "sequence_name"), 
         list(polyA_polyT, n.5mer, dna_shapeR, fimo_summary, deepbind_top, deepsea_top, 
              kmer5, deepsea, 
              deepbind))

# write out summary
write.csv(summary, "summary_even.csv", row.names = F)
