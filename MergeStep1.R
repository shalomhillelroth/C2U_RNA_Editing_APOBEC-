library(data.table)
library(dplyr)
argV <- commandArgs(trailingOnly = TRUE)

bedfiles = list.files(argV[1],
                      pattern="*_sByCDS_wDPStrand.bed*",
                      recursive = T, full.names = T)

bedfilesNames =  sapply(strsplit(basename(bedfiles), "_"),'[[',1)
chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
            "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")
for (chrom in chroms) {
  bedCounts = data.table()
  for (bedf in bedfiles) {
    name = strsplit(basename(bedf), "_")[[1]][1]
    bedCounts <- rbind(bedCounts, fread(bedf)%>% filter(V1 == chrom) %>% mutate(V6 = factor(x = name, levels =bedfilesNames) ))
  }
  
  
  coords <- unique(bedCounts$V2)
  clen <- length(coords)
  print(paste(chrom, ":", clen))
  prev_b <- 1
  for (batch in seq(min(clen,100000), clen , min(clen,100000))) {
    print(paste("Adding Records: ", prev_b, "-", batch))
    fwrite(tidyr::pivot_wider(data=bedCounts %>% filter(V2 %in% coords[prev_b:batch]),
                              id_cols =c ("V1", "V2", "V3", "V4"),
                              names_from = c("V6"),values_from = c("V5"),
                              values_fill = "NA", names_expand = T),
           file=paste0(argV[2], "/allSamples_withDepthPerStrand.bed"), sep="\t",
           col.names = F, quote = F, row.names = F,
           append = T)
    gc()
    prev_b <- batch + 1
  }
  
  if (prev_b < clen) {
    print(paste("Adding Records: ", prev_b, "-", clen))
    fwrite(tidyr::pivot_wider(data=bedCounts %>% filter(V2 %in% coords[prev_b:clen]),
                              id_cols =c ("V1", "V2", "V3", "V4"),
                              names_from = c("V6"),values_from = c("V5"),
                              values_fill = "NA", names_expand = T),
           file=paste0(argV[2], "/allSamples_withDepthPerStrand.bed"), sep="\t",
           col.names = F, quote = F, row.names = F,
           append = T)
    gc()
  }
  
}

