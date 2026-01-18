library(ggplot2)
library(ggthemes)
library(dplyr)
library(data.table)
library(reshape2)

transMM <- function(MM) {
  if(is.na(MM)){return(NA)}
  return(switch (MM,
                 "AG"= "TC",
                 "AC"= "TG",
                 "AT"= "TA",
                 "CA"= "GT",
                 "CG"= "GC",
                 "CT"= "GA",
                 "TC"= "AG",
                 "TG"= "AC",
                 "TA"= "AT",
                 "GT"= "CA",
                 "GC"= "CG",
                 "GA"= "CT",
                 #NA = "NA",
  ))
}


getCanocical <- function(mm, A, C, G, `T`) {
  ref = substr(mm,1,1)
  # if(mm == "AT" || mm == "TA"){
  #   return(A)
  # }
  # if(mm == "CG" || mm == "GC"){
  #   return(C)
  # }
  return(switch (ref,
                 "A" = A,
                 "C" = C,
                 "G" = G,
                 `T` = `T`))
}


getMismatched <- function(mm, A, C, G, `T`) {
  ref = substr(mm,2,2)
  # if(mm == "AT" || mm == "TA"){
  #   return(`T`)
  # }
  # if(mm == "CG" || mm == "GC"){
  #   return(G)
  # }
  return(switch (ref,
                 "A" = A ,
                 "C" = C,
                 "G" = G ,
                 `T` = `T`))
}


argV <- commandArgs(trailingOnly = TRUE)
sitesBED <- argV[1]
DNAPooledBED<- argV[2] 
RNAPooledBED<- argV[3] 
outDir<- argV[4]

### Sites Meta ##############

deNovoSites<-read.table(sitesBED, sep="\t", col.names = c("Chr", "Start", "End", "MM", "Ref"))
deNovoSites <- deNovoSites %>% rowwise() %>%
  mutate(MM = if_else(substr(MM[1],1,1) != Ref,
                       transMM(MM[1]), MM[1]))
deNovoSites <- unique(deNovoSites)


### Pooled WGS ####
pooledWGS <- fread(DNAPooledBED)
pooledWGS <- merge(pooledWGS, deNovoSites,
                   by = c("Chr", "Start", "End"))
pooledWGS <- pooledWGS %>% rowwise() %>%
  dplyr::mutate(
    Canonical = getCanocical(MM[1], A[1], C[1], G[1], `T`[1]), 
    Mismatched = getMismatched(MM[1], A[1], C[1], G[1], `T`[1]),
    HighQCoverage = TotalCoverage - LowQ - N) %>%
  mutate(MMPC = 100*Mismatched/ (Mismatched+Canonical))

pooledWGS <- merge(pooledWGS %>% select(Chr, End, MM, Mismatched, Canonical, MMPC),
                    deNovoSites %>% select(Chr, End, MM), all = T)
pooledWGS <- pooledWGS %>% replace(is.na(.), 0)

dnaEdited<- pooledWGS %>% filter(Mismatched + Canonical >= 100,  MMPC >= 0.1) %>%
  mutate(Coverage = Mismatched + Canonical)
dnaEdited <- rbind(dnaEdited, pooledWGS %>% filter(Mismatched + Canonical < 100) %>%
                     mutate(Coverage = Mismatched + Canonical))
write.table(dnaEdited%>% select(Chr, End, MM, Coverage, MMPC),
            file=paste0(outDir,"/AllTissues.PooledCMPileup.WGSSites.txt"),
            sep = " ", quote = F, col.names = F, row.names = F)
rm(dnaEdited)
rm(pooledWGS)
gc()

### Pooled RNA ####
pooledRNA <- fread(RNAPooledBED)
pooledRNA <- merge(pooledRNA, deNovoSites,
                   by = c("Chr", "Start", "End"))

pooledRNA <- pooledRNA %>% #filter(TotalCoverage - N - LowQ > 10) %>%
  rowwise() %>%
  dplyr::mutate(
    Canonical = getCanocical(MM[1], A[1], C[1], G[1], `T`[1]), 
    Mismatched = getMismatched(MM[1], A[1], C[1], G[1], `T`[1]),
    Coverage= Mismatched+Canonical,
    HighQCoverage = TotalCoverage - LowQ - N) %>%
  mutate(MMPC = 100*Mismatched/ (Mismatched+Canonical)) %>% 
  select(Chr, Start, End, MM, Tissue, Canonical, Mismatched, MMPC, HighQCoverage, Coverage)

tissuesPerSites <- pooledRNA %>% filter(Mismatched + Canonical >= 10) %>%
  group_by(Chr, End, MM) %>%
  summarise(TissuesWithCoverage10 =n())
pooledEditingPerTissue <- pooledRNA %>% 
  filter(Tissue != "Combined",  MMPC >= 1, Mismatched + Canonical >= 10) %>% 
  mutate(Coverage = Mismatched + Canonical)
maxEditing<-pooledEditingPerTissue %>% group_by(Chr, End, MM) %>%summarise(
                                   MMPC= max(MMPC))
pooledEditingPerTissue <- merge(pooledEditingPerTissue, maxEditing)
editedTissuesPerSites <- pooledEditingPerTissue %>% group_by(Chr, End, MM) %>%
  summarise(TissuesWithMaxEditing =n())
editedTissuesPerSites <- merge(editedTissuesPerSites, tissuesPerSites, all=T)
maxEditing<-pooledEditingPerTissue %>% group_by(Chr, End, MM) %>%summarise(
  Coverage= max(Coverage))
pooledEditingPerTissue <- merge(pooledEditingPerTissue, maxEditing)
pooledEditingPerTissue <- pooledEditingPerTissue %>% 
  group_by(Chr, End, MM, Coverage, MMPC) %>% 
  summarise(MaxTissues=paste(Tissue, sep = ";"),
            NMaxTissues = n())

pooledEditingPerTissue <- merge(pooledEditingPerTissue, editedTissuesPerSites, all=T)
write.table(pooledEditingPerTissue %>% select(Chr, End, MM, Coverage, MMPC),
            file=paste0(outDir,"/AllTissues.EditingGt1PC.txt"),
            sep = " ", quote = F, col.names = F, row.names = F)
write.table(pooledEditingPerTissue,
            file=paste0(outDir,"/AllTissues.EditingGt1PC.tab"),
            sep = "\t", quote = F, col.names = T, row.names = F)
rm(pooledEditingPerTissue)
pooledEditingCombined <- pooledRNA %>% 
  filter(Tissue == "Combined",  MMPC >= 1, Mismatched + Canonical >= 100) %>% 
  mutate(Coverage = Mismatched + Canonical)
write.table(pooledEditingCombined %>% select(Chr, End, MM, Coverage, MMPC),
            file=paste0(outDir,"/AllTissues.Combined.EditingGt1PC.txt"),
            sep = " ", quote = F, col.names = F, row.names = F)
write.table(pooledEditingCombined,
            file=paste0(outDir,"/AllTissues.Combined.EditingGt1PC.tab"),
            sep = "\t", quote = F, col.names = T, row.names = F)

gc()

### Zscore ####
zScore <- function(coverage, mmReads, q) {
  p <- mmReads / coverage
  return (sqrt(coverage)*(p-q) / sqrt(q*(1-q)))
}
pooledRNA <- fread(RNAPooledBED)
pooledRNA <- merge(pooledRNA, deNovoSites,
                   by = c("Chr", "Start", "End"))

pooledRNA <- pooledRNA %>% #filter(TotalCoverage - N - LowQ > 10) %>%
  rowwise() %>%
  dplyr::mutate(
    Canonical = getCanocical(MM[1], A[1], C[1], G[1], `T`[1]), 
    Mismatched = getMismatched(MM[1], A[1], C[1], G[1], `T`[1]),
    Coverage= Mismatched+Canonical,
    z_score_0.01 = zScore(coverage = Coverage, mmReads = Mismatched, 0.01),
    z_score_0.005 = zScore(coverage = Coverage, mmReads = Mismatched, 0.005)) %>%
  mutate(MMPC = 100*Mismatched/ (Mismatched+Canonical)) %>% 
  select(Chr, Start, End, MM, Tissue, Canonical, Mismatched, MMPC, Coverage, z_score_0.01, z_score_0.005)

write.table(pooledRNA %>% select(Chr, End, MM, Coverage, MMPC,Tissue, z_score_0.01, z_score_0.005),
            file=paste0(outDir,"/AllTissues.CovGt10.ZScores.tab"),
            sep = "\t", quote = F, col.names = T, row.names = F)

maxScores <- pooledRNA %>% filter(Tissue != "Combined") %>%group_by(Chr, End, MM) %>% 
  filter(z_score_0.01 == max(z_score_0.01))


write.table( maxScores %>%
               select(Chr, End, MM, Tissue, Mismatched, Coverage, z_score_0.01, z_score_0.005),
            file=paste0(outDir,"/AllTissues.MaxZ_scores.Editing.tab"),
            sep = " ", quote = F, col.names = F, row.names = F)


maxScores <- pooledRNA %>%
  filter(!Tissue %in% c("Combined","Cervix_Ectocervix","Cervix_Endocervix", "Fallopian_Tube")) %>%
  group_by(Chr, End, MM) %>% 
  filter(z_score_0.01 == max(z_score_0.01))


write.table( maxScores %>%
               select(Chr, End, MM, Tissue, Mismatched, Coverage, z_score_0.01, z_score_0.005),
             file=paste0(outDir,"/AllTissues.MaxZ_scores.NoSmallTissues.Editing.tab"),
             sep = " ", quote = F, col.names = F, row.names = F)
rm(pooledRNA)
gc()

