library(ggplot2)
library(ggthemes)
library(data.table)
library(dplyr)

transMM <- function(MM) {
  if(is.na(MM)){return(NA)}
  return(switch (MM,
                 "AG"= "AG",
                 "AC"= "AC",
                 "AT"= "AT",
                 "CA"= "CA",
                 "CG"= "CG",
                 "CT"= "CT",
                 "TC"= "AG",
                 "TG"= "AC",
                 "TA"= "AT",
                 "GT"= "CA",
                 "GC"= "CG",
                 "GA"= "CT",
                 #NA = "NA",
  ))
}

compMM <- function(MM) {
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

tissues <- list.dirs("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155", full.names = F,
                     recursive = F)
tissues<- Filter(function(x) !grepl("Step", x, ignore.case = T), tissues)
cdsPostClusterStats <- data.frame()
generalPostClusterStats <- data.frame()
for (Tissue in tissues) {
  if(!file.exists(paste0("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/",
                         Tissue,
                         "/Analysis/all_samples/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.forAnnovar.txt.hg38_multianno.csv"))){
    print(Tissue)
    next
  }
  annovarT <- fread(paste0("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/",
                           Tissue,
                           "/Analysis/all_samples/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.forAnnovar.txt.hg38_multianno.csv"))
  step1NewRunStrands <- unique(fread(paste0("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/",
                                            Tissue, 
                                            "/Clustering/SNPsAsMM/", "allSamples_hg38RefseqCuratedAllExons.NoAltChr.NoAltChr.NoOppStrand.bed")))
  step1NewRunStrands <- step1NewRunStrands %>% mutate(Chr = V1, End =as.numeric(V3), newRunMM = V4, Strand = V6) %>% 
    select(Chr, End, newRunMM, Strand)
  
  postClusterSNPsAsMM <-fread(paste0("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/",
                                     Tissue, 
                                     "/Clustering/SNPsAsMM/", Tissue,
                                     ".AllSites.NoBLAT.ClusterFiltered_100.StrandOblivious.tab"))%>%
    mutate(Chr = V1, End =as.numeric(V3), detectedMM = V4) %>% 
    select(Chr, End, detectedMM) %>% rowwise() %>%mutate(MM6 = transMM(detectedMM[1]))
  
  postClusterSNPsAsMM <- merge(postClusterSNPsAsMM, step1NewRunStrands, all.x=T)%>% rowwise() %>%
    mutate(MMFixed = if_else(Strand == "-", compMM(detectedMM[1]), detectedMM[1]))
  postClusterSNPsAsMM <- merge(postClusterSNPsAsMM, annovarT %>%select(Chr, End, Func.ncbiRefSeq), by = c("Chr", "End"), all.x=T)
  postClusterSitesStatsCDS <- postClusterSNPsAsMM %>% filter(Func.ncbiRefSeq == "exonic")%>%
    group_by(MMFixed) %>% summarise(
      NewRunN = n()
    )%>% filter(!is.na(MMFixed))%>%ungroup() %>%
    mutate(pc = round(100 * NewRunN / sum(NewRunN), digits = 2)) %>%
    mutate(res = paste0(NewRunN, " (", pc, "%)") )%>% select(res, MMFixed) %>%
    tidyr::pivot_wider(names_from = c("MMFixed"),values_from = c("res")) %>%
    mutate(type = "SNPs as MM", Tissue=Tissue)
  cdsPostClusterStats <- rbind(cdsPostClusterStats, postClusterSitesStatsCDS)
  
  postClusterSitesStatsGeneral <- postClusterSNPsAsMM%>%
    group_by(MM6, Func.ncbiRefSeq) %>% summarise(
      NewRunN = n()
    )%>%ungroup() %>%group_by(Func.ncbiRefSeq) %>%
    mutate(pc = round(100 * NewRunN / sum(NewRunN), digits = 2)) %>%
    mutate(res = paste0(NewRunN, " (", pc, "%)") )%>% select(res, MM6) %>%
    tidyr::pivot_wider(names_from = c("MM6"),values_from = c("res")) %>%
    mutate(type = "SNPs as MM", Tissue=Tissue)
  generalPostClusterStats <- rbind(generalPostClusterStats, postClusterSitesStatsGeneral)
  
  postClusterSNPMasked <-fread(paste0("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/",
                                     Tissue, 
                                     "/Clustering/SNPsMasked/", Tissue,
                                     ".AllSites.NoBLAT.ClusterFiltered_100.StrandOblivious.tab"))%>%
    mutate(Chr = V1, End =as.numeric(V3), detectedMM = V4) %>% 
    select(Chr, End, detectedMM) %>% rowwise() %>%mutate(MM6 = transMM(detectedMM[1]))
  
  postClusterSNPMasked <- merge(postClusterSNPMasked, step1NewRunStrands, all.x=T)%>% rowwise() %>%
    mutate(MMFixed = if_else(Strand == "-", compMM(detectedMM[1]), detectedMM[1]))
  postClusterSNPMasked <- merge(postClusterSNPMasked, annovarT %>%select(Chr, End, Func.ncbiRefSeq), by = c("Chr", "End"), all.x=T)
  postClusterSitesStatsCDS <- postClusterSNPMasked %>% filter(Func.ncbiRefSeq == "exonic")%>%
    group_by(MMFixed) %>% summarise(
      NewRunN = n()
    )%>% filter(!is.na(MMFixed))%>%ungroup() %>%
    mutate(pc = round(100 * NewRunN / sum(NewRunN), digits = 2)) %>%
    mutate(res = paste0(NewRunN, " (", pc, "%)") )%>% select(res, MMFixed) %>%
    tidyr::pivot_wider(names_from = c("MMFixed"),values_from = c("res")) %>%
    mutate(type = "SNPs Masked", Tissue=Tissue)
  cdsPostClusterStats <- rbind(cdsPostClusterStats, postClusterSitesStatsCDS)
  postClusterSitesStatsGeneral <- postClusterSNPMasked%>%
    group_by(MM6, Func.ncbiRefSeq) %>% summarise(
      NewRunN = n()
    )%>%ungroup() %>%group_by(Func.ncbiRefSeq) %>%
    mutate(pc = round(100 * NewRunN / sum(NewRunN), digits = 2)) %>%
    mutate(res = paste0(NewRunN, " (", pc, "%)") )%>% select(res, MM6) %>%
    tidyr::pivot_wider(names_from = c("MM6"),values_from = c("res")) %>%
    mutate(type = "SNPs Masked", Tissue=Tissue)
  generalPostClusterStats <- rbind(generalPostClusterStats, postClusterSitesStatsGeneral)
}

allAnnovar <-unique(fread("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/AllTissuesSites.hg38_multianno.csv") %>%
  mutate(MM6 = paste0(Ref, Alt))) %>% rowwise() %>%mutate(MM6 = transMM(MM6[1]))
postClusterSNPsAsMM <-fread("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/AllTissues.AllSites.SNPsAsMM.NoBLAT.ClusterFiltered_100.StrandOblivious.WithMM.bed")%>%
  mutate(Chr = V1, End =as.numeric(V3), detectedMM = V4) %>% 
  select(Chr, End, detectedMM) %>% rowwise() %>%mutate(MM6 = transMM(detectedMM[1])) %>% select(-detectedMM)
postClusterSNPsAsMM <- unique(postClusterSNPsAsMM)
postClusterSNPsAsMM <- unique(merge(postClusterSNPsAsMM,
                             allAnnovar %>%select(Chr, End, Func.ncbiRefSeq, MM6),
                             by = c("Chr", "End", "MM6"), all.x=T))
ccc <- paste(postClusterSNPsAsMMS$Chr, postClusterSNPsAsMMS$End)
postClusterSNPsAsMMC <- postClusterSNPsAsMM %>% filter(paste(Chr, End) %in% ccc) %>%
  group_by(MM6) %>% summarise(n = n())
postClusterSitesStatsGeneral <- postClusterSNPsAsMM%>%
  group_by(MM6, Func.ncbiRefSeq) %>% summarise(
    NewRunN = n()
  )%>%ungroup() %>%group_by(Func.ncbiRefSeq) %>%
  mutate(pc = round(100 * NewRunN / sum(NewRunN), digits = 2)) %>%
  mutate(res = paste0(NewRunN, " (", pc, "%)") )%>% select(res, MM6) %>%
  tidyr::pivot_wider(names_from = c("MM6"),values_from = c("res")) %>%
  mutate(type = "SNPs as MM", Tissue="_Combined")
generalPostClusterStats <- rbind(generalPostClusterStats, postClusterSitesStatsGeneral)


postClusterSNPMasked <-fread("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/AllTissues.AllSites.SNPsMasked.NoBLAT.ClusterFiltered_100.StrandOblivious.WithMM.bed")%>%
  mutate(Chr = V1, End =as.numeric(V3), detectedMM = V4) %>% 
  select(Chr, End, detectedMM) %>% rowwise() %>%mutate(MM6 = transMM(detectedMM[1])) %>%
  select(-detectedMM)
postClusterSNPMasked <- unique(postClusterSNPMasked)
postClusterSNPMasked <- unique(merge(postClusterSNPMasked, 
                                     allAnnovar %>%
                                       select(Chr, End, Func.ncbiRefSeq, MM6),
                                     by = c("Chr", "End", "MM6"), all.x=T))
postClusterSitesStatsGeneral <- postClusterSNPMasked%>%
  group_by(MM6, Func.ncbiRefSeq) %>% summarise(
    NewRunN = n()
  )%>%ungroup() %>%group_by(Func.ncbiRefSeq) %>%
  mutate(pc = round(100 * NewRunN / sum(NewRunN), digits = 2)) %>%
  mutate(res = paste0(NewRunN, " (", pc, "%)") )%>% select(res, MM6) %>%
  tidyr::pivot_wider(names_from = c("MM6"),values_from = c("res")) %>%
  mutate(type = "SNPs Masked", Tissue="_Combined")
generalPostClusterStats <- rbind(generalPostClusterStats, postClusterSitesStatsGeneral)

write.csv(cdsPostClusterStats %>% relocate(Tissue, type ),
          "/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/PostClusteringStats.CDS.csv",
          quote = F, row.names = F)
write.csv(generalPostClusterStats  %>% relocate(Tissue, type),
          "/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/PostClusteringStats.AllRegions.csv",
          quote = F, row.names = F)

library(ggplot2)
library(ggthemes)
library(hrbrthemes)
library(data.table)
library(dplyr)
extrafont::font_import()
extrafont::loadfonts()
extrafont::loadfonts(device = "win")

library(showtext)
font_add("RobotoCondensed", "/home/alu/hillelr/.fonts/RobotoCondensed-Regular.ttf")  # Use the actual file path
showtext_auto() 

postClusterSitesStatsGeneral <- postClusterSNPsAsMM%>%
  group_by(MM6, Func.ncbiRefSeq) %>% summarise(
    NewRunN = n()
  )%>%ungroup() %>%group_by(Func.ncbiRefSeq) %>%
  mutate(pc = round(100 * NewRunN / sum(NewRunN), digits = 2)) %>%
  mutate(res = paste0(NewRunN, " (", pc, "%)") )%>% select(pc, NewRunN, MM6) 

ggplot(data=postClusterSitesStatsGeneral%>% filter(!grepl(";", Func.ncbiRefSeq)))+
  geom_segment( aes(x=MM6, xend=MM6, y=0, yend=NewRunN), color="grey") +
  geom_point( aes(x=MM6, y=NewRunN, color=MM6), size=3 ) +
  geom_text( aes(x=MM6, y=NewRunN,
                 label=paste0(NewRunN , " (", round(pc, digits = 2),"%)")),
             size=3, vjust = 1, hjust = -0.5) +
  
  coord_flip()+
  theme_ipsum(base_family="Arial") +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.spacing = unit(0.1, "lines")
  ) +
  scale_y_continuous(trans="log10")+
  xlab("") +
  ylab("# Sites") +
  ggtitle("Sites Genomic Function")+
  facet_wrap(~Func.ncbiRefSeq, ncol=1)
ggsave("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/PostClusteringStats.AllRegions.FuncLolipop.pdf", height = 14, width = 16)

postClusterSitesStatsGeneral <- postClusterSNPsAsMM%>%
  group_by(MM6, Func.ncbiRefSeq) %>% summarise(
    NewRunN = n()
  )%>%ungroup() %>%group_by(MM6) %>%
  mutate(pc = round(100 * NewRunN / sum(NewRunN), digits = 2)) %>%
  mutate(res = paste0(NewRunN, " (", pc, "%)") )%>% select(pc, NewRunN, MM6, Func.ncbiRefSeq) 


ggplot(data=postClusterSitesStatsGeneral%>% filter(!grepl(";", Func.ncbiRefSeq)))+
  geom_segment( aes(x=Func.ncbiRefSeq, xend=Func.ncbiRefSeq, y=0, yend=pc), color="grey") +
  geom_point( aes(x=Func.ncbiRefSeq, y=pc, color=Func.ncbiRefSeq), size=3 ) +
  geom_text( aes(x=Func.ncbiRefSeq, y=pc,
                 label=paste0(pc, "% (", NewRunN, ")")),
             size=3, vjust = 1, hjust = -0.5) +
  
  coord_flip()+
  theme_ipsum(base_family="Arial") +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.spacing = unit(0.1, "lines")
  ) +
  # scale_y_continuous(trans="log10")+
  xlab("") +
  ylab("# Sites") +
  ggtitle("Sites Genomic Function")+
  facet_wrap(~MM6, ncol=1)
ggsave("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/PostClusteringStats.AllRegions.FuncLolipop.ByMM.pdf", height = 14, width = 16)
#### RMSK ####
allFinalSitesRMSK <- fread("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/AllSitesPostClustering.SNPsAsMM.WithRMSK.bed",
                           col.names = c("Chr", "Start", "End", "MM", "Rep", "Type", "Family"))

allFinalSitesRMSK <- allFinalSitesRMSK %>%rowwise() %>% mutate(Mismatch = transMM(MM[1]))
#decide what to do with <Family>?
allFinalSitesRMSK <- allFinalSitesRMSK %>% mutate(Type = gsub("?", "", Type, fixed = T))
allFinalSitesRMSKStats <- allFinalSitesRMSK %>% 
  group_by(Mismatch, Type) %>%
  summarise(n=n()) %>% mutate(freq = round(100* n / sum(n), digits=2)) %>% 
  arrange(freq, .by_group = T)

ggplot(allFinalSitesRMSKStats) +
  geom_segment( aes(x=Type, xend=Type, y=0, yend=freq), color="grey") +
  geom_point( aes(x=Type, y=freq, color=Mismatch), size=3 ) +
  geom_text( aes(x=Type, y=freq, label=paste0(freq , "% (", n ,")")), size=3, vjust = 1, hjust = -0.5) +
  
  coord_flip()+
  theme_ipsum(base_family="Arial") +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.spacing = unit(0.1, "lines")
  ) +
  theme_bw()+
  xlab("") +
  ylab("# Sites") +
  ggtitle("Sites Within Non-Alu Repeats")+
  facet_grid(Mismatch~., scales = "free_y", space="free")

ggsave("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/PostClusteringStats.AllRegions.RMSKLolipop.pdf", height = 14, width = 16)
allFinalSitesRMSKStats <- allFinalSitesRMSKStats%>% filter(Mismatch %in% c("AG", "CT"))
allFinalSitesRMSKStats$Type <- factor(allFinalSitesRMSKStats$Type, levels=unique(allFinalSitesRMSKStats$Type))

ggplot(allFinalSitesRMSKStats) +
  geom_segment( aes(x=Type, xend=Type, y=0, yend=freq), color="grey") +
  geom_point( aes(x=Type, y=freq, color=Mismatch), size=3 ) +
  geom_text( aes(x=Type, y=freq, label=paste0(freq , "% (", n ,")")), size=3, vjust = 1, hjust = -0.5) +
  
  coord_flip()+
  theme_ipsum(base_family="Arial") +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.spacing = unit(0.1, "lines")
  ) +
  theme_bw()+
  xlab("") +
  ylab("# Sites") +
  ggtitle("Sites Within Non-Alu Repeats")+
  facet_grid(Mismatch~., scales = "free_y", space="free")

ggsave("/private9/Projects/KnownC2T/DeNovoDetect/GTExV8.dbSNP155/PostClusteringStats.AllRegions.AG_CT.pdf", height = 10, width = 12)
