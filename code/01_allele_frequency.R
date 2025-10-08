library(Rsamtools)
library(dplyr)
library(data.table)
library(GenomicAlignments)
library(BuenColors)
library(ggbeeswarm)

# Define SNV of interest
chr="chrX"
pos=20534935
target_gr <- GRanges(chr, IRanges(c(pos-0,pos,pos+0), width=1))
sbp = ScanBamParam(which=target_gr, mapqFilter = 10)

# Function to get allele frequencies split by genotype
process_af <- function(day){
  
  dt <- fread(paste0("../data/hash/mHSPC_2_",day,"_hash_ids_Seurat.csv"), header = FALSE) 
  
  # parse 1 bam into pieces
  one_bam_file <- paste0("../data/bam/LSK.",day,".bam")
  bam_in <- readGAlignments(one_bam_file, param = ScanBamParam(tag = c("CB")))
  
  grab_af <- function(v2_list, what1, what2){
    ss_barcodes <- dt %>% filter(V2 %in% v2_list) %>% pull(V1)
    filterMe <- FilterRules(list(keep_cb=function(x) { data.frame(bam_in)[["CB"]] %in% ss_barcodes }))
    im_bam_file <- paste0("../temp/", day, "_ss.bam")
    filterBam(one_bam_file, im_bam_file, filter=filterMe)
    
    vardf <- pileup(im_bam_file, scanBamParam = sbp, pileupParam = PileupParam(distinguish_strands=FALSE, max_depth = 10000)) %>%
      filter(pos == 20534935)
    data.frame(
      what1, what2, 
      coverage = sum(vardf$count)
    ) %>% mutate(wt_af = (sum(vardf %>% filter(nucleotide == "T") %>% pull(count)))/coverage)
    
  }
  
  rbind(
    rbind(
      grab_af(paste0("R3_Rosa26_", day),"Rosa26", "R3"),
      grab_af(paste0("R3C8_Rosa26_", day),"Rosa26", "R3C8"),
      grab_af(paste0("WT_Rosa26_",day), "Rosa26", "aWT")
    ),
    rbind(
      grab_af(paste0("R3_M41T_", day),"M41T", "R3"),
      grab_af(paste0("R3C8_M41T_", day),"M41T", "R3C8"),
      grab_af(paste0("WT_M41T_",day), "M41T", "aWT")
    )
  ) %>% mutate(day = day)
}

af_df <- rbind(
  process_af("d5"), 
  process_af("d8")
) %>% mutate(mut_af = (1-wt_af)*100) %>%
  mutate(SE = sqrt(mut_af/100*(1-mut_af/100)/coverage)*100)

p1 <- ggplot(af_df, aes(x = interaction(what2, day), y = mut_af, color = what1)) +
  geom_point() +
  scale_color_manual(values = c("dodgerblue3", "lightblue")) +
  labs(x = "Genotype x Day", y = "% M41T editing", color = "gRNA", shape = "genotype") +
  geom_errorbar(aes(ymin=mut_af-SE, ymax=mut_af+SE), width=.2) + pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") 
cowplot::ggsave2(p1, file = "../figures/allele_frequency.pdf", width = 2.5, height = 2)

