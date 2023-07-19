# Reformating LC M005 proteomics peptide data
# Burns et al. (2023) : Proteomic changes induced by longevity promoting interventions in mice 
# Adam Burns - 2023-07-06

## Load required packages
library(tidyverse)

## Loop through tissues
for(x in c('Liver', 'Kidney', 'Gastroc')){
  
  ## Import data
  dat.0 <- read.csv("Mouse_2021.10_Liver_peptide_unnormalized_proteotypic-Phase2.csv", header=TRUE, quote="")
  
  ## Combine charge states
  ## & filter degenerate peptides and singleton
  dat.0 <- dat.0 %>%
    separate(col = Peptide, into = c('Peptide', 'Charge'), sep = "\\.", extra="merge", fill="right") %>%
    select(-Charge)
    
  # Sum across charge states
  dat.0[is.na(dat.0)] <- 0 
    
  dat.0 <- aggregate(. ~ Peptide + Protein, data = dat.0, sum)
  
  # Replace 0 with NA
  dat.0[dat.0 == 0] <- NA 
  
  
  ## Convert to long
  dat.1 <- as.data.frame(dat.0) %>%
    pivot_longer(cols = colnames(dat.0)[!(colnames(dat.0) %in% c('Peptide', 'Protein'))], names_to = "SampleID", values_to = "y") %>%
    filter(!(is.na(y))) %>%
    mutate(tissue = x)
  
  
  ## Filter degenerate peptides (match to multiple proteins)
  ## Filter singleton (only one peptide per protein) peptides
  ## Filter proteins detected in less than 4 samples
  dat.2 <- dat.1 %>%
    group_by(Peptide) %>%
    mutate(protperpep =  length(unique(Protein))) %>%
    ungroup() %>%
    filter(protperpep == 1) %>%
    group_by(Protein) %>%
    mutate(pepperprot = length(unique(Peptide))) %>%
    mutate(sampperprot = length(unique(SampleID))) %>%
    ungroup() %>%
    filter(protperpep == 1, pepperprot > 1, sampperprot >= 4)
  
  # Assign
  assign(paste0('dat.', x), dat.2)

}


## Combine tissues
dat.3 <- bind_rows(dat.Liver, dat.Kidney, dat.Gastroc)


## Create peptide key
peptide.key <- as.data.frame(cbind(unique(dat.3$Peptide), paste0("p", seq(1:length(unique(dat.3$Peptide))))))
colnames(peptide.key) <- c("Peptide.seq", "PeptideID")


## Add peptide ID
dat.4 <- dat.3 %>%
  left_join(peptide.key, by=c("Peptide" = "Peptide.seq")) %>%
  select(PeptideID, Peptide, Protein, pepperprot, SampleID, tissue, y) %>%
  distinct()


## Export/Save file
write.csv(dat.4, file="lcproteomics_m005_formatted.csv", quote=FALSE, row.names=FALSE)

