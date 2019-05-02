library(dplyr)
library(tidyr)
library(data.table)
library(dendextend)
library(circlize)
library(seqinr)
library(DECIPHER)
library(ape)

## Load timetree, then relabel with underscore to match names in habitat data

fulltimetree <- ReadDendrogram(file = "actinopt_12k_treePL.tre")
fulltimetree_hclust <- as.hclust(fulltimetree)
fulltimetree_hclust$labels <- fulltimetree_hclust$labels %>% 
  base::gsub(" ", "_", .) %>% #replace spaces in labels with _
  base::gsub("'", "", .)
fulltimetree_dend <- as.dendrogram(fulltimetree_hclust)

## Load habitat data, remove duplicates

allhab <- fread(file = "Enbody_fish_species_habitat_combined_formatted.txt",
                col.names = c("spp", "marine", "fresh", "brackish")) %>% 
  distinct()

length(allhab$spp) == length(unique(allhab$spp)) # checks out

## Load allele alignment file

raboskyalleles <- fread(file = "final_alignment.selection.pro.213_261.tsv", header = F, col.names = c("spp", "allele")) %>% 
  extract(allele, into=c("r213", "r261"), '(.)(.)') %>%
  filter(r261 == "F" | r261 == "Y")

## Find intersection between habitat species, alignment species, and phylogeny species and remove problematic alignments then no hab info
sppToDrop_badAln <- scan("spp_to_drop_badAln.txt", character())
phyloWithHab <- allhab %>% 
  filter(spp %in% fulltimetree_hclust$labels) %>% 
  filter(spp %in% raboskyalleles$spp) %>% 
  filter(!spp %in% sppToDrop_badAln) %>%
  filter(!(marine == "F" & fresh == "F" & brackish == "F"))
  
## Export reduced species list for dendrogram pruning

write.table(phyloWithHab$spp, file = "spp_to_keep.txt", col.names = F, row.names = F, quote = F)

## Manually trim off root with nano and load trimmed timetree

trimmredTree <- ReadDendrogram(file = "timetree_trimmed.tre")
trimmredTree_hclust <- as.hclust(trimmredTree)
trimmredTree_hclust$labels <- trimmredTree_hclust$labels %>%
  base::gsub(" ", "_", .) %>% #replace spaces in labels with _
  base::gsub("'", "", .)
trimmredTree_dend <- as.dendrogram(trimmredTree_hclust)

rabo_allele_hab <- tibble(spp = trimmredTree_hclust$labels) %>% 
  inner_join(raboskyalleles) %>%
  inner_join(phyloWithHab) %>% 
  mutate(ismarine = if_else(marine == "T" & fresh == "F" & brackish == "F", TRUE, FALSE)) %>% 
  mutate(habcol = if_else(ismarine == FALSE & r261 == "F", "#FF000055",
                          if_else(ismarine == FALSE & r261 == "Y", "#FF0000FF",
                                  if_else(ismarine == TRUE & r261 == "Y", "#0000FF55","#0000FFFF"))))

raboredlabs <- rabo_allele_hab %>% filter(r261 == "Y") %>% .$spp # Get red labels
rabodendpruned <- trimmredTree_dend %>%  #Convert back to dendrogram
  color_branches(k = 1, col = 4) %>%
  set("by_labels_branches_col", value = c(raboredlabs), type = "all")
  

nspprabo2 <- length(trimmredTree_hclust$labels)
rabodendpruned_height <- attr(rabodendpruned, "height")
rabodendpruned_habcol <- tibble(spp = trimmredTree_hclust$labels) %>% 
  mutate(x = row_number()) %>% 
  filter(spp %in% rabo_allele_hab$spp)

rabo_order <- order.hclust(trimmredTree_hclust)
rabo_allele_hab <- rabo_allele_hab[c(rabo_order),] %>% 
  mutate(tiporder = row_number())
  

circos.par(cell.padding = c(0, 0, 0, 0),
           start.degree = -55,
           gap.after = 10,
           track.height = 0.075)
circos.initialize(factors = "a", xlim = c(0, nspprabo2)) # only one sector

circos.track(ylim = c(0, 1),
             bg.border = NA,
             panel.fun = function(x, y) {
               for(i in 1:nspprabo2) {
                 # circos.points(x = 2013, y = .5)
                 circos.rect(i-.5,0,i+.5,1,
                             col = rabo_allele_hab$habcol[i], border = "#00000000")
               }
               # for(j in seq(from = 5, to = nspprabo2, by = 10)) {
               #   circos.text(x = j, y = 1, labels = j, facing = "clockwise")
               # }
             })

circos.track(ylim = c(0, rabodendpruned_height), bg.border = NA, 
             track.height = .85, panel.fun = function(x, y) {
               circos.dendrogram(rabodendpruned)
             })

circos.clear()

## Overlap of species between rabo set and original
newspp <- raboskyalleles$spp[aa_hab_df_habcol$spp %in% raboskyalleles$spp]
temp <- setdiff(raboskyalleles$spp, aa_hab_df_habcol$spp)
temp2 <- intersect(raboskyalleles$spp, aa_hab_df_habcol$spp)
temp3 <- setdiff(aa_hab_df_habcol$spp, rabotre$tip.label)
temp4 <- setdiff()
temp <- 

# From the original phylogeny of 1054 species with habitat info, 178 do not appear in the big phylogeny but 
# 17 of those are because the big phylogeny has sub species names that don't match exactly but can be reconciled 
# (i.e. Acanthopagrus_schlegelii_schlegelii in big phylo but Acanthopagrus_schlegelii in the original one). 
# That leaves 161 species that we have data for that are not present in the big phylogeny, or 893 which overlap.

temp5 <- rabotre$tip.label[grep(paste(aa_hab_df_habcol$spp, collapse = "|"), rabotre$tip.label)]

# 876 of the original phylogeny tip labels appear in the big phylogeny. Maybe a few more due to subspp

temp6 <- raboskyalleles$spp[grep(paste(aa_hab_df_habcol$spp, collapse = "|"), raboskyalleles$spp, invert = T)]

# 1666 species appear in the big phylogeny with RH1 codon 261 but not original and so dont have habitat

# Test for significance of various alleles using the pruned phylogeny


