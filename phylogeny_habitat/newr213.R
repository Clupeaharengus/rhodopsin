library(tidyr)
library(dplyr)
library(Biostrings)
library(ggplot2)

## Read in FASTA alignment of robosppecies

rhoAAaln <- readAAMultipleAlignment("final_alignment.translated.fullrhodopsin.fasta", format = "FASTA")
rhoAAss <- AAStringSet(rhoAAaln)  #Convert to string set
rhoAAss <- rhoAAss[rabodendpruned_habcol$spp]

r213 <- subseq(rhoAAss, 213, width = 1)

writeXStringSet(r213, file = "r213.fasta")

r213.df <- r213 %>% 
  as.data.frame() %>% 
  add_rownames()
colnames(r213.df) <- c("spp","r213")
new213hab.df <- rabo_allele_hab %>% 
  select(-r213) %>% 
  left_join(r213.df)

pos213_aacount_261Y <- new213hab.df %>% filter(r261 == "Y") %>% group_by(r213) %>% tally() %>% mutate(Y261 = n) %>% 
  select(-n)

pos213_aacount_261F <- new213hab.df %>% filter(r261 == "F") %>% group_by(r213) %>% tally() %>% mutate(F261 = n) %>% 
  select(-n)

full_join(pos213_aacount_261F, pos213_aacount_261Y) %>% 
  write.table(file = "pos213_pos261YF.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

##

pos213_aacounts_marine <- marineFisherscore.df.gather %>% filter(aapos == 213)
pos213_aacounts_nonmarine <- nonmarineFisherscore.df.gather %>% filter(aapos == 213)
pos213_aacounts <- full_join(pos213_aacounts_marine, pos213_aacounts_nonmarine)
pos213_aacounts %>% select(aa, marine_count, non_marine_count) %>% 
  write.table(file = "pos213_aacounts.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

pos214_aacounts_marine <- marineFisherscore.df.gather %>% filter(aapos == 214)
pos214_aacounts_nonmarine <- nonmarineFisherscore.df.gather %>% filter(aapos == 214)
pos214_aacounts <- full_join(pos214_aacounts_marine, pos214_aacounts_nonmarine)
pos214_aacounts %>% select(aa, marine_count, non_marine_count) %>% 
  write.table(file = "pos214_aacounts.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

## Statistical tests for table1

new213hab.df %>% group_by(r261, ismarine) %>% tally()
fisher.test(x = c(15,145), y = c(61,116))

f_scratch <- matrix(c(15,145,61,116), nrow = 2) %>% fisher.test()
f_scratch$p.value

marineFisherscore.df[213,] %>% View()
nonmarineFisherscore.df[213,] %>% View()
