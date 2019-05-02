library(tidyr)
library(data.table)
library(ggplot2)
library(dplyr)

## Calculate allele frequency significance from rabodend phylogeny
## Load final species list / alleles / habitat
rabo_allele_hab <- read.table(file = "rabo_allele_hab.tsv", header = T)

## Table 1 - Allele vs habitat sifnificance test
rabo_allele_hab %>% filter(r261 == "Y" | r261 == "F") %>% select(r261,ismarine) %>% group_by(r261, ismarine) %>% tally()
fisher261 <- fisher.test(tibble(marine = c(875,10), nonmarine =c(815,356)))
fisher261$p.value

rabo_allele_hab %>% filter(r213 == "I" | r213 == "T") %>% select(r213,ismarine) %>% group_by(r213, ismarine) %>% tally()
fisher213 <- fisher.test(tibble(marine = c(573,182), nonmarine =c(1033,77)))
fisher213$p.value

## Table S2 - Co-occurance of 213 alleles with 261
rabo_allele_hab %>% filter(r261 == "Y" | r261 == "F") %>% select(r261,r213) %>% group_by(r261,r213) %>% summarise(n = n())


## Deprecated
fishall <- fread(file = "allfish_alleles.tsv", header = F, col.names = c("accession","spp", "allele")) %>% 
  extract(allele, into=c("r213", "r261"), '(.)(.)')

fish_u <- fishall %>% 
  distinct(spp, r213, r261, .keep_all = T)

r261Y <- fish_u %>%
  filter(r261 == "Y") %>% 
  count(r213) %>% 
  mutate(r261_Y =n) %>% 
  select(-n)

r261F <- fish_u %>%
  filter(r261 == "F") %>% 
  count(r213) %>% 
  mutate(r261_F =n) %>% 
  select(-n)

all_tab <- full_join(r261F, r261Y, by = "r213")
all_tab[is.na(all_tab)] <- 0
summarise(all_tab, tot = sum(r261_Y)+sum(r261_F)) #Sanity check 1146 OK

all_tab_Xsq <- all_tab %>% 
  filter(r213 == "I" | r213 == "T") %>% 
  select(-r213) %>% 
  chisq.test()

aa.glm <- glm(cbind(r261_F, r261_Y) ~ r213, data = all_tab, family = binomial())
anova(aa.glm)
lsmeans(aa.glm, pairwise ~ r213)

## Using the fish_hab generated in the cladogram.R script compare frequencies of habitat with amino acid identity
hab_261 <-
  fish_hab %>% 
  mutate(habitat = ifelse(marine == "T" & freshwater == "T" & brackish == "T", "marine;freshwater;brackish",
                          ifelse(marine == "F" & freshwater == "T" & brackish == "T", "freshwater;brackish",
                                 ifelse(marine == "T" & freshwater == "F" & brackish == "T", "marine;brackish",
                                        ifelse(marine == "F" & freshwater == "F" & brackish == "T", "brackish",
                                               ifelse(marine == "T" & freshwater == "F" & brackish == "F", "marine",
                                                      ifelse(marine == "F" & freshwater == "T" & brackish == "F", "freshwater", NA))))))) %>% 
  group_by(r261, habitat) %>% 
  tally() %>% 
  spread(r261,n)

write.table(hab_261, file = "hab261.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
                
         