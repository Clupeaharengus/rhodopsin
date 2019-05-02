## RH1 dAF and chi2 calcs using refined population sets

#Load Libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(viridis)

#Define Functions
chiFunc <- function(x,y) {
  a1 <- x[1]
  a2 <- y[1]
  b1 <- x[2]
  b2 <- y[2]
  sumall <- a1+a2+b1+b2
  expectA1 <- (a1 + a2)*(a1 + b1)/sumall
  expectA2 <- (a1 + a2)*(a2 + b2)/sumall
  expectB1 <- (b1 + b2)*(a1 + b1)/sumall
  expectB2 <- (b1 + b2)*(a2 + b2)/sumall
  chi <- (a1 - expectA1)^2 / expectA1 + 
    (a2 - expectA2)^2 / expectA2 + 
    (b1 - expectB1)^2 / expectB1 + 
    (b2 - expectB2)^2 / expectB1
  return(chi)
}

dafFunc <- function(ref,alt) {
  ref_A <- ref[1]
  alt_A <- alt[1]
  ref_B <- ref[2]
  alt_B <- alt[2]
  refAFA <- ref_A/(ref_A+alt_A)
  altAFA <- alt_A/(ref_A+alt_A)
  refAFB <- ref_B/(ref_B+alt_B)
  altAFB <- alt_B/(ref_B+alt_B)
  refdeltaAF <- abs(refAFA - refAFB)
  return(refdeltaAF)
}

chi_daf_test <- function(x) {
  a <- x %>% 
    group_by(CHROM, POS) %>%
    mutate(refsum = sum(ref), altsum = sum(alt)) %>% 
    filter(refsum > 0 & altsum > 0) %>% 
    select(-altsum, -refsum) %>%
    summarise(chisq = chiFunc(ref,alt)) %>% 
    mutate(logpval = -log(pchisq(chisq,1,lower.tail=FALSE), base = 10))
  b <- x %>% 
    group_by(CHROM, POS) %>% 
    summarise(daf = dafFunc(ref,alt))
  y <- full_join(x, b) 
  z <- full_join(y,a)
  return(z)
}

groupsum <- function(dfin,group1,group2) {
  dfin %>%
    mutate(group = ifelse(grepl(paste(group1,collapse="|"), pop), "group1",
                          ifelse(grepl(paste(group2,collapse="|"), pop), "group2", NA))) %>% 
    group_by(group, CHROM, POS) %>% 
    summarise(ref = sum(ref), alt = sum(alt)) %>%  
    ungroup() %>% 
    arrange(CHROM, POS)
}

## Read in pools
allpools <- fread("rh1pools_BGI_s14.AD", header = T)

rh1poolnames <- fread("rh1_samplelist.txt", skip = 2, header = F, col.names = "pop")
rh1poolnamesvec <- scan("rh1_samplelist.txt", skip = 2, what = "character")

##Separate alleles
alleles_rh1 <- allpools %>%
  select(one_of("CHROM", "POS", rh1poolnamesvec)) %>% 
  gather("pop", "AD", -"CHROM", -"POS") %>% 
  separate(AD, c("ref","alt"), sep = ",", convert = T)

##Filter based on allele depth then sites
alleles_rh1_filt <- alleles_rh1 %>% 
  filter((ref+alt >= 10) & (ref + alt <= 100))
numberpops <- alleles_rh1_filt %>% 
  group_by(CHROM, POS) %>% 
  summarise(ns = n()) %>% 
  ungroup() %>% 
  arrange(-ns) %>% 
  slice(1) %>% 
  .[['ns']]

alleles_rh1_filt <- alleles_rh1_filt %>% 
  group_by(CHROM, POS) %>% 
  filter(n() > .7*numberpops) %>% 
  ungroup()

medcov <- alleles_rh1_filt %>% 
  summarise(med=median(ref+alt)) %>% 
  .[['med']]# 40

alleles_rh1_filt_downsamp <- alleles_rh1_filt %>%
  mutate(sampfactor = (ref+alt)/medcov) %>%
  mutate(ref = if_else(sampfactor > 1, ref/sampfactor, as.double(ref)),
         alt = if_else(sampfactor > 1, alt/sampfactor, as.double(alt)))

## Get order of populations using dAF 
secchi <- read.csv("secchi_rough_approx_rho_GENSINC_pool_long-3_190130.csv")
dAF_pops <- alleles_rh1_filt_downsamp %>%
  filter(POS == 11640748) %>% 
  mutate(F261Y = ref/(ref+alt)) %>% 
  select(pop, F261Y)
dAF_pops2 <- alleles_rh1_filt_downsamp %>%
  filter(POS == 11640892) %>% 
  mutate(I213T = 1-(ref/(ref+alt))) %>% 
  select(pop, I213T)
# pool_secchi <- rh1poolnames %>% 
#   left_join(secchi) %>% 
#   select(-AF) %>% 
#   left_join(dAF_pops) %>% 
#   left_join(dAF_pops2)
pool_secchi <- fread("rh1_poolAF_secchi_sal_abs_20190204.csv") %>% 
  mutate(pop = factor(pop,rh1_AF_order))
rh1_AF_order <- dAF_pops %>% 
  arrange(F261Y) %>% 
  .[['pop']]
gathered_pool_secchi <- pool_secchi %>% 
  select(pop,Salinity,secchi,Absorbance_412nm,lon,lat,F261Y,I213T) %>% 
  gather(site,dAF,7:8) %>%
  mutate(site = factor(site)) %>% 
  mutate(pop = factor(pop,rh1_AF_order))

## Generate dAF heatmap
p1 <- gathered_pool_secchi %>% 
  ggplot(aes(y = site, x = pop, fill = dAF)) +
  geom_tile(color='White', size = 0.1) +
  scale_fill_viridis(name="Baltic AF", option = "viridis", limits = c(0,1)) +
  scale_x_discrete(expand = c(0,0),
                   labels = 1:38) +
  scale_y_discrete(expand = c(0,0),
                   limits=c("I213T","F261Y")) +
  theme_minimal() +
  theme(text = element_text(size=24)) +
  theme(axis.title.x=element_blank(),
        # axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        # legend.position="none",
        axis.text.y=element_blank(),
        axis.text.x = element_text(color = "#000000"),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_line(size = 1),
        legend.title = element_blank(),
        legend.text = element_blank())

## Make salinity and secchi charts
# F261Y_AF_secchi_order <- pool_secchi %>% 
#   arrange(-F261Y, secchi) %>% 
#   .[['pop']]
pAbs <- pool_secchi %>%
  ggplot(aes(x = pop, y = Absorbance_412nm)) +
  # geom_step(group = 1, size = .75) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey70"),
        axis.text.y=element_text(size = 20, color = "#000000"),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90),
        axis.title.y=element_blank()) 
  
pSalinity <- pool_secchi %>%
  ggplot(aes(x = pop, y = Salinity)) +
  # geom_step(group = 1, size = .75) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey70"),
        # axis.text.y=element_text(size = 20, color = "#000000"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90),
        axis.title.y=element_blank()) 

dummypop <- pool_secchi %>% 
  filter(pop == "HGS15_NSSH_Atlantic_Spring") %>% 
  mutate(pop = "dummy")
secchi_sal_dummypop <- dummypop %>%
  bind_rows(pool_secchi)
dummy_order <- secchi_sal_dummypop %>% 
  arrange(-F261Y, secchi) %>% 
  .[['pop']]

p2 <- secchi_sal_dummypop %>%
  mutate(pop = factor(pop,dummy_order)) %>% 
  ggplot(aes(x = pop, y = secchi)) +
  geom_step(group = 1, size = .75) +
  theme_minimal() +
  theme(axis.title.y=element_blank(),
        panel.grid.major.x = element_line(color = "grey70"),
        panel.grid.minor.x = element_line(color = "grey80"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  coord_flip() +
  scale_y_continuous(breaks = c(3,6,9),
                     limits = c(3,9.2))

dummy_order <- secchi_sal_dummypop %>% 
  arrange(-F261Y, Salinity) %>% 
  .[['pop']]

p3 <- secchi_sal_dummypop %>%
  mutate(pop = factor(pop,dummy_order)) %>% 
  ggplot(aes(x = pop, y = Salinity)) +
  geom_step(group = 1, size = .75) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(color = "grey70"),
        panel.grid.minor.x = element_line(color = "grey80"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank()) +
  coord_flip() +
  scale_y_continuous(breaks = c(0,20,40),
                     limits = c(0,40))


## Define Groups
group1 <- pool_secchi %>% 
  filter(F261Y > 0 & I213T > 0) %>% #Filter to compare Atlantic vs Central Baltic
  filter(F261Y > 0) %>% #Filter based on salinity
  .[['pop']]
group2 <- pool_secchi %>% 
  filter(F261Y == 0) %>% 
  .[['pop']]

## Calculate allele depth sums for groups
groupsums <- groupsum(alleles_rh1_filt_downsamp, group1, group2) %>% 
  filter(!is.na(group))

#Calculate dAF and chi2 values
dAFchi <- groupsums %>%
  chi_daf_test()

#Plot chi
chitest <- dAFchi
# pdf("rh1_chi2.pdf", useDingbats=FALSE, width = 2.375, height = 1.58)
chitest %>% 
  # filter(POS > low2 & POS < hi2) %>% 
  ggplot() + 
  geom_point(aes(x = POS, y = logpval), alpha = .25) + 
  geom_vline(aes(xintercept = 11640748), color = "red", size = 1) +
  theme_bw() +
  theme(axis.title = element_blank(),
        # axis.text = element_blank(),
        legend.position="none") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,175))
# dev.off()

## Fixed SNP figure

fixedAtlantic <- function(ref, alt) {
  ref_1 <- ref[1]
  alt_1 <- alt[1]
  ref_2 <- ref[2]
  alt_2 <- alt[2]
  AtlFixedBool <- if_else((ref_1 > 0 && alt_1 > 0) && xor(ref_2 == 0,alt_2 == 0), T, 
                          if_else((ref_1 > 0 && alt_1 == 0) && (ref_2 == 0 && alt_2 > 0), T,
                                  if_else((ref_1 == 0 && alt_1 > 0) && (ref_2 > 0 && alt_2 == 0), T, F)))
  return(AtlFixedBool)
}

AtlanticFixedPos <- dAFchi %>% 
  group_by(POS) %>%
  summarise(Atlantic_fixed = fixedAtlantic(ref, alt))
AtlanticFixedPos <- AtlanticFixedPos %>% filter(Atlantic_fixed == T)

# Plot daf chi with color coded atlatic fixed alleles
chitest_atlfixed <- dAFchi %>% 
  right_join(AtlanticFixedPos)
low2 <- 9000000
hi2 <- 15000000
# pdf("rh1_dAF.pdf", useDingbats=FALSE, width = 2.375, height = 1.58)

chitest_atlfixed %>%
  filter(POS > low2 & POS < hi2) %>% 
  ggplot() + 
  geom_point(aes(x = POS, y = daf),
             alpha = .25, 
             na.rm = T) +
  geom_point(aes(x = POS, y = daf),
             color = "red",
             size = 3,
             na.rm = T,
             data = subset(chitest_atlfixed,POS == 11640748)) +
  geom_point(aes(x = POS, y = daf),
             color = "red",
             # shape = 17,
             size = 3,
             na.rm = T,
             data = subset(chitest_atlfixed,POS == 11640892)) +
  # geom_vline(aes(xintercept = 11640748), color = "red", size = 1, alpha = .5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        # axis.text = element_blank(),
        legend.position="none") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,1))
  

# dev.off()

chitest_atlfixed %>% 
  filter(POS > low2 & POS < hi2) %>% 
  ggplot() + 
  facet_wrap(~CHROM, scales = "free_x") +
  geom_point(aes(x = POS, y = logpval)) + 
  geom_vline(data = subset(chitest, CHROM == "Clupeapallasi.scafSeq.final_ovlk_hic_scaffold_14"),
             aes(xintercept = 11641000), color = "red", size = 2) +
  theme_minimal()
