#Rhodopsin locus in forskarhjälpen
load("~/Projects/Herring/data/SNP_chip/Forskarhjalpen/SNP_map.Rdata")
fh_geno <- read.delim("~/Projects/Herring/data/SNP_chip/Forskarhjalpen/Calls_Final.txt", comment.char = "#",row.names = 1, stringsAsFactors = F)
fh_geno_all_original_snps <- fh_geno[rownames(fh_geno) %in% SNP_map[,"SNP_name"],]

Rhodopsin_SNP_map_raw<- SNP_map[SNP_map[,"scaffold"] == 341 & SNP_map[,"pos"] > 8e5 & SNP_map[,"pos"] < 12e5,]
Rhodopsin_SNP_map <- Rhodopsin_SNP_map_raw[order(Rhodopsin_SNP_map[,"pos"]),]

all_ind_names <- names(fh_geno_all_original_snps)
school_filter <- "PA234_([A-Z]+)[0-9]+_.+"
school_codes <- unique(sub(school_filter, "\\1", all_ind_names[grep(school_filter, all_ind_names)]))
School_loc <- read.table("~/Projects/Herring/data/SNP_chip/Forskarhjalpen/README.Schools.Table2.txt", stringsAsFactors=F, sep = "\t", header = T)

Rho_by_school <- Rhodopsin_SNP_map

for(sc in school_codes){
	Rho_by_school[,sc] <- sc_freq_v2(sc, marker_map=Rhodopsin_SNP_map)
}

tmp_rho_vec <- Rho_by_school[Rho_by_school[,"SNP_name"] == "AX-99089746",-(1:4)]
tmp_rho_vec <- tmp_rho_vec[match(School_loc[,"Kod"], names(Rho_by_school)[-(1:4)])]
School_loc[,"AX-99089746_freq"] <- t(tmp_rho_vec)


pdf(file = "~/Projects/Herring/doc/Rhodopsin/Rhodopsin_FH_heatmap.pdf", height = 10, width = 10)
par(xpd = NA)
cr_palette <- colorRampPalette(c("steelblue", "darkorchid3", "firebrick3"))
heatmap(as.matrix(Rho_by_school[,-(1:4)]), Rowv = NA, col = cr_palette(100), na.rm = T, scale = "none", cexRow = 0.15, labRow = Rho_by_school[,"SNP_name"])
#legend(y = par("usr")[4]*1.05, x= par("usr")[1], legend = c("Aut/Aut", "Aut/Spr", "Spr/Spr"), col = c("steelblue", "darkorchid3", "firebrick3"), pch = 15, pt.cex = 2)
dev.off()

require(maptools)
require(raster)
wrld_map <- readShapePoly("~/Projects/Herring/data/SNP_chip/SEC16B/10m-admin-0-countries/10m_admin_0_countries.shp")
cr <- colorRamp(c("blue", "red"))
pdf("~/Projects/Herring/doc/Rhodopsin/School_map.pdf", width = 8, height = 12)
plot(x = 15, y = 57, type = "n", xlim = c(5,22), ylim = c(54,66), xlab = "Longitude", ylab = "Latitude", main = "")
plot(wrld_map, add = T, col = "grey90", border = "grey70")
points(x= School_loc[,"Long"], y = School_loc[,"Lat"], col = rgb(cr(School_loc[,"AX-99089746_freq"]), maxColorValue=255), pch = 20, cex = 4)
text(x= School_loc[,"Long"], y = School_loc[,"Lat"], labels= School_loc[,"Kod"], pos = 3, cex = 1)
#text(x= s44_by_school[,"Long"], y = s44_by_school[,"Lat"], labels= paste("HWE: ", signif(-log10(s44_by_school[,"HWE_fisher"]), 2), sep  =""),pos = 1, cex = 0.4)
dev.off()



#LD calculations
#Re-seq data - moderately useful
#load("~/Projects/Herring/data/genomic_data/new_individual_SNP_set/Herring_67_HapDist.RData")
#load("~/Projects/Herring/data/v2.0.2_genotypes/salinity_selection_contrast.RData")
#h67_LD_data_v2.0.2 <- LD_GenABLE_export(herring_67$geno, sal_sel_df_Chv2)
#h67_LD_data_v2.0.2_clean <- h67_LD_data_v2.0.2
#h67_LD_data_v2.0.2_clean$GenABEL <- h67_LD_data_v2.0.2$GenABEL[!is.na(h67_LD_data_v2.0.2$ordered_df$HiC_chr),]
#h67_LD_data_v2.0.2_clean$ordered_df <- h67_LD_data_v2.0.2$ordered_df[!is.na(h67_LD_data_v2.0.2$ordered_df$HiC_chr),]
#rho_1 <- which(h67_LD_data_v2.0.2_clean$ordered_df[,"scaffold"] == "scaffold341" & h67_LD_data_v2.0.2_clean$ordered_df[,"pos"] == 1098026)
#rho_2 <- which(h67_LD_data_v2.0.2_clean$ordered_df[,"scaffold"] == "scaffold341" & h67_LD_data_v2.0.2_clean$ordered_df[,"pos"] == 1098170)
#h67_LD_data_v2.0.2_clean$ordered_df[c(rho_1, rho_2),]
#matchPattern("T", h67_LD_data_v2.0.2_clean$ordered_df[rho_1,"genotypes"])

#atl_rho_GenABEL <- h67_LD_data_v2.0.2_clean$GenABEL[c(1:30, 42:67),]
#balt_rho_GenABEL <- h67_LD_data_v2.0.2_clean$GenABEL[31:41,]

#balt_rho_LD <- chr_LD_calculation(4, geno_df = h67_LD_data_v2.0.2_clean$ordered_df , genabel_data = balt_rho_GenABEL, draw_plot = T, plot_prefix = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Baltic_rho_type_LD_chr")
#atl_rho_LD <- chr_LD_calculation(4, geno_df = h67_LD_data_v2.0.2_clean$ordered_df , genabel_data = atl_rho_GenABEL, draw_plot = T, plot_prefix = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Atlantic_rho_type_LD_chr")

load("~/Projects/Herring/data/HiC_assemblies/Version_2_release/satsuma_liftover/Ch.v2.0.2_v_Ilumina_satsuma.RData")
load("~/Projects/Herring/data/SNP_chip/Forskarhjalpen/fh_LD_all_original_SNPs_data.Rdata")

#Forskarhjälpen - more useful
Rhodopsin_SNP_map[Rhodopsin_SNP_map$pos == 1098026,]
Rhodopsin_SNP_map[Rhodopsin_SNP_map$pos == 1098170,]
rho_1_gt <- as.numeric(fh_LD_all_original_SNPs_data[,"AX-99134661"]@gtdata)
rho_1_gt[is.na(rho_1_gt)] <- 99
rho_2_gt <- as.numeric(fh_LD_all_original_SNPs_data[,"AX-99145279"]@gtdata)
rho_2_gt[is.na(rho_2_gt)] <- 99

double_baltic <- rho_1_gt == 2 & rho_2_gt == 2
double_mix <- rho_1_gt == 2 & rho_2_gt == 0
double_atlantic <- rho_1_gt == 0 & rho_2_gt == 0



#Need liftover like this
target_scaff_vec <- c(341, 634) #Contains the Rhodopsin locus
LD_marker_map <- as.data.frame(cbind(as.character(fh_LD_all_original_SNPs_data[,fh_LD_all_original_SNPs_data@gtdata@chromosome %in% target_scaff_vec]@gtdata@chromosome), fh_LD_all_original_SNPs_data[,fh_LD_all_original_SNPs_data@gtdata@chromosome %in% target_scaff_vec]@gtdata@map), stringsAsFactors = F)
names(LD_marker_map) <- c("chr", "pos")
LD_marker_map[, "chr"] <- paste("scaffold", LD_marker_map[, "chr"], sep = "")
LD_marker_map[, "pos"] <- as.numeric(LD_marker_map[, "pos"])
LD_marker_map[, "LD_idx"] <- 1:dim(LD_marker_map)[1]
LD_lo_df <- HiC_liftover(scaffold_data=LD_marker_map, liftover_df=Ch_v2.0.2_v_Ilu_satsuma, chr_size_df=Ch_v2.0.2_sizes, lo_cols="LD_idx")

### Implemented in function below
#Then LD like this
#tmp_double_balt_LD_data <- fh_LD_all_original_SNPs_data[double_baltic,fh_LD_all_original_SNPs_data@gtdata@chromosome %in% target_scaff_vec]
#double_balt_LD <- r2fast(tmp_double_balt_LD_data)
#double_balt_LD[lower.tri(double_balt_LD)] <- t(double_balt_LD)[lower.tri(double_balt_LD)]
#rm(tmp_double_balt_LD_data)
###

#And finally re-ordering like this
#image(x = LD_lo_df[,"SNP_HiC_pos"], y = LD_lo_df[,"SNP_HiC_pos"], z = 1-double_balt_LD[LD_lo_df[,"LD_idx"], LD_lo_df[,"LD_idx"]], main = "Rhodopsin LD; \"Double Baltic\" individuals")

#marker_subset_filter <- (fh_LD_all_original_SNPs_data@gtdata@chromosome == 341 & fh_LD_all_original_SNPs_data@gtdata@map < 2e6)
marker_subset_filter <- (fh_LD_all_original_SNPs_data@gtdata@chromosome %in% target_scaff_vec)
GenABEL_subset <- fh_LD_all_original_SNPs_data[,marker_subset_filter]
GenABEL_subset_ordered <- GenABEL_subset[,LD_lo_df$LD_idx]
GenABEL_subset_ordered@gtdata@map <- LD_lo_df$SNP_HiC_pos
GenABEL_subset_ordered@gtdata@chromosome <- factor(4)
GenABEL_subset_ordered <- GenABEL_subset_ordered[,GenABEL_subset_ordered@gtdata@map > 9e6 & GenABEL_subset_ordered@gtdata@map < 13e6]

all_rho_LD <- plot_subset_LD(marker_lo_map = LD_lo_df, ts_vec = target_scaff_vec, png_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/All_inds_FH_rhodopsin_LD.png", main_title = "Rhodopsin LD; All individuals", focal_pos = LD_lo_df["AX-99134661",]$SNP_HiC_pos)
rho_ld_decay_all <- plot.LD.decay(data = GenABEL_subset_ordered, dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/All_inds_FH_rhodopsin_LD_decay.pdf")

double_balt_rho_LD <- plot_subset_LD(marker_lo_map = LD_lo_df, ind_list = double_baltic, ts_vec = target_scaff_vec, png_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Double_baltic_inds_FH_rhodopsin_LD.png", main_title = "Rhodopsin LD; Double baltic individuals", focal_pos = LD_lo_df["AX-99134661",]$SNP_HiC_pos)
rho_ld_decay_double_balt <- plot.LD.decay(data = GenABEL_subset_ordered[double_baltic,], dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Double_baltic_inds_FH_rhodopsin_LD_decay.pdf")

double_atl_rho_LD <- plot_subset_LD(marker_lo_map = LD_lo_df, ind_list = double_atlantic, ts_vec = target_scaff_vec, png_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Double_atlantic_inds_FH_rhodopsin_LD.png", main_title = "Rhodopsin LD; Double atlantic individuals", focal_pos = LD_lo_df["AX-99134661",]$SNP_HiC_pos)
rho_ld_decay_double_balt <- plot.LD.decay(data = GenABEL_subset_ordered[double_atlantic,], dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Double_atlantic_inds_FH_rhodopsin_LD_decay.pdf")

double_mix_rho_LD <- plot_subset_LD(marker_lo_map = LD_lo_df, ind_list = double_mix, ts_vec = target_scaff_vec, png_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Double_mix_inds_FH_rhodopsin_LD.png", main_title = "Rhodopsin LD; Double mix individuals", focal_pos = LD_lo_df["AX-99134661",]$SNP_HiC_pos)
rho_ld_decay_double_mix <- plot.LD.decay(data = GenABEL_subset_ordered[double_mix,], dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Double_mix_inds_FH_rhodopsin_LD_decay.pdf")

rho_1_balt_LD <- plot_subset_LD(marker_lo_map = LD_lo_df, ind_list = (rho_1_gt == 2), ts_vec = target_scaff_vec, png_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho1_baltic_inds_FH_rhodopsin_LD.png", main_title = "Rhodopsin LD; Rho1 baltic individuals", focal_pos = LD_lo_df["AX-99134661",]$SNP_HiC_pos)
rho_1_balt_ld_decay <- plot.LD.decay(data = GenABEL_subset_ordered[(rho_1_gt == 2),], dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho1_baltic_inds_FH_rhodopsin_LD_decay.pdf")

rho_2_balt_LD <- plot_subset_LD(marker_lo_map = LD_lo_df, ind_list = (rho_2_gt == 2), ts_vec = target_scaff_vec, png_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho2_baltic_inds_FH_rhodopsin_LD.png", main_title = "Rhodopsin LD; Rho2 baltic individuals", focal_pos = LD_lo_df["AX-99134661",]$SNP_HiC_pos)
rho_2_balt_ld_decay <- plot.LD.decay(data = GenABEL_subset_ordered[(rho_2_gt == 2),], dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho2_baltic_inds_FH_rhodopsin_LD_decay.pdf")

rho_1_atl_LD <- plot_subset_LD(marker_lo_map = LD_lo_df, ind_list = (rho_1_gt == 0), ts_vec = target_scaff_vec, png_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho1_atlantic_inds_FH_rhodopsin_LD.png", main_title = "Rhodopsin LD; Rho1 atlantic individuals", focal_pos = LD_lo_df["AX-99134661",]$SNP_HiC_pos)
rho_1_atl_ld_decay <- plot.LD.decay(data = GenABEL_subset_ordered[(rho_1_gt == 0),], dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho1_atlantic_inds_FH_rhodopsin_LD_decay.pdf")

rho_2_atl_LD <- plot_subset_LD(marker_lo_map = LD_lo_df, ind_list = (rho_2_gt == 0), ts_vec = target_scaff_vec, png_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho2_atlantic_inds_FH_rhodopsin_LD.png", main_title = "Rhodopsin LD; Rho2 atlantic individuals", focal_pos = LD_lo_df["AX-99134661",]$SNP_HiC_pos)
rho_2_atl_ld_decay <- plot.LD.decay(data = GenABEL_subset_ordered[(rho_2_gt == 0),], dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho2_atlantic_inds_FH_rhodopsin_LD_decay.pdf")

#Defining the extent of the sweep
pdf("~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rhodposin_sweep_extent.pdf", width = 10)
  plot(y = colSums(all_rho_LD[LD_lo_df$LD_idx, LD_lo_df$LD_idx], na.rm = T), x = LD_lo_df$SNP_HiC_pos, pch = 16, ylab = "Sum r^2", xlab = "Position")
  abline(v = c(9.3e6, 12.5e6), col = "red", lwd = 2)
  rho_pos <- 11217453 # Location of presumably causative SNP Rho 1 - see LD_lo_df["AX-99134661",]
  abline(v = rho_pos, col = "grey50", lwd = 2)
dev.off()

#LD decay from the central SNPs
#Distances
#marker_dist <- as.matrix(dist(as.matrix(LD_lo_df$SNP_HiC_pos)))
#rho1_dist <- marker_dist[,which(rownames(LD_lo_df) == "AX-99134661")]
#rho1_dist[1:length(rho1_dist) < which(rownames(LD_lo_df) == "AX-99134661")] <- -rho1_dist[1:length(rho1_dist) < which(rownames(LD_lo_df) == "AX-99134661")] 
#all_rho_LD_phys_order <- all_rho_LD[LD_lo_df$LD_idx, LD_lo_df$LD_idx]
#rho1_nhood <- abs(rho1_dist) < 1e4
#rho1_bps <- (2^(1:11)) * 1e3
#rho1_bps <- c(min(rho1_dist), -rev(rho1_bps), rho1_bps, max(rho1_dist))
#rho1_idx <- cut(rho1_dist, rho1_bps, include.lowest = T)
#rho1_dist_mean_LD <- sapply(levels(rho1_idx), FUN = function(X){mean(all_rho_LD_phys_order[rho1_idx == X,rho1_nhood])})
#plot(x = (rho1_bps[-1] + rho1_bps[-length(rho1_bps)])/2, y = rho1_dist_mean_LD, xlim = c(-2e6, 2e6), type = "b", col = "red")

plot_focal_LD_decay("AX-99134661", rho_1_balt_LD, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho1_baltic_focal_LD_decay.pdf", main = "Homzygous \"Baltic\" at Rho 1", ylim = c(0,1), xlim = c(-2e6, 2e6))
plot_focal_LD_decay("AX-99134661", rho_2_balt_LD, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho2_baltic_focal_LD_decay.pdf", main = "Homzygous \"Baltic\" at Rho 2", ylim = c(0,1), xlim = c(-2e6, 2e6))
plot_focal_LD_decay("AX-99134661", rho_1_atl_LD, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho1_atlantic_focal_LD_decay.pdf", main = "Homzygous \"Atlantic\" at Rho 1", ylim = c(0,1), xlim = c(-2e6, 2e6))
plot_focal_LD_decay("AX-99134661", rho_2_atl_LD, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho2_atlantic_focal_LD_decay.pdf", main = "Homzygous \"Atlantic\" at Rho 2", ylim = c(0,1), xlim = c(-2e6, 2e6))
plot_focal_LD_decay("AX-99134661", all_rho_LD, pdf_file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/All_focal_LD_decay.pdf", main = "All individuals")

#Heterozygosity
pdf("~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rhodposin_locus_heterozygosity.pdf", width = 10)
Double_baltic_GTs <- as.numeric(GenABEL_subset_ordered[double_baltic,]@gtdata)
plot(y = colSums(Double_baltic_GTs == 1, na.rm = T)/dim(Double_baltic_GTs)[1], x = GenABEL_subset_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "Homozygous Baltic at Rho 1 & 2")
Double_atlantic_GTs <- as.numeric(GenABEL_subset_ordered[double_atlantic,]@gtdata)
plot(y = colSums(Double_atlantic_GTs == 1, na.rm = T)/dim(Double_atlantic_GTs)[1], x = GenABEL_subset_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "Homozygous Atlantic at Rho 1 & 2")
All_GTs <- as.numeric(GenABEL_subset_ordered@gtdata)
plot(y = colSums(All_GTs == 1, na.rm = T)/dim(All_GTs)[1], x = GenABEL_subset_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "All individuals")
Rho2_baltic_GTs <- as.numeric(GenABEL_subset_ordered[(rho_2_gt == 2),]@gtdata)
plot(y = colSums(Rho2_baltic_GTs == 1, na.rm = T)/dim(Rho2_baltic_GTs)[1], x = GenABEL_subset_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "Homozygous Baltic at Rho 2")
Rho1_baltic_GTs <- as.numeric(GenABEL_subset_ordered[(rho_1_gt == 2),]@gtdata)
plot(y = colSums(Rho1_baltic_GTs == 1, na.rm = T)/dim(Rho1_baltic_GTs)[1], x = GenABEL_subset_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "Homozygous Baltic at Rho 1")
Rho2_atlantic_GTs <- as.numeric(GenABEL_subset_ordered[(rho_2_gt == 0),]@gtdata)
plot(y = colSums(Rho2_atlantic_GTs == 1, na.rm = T)/dim(Rho2_atlantic_GTs)[1], x = GenABEL_subset_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "Homozygous Atlantic at Rho 2")
Rho1_atlantic_GTs <- as.numeric(GenABEL_subset_ordered[(rho_1_gt == 0),]@gtdata)
plot(y = colSums(Rho1_atlantic_GTs == 1)/dim(Rho1_atlantic_GTs)[1], x = GenABEL_subset_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "Homozygous Atlantic at Rho 1")
dev.off()





#LD decay-to-age according to Risch et al 1995
Double_baltic_freqs <- colSums(Double_baltic_GTs, na.rm = T)/(dim(Double_baltic_GTs)[1]*2)
Rho1_baltic_freqs <- colSums(Rho1_baltic_GTs, na.rm = T)/(dim(Rho1_baltic_GTs)[1]*2)
Double_atlantic_freqs <- colSums(Double_atlantic_GTs, na.rm = T)/(dim(Double_atlantic_GTs)[1]*2)
baltic_site_filter <- Double_atlantic_freqs < 0.05
plot( x = GenABEL_subset_ordered@gtdata@map[baltic_site_filter], y = Double_baltic_freqs[baltic_site_filter], pch  = 16, col = "red") #Isolating markers that are missing/nearly missing from the atalantic background
rec_rate <-  (56.0 - 49.7)/ (14.36 - 8.24) #Recombination rate in cM/Mb,
#phys_dist <-  abs(GenABEL_subset_ordered@gtdata@map[baltic_site_filter] -LD_lo_df["AX-99134661",]$SNP_HiC_pos) # Rho1 version
phys_dist <-  abs(GenABEL_subset_ordered@gtdata@map[baltic_site_filter] -LD_lo_df["AX-99145279",]$SNP_HiC_pos) # Rho2 version
theta_vec  <- phys_dist * (rec_rate/1e6/100)
#theta_vec[theta_vec < 0.001] <- NA #Removing markers very close to the target markers, for numerical stability

#delta_vec <- Double_baltic_freqs[baltic_site_filter] # Both mutations
#pdf("~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho2_age_estimation.pdf")

delta_vec <- Rho1_baltic_freqs[baltic_site_filter] #Rho1 only
pdf("~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho1_age_estimation.pdf")

delta_vec[delta_vec < 0.2] <- NA #Removing presumably new, derived variants

par(mar = c(5,4,3,4))
#plot( x = GenABEL_subset_ordered@gtdata@map[baltic_site_filter], y = delta_vec, pch  = 16, col = "red")
g_vec <- (log10(delta_vec)/log10(1-theta_vec))^2
#g_vec <- (log10(delta_vec[-length(delta_vec)])/log10(1-diff(theta_vec))) * (log10(delta_vec[-1])/log10(1-diff(theta_vec))) # all neighbouring pairs
#plot( x = GenABEL_subset_ordered@gtdata@map[baltic_site_filter], y = g_vec, pch  = 16, col = "red", ylim = c(0, 5e4), xlim = c(1.11e7, 1.18e7)) #Right side
plot( x = GenABEL_subset_ordered@gtdata@map[baltic_site_filter], y = g_vec, pch  = 16, col = "red", ylim = c(0, 1e5), ylab = "Generations", xlab = "Position")
axis(4, at = seq(from = 0, to = 1e5, by = 1e4), labels = seq(from = 0, to = 1e5, by = 1e4)/1e5)
mtext("Major allele frequency", 4,2, las = 0)
points( x = GenABEL_subset_ordered@gtdata@map[baltic_site_filter], y = delta_vec*1e5, pch  = 17, col = "darkorchid")
#points( x = Ch_v2_recombination_profile[Ch_v2_recombination_profile[,1] == 4, "BIN_START"], y = Ch_v2_recombination_profile[Ch_v2_recombination_profile[,1] == 4, "RR"]*5e3, pch  = 18, col = "grey50")
summary(g_vec[ GenABEL_subset_ordered@gtdata@map[baltic_site_filter] > 1.11e7 & GenABEL_subset_ordered@gtdata@map[baltic_site_filter] < 1.18e7])
abline(h = summary(g_vec[ GenABEL_subset_ordered@gtdata@map[baltic_site_filter] > 1.11e7 & GenABEL_subset_ordered@gtdata@map[baltic_site_filter] < 1.18e7])[3], col = "grey50")
#summary(g_vec[ GenABEL_subset_ordered@gtdata@map[baltic_site_filter] > 1.06e7 & GenABEL_subset_ordered@gtdata@map[baltic_site_filter] < 1.13e7])
abline(v = rho_pos, col = "grey50", lwd = 2)
dev.off()
#plot( x = GenABEL_subset_ordered@gtdata@map[baltic_site_filter][-1], y = g_vec, pch  = 16, col = "red")


#Table across schools/locations for the 2 key mutations, using GenABLE object for consistency
Rho_by_school_GenABEL <- School_loc
Rho_by_school_GenABEL <- Rho_by_school_GenABEL[-15,] #Removing duplicate entry for Lommarskolan
Rho_by_school_GenABEL[,"Rho1_freq"] <- NA
Rho_by_school_GenABEL[,"Rho2_freq"] <- NA
Rho_by_school_GenABEL[,"r^2"] <- NA
Rho_by_school_GenABEL[,"D'"] <- NA
for(sc in school_codes){
  sc_inds <- grep(sc, names(fh_geno_all_original_snps))
  GenA_subset <- fh_LD_all_original_SNPs_data[sc_inds,c("AX-99134661", "AX-99145279")]
  GenA_subset_GTs <- as.numeric(GenA_subset@gtdata)
  Rho_by_school_GenABEL[Rho_by_school_GenABEL[,"Kod"] == sc, c("Rho1_freq", "Rho2_freq")] <- colSums(GenA_subset_GTs, na.rm = T)/(dim(GenA_subset_GTs)[1]*2)
  Rho_by_school_GenABEL[Rho_by_school_GenABEL[,"Kod"] == sc, c("r^2")] <- r2fast(GenA_subset)[1,2]
  Rho_by_school_GenABEL[Rho_by_school_GenABEL[,"Kod"] == sc, c("D'")] <- dprfast(GenA_subset)[1,2]
}

school_n_ind <- array()
for(sc in school_codes){
  sc_inds <- grep(sc, names(fh_geno_all_original_snps))
  Rho_by_school_GenABEL[Rho_by_school_GenABEL[,"Kod"] == sc,"n_ind"]  <- length(sc_inds)
}



Rho_by_school_GenABEL[,9:12] <- round(Rho_by_school_GenABEL[,9:12], 2)
Rho_by_school_GenABEL[,"IY_freq"] <- Rho_by_school_GenABEL[,"Rho1_freq"] -  Rho_by_school_GenABEL[,"Rho2_freq"]
Rho_by_school_GenABEL[,"TY_freq"] <- Rho_by_school_GenABEL[,"Rho2_freq"]
Rho_by_school_GenABEL[,"IF_freq"] <- 1 - Rho_by_school_GenABEL[,"Rho1_freq"]

school_order_vec <- Rho_by_school_GenABEL$Lat
school_order_vec[is.na(school_order_vec)] <- -99
school_order_vec[which(Rho_by_school_GenABEL$Long > 13)] <- -school_order_vec[which(Rho_by_school_GenABEL$Long > 13)]
school_idx <- order(school_order_vec, decreasing = T)
write.table(Rho_by_school_GenABEL[school_idx,c(1,5,6,15,13,14,11:12)], file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho_freq_by_school.txt", quote = F, row.names = F, sep = "\t")
write.table(Rho_by_school_GenABEL[school_idx,c("Plats","n_ind")], file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/N_ind_by_school.txt", quote = F, row.names = F, sep = "\t")
write.table(Rho_by_school_GenABEL[school_idx,c("Plats","Datum")], file = "~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/date_by_school.txt", quote = F, row.names = F, sep = "\t")

#Extracting individuals from sampling locations in the Atlantic
atl_school_inds <- integer()
for(sc in Rho_by_school_GenABEL[school_idx,"Kod"][1:7]){
  sc_inds <- grep(sc, names(fh_geno_all_original_snps))
  atl_school_inds <- c(atl_school_inds, sc_inds)
}

#Checking TSHR to find autumn spawners
SNP_map[SNP_map$scaffold == 1420 & SNP_map$pos %in% c(133030, 108994, 109295,100524,137485),]
#SNP_name scaffold    pos                 raw_name
#52869 AX-99110137     1420 133030 AX-99110137_scaffold1420
TSHR_by_school_GenABEL <- School_loc
TSHR_by_school_GenABEL <- TSHR_by_school_GenABEL[-15,] #Removing duplicate entry for Lommarskolan
TSHR_by_school_GenABEL[,"AX_99110137_freq"] <- NA
TSHR_by_school_GenABEL[,"AX_99030516_freq"] <- NA
TSHR_by_school_GenABEL[,"AX_99034443_freq"] <- NA
TSHR_by_school_GenABEL[,"AX_99136256_freq"] <- NA
TSHR_by_school_GenABEL[,"AX_99143205_freq"] <- NA

for(sc in school_codes){
  sc_inds <- grep(sc, names(fh_geno_all_original_snps))
  GenA_subset <- fh_LD_all_original_SNPs_data[sc_inds,c("AX-99110137", "AX-99030516", "AX-99034443", "AX-99136256", "AX-99143205")]
  GenA_subset_GTs <- as.numeric(GenA_subset@gtdata)
  TSHR_by_school_GenABEL[TSHR_by_school_GenABEL[,"Kod"] == sc, c("AX_99110137_freq","AX_99030516_freq", "AX_99034443_freq", "AX_99136256_freq", "AX_99143205_freq")] <- colSums(GenA_subset_GTs, na.rm = T)/(dim(GenA_subset_GTs)[1]*2)
}
TSHR_by_school_reorder <- TSHR_by_school_GenABEL[school_idx,c(1:2, 5:6,9,10,11,13)]
rownames(TSHR_by_school_reorder) <- 1:dim(TSHR_by_school_reorder)[1]


pdf("~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rhodposin_locus_heterozygosity_atl_schools.pdf", width = 10)
atl_school_GTs <- as.numeric(GenABEL_subset_ordered[atl_school_inds,]@gtdata)
plot(y = colSums(atl_school_GTs == 1, na.rm = T)/dim(atl_school_GTs)[1], x = GenABEL_subset_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "Samples from Atalntic locations")
dev.off()

# Support functions
sc_freq_v2 <- function(sc, geno_data = fh_geno_all_original_snps, marker_map){
  sc_inds <- grep(sc, names(geno_data))
  #ts_map <- marker_map[marker_map[,"scaffold"] == target_scaff,]
  #map_order <- order(ts_map[,"pos"])
  SNP_genos <- geno_data[marker_map[,"SNP_name"], sc_inds]
  #ts_fh_SNPs <- ts_fh_SNPs[map_order,]
  SNP_genos[SNP_genos == -1] <- NA
  sc_alt_freq <- rowSums(SNP_genos, na.rm=T) / (rowSums(!is.na(SNP_genos))*2)
  return(sc_alt_freq)
}

plot_subset_LD <- function(marker_lo_map, ts_vec,  ind_list = NULL, snp_data = fh_LD_all_original_SNPs_data, png_file, main_title = "", focal_pos = NULL){
  if(!is.null(ind_list)){
    tmp_LD_data <- snp_data[ind_list,snp_data@gtdata@chromosome %in% ts_vec]
  } else {
    tmp_LD_data <- snp_data[,snp_data@gtdata@chromosome %in% ts_vec]
  }
  sub_LD <- r2fast(tmp_LD_data)
  sub_LD[lower.tri(sub_LD)] <- t(sub_LD)[lower.tri(sub_LD)]
  png(png_file, width = 1000, height = 1000)
  image(x = marker_lo_map[,"SNP_HiC_pos"], y = marker_lo_map[,"SNP_HiC_pos"], z = 1-sub_LD[marker_lo_map[,"LD_idx"], marker_lo_map[,"LD_idx"]], main = main_title, xlab = "Position", ylab = "Position")
  if(!is.null(focal_pos)){
    abline(v = focal_pos)
    abline(h = focal_pos)
  }
  dev.off()
  return(sub_LD)
}

plot_focal_LD_decay <- function(focal_marker, LD_matrix, LD_map = LD_lo_df, pdf_file = NULL, ...){
  marker_dist <- as.matrix(dist(as.matrix(LD_map$SNP_HiC_pos)))
  focal_dist <- marker_dist[,which(rownames(LD_map) == focal_marker)]
  focal_dist[1:length(focal_dist) < which(rownames(LD_map) == focal_marker)] <- -focal_dist[1:length(focal_dist) < which(rownames(LD_lo_df) == focal_marker)] 
  LD_phys_order <- LD_matrix[LD_map$LD_idx, LD_map$LD_idx]
  focal_nhood <- abs(focal_dist) < 1e4
  focal_bps <- (2^(1:11)) * 1e3
  focal_bps <- c(min(focal_dist), -rev(focal_bps), focal_bps, max(focal_dist))
  focal_idx <- cut(focal_dist, focal_bps, include.lowest = T)
  focal_dist_mean_LD <- sapply(levels(focal_idx), FUN = function(X){mean(LD_phys_order[focal_idx == X,focal_nhood])})
  if(!is.null(pdf_file)) pdf(pdf_file)
  plot(x = (focal_bps[-1] + focal_bps[-length(focal_bps)])/2, y = focal_dist_mean_LD, type = "b", col = "red", xlab = "Distance from focal SNP", ylab = "Mean LD (r^2)", ...)
  if(!is.null(pdf_file)) dev.off()
}

##' @title Plots LD decay pattern for a given data
##' 
##' @description Function plots the LD decay with distance between markers for a given data. 
##' LD is measured using r2 statistics.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param data a  \code{\link[GenABEL]{gwaa.data-class}} object 
##' @param N number of bins used to partition the distance between markers
##' @param dmin minimal distance to consider
##' @param maximal distance to consider
##' @param ... the remaining parameters that pass to the plot function
##' @return NULL
##' @examples \dontrun{
##' dat <- dat[,dat@@gtdata@@chromosome == 3]
##' plot.LD.decay(data=dat, N=200, dmin=0, dmax=1e6, main="Chr3")
##' }
##' @keywords LD, LD decay, plot, r2
##' @export plot.LD.decay
plot.LD.decay <- function(data, N=200, dmin=NA, dmax=NA, pdf_file = NULL, dm = NULL, ...){
  m <- as.matrix(data@gtdata@map)
  if(is.null(dm)) dm <- as.matrix(dist(m, diag=T))
  dm <- dm[upper.tri(dm)]
  r2m <- r2fast(data)
  r2m <- r2m[upper.tri(r2m)]
  # filtering
  if (!is.na(dmin) & !is.na(dmax)) {
    toKeep <- which(dm < dmax & dm > dmin) 
    dm <- dm[toKeep]
    r2m <- r2m[toKeep]
  }
  bpts <- pretty(dm, n=N)
  INDEX <- cut(dm, bpts, include.lowest = T)
  pts <- tapply(r2m, INDEX, mean)
  pts2 <- tapply(r2m, INDEX, median)
  dist <- bpts[1:length(pts)]
  if(!is.null(pdf_file)) pdf(pdf_file)
  plot(dist, pts, cex=.5, col='olivedrab', xlab="distance in bp", ylab="r2", type='n', las=1, bty='n', ...)
  grid()
  abline(v=seq(1,max(dist), by=1e6), lty=3, col="grey")
  points(dist, pts, type='l', cex=.5, col='olivedrab')
  points(dist, pts2, type='l', cex=.5, col='tomato')
  legend("topright", bty='n', legend=c("mean", "median"), col=c("olivedrab","tomato"), lty=1)
  if(!is.null(pdf_file)) dev.off()
  return(list(d = dist, p1 = pts, p2 = pts2))
}

HiC_liftover <- function(scaffold_data, liftover_df, chr_size_df, lo_cols){
  require(GenomicRanges)
  #Selection_GR <- GRanges(seqnames = selection_data[,1],ranges = IRanges(start = selection_data[,2], end = selection_data[,2]))
  scaffold_GR <- GRanges(seqnames = scaffold_data[,"chr"],ranges = IRanges(start = scaffold_data[,"pos"], end = scaffold_data[,"pos"]))
  #BGI_v_Ilu_df
  liftover_GR <- GRanges(seqnames= liftover_df[,"Ilu_seqnames"],ranges=IRanges(start = as.numeric(liftover_df[,"Ilu_start"]), end = as.numeric(liftover_df[,"Ilu_end"])))
  
  #Matching SNP positions with entries in the BGI to Ilu satsuma alignment
  scaffold_matches <- findOverlaps(query = scaffold_GR, subject = liftover_GR, type = "any" )
  if(class(lo_cols) == "numeric") lo_cols <- names(scaffold_data)[lo_cols]
  lo_df <- cbind(scaffold_data[scaffold_matches@from,c("chr","pos",lo_cols)], liftover_df[scaffold_matches@to,])
  
  #Adjusting positions within each matched interval
  target_SNPs <-  sign(lo_df$direction_est) == 1
  lo_df[target_SNPs,"SNP_HiC_pos"] <- lo_df[target_SNPs,"start"] + (lo_df[target_SNPs,"pos"] - lo_df[target_SNPs,"Ilu_start"])
  target_SNPs <-  sign(lo_df$direction_est) != 1
  lo_df[target_SNPs,"SNP_HiC_pos"] <- lo_df[target_SNPs,"start"] + (lo_df[target_SNPs,"Ilu_end"] - lo_df[target_SNPs,"pos"])
  
  #Add a cumulative position
  #HiC_chr_vec <- chromosome_list[order(as.numeric(sub("hic_scaffold_","",chromosome_list)))]
  pos_adj_vec <- c(0,cumsum(chr_size_df[,2]))
  pos_adj_df <- cbind(chr_size_df,pos_adj_vec[-length(pos_adj_vec)])
  lo_df[,"SNP_cumulative_pos"] <- lo_df[,"SNP_HiC_pos"] + pos_adj_df[match(lo_df[,"seqnames"], pos_adj_df[,1]),3]
  chr_col_vec <- c("grey30", "grey70")[match(lo_df[,"seqnames"], pos_adj_df[,1]) %% 2 + 1]
  lo_df[,"col"] <- chr_col_vec
  
  lo_df <- lo_df[!is.na(lo_df[,"SNP_cumulative_pos"]),]
  
  lo_df <- lo_df[order(lo_df[,"SNP_cumulative_pos"]),]
  return(lo_df)
}





