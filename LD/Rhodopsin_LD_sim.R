#Rhodopsin LD breakdown simulation

plot(y = colSums(all_rho_LD[LD_lo_df$LD_idx, LD_lo_df$LD_idx], na.rm = T), x =LD_lo_df$SNP_HiC_pos, pch = 16, xlab = "Position", ylab = "Sum of LD across all individuals") # Deducing the sweep size
abline(v = c(9.3e6, 12.5e6), col = "red")
text(x = c(9.3e6, 12.5e6), y = 100, labels = c("9.3 Mb", "12.5 Mb"), pos = c(2,4))

#Analytical version
region_phys_size <- 3.2 #Size of the sweep reion in Mb
rec_rate <-  (56.0 - 49.7)/ (14.36 - 8.24) #Recombination rate in cM/Mb
p_rec_0 <- (1-rec_rate/100)^region_phys_size
plot(x = 0:100, y = p_rec_0^(0:100), pch = 16, ylab = "Probability of intact block", xlab = "Generations after mutation", main = "(p_zero_recombinations)^generations") #Probability of zero recombiations
abline(h  = 0.5, col = "red")
abline(v  = max(which(p_rec_0^(0:100) >= 0.5)), col = "blue")
text(x= 22, y = 0.8, labels = paste("last value >= 0.5:", max(which(p_rec_0^(0:100) >= 0.5)), "generations"), pos = 4)
      
#Simulation version
phys_size <- 32267647 #Size of chromsome 4 in the Ch_v2.0.2 assembly in bp
rho_pos <- 11217453 # Location of presumably causative SNP Rho 1 - see LD_lo_df["AX-99134661",]

#Generating a data frame with interval-based recombination rate  estimates
load("Ch_v2.0.2_linkage.RData")
chr4_linkage_map <- Ch_v2_linkage_map[Ch_v2_linkage_map$chr == 4,]

#Evening out the rates across each "step" of the map
chr4_linkage_map <- cbind(chr4_linkage_map[which(diff(chr4_linkage_map[,"all_map"]) > 0) +1 ,], diff(chr4_linkage_map[,"all_map"])[which(diff(chr4_linkage_map[,"all_map"]) > 0)])
chr4_linkage_map$HiC_pos[dim(chr4_linkage_map)[1]] <- phys_size #Stretching the fianl interval to the end of the chromosome
chr4_linkage_map <- chr4_linkage_map[,c(2,5,7)]
names(chr4_linkage_map)[3] <- "interval_dist"
chr4_linkage_map$interval_size <- diff(c(0,chr4_linkage_map$HiC_pos))

n_runs <- 1000
n_generations <- 200
block_size_df <- data.frame(dummy = 1:n_generations)
for (i in 1:n_runs) block_size_df[,i] <- simulate_linkage_breakdown(sim_linkage_map = chr4_linkage_map, focus_pos = rho_pos,n_gen = n_generations, draw_plot = F)$size

pdf("~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho_simulation_summary.pdf")
plot(x = 1:n_generations, y = rowSums(block_size_df)/n_runs/1e6, xlab = "Generations after mutation", ylab = "Average block size", pch = 16, col = "firebrick", main = "Results from 1000 simulations", type  ="n")
points(x = 1:n_generations, y = (apply(block_size_df,1, "summary")/1e6)[2,], pch = 16, col = "grey50")
points(x = 1:n_generations, y = (apply(block_size_df,1, "summary")/1e6)[3,], pch = 17, col = "darkorchid")
points(x = 1:n_generations, y = (apply(block_size_df,1, "summary")/1e6)[5,], pch = 18, col = "grey50")

abline(h  = 3.2, lwd  = 2, lty = "dashed")
#abline(v  = max(which(rowSums(block_size_df)/n_runs/1e6 >= 3.2)), col = "firebrick")
abline(v  = max(which(apply(block_size_df,1, "median")/1e6 >= 3.2)), col = "darkorchid")
abline(v  = max(which((apply(block_size_df,1, "summary")/1e6)[2,] >= 3.2)), col = "grey50")
abline(v  = max(which((apply(block_size_df,1, "summary")/1e6)[5,] >= 3.2)), col = "grey50")

#text(x= 36, y = 20., labels = paste("last value >= 3.2:", max(which(rowSums(block_size_df)/n_runs/1e6 >= 3.2)), "gen."), pos = 4)
#text(x= 22, y =25., labels = paste("last value >= 3.2:", max(which(apply(block_size_df,1, "median")/1e6 >= 3.2)), "gen."), pos = 4)
legend(x="topright",  legend = c("First quartile", "Median", "Third quartile"), pch = c(16, 17, 18), col = c("grey50", "darkorchid", "grey50"))
dev.off()

pdf("~/Projects/Herring/doc/Rhodopsin/LD_v2.0.2/Rho_simulation_example.pdf")
simulate_linkage_breakdown(sim_linkage_map = chr4_linkage_map, focus_pos = rho_pos,n_gen = n_generations, draw_plot = T)
#segments(x0 = c(0, chr4_linkage_map$HiC_pos[-length(chr4_linkage_map$HiC_pos)]), x1 = chr4_linkage_map$HiC_pos, y0 = (chr4_linkage_map$interval_dist/chr4_linkage_map$interval_size)*1e6*8, col = "grey50", lwd = 2, lty  = "dashed")
lines(x = c(0, chr4_linkage_map$HiC_pos), y = c(0, chr4_linkage_map$all_map), col = "grey50", lwd = 2, lty  = "dashed")
axis(4)
dev.off()

#Testing some s-values according to Fisher (delta_p = p*s)
p0 <- 1/1e6
s_array <- seq(from = 0, to = 1, by = 0.05)
s_screen_df <- data.frame(gen = 1:n_generations) 
for (s in s_array){
  #s <- 0.4
  p_vec <- numeric(n_generations)
  p_vec[1] <- p0
  for(i in 2:n_generations){
    w_avg <- p_vec[i-1] * (1+s) + (1-p_vec[i-1])
    s_i <- ((1+s) - w_avg)/w_avg
    p_vec[i] <- p_vec[i-1] + s_i*p_vec[i -1]
  }
  s_screen_df[,paste("s_", s , sep  ="")] <- p_vec
}
plot(x =0 , y = 0, xlim = c(0,200), ylim = c(0,1), xlab = "Generation", ylab = "Frequency", type = "n")
for(i in 2:dim(s_screen_df)[2]){
  points(y = s_screen_df[,i], x = 1:n_generations, pch = i %% 7, col = i %% 8, cex = 0.5, type = "b")
}
abline(h = 0.5)
abline(v  = max(which(apply(block_size_df,1, "median")/1e6 >= 3.2)), col = "darkorchid")
abline(v  = max(which((apply(block_size_df,1, "summary")/1e6)[2,] >= 3.2)), col = "grey50")
abline(v  = max(which((apply(block_size_df,1, "summary")/1e6)[5,] >= 3.2)), col = "grey50")
points(y = s_screen_df[,7], x = 1:n_generations, pch = 16, col = "steelblue", cex = 1, type = "b", lwd = 2.5)

#Support functions
simulate_linkage_breakdown <- function(sim_linkage_map, focus_pos, n_gen = 100 , draw_plot = F){
  rec_df <- data.frame(rec_pos = 1:n_gen, lower_bp = NA, upper_bp = NA, size = NA)
  lower_bp <- 0
  upper_bp <- max(sim_linkage_map$HiC_pos)
  
  rec_vec <- runif(n=n_gen, max = max(c(sim_linkage_map$all_map, 100)))
  for(i in 1:n_gen){
    if(rec_vec[i] <= max(sim_linkage_map$all_map)){
      rec_interval <- min(which(sim_linkage_map$all_map > rec_vec[i]))
      rec_pos <- sim_linkage_map$HiC_pos[rec_interval] - sample(0:sim_linkage_map$interval_size[rec_interval], 1)
      if(rec_pos < focus_pos & rec_pos > lower_bp) lower_bp <- rec_pos
      if(rec_pos > focus_pos & rec_pos < upper_bp) upper_bp <- rec_pos
    } else{
      rec_pos <- NA
    }
    rec_df[i,] <- c(rec_pos, lower_bp, upper_bp, upper_bp-lower_bp)
  }
  if(draw_plot){
    plot(x = c(0, phys_size), y = c(0, n_gen), xlab = "Position", ylab = "Generation", type = "n", axes = F)
    axis(1)
    axis(2, at = seq(from = 0, to = n_gen, by = 10), labels = rev(seq(from = 0, to = n_gen, by = 10)))
    segments(x0 = rec_df[,"lower_bp"], x1 = rec_df[,"upper_bp"], y0 = rev(1:n_gen), col = "blue", lwd = 2)
    points(x = rec_df[,"rec_pos"], y = rev(1:n_gen), pch = 4, col = "red")
  }
  return(invisible(rec_df))
}




