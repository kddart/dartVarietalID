library(dartVarietalID)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(dartR.base)
rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
                         pattern = ".rds",
                         full.names = T)
crop_names <- lapply(basename(rds_files),function(x){
  gsub(
    "_DAP.rds$",
    "",x)
})
crop_names <- unlist(crop_names)

diff_list <- as.list(1:10)
names(diff_list) <- crop_names

het_list <- as.list(1:10)
names(het_list) <- crop_names

for(i in 1:length(rds_files)){
t1 <- readRDS(rds_files[i])
test_pop_ref <- t1$gl.references
res_dif <- gl.report.fstat(test_pop_ref)
res_dif2 <- res_dif$Stat_matrices$Fstp
res_dif2[upper.tri(res_dif2)] <- NA
res_dif2 <- as.data.frame(as.table(as.matrix(res_dif2)))
res_dif2 <- as.data.frame(res_dif2$Freq)
colnames(res_dif2) <- "diff"
diff_list[[i]] <- res_dif2


res_het <- gl.report.heterozygosity(test_pop_ref,plot.display = F)
res_het2 <- as.data.frame(res_het$uHe)
colnames(res_het2) <- "het"
het_list[[i]] <- res_het2

}


corr_list <- as.list(1:10)
names(corr_list) <- crop_names
for(i in 1:length(rds_files)){

  t1 <- readRDS(rds_files[i])
  res_corr <- as.data.frame(t1$res_summary$corr)
  colnames(res_corr) <- "corr"
  res_corr$crop <- crop_names[i]
  corr_list[[i]] <- res_corr


}
saveRDS(corr_list,"corr_list.rds")



saveRDS(het_list,"het_list.rds")
saveRDS(diff_list,"diff_list.rds")

