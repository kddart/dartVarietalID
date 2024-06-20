library(dartVarietalID)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(patchwork)
library(data.table)

counts.file <- "/Users/mijangos/davis/Report_DCob22-7681_Counts.csv"
info.file <- "/Users/mijangos/davis/Info_file_davis_DCob22-7681.csv"

x <- runSampleAnalysis(
  counts.file = counts.file,
  info.file = info.file ,
  ncores = parallel::detectCores(),
  na.perc.threshold = 80
)

filename <- gsub("_counts.csv$", "_DAP", basename(counts.file))
filename.rds <- paste0(filename, "_DAP.rds")
saveRDS(x, filename.rds)
DAP_write_excel(x, info.n = info.file, filename)
#### HEATMAP REFERENCES
res.heatmap <- run.heatmap(rds.file=filename.rds)
saveRDS(res.heatmap,paste0(filename,"_heatmap.rds"))
#### CALCULATE PURITY
res.purity <- run.purity(
  rds.file = filename.rds,
  counts.file = counts.file,
  info.file = info.file,
  n.ref = 1,
  ncores = 12
)
saveRDS(res.purity,paste0(filename,"_TargetID_purity.rds"))
# res.purity <- readRDS(paste0(filename,"_TargetID_purity.rds"))
#### CALCULATE DIFFERENTIATION
res.diff <- run.diff(rds.file = filename.rds,
                     n.ref = 1)
saveRDS(res.diff,paste0(filename,"_diff.rds"))
# res.diff <- readRDS(paste0(filename,"_TargetID_diff.rds"))

#### CALCULATE CORRELATION
res.corr <- run.correlation(
  rds.file = filename.rds,
  counts.file = counts.file ,
  info.file = info.file,
  n.ref = 1
)
saveRDS(res.corr,paste0(filename,"_corr.rds"))
# res.corr <- readRDS(paste0(filename,"_TargetID_corr.rds"))
#### CALCULATE OVERLAP
res.overlap <- run.overlap(
  rds.file = filename.rds,
  counts.file = counts.file,
  info.file = info.file,
  n.ref = 1,
  ncores = 4,
  parallel = F
)

saveRDS(res.overlap,paste0(filename,"_overlap.rds"))
#### PLOTTING RESULTS
res.plot <- plot.stats(
  rds.file = filename.rds,
  corr.data = res.corr,
  diff.data = res.diff,
  overlap.data = NULL,
  purity.data = NULL,
  n.ref = 1,
  title_plot = "DCob22-7681"
)

# p_fin <- p_list[[1]]/p_list[[2]]/p_list[[3]]/p_list[[4]]/p_list[[5]]/p_list[[6]]/p_list[[7]]/p_list[[8]]/p_list[[9]]/p_list[[10]]
ggsave(paste0(filename,"DCob22-7681_stats.png"), plot = res.plot, width = 10, height = 5, units = "in", dpi="retina", bg = "transparent",limitsize = FALSE)

ggplot(res_sum,aes(x=Probability_corr_scaled,y=Het)) +
  geom_point()+
  xlab("Probability score") +
  ylab("Heterozygosity") +
  theme_bw()

