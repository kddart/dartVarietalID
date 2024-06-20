library(ggplot2)
library(ggpointdensity)
library(viridis)
library(patchwork)
library(data.table)

rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
                         pattern = "_DAP_240514.rds",
                         full.names = T)
crop_names <- lapply(basename(rds_files),function(x){
  gsub(
    "_DAP_240514.rds$",
    "",x)
})
crop_names <- unlist(crop_names)

diff_list <- readRDS("fst_list_Neis_2.rds")
corr_list <- readRDS("corr_res.rds")
# overlap_list <- readRDS("overlap_res.rds")
# purity_list <- readRDS("purity_list.rds")
# purity_list[[7]] <- rep(0,length(overlap_list[[7]]))

p_list <- as.list(1:10)

for(i in c(1:10)){

  t1 <- readRDS(rds_files[i])
  # res_sum <- t1$res_summary
  res_sum <-  t1$res_full
  names_sam <- names(res_sum)
  res_sum <- lapply(1:length(res_sum), function(x){
    tmp <- res_sum[[x]]
    tmp <- tmp[1:10,]
    tmp$TargetID.sample <- names_sam[x]
    tmp$TargetID.reference <- tmp$variety
    return(tmp)
  })
  res_sum <- rbindlist(res_sum)

  res_sum$fst <- diff_list[[i]]
  res_sum$corr <- corr_list[[i]]

  # first_match <- seq(1, nrow(res_sum),10)
  # res_sum <- res_sum[first_match,]


  # res_sum$overlap <- overlap_list[[i]]
  # purity_list2 <- purity_list[[i]]
  # purity_list2 <- purity_list2[which(purity_list2 < 1)]
   res_sum$purity <- purity_list[[i]]
  # res_sum[which(res_sum$purity>99),"purity"] <- NA

  # p1 <-   ggplot(res_sum,aes(x= Probability.reference,
  #                            y= overlap/100))+
  #   geom_pointdensity()+
  #   scale_color_viridis() +
  #   geom_smooth()  +
  #   labs(title=crop_names[i]) +
  #   xlim(0.60,1)+
  #   ylim(0,1) +
  #   xlab("Probability") +
  #   ylab("Overlap") +
  #   theme_bw() +
  #   theme(plot.title = element_text(size=14)) +
  #   theme(legend.position="none")

  p2 <-   ggplot(res_sum,aes(x= Probability_corr_scaled,
                             y= corr))+
    geom_pointdensity()+
    scale_color_viridis() +
    geom_smooth()  +
    xlim(0.50,1)+
    ylim(0,1) +
    xlab("Probability_corr_scaled") +
    ylab("Correlation") +
    theme_bw() +
    theme(legend.position="none") +
    labs(title=crop_names[i])

  p3 <- ggplot(res_sum,aes(x= Probability_corr_scaled,
                           y= fst))+
    geom_pointdensity()+
    scale_color_viridis() +
    geom_smooth()  +
    xlim(0.50,1)+
    ylim(0,0.5) +
    xlab("Probability_corr_scaled") +
    ylab("Nei's distance") +
    theme_bw() +
    theme(legend.position="none")

  p4 <- ggplot(res_sum,aes(x= Probability_corr_scaled,
                           y= purity))+
    geom_pointdensity()+
    scale_color_viridis() +
    geom_smooth()  +
    xlim(0.50,1)+
    ylim(0.7,1.45) +
    xlab("Probability") +
    ylab("Purity score") +
    theme_bw()

  p5 <- p2+p3+p4

  p_list[[i]] <- p5
  # print(p5)

}

p_fin <- p_list[[1]]/p_list[[2]]/p_list[[3]]/p_list[[4]]/p_list[[5]]/p_list[[6]]/p_list[[7]]/p_list[[8]]/p_list[[9]]/p_list[[10]]
ggsave("fin_stats_7.png",  width = 10, height = 50, units = "in", dpi="retina", bg = "transparent",limitsize = FALSE)
