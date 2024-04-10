library(dartVarietalID)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(patchwork)
library(data.table)
library(ggridges)
rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
                         pattern = ".rds",
                         full.names = T)
crop_names <- lapply(basename(rds_files),function(x){
  gsub(
    "_DAP.rds$",
    "",x)
})
crop_names <- unlist(crop_names)

het_list <- readRDS("het_list.rds")
het_list2 <- lapply(1:length(het_list),function(x){
  tmp <- het_list[[x]]
  tmp$crop <- crop_names[x]
  return(tmp)
})
het_list2 <- rbindlist(het_list2)
het_list2$stat <- "Diversity"
colnames(het_list2) <- c("val","crop","stat")


diff_list <- readRDS("diff_list.rds")
diff_list2 <- lapply(1:length(diff_list),function(x){
  tmp <- diff_list[[x]]
  tmp$crop <- crop_names[x]
  return(tmp)
})
diff_list2 <- rbindlist(diff_list2)
diff_list2$stat <- "Differentiation"
colnames(diff_list2) <- c("val","crop","stat")

corr_list <- readRDS("corr_list.rds")
corr_list2 <- rbindlist(corr_list)
corr_list2$stat <- "Correlation"
colnames(corr_list2) <- c("val","crop","stat")

purity_list <-readRDS("purity_list.rds")
purity_list <- lapply(1:length(purity_list),function(x){
  tmp <- as.data.frame(purity_list[[x]])
  tmp$crop <- crop_names[x]
  tmp$stat <- "Purity"
  colnames(tmp) <- c("val","crop","stat")
  return(tmp)
})
purity_list2 <- rbindlist(purity_list)
purity_list2 <- purity_list2[which(purity_list2$val < 1),]


fin_stats <- rbind(corr_list2,diff_list2,het_list2)
p_list <- as.list(1:10)

for(i in 1:length(rds_files)){

  t1 <- readRDS(rds_files[i])
  res_sum <- t1$res_summary

  r2_prob_corr <- round(summary(lm(res_sum$Probability.reference ~ res_sum$overlap))$r.squared,2)

p1 <-   ggplot(res_sum,aes(x= Probability.reference,
                         y= overlap/100))+
    geom_pointdensity()+
    scale_color_viridis() +
    geom_smooth()  +
    # labs(title=crop_names[[i]]) +
    # geom_text(aes( x=0.55,y=0.9,
    #           label = paste("R^2: ", r2_prob_corr,sep="")),
    #           show.legend=F,
    #           position = position_dodge(width=0.9),  size=10)+
    xlim(0.50,1)+
    ylim(0,1) +
  xlab("Probability") +
  ylab("Overlap") +
  theme_bw() +
  theme(plot.title = element_text(size=14)) +
  theme(legend.position="none")

p2 <-   ggplot(res_sum,aes(x= Probability.reference,
                           y= corr))+
  geom_pointdensity()+
  scale_color_viridis() +
  geom_smooth()  +
  # geom_text(aes( x=0.55,y=0.9,
  #                label = paste("R^2: ", r2_prob_corr,sep="")),
  #           show.legend=F,
  #          size=4
  #           )+
  xlim(0.50,1)+
  ylim(0,1) +
  xlab("Probability") +
  ylab("Correlation") +
  theme_bw() +
  theme(legend.position="none")

res_sum2 <- t1$res_summary
res_sum2$fst <-  fst_list[[i]]
res_sum_fin[[i]] <- res_sum2

p3 <- ggplot(res_sum2,aes(x= Probability.reference/100,
                          y= fst))+
  geom_pointdensity()+
  scale_color_viridis() +
  geom_smooth()  +

  xlim(0.50,1)+
  ylim(0,1) +
  xlab("Probability") +
  ylab("FST") +
  theme_bw()


# p3 <-   ggplot(res_sum,aes(y= overlap/100,
#                            x= corr))+
#   geom_pointdensity()+
#   scale_color_viridis() +
#   geom_smooth()  +
#   xlim(0,1)+
#   ylim(0,1) +
#   xlab("Overlap") +
#   ylab("Correlation") +
#   theme_bw() +
#   theme(legend.position="none")

p4 <- p1+p2+p3

p_list[[i]] <- p4

}

p_fin <- p_list[[1]]/p_list[[2]]/p_list[[3]]/p_list[[4]]/p_list[[5]]/p_list[[6]]/p_list[[7]]/p_list[[8]]/p_list[[9]]/p_list[[10]]
p_fin
ggsave("fin_stats.pdf",  width = 15, height = 50, units = "in", dpi="retina", bg = "transparent",limitsize = FALSE)

den_corr <- ggplot(corr_list2,aes(x= corr,fill=crop,y=crop))+
  geom_density_ridges()+
  theme_ridges(grid = T,
               center_axis_labels = TRUE,font_size = 18,font_family="Helvetica")+
  theme(legend.position = "none",strip.background = element_rect(colour="black",fill="white"))+
  labs(y = "", x = "Correlation")
den_corr
ggsave("correlation_distributions.pdf",  width = 8, height = 8, units = "in", dpi="retina", bg = "transparent"  )

den_diff <- ggplot(diff_list2,aes(x= diff,fill=crop,y=crop))+
  geom_density_ridges()+
  theme_ridges(grid = T, center_axis_labels = TRUE,font_size = 18,font_family="Helvetica")+
  theme(legend.position = "none",strip.background = element_rect(colour="black",fill="white"))+
  labs(y = "", x = "Differentiation")
den_diff
ggsave("diff_distributions.pdf",  width = 8, height = 8, units = "in", dpi="retina", bg = "transparent"  )

den_het <- ggplot(het_list2,aes(x= val,fill=crop,y=crop))+
  geom_density_ridges()+
  theme_ridges(grid = T, center_axis_labels = TRUE,font_size = 18,font_family="Helvetica")+
  theme(legend.position = "none",strip.background = element_rect(colour="black",fill="white"))+
  labs(y = "", x = "Diversity")
den_het
ggsave("het_distributions.pdf",  width = 8, height = 8, units = "in", dpi="retina", bg = "transparent"  )

den_pur <- ggplot(purity_list2,aes(x= val,fill=crop,y=crop))+
  geom_density_ridges()+
  theme_ridges(grid = T, center_axis_labels = TRUE,font_size = 18,font_family="Helvetica")+
  theme(legend.position = "none",strip.background = element_rect(colour="black",fill="white"))+
  labs(y = "", x = "Diversity")
den_pur
ggsave("het_distributions.pdf",  width = 8, height = 8, units = "in", dpi="retina", bg = "transparent"  )

den_fin <- ggplot(fin_stats,aes(x= val,fill=crop,y=crop))+
  geom_density_ridges()+
  theme_ridges(grid = T, center_axis_labels = TRUE,font_size = 18,font_family="Helvetica")+
  theme(legend.position = "none",strip.background = element_rect(colour="black",fill="white"))+
  labs(y = "", x = "") +
  facet_grid(~stat, scales = "free")
den_fin

ggsave("distributions.pdf",  width = 16, height = 8, units = "in", dpi="retina", bg = "transparent"  )




den_corr + den_diff + den_het

