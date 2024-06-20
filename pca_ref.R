library(dartR)

rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
                         pattern = "self.rds",
                         full.names = T)
pcoa_list <- as.data.frame(matrix(nrow = 10,ncol = 11))
colnames(pcoa_list) <- c("Dataset",paste0("Axis_",1:10))

crop_names <- lapply(basename(rds_files),function(x){
  gsub(
    "_DAP_self.rds$",
    "",x)
})

for(i in 1:length(rds_files)){
  print(crop_names[[i]])

  t1 <- readRDS(rds_files[i])
  t2 <- t1$gl.references
  class(t2) <- "genlight"
  pcoa <- gl.pcoa(t2,
                  nfactors = 10)

  s <- sum(pcoa$eig[pcoa$eig >= 0])
  e <- round(pcoa$eig * 100 / s, 1)

  pcoa_list[i,] <- c(crop_names[[i]], e[1:10])
  saveRDS(pcoa,paste0(crop_names[[i]],"_pca.rds"))

}

write.csv(pcoa_list,"pcoa_list.csv")
library(reshape2)
pcoa_list2 <- read.csv("pcoa_list.csv")
pcoa_list2 <- pcoa_list2[,-1]
colnames(pcoa_list2) <- c("Dataset",1:10)
pcoa_list2$mean <- rowMeans(pcoa_list2[,2:11])

pcoa_list3 <- melt(pcoa_list2[,1:11] ,id.vars = c("Dataset"))
# str(pcoa_list2)
pcoa_list3$variable <- as.numeric(pcoa_list3$variable)
pcoa_list3$value <- as.numeric(pcoa_list3$value)

library(dartR.base)
ggplot(pcoa_list3,aes(x = variable, y = value, color = Dataset)) +
  geom_hline(data =pcoa_list2, aes(yintercept = mean),size =1 )+
  geom_line() +
  geom_point() +
  facet_wrap(~Dataset,scales = "free_y") +
  theme(legend.position = "none") +
  labs(x = "Dimensions"  ,y= "eigenvalues (proportional to variance explained)"  ) +
  scale_x_continuous(breaks = seq(1,10, 1))




