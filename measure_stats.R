library(data.table)
library(ggplot2)
library(patchwork)
library(ggpointdensity)
library(viridis)

 t1 <- sam_pops_sep[[4]]
 t1 <- as.matrix(t1)
 t2 <- t1[,1:10]
geno_table <- apply(t2,2,table, useNA ="always")
geno_prop <- Reduce(rbind,geno_table)

cdata<- ref_sam$counts
cdata <- cdata[1:20,1:10]
nsnp <- nrow(cdata) / 2
npop <- ncol(cdata)
esl <- seq(2, nrow(cdata), 2)
# calculating allele frequencies
allele_1 <- cdata[esl,]
allele_2 <- cdata[-esl,]
p_freq <- allele_1 / (allele_1 + allele_2)


het_ref <- lapply(ref_pops_sep,gl.Ho)
het_ref <- lapply(het_ref,mean,na.rm =T)
het_ref_df <- data.frame(variety=names(het_ref),
                         het= unlist(unname(het_ref)))

merge_het <- lapply(1:27,function(x){
  tmp <- merge(res_tmp[[x]],het_ref_df,by="variety")
  tmp <- tmp[order(tmp$Probability,decreasing = T),]
  tmp <- tmp[1,]
  return(tmp)
})

rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
                         pattern = "_DAP_ref_mac10.rds",
                         full.names = T)
fin_list <- as.list(1:10)

crop_names <- lapply(basename(rds_files),function(x){
  gsub(
    "_DAP_ref_mac10$",
    "",x)
})
for(i in 1:length(rds_files)){

  t1 <- readRDS(rds_files[i])
  full <- t1$res_full
  res_sum <- t1$res_summary
  res_sum <- res_sum[,c("TargetID.sample",
                        "TargetID.reference",
                        "Probability",
                        "Probability_corr",
                        "Probability_scaled",
                        "Probability_corr_scaled",
                        "diversity")]
  colnames(res_sum) <- c("IDsam",
                         "IDref",
                         "Probability",
                         "Probability_corr",
                         "Probability_scaled",
                         "Probability_corr_scaled",
                         "Het")
  res_sum$set <- "match"
  res_sum$pair <- paste0(res_sum$IDsam,"_",res_sum$IDref)
  full_names <- names(full)
  full <- lapply(1:length(full),function(x){
    tmp_full <- full[[x]]
    tmp_full$IDsam <- full_names[x]
    return(tmp_full)
  })
  full_bind <- rbindlist(full)
  full_bind <- full_bind[,c("IDsam",
                            "TargetID",
                            "Probability",
                            "Probability_corr",
                            "Probability_scaled",
                            "Probability_corr_scaled",
                            "Het")]
  colnames(full_bind) <- c("IDsam",
                           "IDref",
                           "Probability",
                           "Probability_corr",
                           "Probability_scaled",
                           "Probability_corr_scaled",
                           "Het")
  full_bind$set <- "all"
  full_bind$pair <- paste0(full_bind$IDsam,"_",full_bind$IDref)

  # full_bind_2 <- rbind(res_sum,full_bind)

  full_bind$dataset <- crop_names[[i]]

  fin_list[[i]] <- full_bind

}

fin_list2 <- rbindlist(fin_list)
fin_list2$type <- "other"
fin_list2[which(fin_list2$IDsam == fin_list2$IDref),"type"] <- "self"
fin_list2$type <- as.factor(fin_list2$type)
fin_list2$dataset <- as.factor(fin_list2$dataset)

################################################################
################################################################
################################################################
Probability <- ggplot(fin_list2) +
  geom_point(data = fin_list2[fin_list2$type == "other",],
             aes(x= Het, y =Probability),color= "gray") +
  geom_smooth(
    aes(x = Het, y =Probability),
    color = "darkseagreen4",
    data = fin_list2[fin_list2$type == "other",],
    se = FALSE
  ) +
  geom_point(data = fin_list2[fin_list2$type == "self",],
             aes(x= Het, y =Probability),color= "deepskyblue") +
  geom_smooth(
    aes(x = Het, y =Probability),
    color = "deepskyblue4",
    data = fin_list2[fin_list2$type == "self",],
    se = FALSE
  )  +
  facet_wrap(~dataset) +
  theme_bw()

ggsave("Probability_variety_self_by_ID_80NA_nodup.png",  width = 15, height = 15, units = "in")

################################################################
################################################################
################################################################
Probability_corr <- ggplot(fin_list2) +
  geom_point(data = fin_list2[fin_list2$type == "other",],
             aes(x= Het, y =Probability_corr),color= "gray") +
  geom_smooth(
    aes(x = Het, y =Probability_corr),
    color = "darkseagreen4",
    data = fin_list2[fin_list2$type == "other",],
    se = FALSE
  ) +
  geom_point(data = fin_list2[fin_list2$type == "self",],
             aes(x= Het, y =Probability_corr),color= "deepskyblue") +
  geom_smooth(
    aes(x = Het, y =Probability_corr),
    color = "deepskyblue4",
    data = fin_list2[fin_list2$type == "self",],
    se = FALSE
  )  +
  facet_wrap(~dataset) +
  theme_bw()

ggsave("Probability_corr_variety_self_by_ID_80NA_nodup.png",  width = 15, height = 15, units = "in")
################################################################
################################################################
################################################################
Probability_corr_scaled <- ggplot(fin_list2) +
  geom_point(data = fin_list2[fin_list2$type == "other",],
             aes(x= Het, y =Probability_corr_scaled),color= "gray") +
  geom_smooth(
    aes(x = Het, y =Probability_corr_scaled),
    color = "darkseagreen4",
    data = fin_list2[fin_list2$type == "other",],
    se = FALSE
  ) +
  geom_point(data = fin_list2[fin_list2$type == "self",],
             aes(x= Het, y =Probability_corr_scaled),color= "deepskyblue") +
  geom_smooth(
    aes(x = Het, y =Probability_corr_scaled),
    color = "deepskyblue4",
    data = fin_list2[fin_list2$type == "self",],
    se = FALSE
  )  +
  facet_wrap(~dataset) +
  theme_bw()

ggsave("Probability_corr_scaled_self_by_ID_80NA_nodup_2.png",  width = 15, height = 15, units = "in")
################################################################
################################################################
################################################################

Probability_corr_scaled <- ggplot(fin_list2) +
  # geom_density()+
  # facet_wrap(~dataset)

  geom_point(data = fin_list2[fin_list2$type == "other",],
             aes(x= Het, y =Probability_corr_scaled),color= "gray") +
  geom_smooth(
    aes(x = Het, y =Probability_corr_scaled),
    color = "darkseagreen4",
    data = fin_list2[fin_list2$type == "other",],
    se = FALSE
  ) +
  geom_point(data = fin_list2[fin_list2$type == "self",],
             aes(x= Het, y =Probability_corr_scaled),color= "deepskyblue") +
  geom_smooth(
    aes(x = Het, y =Probability_corr_scaled),
    color = "deepskyblue4",
    data = fin_list2[fin_list2$type == "self",],
    se = FALSE
  )  +
  # facet_wrap(~dataset) +
  theme_bw()

ggsave("Probability_corr_scaled_ALL_self_by_ID_80NA_nodup.png",  width = 10, height = 10, units = "in")
################################################################
################################################################
################################################################

true_pos_tmp <- fin_list2[fin_list2$type == "self",]
true_pos_tmp <- true_pos_tmp[which(true_pos_tmp$Probability_corr_scaled>0.9),]

total_tests_self <- nrow(true_pos_tmp)
true_pos <- nrow(true_pos_tmp[which(true_pos_tmp$Probability_corr_scaled >= 0.95),])
false_neg <- nrow(true_pos_tmp[which(true_pos_tmp$Probability_corr_scaled < 0.95),])

sensitivity <- true_pos /  (true_pos+false_neg)

true_neg_tmp <- fin_list2[fin_list2$type == "other",]
total_tests_other <- nrow(true_neg_tmp)
total_tests <- total_tests_self + total_tests_other

true_neg_tmp <- true_neg_tmp[which(true_neg_tmp$Probability_corr_scaled>0.9),]

# true_neg_tmp2 <- which(true_neg_tmp$pair %in% true_pos_tmp$pair)
# true_neg_tmp3 <- true_neg_tmp[-true_neg_tmp2,]

true_neg <- nrow(true_neg_tmp[which(true_neg_tmp$Probability_corr_scaled < 0.95),])
false_pos <- nrow(true_neg_tmp[which(true_neg_tmp$Probability_corr_scaled >= 0.95),])

specificity <- true_neg / (true_neg + false_pos)

false_pos_rate <- false_pos /  (true_neg + false_pos )
false_neg_rate <- false_neg / (true_pos + false_neg)

total_error_rate <- (false_pos + false_neg)/total_tests

accuracy <- (true_pos + true_neg) / (true_pos + true_neg + false_pos + false_neg)

precision <- true_pos / (true_pos + false_pos)



t2 <- t1$gl.references
class(t2) <- "genlight"
pcoa <- gl.pcoa(t2,
                nfactors = 5)

s <- sum(pcoa$eig[pcoa$eig >= 0])
e <- round(pcoa$eig * 100 / s, 1)



merge_het2 <- rbindlist(sum_list)
# merge_het2$het <- as.factor(round(merge_het2$het,2))

p1 <- ggplot(merge_het2,aes(x= diversity, y =Probability )) +
  geom_point() +
  geom_smooth() +
  ylim(50,110)

p2 <- ggplot(merge_het2,aes(x= diversity, y =Probability_scaled * 100 )) +
  geom_point() +
  geom_smooth()+
  ylim(50,110)

p3 <- ggplot(merge_het2,aes(x= diversity, y =Probability_corr  )) +
  geom_point() +
  geom_smooth()+
  ylim(50,110)

p4 <- ggplot(merge_het2,aes(x= diversity, y =Probability_corr_scaled *100 )) +
  geom_point() +
  geom_smooth()+
  ylim(50,110)


p1 + p2 + p3 + p4

p5 <- ggplot(merge_het2,aes(x= Probability_corr  )) +
  geom_density()

p6 <- ggplot(merge_het2,aes(x =Probability )) +
  geom_density()



res1 <- res_sum$Probability_corr
res_tmp_range3 <- range(res1)
res2 <- scales::rescale(res1,
                to = c(0.5,1),
                from = res_tmp_range3 )

plot(res_sum$Probability_corr_scaled


