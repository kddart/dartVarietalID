rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
                         pattern = ".rds",
                         full.names = T)
crop_names <- lapply(basename(rds_files),function(x){
  gsub(
    "_DAP.rds$",
    "",x)
})
crop_names <- unlist(crop_names)
fst_list <- as.list(1:10)
names(fst_list) <- crop_names
for(i in 1:length(rds_files)){
  print(crop_names[i])
  t1 <- readRDS(rds_files[i])
  TargetID.sample <- t1$res_summary$TargetID.sample
  res_tmp <- t1$res_full
  test_pop_ref <- t1$gl.references
  test_pop_sam <- t1$gl.samples

  res_sum <- t1$res_summary

  # Separating populations
  sam_pops_sep <- seppop(test_pop_sam)
  ref_pops_sep <- seppop(test_pop_ref)

  dif_res_f <- NULL
  for(y in 1:length(sam_pops_sep)){

    sam <- names(sam_pops_sep)[y]
    ref <- names(ref_pops_sep)[which(names(ref_pops_sep)==res_sum[y,"variety.reference"])]
    sam_gl <- sam_pops_sep[[which(names(sam_pops_sep)==sam)]]
    ref_gl <- ref_pops_sep[[which(names(ref_pops_sep)==ref)]]

    sam_ref <- rbind(sam_gl,ref_gl)
    pop(sam_ref) <- as.factor(pop(sam_ref))

    dif_res <- dartR::utils.basic.stats(sam_ref)
    dif_res <- dif_res$overall
    dif_res <- dif_res["Fstp"]
    dif_res_f <- c(dif_res_f,dif_res)
  }
  fst_list[[i]] <- unname(dif_res_f)

}

saveRDS(fst_list,"fst_list.rds")

res_sum_fin <- as.list(1:10)
names(res_sum_fin) <- crop_names
for(i in 1:length(rds_files)){
  t1 <- readRDS(rds_files[i])
  res_sum2 <- t1$res_summary
  res_sum2$fst <-  fst_list[[i]]
  res_sum_fin[[i]] <- res_sum2

  p1 <-   ggplot(res_sum2,aes(x= Probability.reference/100,
                               y= fst))+
    geom_pointdensity()+
    scale_color_viridis() +
    geom_smooth()  +
    labs(title=crop_names[[i]]) +
    xlim(0.50,1)+
    ylim(0,1) +
    xlab("Probability") +
    ylab("fst") +
    theme_bw() +
    theme(plot.title = element_text(size=14)) +
    theme(legend.position="none")

  print(p1)

}







