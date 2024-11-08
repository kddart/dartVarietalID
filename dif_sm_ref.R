library(dartR)
library(data.table)

shannon <- function(x) {
  x <- x[x > 0]
  p <- x / sum(x)
  out <- -sum(p * log(p))
  return(out)
}

MI <- function(x) {
  mat <- as.matrix(x)

  n <- sum(mat)
  prob.hat <- mat / n
  px.hat <- apply(prob.hat, 1, sum)
  py.hat <- apply(prob.hat, 2, sum)
  I.hat <- shannon(px.hat) + shannon(py.hat) - shannon(prob.hat)
  # MLE of Mutual Information!
  return(I.hat)
}

mutual_information <- function(mat) {
  EstMLEFun <- function(mat) {
    # MLE
    entropyFun <- function(p) {
      p <- p[p > 0]
      out <- -sum(p * log(p))
      return(out)
    }
    n <- sum(mat)
    prob.hat <- mat / n
    px.hat <- apply(prob.hat, 1, sum)
    py.hat <- apply(prob.hat, 2, sum)
    I.hat <-
      entropyFun(px.hat) + entropyFun(py.hat) - entropyFun(prob.hat)
    # MLE of Mutual Information!
    return(I.hat)
  }
  mydata <- as.matrix(mat)
  est <- EstMLEFun(mydata)
  return(est)
}

rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
                         pattern = "_DAP_240514.rds",
                         full.names = T)
crop_names <- lapply(basename(rds_files), function(x) {
  gsub("_DAP_240514.rds$",
       "", x)
})
crop_names <- unlist(crop_names)
fst_list <- as.list(1:10)
names(fst_list) <- crop_names
for (i in 1:length(rds_files)) {
  print(crop_names[i])
  t1 <- readRDS(rds_files[i])
  # TargetID.sample <- t1$res_summary$TargetID.sample
  res_tmp <- t1$res_full
  names_sam <- names(res_tmp)
  res_tmp <- lapply(1:length(res_tmp), function(x) {
    tmp <- res_tmp[[x]]
    tmp$TargetID.sample <- names_sam[x]
    tmp$TargetID.reference <- tmp$variety
    return(tmp)
  })
  res_summary <- res_tmp
  test_pop_ref <- t1$gl.references
  test_pop_sam <- t1$gl.samples
  # res_sum <- t1$res_summary
  res_sum <- res_summary
  # Separating populations
  sam_pops_sep <- seppop(test_pop_sam)
  ref_pops_sep <- seppop(test_pop_ref)
  dif_res_fin <- NULL
  for (y in 1:length(sam_pops_sep)) {
    sam <- names(sam_pops_sep)[y]
    sam_gl <- sam_pops_sep[[which(names(sam_pops_sep) == sam)]]
    class(sam_gl) <- "genlight"
    dif_res_tmp <- NULL
    for (z in 1:10) {
      # ref <- names(ref_pops_sep)[which(names(ref_pops_sep)==res_sum[[y]][x,"variety"])]
      # ref_gl <- ref_pops_sep[[which(names(ref_pops_sep)ref)]]
      # names_ref <- res_sum[[y]][1:10,"TargetID"]
      names_ref <- res_sum[[y]][z, "TargetID"]
      # names_ref <- res_sum[[y]][1,"variety"]
      # ref_gl_tmp <- lapply(1:length(names_ref), function(x) {
      #   tmp_ref <- which(names(ref_pops_sep) == names_ref[x])
      #   tmp <- ref_pops_sep[[tmp_ref]]
      #   return(tmp)
      # })
      # # ref_gl <- Reduce(rbind,ref_gl_tmp)
      # ref_gl <- ref_gl_tmp[[1]]
      ref_gl <- ref_pops_sep[which(names(ref_pops_sep)==names_ref)][[1]]
      class(ref_gl) <- "genlight"
      sam_ref <- rbind(sam_gl, ref_gl)
      # pop(sam_ref) <- as.factor(pop(sam_ref))
      sam_ref$other$loc.metrics <- sam_gl$other$loc.metrics
      sam_ref@other$loc.metrics.flags$monomorphs <- FALSE
      dif_res <- dartR::gl.dist.pop(
        sam_ref,
        method = "nei",
        plot.out = FALSE,
        verbose = 0
      )
      # dartR::gl.propShared()
      # dif_res <- unname(as.matrix(dif_res)[-1,1])
      dif_res <- unname(as.numeric(dif_res))
      # sam_ref$other$loc.metrics <- ref_gl$other$loc.metrics
      # sam_ref@other$loc.metrics.flags$monomorphs <- FALSE
      # sam_ref2 <- gl.filter.callrate(sam_ref,threshold = 1,verbose = 0, plot.out = F)
      # dif_res <- mutual_information(as.matrix(sam_ref2))
      # dif_res <- dartR::utils.basic.stats(sam_ref)
      # dif_res <- dif_res$overall
      # dif_res <- dif_res["Gst_H"]
      dif_res_tmp <- c(dif_res_tmp, dif_res)
    }
    dif_res_fin <- c(dif_res_fin,dif_res_tmp)
  }
  fst_list[[i]] <- dif_res_fin
}

saveRDS(fst_list, "fst_list_Neis_2.rds")

res_sum_fin <- as.list(1:10)
names(res_sum_fin) <- crop_names
for (i in 1:length(rds_files)) {
  t1 <- readRDS(rds_files[i])
  res_sum2 <- t1$res_summary
  res_sum2$fst <-  fst_list[[i]]
  res_sum_fin[[i]] <- res_sum2

  p1 <-   ggplot(res_sum2, aes(x = Probability.reference / 100,
                               y = fst)) +
    geom_pointdensity() +
    scale_color_viridis() +
    geom_smooth()  +
    labs(title = crop_names[[i]]) +
    xlim(0.50, 1) +
    ylim(0, 1) +
    xlab("Probability") +
    ylab("fst") +
    theme_bw() +
    theme(plot.title = element_text(size = 14)) +
    theme(legend.position = "none")

  print(p1)

}
