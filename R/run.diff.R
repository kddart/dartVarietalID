run.diff <- function(rds.file,
                     n.ref = 10){

crop_name <- gsub(".rds$","", basename(rds.file))

  print(crop_name)
  t1 <- readRDS(rds.file)
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
  pop(test_pop_ref) <- test_pop_ref$other$ind.metrics$TargetID
  ref_pops_sep <- seppop(test_pop_ref)
  dif_res_fin <- NULL
  for (y in 1:length(sam_pops_sep)) {
    sam <- names(sam_pops_sep)[y]
    sam_gl <- sam_pops_sep[[which(names(sam_pops_sep) == sam)]]
    class(sam_gl) <- "genlight"
    dif_res_tmp <- NULL
    for (z in 1:n.ref) {
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
  fst_list <- dif_res_fin
  return(fst_list)
}
