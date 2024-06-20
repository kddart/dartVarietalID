run.correlation <- function(rds.file,
                            counts.file,
                            info.file,
                            n.ref = 10){

  # read in counts file and info file
  ref_sam <- read.dart.counts(counts.file = counts.file,
                              info.file = info.file)
  t1 <- readRDS(rds.file)
  res_summary <-  t1$res_full
  names_sam <- names(res_summary)
  res_summary <- lapply(1:length(res_summary), function(x){
    tmp <- res_summary[[x]]
    tmp <- tmp[1:n.ref,]
    tmp$TargetID.sample <- names_sam[x]
    tmp$TargetID.reference <- tmp$TargetID
    return(tmp)
  })
  res_summary <- data.table::rbindlist(res_summary)

  counts <- ref_sam$counts
  cor_df <- res_summary[, c("TargetID.sample", "TargetID.reference")]
  r1 <- apply(cor_df, 1, function(x) {
    summary(lm(counts[, x[1]] ~
                 counts[, x[2]]))$r.squared
  })

  return(r1)

}
