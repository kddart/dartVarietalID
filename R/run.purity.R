
run.purity <- function(rds.file,
                       counts.file,
                       info.file,
                       n.ref = 10,
                       ncores = 12) {

  t1 <- readRDS(rds.file)
  res_summary <-  t1$res_summary

  res_sum <-  t1$res_full
  names_sam <- names(res_sum)
  res_sum2 <- lapply(1:length(res_sum), function(x) {
    tmp <- res_sum[[x]]
    tmp <- tmp[1:n.ref, ]
    tmp$TargetID.sample <- names_sam[x]
    tmp$TargetID.reference <- tmp$TargetID
    return(tmp)
  })

    res_sum3 <- data.table::rbindlist(res_sum2)

    genotypic_counts <- dartVarietalID::ds14.genotypic(
      dartVarietalID::ds14.read(counts.file))
    infoFile <- dartVarietalID::readTargetInfoFile(file = info.file)
    assigned_test_reference <- res_sum3$TargetID
    names(assigned_test_reference) <- res_sum3$TargetID.sample
    res_purity <- dartVarietalID::calculatePurity(genotypic_counts,
                                  infoFile,
                                  assigned_test_reference,
                                  ncores)
    res_purity2 <- unlist(unname(res_purity))
    return(res_purity2)

  }
