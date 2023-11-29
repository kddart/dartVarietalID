calculatePurity <- function(genotypic_counts,
                            infoFile,
                            assigned_test_reference,
                            ncores) {
  infoFile_references <- infoFile$getReferences()
  refsWithTargetID <- data.frame(TargetID = rownames(infoFile_references),
                                 RefType = infoFile_references$RefType)

  varietiesGenotypes <- extractVarietyGenotypes(genotypic_counts,
                                                refsWithTargetID)
  varietiesGenotypesWithStats <- summaryStatsVarietyGenotypes(varietiesGenotypes,
                                                              n.cores = ncores)
  varietiesDiscretised <- lapply(varietiesGenotypesWithStats,
                                 discretiseVarietyMeanAndNZThresh)
  names(varietiesDiscretised) <-
    sapply(varietiesDiscretised, function(x)
      x$name)
  unknown_genotypes <- extractVarietyGenotypesAsUnknown(genotypic_counts,
                                                        names(assigned_test_reference))
  unknownSumList <- summaryStatsVarietyGenotypes(unknown_genotypes,
                                                 n.cores = ncores)
  do.call(rbind, lapply(1:length(assigned_test_reference), function(i) {
    result <- profileSampleAgainstVariety(unknownSumList[[i]],
                                          varietiesDiscretised[[assigned_test_reference[[i]]]])

    data.frame(
      absent_score = result$absent_score,
      present_score = result$present_score,
      purityPercent = result$purityPercent
    )
  }))
}

extractVarietyGenotypesAsUnknown <- function(genotypes,
                                             targets) {
  return(lapply(targets, function(x) {
    g <- matrix(genotypes[, as.character(x)], ncol = 1)

    colnames(g) <- x

    rownames(g) <- rownames(genotypes)

    obj <- {

    }
    obj$name <- "UNKNOWN"

    obj$genotypes <- g

    return(obj)
  }))
}

profileSampleAgainstVariety <- function(sample,
                                        discretise_variety) {
    obj <- {}

    obj$sample$n_total <- dim(sample$genotypes)[1]

    if (class(discretise_variety$present$genotypes)[1] == "list") {
      obj$variety$n_present <-
        length(discretise_variety$present$genotypes)
      obj$variety$n_absent <-
        length(discretise_variety$absent$genotypes)
      obj$variety$n_total <-
        obj$variety$n_present + obj$variety$n_absent

    } else {
      obj$variety$n_present <-
        dim(discretise_variety$present$genotypes)[1]
      obj$variety$n_absent <-
        dim(discretise_variety$absent$genotypes)[1]
      obj$variety$n_total <-
        obj$variety$n_present + obj$variety$n_absent

    }

    vrownamespresent <-
      rownames(discretise_variety$present$genotypes)
    vrownamesabsent <- rownames(discretise_variety$absent$genotypes)

    sampleFilteredByPresent <-
      sample$genotypes[rownames(sample$genotypes) %in%
                         vrownamespresent, ]
    sampleFilteredByAbsent <-
      sample$genotypes[rownames(sample$genotypes) %in%
                         vrownamesabsent, ]

    obj$present$filtered <- length(sampleFilteredByPresent)
    obj$present$matched <- sum(sampleFilteredByPresent > 0)
    obj$present$mismatched <-
      obj$present$filtered - obj$present$matched

    obj$absent$filtered <- length(sampleFilteredByAbsent)
    obj$absent$matched <- sum(sampleFilteredByAbsent <= 1)
    obj$absent$mismatched <-
      obj$absent$filtered - obj$absent$matched

    if (obj$present$matched + obj$present$mismatched == 0) {
      obj$present_score <- 1
    } else {
      obj$present_score <- obj$present$matched / (obj$present$matched +
                                                    obj$present$mismatched)
    }

    if (obj$absent$matched + obj$absent$mismatched == 0) {
      obj$absent_score <- 1
    } else {
      obj$absent_score <- obj$absent$matched / (obj$absent$matched +
                                                  obj$absent$mismatched)
    }

    obj$match <- obj$present_score + obj$absent_score
    obj$purityPercent <- obj$absent_score * 100

    return(obj)

  }

discretiseVarietyMeanAndNZThresh <- function(variety,
                                             presenceMeanCountThreshold = 5,
                                             absentThresholdMean = 0.1) {
  gs <- cbind(variety$genotypes, variety$stats)
  gs_sorted <- gs[order(unlist(gs[, "mean"]), decreasing = TRUE), ]
  presentfilt <- gs_sorted[, "mean"] >= presenceMeanCountThreshold
  obj <- {
  }
  obj$name <- variety$name
  # Extract present
  obj$present <- {
  }
  obj$present$genotypes <-
    gs_sorted[presentfilt, 1:ncol(variety$genotypes),
              drop = FALSE]
  if (is.null(dim(obj$present$genotypes))) {
    obj$present$genotypes <- as.data.frame(t(obj$present$genotypes),
                                           row.names = rownames(gs_sorted)[presentfilt])
  } else {
    obj$present$genotypes <- as.matrix(obj$present$genotypes)
  }
  stats <- gs_sorted[presentfilt, (ncol(variety$genotypes) +
                                     1):ncol(gs_sorted), drop = F]
  if (ncol(variety$genotypes) == 1) {
    obj$present$stats <- as.matrix(t(data.frame(stats)))
  } else {
    obj$present$stats <- as.matrix(stats)
  }
  if (is.null(dim(obj$present$stats))) {
    obj$present$stats <- as.data.frame(t(obj$present$stats),
                                       row.names = rownames(gs_sorted)[presentfilt])
  } else {
    obj$present$stats <- as.matrix(obj$present$stats)
  }
  # Extract absent
  obj$absent <- {
  }
  absentfilt <- which(TRUE == apply(as.data.frame(gs_sorted[,
                                                            colnames(variety$genotypes)]), 1, function(x) {
                                                              a <- all(x <= 1) && mean(unlist(x)) < absentThresholdMean

                                                              return(a)
                                                            }))

  obj$absent$genotypes <-
    gs_sorted[absentfilt, 1:ncol(variety$genotypes),
              drop = FALSE]
  if (is.null(dim(obj$absent$genotypes))) {
    obj$absent$genotypes <- as.data.frame(t(obj$absent$genotypes),
                                          row.names = rownames(gs_sorted)[absentfilt])
  } else {
    obj$absent$genotypes <- as.matrix(obj$absent$genotypes)
  }
  obj$absent$stats <-
    gs_sorted[absentfilt, (ncol(variety$genotypes) +
                             1):ncol(gs_sorted), drop = FALSE]
  if (is.null(dim(obj$absent$stats))) {
    obj$absent$stats <- as.data.frame(t(obj$absent$stats),
                                      row.names = rownames(gs_sorted)[absentfilt])
  } else {
    obj$absent$stats <- as.matrix(obj$absent$stats)
  }

  return(obj)
}

summaryStatsVarietyGenotypes <- function(varietiesGenotypes,
                                         n.cores) {
  if (n.cores == 1) {
    d <- lapply(
      X = varietiesGenotypes,
      FUN = function(v) {
        v$stats <- t(apply(as.data.frame(v$genotypes), 1,
                           calcStatsMarker))
        return(v)
      }
    )

  } else {
    # if unix
    if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
      d <- mclapply(
        X = varietiesGenotypes,
        FUN = function(v) {
          v$stats <- t(apply(as.data.frame(v$genotypes),
                             1, calcStatsMarker))
          return(v)
        },
        mc.cores = n.cores
      )
    }

    ## if windows
    if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
      d <- lapply(
        X = varietiesGenotypes,
        FUN = function(v) {
          v$stats <- t(apply(as.data.frame(v$genotypes),
                             1, calcStatsMarker))
          return(v)
        }
      )
    }
  }

  return(d)

}

calcStatsMarker <- function(m) {
  return(c(
    mean = mean(m),
    sd = sd(m),
    nz = sum(m == 0),
    nnz = sum(m != 0),
    sum = sum(m),
    sum_nnz = sum(m[m != 0])
  ))

}

extractVarietyGenotypes <- function(genotypes,
                                    meta) {
  varieties <- unique(as.character(meta$RefType))
  varieties <- varieties[varieties != ""]
  meta_grouped <- lapply(varieties, function(v) {
    return(meta[meta$RefType == v, ])
  })
  return(lapply(meta_grouped, function(g) {
    obj <- {
    }
    obj$name <- as.character(g$RefType[1])
    obj$genotypes <- sapply(g$TargetID, function(t) {
      g <- t(as.matrix(genotypes[, as.character(t)]))
      if (is.null(g)) {
        warning(
          paste(
            "Missing sample",
            t,
            "in genotypes specified in meta (Variety=",
            g$RefType[1]
          )
        )
      }
      return(g)
    })
    colnames(obj$genotypes) <- g$TargetID
    rownames(obj$genotypes) <- rownames(genotypes)
    return(obj)
  }))
}
