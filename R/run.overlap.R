

run.overlap <- function(rds.file,
                        counts.file,
                        info.file,
                        n.ref = 10,
                        ncores = 10,
                        parallel = TRUE){

n.varieties <- n.ref
  print(rds.file)
  t1 <- readRDS(rds.file)
  res_tmp <- t1$res_full
  test_pop_ref <- t1$gl.references
  test_pop_sam <- t1$gl.samples
  pop(test_pop_ref) <- test_pop_ref$other$ind.metrics$TargetID
  # Separating populations
  sam_pops_sep <- seppop(test_pop_sam)
  ref_pops_sep <- seppop(test_pop_ref)

  #identify samples with all missing data
  NAs <-
    lapply(sam_pops_sep, function(x) {
      sum(sapply(x@gen, function(e) {
        length(e@NA.posi)
      }))
    })
  #total number of genotypes
  total_geno <- nInd(sam_pops_sep[[1]]) * nLoc(sam_pops_sep[[1]])
  all_NAs <-  which(NAs == total_geno)

  # remove sample if all data is missing
  if (length(all_NAs) > 0) {
    sam_pops_sep <- sam_pops_sep[-(all_NAs)]
  }
  print(paste(length(sam_pops_sep), "samples"))
  test_gl <- lapply(1:length(sam_pops_sep), function(z) {
    sam_name <- names(sam_pops_sep[z])
    # get the sample to test
    sam_gl <- sam_pops_sep[[z]]
    # get the closest references
    # ref_pops <- res_tmp[[sam_name]][1:n.varieties, "variety"]
    ref_pops <- res_tmp[[sam_name]][1:n.varieties, "TargetID"]

    ref_gl <-
      dartR::gl.keep.pop(test_pop_ref, pop.list = ref_pops, verbose = 0)
    tog <- rbind(ref_gl, sam_gl)
    tog <- dartR::gl.filter.monomorphs(tog, verbose = 0)

    return(list(
      sam = sam_name,
      ref = ref_pops[1],
      tog_gl = tog
    ))
  })

  if(parallel){
    res <- parallel::mclapply(X = 1:length(test_gl),
                              FUN = tryCatch({
                                function(x) {
                                  overlap_proportion(
                                    tog_gl = test_gl[[x]][[3]],
                                    test.sample = test_gl[[x]][[1]],
                                    test.ref = test_gl[[x]][[2]],
                                    plot = FALSE
                                  )
                                }
                              },
                              # error handling
                              error = function(e) {
                                return(NULL)
                              },
                              # warning handling
                              warning = function(w) {
                                return(NULL)
                              }),
                              mc.cores = ncores)

  }else{

    res <- lapply(1:length(test_gl), function(x) {
      tmp <- tryCatch({
        overlap_proportion(
          tog_gl = test_gl[[x]][[3]],
          test.sample = test_gl[[x]][[1]],
          test.ref = test_gl[[x]][[2]],
          # test.sample = TargetID.sample[x],
          # full.report = res_tmp,
          # ref = test_pop_ref,
          # sam = test_pop_sam,
          # n.varieties=10,
          plot = F
        )
      },
      # error handling
      error = function(e) {
        return(NULL)
      },
      # warning handling
      warning = function(w) {
        return(NULL)
      })

      return(tmp)

    })

  }

  res <- lapply(res, function(x) {
    if (is.null(x)) {
      x <- matrix(nrow = 1, ncol = 2)
      return(x)
    } else{
      return(x)
    }
  })

  # res4 <- unlist(lapply(res, "[", 2))
  # res2[na.res] <- 0
  res2 <- lapply(res, "[", 2)
  res2 <- unlist(res2)
  overlap_res <- as.numeric(res2)

  return(overlap_res)

}

