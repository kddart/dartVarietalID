#' @name runSampleAnalysis
#' @title Main script
#' @description
#' Main script
#' @param counts.file [required].
#' @param info.file [required].
#' @param ncores [default parallel::detectCores() - 1].
#' @param pop.size [default 10].
#' @details
#' Main script
#' @return A list with the following elements:
#'
#' res_summary. A data frame containing the best reference match for each sample
#' and its probability.
#'
#' res_full. A list of data frames with the probability match for each sample in
#'  the dataset.
#'
#' gl.references. A genlight object containing the references.
#'
#' gl.samples. A genlight object containing the samples.
#' @author Luis Mijangos
#' @examples
#'
#' filename <- system.file('extdata','SEQ_SNPs_counts_0_Target.csv', package='dartVarietalID')
#' info <- system.file('extdata','InfoFile_corrected.csv', package='dartVarietalID')
#' ID_res <- runSampleAnalysis(counts.file = filename, info.file = info,
#' ncores = 1, pop.size = 10)
#'
#' @import dplyr
#' @import dartR
#' @import stats
#' @import parallel
#' @import data.table
#' @import Rcpp
#' @import shinyWidgets
#' @import shiny
#' @import shinyjs
#' @import tableHTML
#' @import colorspace
#' @import adegenet
#' @import utils
#' @import plotly
#' @import methods
#' @import semantic.dashboard
#' @import shiny.semantic
#' @import DT
#' @export

runSampleAnalysis <- function(counts.file,
                                  info.file,
                                  ncores = parallel::detectCores() - 1,
                                  pop.size = 10) {

  # read in counts file and info file
  ref_sam <- read.dart.counts(counts.file = counts.file,
                              info.file = info.file)
  # converting to genotypes
  ref_sam_pops <- counts2geno(count.data = ref_sam,
                              pop.size = pop.size)
  # separating references from samples
  pop(ref_sam_pops) <- ref_sam_pops$other$ind.metrics$SampleType
  test_pop_ref <- dartR::gl.keep.pop(ref_sam_pops,
                                     pop.list = "reference",
                                     verbose = 0)
  test_pop_sam <- dartR::gl.keep.pop(ref_sam_pops,
                                     pop.list = "sample",
                                     verbose = 0)
  pop(test_pop_ref) <- test_pop_ref$other$ind.metrics$RefType
  pop(test_pop_sam) <- test_pop_sam$other$ind.metrics$TargetID

  # Separating populations
  sam_pops_sep <- seppop(test_pop_sam)
  ref_pops_sep <- seppop(test_pop_ref)

  # selecting the representative individual from the sample using PCA
  top_ind <- as.list(1:length(sam_pops_sep))
  for (y in 1:length(sam_pops_sep)) {
    pop_test_hold <- sam_pops_sep[[y]]
    # removing missing data for PCA
    pop_test <- dartR::gl.filter.callrate(pop_test_hold,
                                          threshold = 1,
                                          verbose = 0)
    # test whether all individuals are the same (when all loci are monomorphic)
    pop_test_mat <- as.matrix(pop_test)
    test_var <- sum(apply(pop_test_mat, 2, function(x) {
      var(x) != 0
    }), na.rm = TRUE)
    # if individuals are different
    if (test_var > 0) {
      pcoa <- adegenet::glPca(pop_test,
                    nf = 3,
                    loadings = FALSE)

      pcoa_scores <- pcoa$scores
      means <- colMeans(pcoa_scores)
      covariance <- stats::cov(pcoa_scores)
      D <- stats::mahalanobis(pcoa_scores, means, covariance, toll = 1e-20)
      top_ind[[y]] <- pop_test_hold[which.min(D),]
      # if individuals are the same, get the first individual
    } else{
      top_ind[[y]] <- pop_test_hold[1,]
    }
  }

  # assigning reference to samples
  if(ncores == 1){
    res_tmp <- lapply(X = top_ind,
                      FUN = dart.assignment,
                      ref = ref_pops_sep)
  }else{

    # if unix
    if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
      res_tmp <- parallel::mclapply(X = top_ind,
                                    FUN = dart.assignment,
                                    ref = ref_pops_sep,
                                    mc.cores = ncores)
    }

    ## if windows
    if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {

      cl <- parallel::makePSOCKcluster(rep("localhost",
                                           ncores))
      res_tmp <- parallel::parLapply(cl = cl,
                                     X = top_ind,
                                     fun = dart.assignment,
                                     ref = ref_pops_sep)
      stopCluster(cl)
    }
  }

  # summary results dataframe
  TargetID.sample <- unlist(lapply(top_ind,function(x){
x$other$ind.metrics$TargetID
  }))
  Genotype.sample <- unlist(lapply(top_ind,function(x){
    x$other$ind.metrics$Genotype
  }))
  res_tmp2 <- lapply(res_tmp,"[",1,)
  res_tmp3 <- data.table::rbindlist(res_tmp2)

  TargetID.reference <- res_tmp3$TargetID
  Genotype.reference <- res_tmp3$Genotype
  RefType.reference <- res_tmp3$RefType
  Probability.reference <- res_tmp3$Probability
  NA.percentage <- round((1 - res_tmp3$NumLoci / nLoc(ref_sam_pops)) * 100,2)

  res_summary <- data.frame(TargetID.sample = TargetID.sample,
                            Genotype.sample = Genotype.sample,
                            TargetID.reference = TargetID.reference,
                            Genotype.reference = Genotype.reference,
                            RefType.reference = RefType.reference,
                            NA.percentage = NA.percentage,
                            Probability.reference = Probability.reference
  )

  names(res_tmp) <- TargetID.sample

  # Calculating purity
  genotypic_counts <- ds14.genotypic(ds14.read(counts.file))
  infoFile <- readTargetInfoFile(info.file)
  assigned_test_reference <- res_summary$RefType.reference
  names(assigned_test_reference) <- res_summary$TargetID.sample
  res_purity <- calculatePurity(genotypic_counts,
                                infoFile,
                                assigned_test_reference,
                                ncores)

  res_summary <- cbind(res_summary,res_purity)
  # Setting NAs to samples with more than 50% of missing data
  col_NAs <- c("Probability.reference",
               "absent_score",
               "present_score",
               "purityPercent")
  res_summary[which(res_summary$NA.percentage>50),col_NAs ] <- NA

  return(list(
    res_summary = res_summary,
    res_full = res_tmp,
    gl.references = test_pop_ref,
    gl.samples = test_pop_sam
  ))

}
