#' @name runSampleAnalysis
#' @title Main script
#' @description
#' Main script
#' @param counts.file File with count data [required].
#' @param info.file File with information of samples and references[required].
#' @param ncores Number of cores to run analysis
#' [default parallel::detectCores() - 1].
#' @param pop.size Number individuals to simulate [default 10].
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
#' filename <- system.file('extdata','SEQ_SNPs_counts_0_Target.csv',
#' package='dartVarietalID')
#' info <- system.file('extdata','InfoFile_corrected.csv',
#' package='dartVarietalID')
#' ID_res <- runSampleAnalysis(counts.file = filename, info.file = info,
#' ncores = 1, pop.size = 10)
#'
#' @rawNamespace import(data.table, except = c(set,first,last,between))
#' @rawNamespace import(colorspace, except = c(plot,show))
#' @rawNamespace import(graphics, except = c(layout,box,grid))
#' @rawNamespace import(DT, except = c(dataTableOutput,renderDataTable))
#' @rawNamespace import(semantic.dashboard, except = c(column,icon,
#' dropdown_menu,menu_item))
#' @rawNamespace import(shiny, except = c(textInput,showNotification,
#' incProgress,modalDialog,removeModal,sliderInput,selectInput,setProgress,
#' fileInput,withProgress,verticalLayout,numericInput,actionButton,updateSliderInput,
#' updateActionButton,removeNotification,textAreaInput,dateInput,Progress,
#' checkboxInput,splitLayout,icon,updateSelectInput,flowLayout,runExample))
#' @rawNamespace import(adegenet, except = c(plot))
#' @rawNamespace import(dplyr, except = c(lag,filter))
#' @rawNamespace import(dendextend, except = c(cutree))
#' @rawNamespace import(Rcpp, except = c(.DollarNames,prompt,show))
#' @rawNamespace import(plotly, except = c(filter))
#' @rawNamespace import(shiny.semantic, except = c(toggle,menu))
#' @rawNamespace import(shinyWidgets, except = c(alert))
#' @rawNamespace import(methods, except = c(removeClass,show))


#' @import dartR
#' @import stats
#' @import parallel
#' @import utils
#' @import shinyjs
#' @import tableHTML
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
      D <- stats::mahalanobis(x = pcoa_scores,
                              center = means,
                              cov = covariance,
                              toll = 1e-20)
      top_ind[[y]] <- pop_test_hold[which.min(D),]
      # if individuals are the same, get the first individual
    } else{
      top_ind[[y]] <- pop_test_hold[1,]
    }
  }


  # if unix
  if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {

    system.time(
    res <- dart.assignment(ref = ref_pops_sep,unknown = top_ind[[1]])
   )

    time_1 <- system.time(
    res_tmp <- parallel::mclapply(X = top_ind,
                                  FUN = dart.assignment,
                                  ref = ref_pops_sep,
                                  mc.cores = ncores
                                  )
    )
  }

  ## if windows
  if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {

    res_tmp <- lapply(X = top_ind,
                      FUN = dart.assignment,
                      ref = ref_pops_sep)

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
  infoFile <- readTargetInfoFile(file = info.file)
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
