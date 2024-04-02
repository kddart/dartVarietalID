#' @name runSampleAnalysis
#' @title Main script
#' @description
#' Main script
#' @param counts.file File with count data [required].
#' @param info.file File with information of samples and references[required].
#' @param ncores Number of cores to run analysis
#' [default parallel::detectCores() - 1].
#' @param pop.size Number individuals to simulate [default 10].
#' @param dis.mat Whether to create and save to wd a genetic distance plot of
#'  the references [default TRUE].
#' @param plot.ref [default TRUE].
#' @param gen_dif [default TRUE].
#' @param purity [default TRUE].
#' @param overlap [default TRUE].
#' @param correlation [default TRUE].
#' @param na.perc.threshold Threshold for missing data to remove references
#' and samples [default 50].
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
#' @import dendextend
#' @import rgl
#' @import SIBER
#' @export

runSampleAnalysis <- function(counts.file,
                              info.file,
                              ncores = parallel::detectCores() -1,
                              pop.size = 10,
                              dis.mat = TRUE,
                              plot.ref = TRUE,
                              gen_dif = TRUE,
                              purity = TRUE,
                              overlap = TRUE,
                              correlation =TRUE,
                              na.perc.threshold = 50) {

  # read in counts file and info file
  ref_sam <- read.dart.counts(counts.file = counts.file,
                              info.file = info.file)
  # converting to genotypes
  ref_sam_pops <- counts2geno(count.data = ref_sam,
                              pop.size = pop.size)

  # separating references from samples
  pop(ref_sam_pops) <- ref_sam_pops$other$ind.metrics$reference
  test_pop_ref <- dartR::gl.keep.pop(ref_sam_pops,
                                     pop.list = "reference",
                                     verbose = 0)
  pop(test_pop_ref) <- test_pop_ref$other$ind.metrics$TargetID
  # removing references with more than the threshold of missing data
  test_pop_ref_NA <- seppop(test_pop_ref)
  na_check_ref <- unlist(lapply(test_pop_ref_NA, na_check))
  pop_drop_ref <- which(na_check_ref > na.perc.threshold)
  if(length(pop_drop_ref)>0){
    test_pop_ref <- dartR::gl.drop.pop(test_pop_ref,
                                pop.list = popNames(test_pop_ref)[pop_drop_ref])
    message(length(pop_drop_ref)," references with more than ",na.perc.threshold,
    " percentage of missing data were removed: ", paste(names(pop_drop_ref)," "))
  }

  # test_pop_ref$other$ind.metrics$variety <- trimws(test_pop_ref$other$ind.metrics$variety, which = "both")
  # test_pop_ref$other$ind.metrics$variety <- gsub(" ","_",test_pop_ref$other$ind.metrics$variety)

  # pop(test_pop_ref) <- test_pop_ref$other$ind.metrics$variety
  pop(test_pop_ref) <- test_pop_ref$other$ind.metrics$variety

  test_pop_sam <- dartR::gl.keep.pop(ref_sam_pops,
                                     pop.list = "sample",
                                     verbose = 0)
  pop(test_pop_sam) <- test_pop_sam$other$ind.metrics$TargetID
  # removing samples with more than the threshold of missing data
  test_pop_sam_NA <- seppop(test_pop_sam)
  na_check_sam <- unlist(lapply(test_pop_sam_NA, na_check))
  pop_drop_sam <- which(na_check_sam > na.perc.threshold)
  if(length(pop_drop_sam)>0){
    test_pop_sam <- dartR::gl.drop.pop(test_pop_sam,
                                pop.list = popNames(test_pop_sam)[pop_drop_sam])
    message(length(pop_drop_sam)," samples with more than ",na.perc.threshold,
            " percentage of missing data were removed: ", paste(names(pop_drop_sam)," "))
  }

  if (dis.mat || plot.ref) {
    test_pop_ref_2 <- test_pop_ref
    pop(test_pop_ref_2) <-
      paste0(
        test_pop_ref_2$other$ind.metrics$TargetID,
        "_",
        test_pop_ref_2$other$ind.metrics$variety
      )

    test_pop_ref_2$other$ind.metrics$variety <- as.factor(test_pop_ref_2$other$ind.metrics$variety)
    t1 <- dartR::gl.dist.pop(test_pop_ref_2, method = "nei",
                             plot.out = FALSE,
                             verbose = 0)
    t1 <- as.matrix(t1)

	if (plot.ref) {

		colors_pops <-
		  polychrome(length(levels(
			test_pop_ref_2$other$ind.metrics$variety
		  )))
		names(colors_pops) <-
		  as.character(levels(test_pop_ref_2$other$ind.metrics$variety))

		df_colors_temp_1 <-
		  as.data.frame(cbind(
			as.character(pop(test_pop_ref_2)),
			as.character(test_pop_ref_2$other$ind.metrics$variety)
		  ))
		df_colors_temp_1 <- unique(df_colors_temp_1)
		df_colors_temp_1$order <- 1:nPop(test_pop_ref_2)
		colnames(df_colors_temp_1) <- c("ind", "pop", "order")

		df_colors_temp_2 <- as.data.frame(cbind(names(colors_pops), colors_pops))
		colnames(df_colors_temp_2) <- c("pop", "color")
		df_colors <- merge(df_colors_temp_1, df_colors_temp_2, by = "pop")
		df_colors$order <- as.numeric(df_colors$order)
		df_colors <- df_colors[order(df_colors$order),]

		df_colors_2 <- merge(data.frame(ind=colnames(t1)),
										df_colors,by="ind" )

		palette.divergent <- dartR.base::gl.colors("div")

		pdf(
		  paste0(strsplit(basename(counts.file), "_")[[1]][1], "_ref_distance.pdf"),
		  width = nPop(test_pop_ref_2) / 5,
		  height = nPop(test_pop_ref_2) / 5
		)
		heatmap.3(
		  t1,
		  margins = c(10, 10),
		  ColSideColors = df_colors_2$color,
		  RowSideColors = df_colors_2$color,
		  sepcolor = "black",
		  dendrogram = "column",
		  trace = "none",
		  col = viridis::turbo(255),
		  colRow = df_colors_2$color,
		  colCol = df_colors_2$color,
		  density.info = "none",
		  reorderfun = function(d, w){reorder(d, w, agglo.FUN = mean, na.rm = TRUE)},
		  main = "Genetic distance (Nei's distance) between references",
		   na.rm = TRUE,
		  na.color = "grey"
		)
		# Close device
		dev.off()
	}
  }

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
  # if unix
  if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
    top_ind <- parallel::mclapply(X = sam_pops_sep,
                                FUN = rep_ind,
                                mc.cores = ncores)
  }

  ## if windows
  if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
    top_ind <- lapply(sam_pops_sep,rep_ind)
  }

  # if unix
  if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {

    res_tmp <- parallel::mclapply(X = top_ind,
                                  FUN = dart.assignment,
                                  ref = ref_pops_sep,
                                  mc.cores = ncores)

  }

  ## if windows
  if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {

    res_tmp <- lapply(X = top_ind,
                      FUN = dart.assignment,
                      ref = ref_pops_sep)

  }

  if(gen_dif){

    res_dif <- parallel::mclapply(X = sam_pops_sep,
                                  FUN = dart.differentiation,
                                  ref = ref_pops_sep,
                                  mc.cores = ncores)
    res_tmp2 <- lapply(1:length(res_dif),function(x){
    tmp1 <- merge(res_tmp[[x]],res_dif[[x]],by="variety")
    tmp1 <- tmp1[order(tmp1$Probability,decreasing = TRUE),]
    })
    res_tmp <- res_tmp2
    }

  # summary results dataframe
  TargetID.sample <- unlist(lapply(top_ind,function(x){
    x$other$ind.metrics$TargetID
  }))
  sample.sample <- unlist(lapply(top_ind,function(x){
    x$other$ind.metrics$sample
  }))
  res_tmp2 <- lapply(res_tmp,"[",1,)
  res_tmp3 <- data.table::rbindlist(res_tmp2)

  TargetID.reference <- res_tmp3$TargetID
  sample.reference <- res_tmp3$sample
  variety.reference <- res_tmp3$variety
  Probability.reference <- res_tmp3$Probability
  NA.percentage <- round((1 - res_tmp3$NumLoci / nLoc(ref_sam_pops)) * 100,2)

  if(gen_dif){
  res_summary <- data.frame(TargetID.sample = TargetID.sample,
                            sample.sample = sample.sample,
                            TargetID.reference = TargetID.reference,
                            sample.reference = sample.reference,
                            variety.reference = variety.reference,
                            NA.percentage = NA.percentage,
                            Probability.reference = Probability.reference,
                            Fst = res_tmp3$Fst,
                            Fstp = res_tmp3$Fstp,
                            Dest = res_tmp3$Dest,
                            Gst_H = res_tmp3$Gst_H
  )
  }else{
    res_summary <- data.frame(TargetID.sample = TargetID.sample,
                              sample.sample = sample.sample,
                              TargetID.reference = TargetID.reference,
                              sample.reference = sample.reference,
                              variety.reference = variety.reference,
                              NA.percentage = NA.percentage,
                              Probability.reference = Probability.reference
    )

  }

  names(res_tmp) <- TargetID.sample

  if(purity){
  # Calculating purity
  genotypic_counts <- ds14.genotypic(ds14.read(counts.file))
  infoFile <- readTargetInfoFile(file = info.file)
  assigned_test_reference <- res_summary$variety.reference
  names(assigned_test_reference) <- res_summary$TargetID.sample
  res_purity <- calculatePurity(genotypic_counts,
                                infoFile,
                                assigned_test_reference,
                                ncores)
  res_summary <- cbind(res_summary,res_purity)
  # Setting NAs to samples with more than 50% of missing data
  col_NAs <- c("Probability.reference",
               "purityPercent")

  }else{

    col_NAs <- c("Probability.reference")
  }

  if(overlap){

    # # if unix
    # if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
    #
    #   res <- parallel::mclapply(X = 1:length(TargetID.sample),
    #                             FUN =
    #                             function(x){
    #                               overlap_proportion(
    #                               test.sample = unname(TargetID.sample)[[x]],
    #                               full.report = res_tmp,
    #                               ref = test_pop_ref,
    #                               sam = test_pop_sam,
    #                               n.varieties=10,
    #                               plot = FALSE)
    #                             }
    #                   ,
    #                   mc.cores = ncores)
    #
    # }
    #
    # ## if windows
    # if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {

      res <- lapply(1:length(TargetID.sample),function(x){

        tmp <-  overlap_proportion(test.sample= names(TargetID.sample)[x],
                                   full.report = res_tmp,
                                   ref = test_pop_ref,
                                   sam = test_pop_sam,
                                   n.varieties=10,
                                   plot = FALSE)

        return(tmp)

      })

    # }

    res2 <- as.data.frame(Reduce(rbind,res))
    colnames(res2) <- c("id","overlap")
    res_summary <- cbind(res_summary,res2)
    res_summary$overlap <- as.numeric(res_summary$overlap)
  }

  if(correlation){

    counts <- ref_sam$counts
    cor_df <- res_summary[,c("TargetID.sample","TargetID.reference")]
    r1 <-apply(cor_df,1,function(x){
      summary(lm(counts[,x[1]]~
                   counts[,x[2]]))$r.squared
    })

    res_summary$corr <- r1


  }

  # Setting NAs to samples with more than 50% of missing data
  res_summary[which(res_summary$NA.percentage>50),col_NAs ] <- NA
  res_summary$NA.percentage <- round(res_summary$NA.percentage, 2)
  res_summary$Probability.reference <- round(res_summary$Probability.reference, 2)
  # res_summary$purityPercent <- round(res_summary$purityPercent, 2)

  return(list(
    res_summary = res_summary,
    res_full = res_tmp,
    gl.references = test_pop_ref,
    gl.samples = test_pop_sam,
	ref_distance = t1
  ))

}
