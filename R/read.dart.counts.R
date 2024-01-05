#' @name read.dart.counts
#' @title Reads count data and info file
#' @param counts.file csv file containing the count data [required].
#' @param info.file csv file containing samples information [required].
#' @return A list with the following elements: counts, loc_metrics, ind_metrics
#' @export

read.dart.counts <- function(counts.file,
                             info.file
                             ) {

  getLastMarkerMetaDataField <- function(filepath) {
    top <- read.csv(
      filepath,
      header = FALSE,
      nrows = 20,
      stringsAsFactors = FALSE
    )

    last_metric <-
      top[dplyr::last(which(top[, 1] == "*")) + 1, dplyr::last(which(top[1,] == "*"))]
    return(last_metric)
  }

  # Read in info file
<<<<<<< HEAD
  ind_metrics <-
    read.csv(
      info.file,
      na.strings = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

=======
  if (inherits(info.file, "data.frame")) {
	  ind_metrics <- info.file  
  } else {
	  ind_metrics <-
		read.csv(
		  info.file,
		  na.strings = "",
		  check.names = FALSE,
		  stringsAsFactors = TRUE
		)
  }	
  
>>>>>>> e65c59229423248f3d8467a05bcb737e17e30e29
  # check whether the TargetID column is present
  id.col <- match("TargetID", names(ind_metrics))
  if (is.na(id.col)) {
    stop("Fatal Error: There is no 'TargetID' column in the info file\n")
  }
  # check whether the SampleType column is present
  type.col <- match("SampleType", names(ind_metrics))
  if (is.na(type.col)) {
    stop("Fatal Error: There is no 'SampleType' column in the info file\n")
  }
  # check whether the Genotype column is present
  gen.col <- match("Genotype", names(ind_metrics))
  if (is.na(gen.col)) {
    stop("Fatal Error: There is no 'Genotype' column in the info file\n")
  }
  # check whether the SampleType column is present
  ref.col <- match("RefType", names(ind_metrics))
  if (is.na(ref.col)) {
    stop("Fatal Error: There is no 'RefType' column in the info file\n")
  }

  # ind_metrics$SampleType <- as.character(ind_metrics$SampleType)
  # ind_metrics[!is.na(ind_metrics$SampleType),"SampleType"] <- "sample"
  ind_metrics[which(!is.na(ind_metrics$SampleType)),"SampleType"] <- "sample"
  ind_metrics[which(is.na(ind_metrics$SampleType)),"SampleType"] <- "reference"

  # Read in headings counts file
  tdummy <-
    read.csv(
      counts.file,
      na.strings = "-",
      check.names = FALSE,
      nrows = 20,
      header = FALSE,
      stringsAsFactors = TRUE
    )

  # get first column and first row
  colskip <- sum(tdummy[1,] == "*") + 1
  topskip <- sum(tdummy[, 1] == "*")

  # Read in counts file
  snpraw <-
    read.csv(
      counts.file,
      na.strings = "-",
      skip = topskip,
      check.names = FALSE,
      stringsAsFactors = TRUE
    )

  # pairing IDs between info file and counts file
  counts_names <- unname(unlist(tdummy[topskip + 1, colskip:ncol(tdummy)]))
  # convert to character
  counts_names <- as.character(counts_names)
  ind_metrics$TargetID <- as.character(ind_metrics$TargetID)
  # trim spaces from id names
  counts_names <- trimws(counts_names, which = "both")
  counts_names <- gsub(" ", "_", counts_names)
  ind_metrics$TargetID <- trimws(ind_metrics$TargetID, which = "both")
  ind_metrics$TargetID <- gsub(" ", "_", ind_metrics$TargetID)
  merge_names <- dplyr::left_join(
    data.frame(TargetID = counts_names),
    ind_metrics,
    by = "TargetID",
    relationship = "many-to-many"
  )
  merge_names <- merge_names[stats::complete.cases(merge_names$SampleType), ]

  # get loc metrics
  lastmetric <- getLastMarkerMetaDataField(counts.file)
  lmet <- which(colnames(snpraw) == lastmetric)
  covmetrics <- snpraw[, 1:lmet]

  # find unique id for loci
  if ("AlleleID" %in% colnames(covmetrics)) {
    covmetrics$clone <-
      sub("\\|.*", "", covmetrics$AlleleID, perl = TRUE)
    spp <- sub(".+-+(\\d{1,3}):.+", "\\1", covmetrics$AlleleID)
    covmetrics$uid <- paste(covmetrics$clone, spp, sep = "-")
  }
  if ("MarkerName" %in% colnames(covmetrics)) {
    covmetrics$uid <- covmetrics$MarkerName
  }

  # getting samples from counts file based on info file
  datas <- snpraw[, (lmet + 1):ncol(snpraw)]

  keep_names <- which(colnames(datas) %in% make.unique(merge_names$TargetID))
  datas <- datas[, keep_names]
  merge_names$TargetID <- make.unique(merge_names$TargetID)
  merge_names <- dplyr::left_join(data.frame(TargetID = colnames(datas)),
                           merge_names,
                           by = "TargetID",
                           relationship = "many-to-many")
  merge_names <- merge_names[stats::complete.cases(merge_names$SampleType), ]

  if(sum(merge_names$TargetID != colnames(datas))){
    stop(cat(
      "TargetID and ID row in counts file do not match."
    ))
  }

  # Provide a summary of the data
  n_ref_sam <- table(merge_names$SampleType)
  cat("\nSummary of the dataset\n")
  cat("  No. of loci:", nrow(datas) / 2, "\n")
  cat("  No. of references:",  n_ref_sam["reference"], "\n")
  cat("  No. of samples:", n_ref_sam["sample"] , "\n")

  return(list(
    counts = datas,
    loc_metrics = covmetrics,
    ind_metrics = merge_names
  ))
}
