
DAP_write_excel <- function(x, info.n, filename) {
  info <- read.csv(info.n)
  nms <- names(x$res_full)
  resfull = do.call(rbind, lapply(
    1:length(x$res_full),
    \(i) data.frame(field_TID = nms[i], x$res_full[[i]])
  ))
  colnames(resfull)[1:4] <-
    c("field_Tid", "ref_Tid", "ref_id", "variety")
  info <- info[, c("TargetID", "sample")]
  colnames(info)[2] <- "field_id"
  full <- merge(resfull, info, by = 1, all.x = TRUE)
  full$var_rank <-
    with(full, stats::ave(Probability, field_id, FUN = \(x) rank(1000 - x, ties.method =
                                                                   "min")))
  full$Probability <- full$Probability / 100
  x[[2]] <- full
  colnames(x$res_summary)[1:7] <-
    c("field_Tid",
      "field_id",
      "ref_Tid",
      "ref_id",
      "variety",
      "NA.perc",
      "Probability")
  x$res_summary$Probability <- x$res_summary$Probability / 100
  x[[3]] <- as.data.frame(x$ref_distance, check.names = FALSE)
  writexl::write_xlsx(x[1:3], paste0(filename, ".xlsx"))
}

library(dartVarietalID)

host <- system("hostname", TRUE)
if (host == "LAPTOP-IVSPBGCA") {
  path <-
    "G:/.shortcut-targets-by-id/1mfeEftF_LgRcxOT98CBIaBbYN4ZHkBr_/share/image"
} else if (grepl("farm.hpc.ucdavis.edu", host)) {
  path <- "."
} else if (host == "DESKTOP-M2BA7AA") {
  path <- "google drive path"
} else if (host == "Luiss-MBP.diversityarrays.com") {
  path <- "/Users/mijangos"
} else if (host == "Luiss-MacBook-Pro.local") {
  path <- "/Users/mijangos"
}

setwd(path)

ff <-
  list.files("DAP_input",
             pattern = "Counts.csv$",
             recursive = TRUE,
             full = TRUE)

for (i in 7:length(ff)) {
  counts.file <- ff[i]
  ordnr <- gsub("_Counts.csv", "", basename(counts.file))
  filename <- file.path(gsub("input", "output/DAP", counts.file))
  filename <- gsub("_Counts.csv$", "_DAP", filename)
  print(filename)

  #	if (file.exists(paste0(filename, ".rds"))) next
  dir.create(dirname(filename), FALSE, TRUE)

  # ojo: using the fixed references
  info.file <- file.path("DAP_input",
                         gsub(
                           "Counts.csv$",
                           "variety-info-refined.csv",
                           basename(counts.file)
                         ))

  x <- runSampleAnalysis(
    counts.file = counts.file,
    info.file = info.file ,
    ncores = parallel::detectCores()-2,
    pop.size = 10,
    dis.mat = FALSE,
    plot.ref = FALSE,
    gen_dif = FALSE,
    purity = FALSE,
    overlap = FALSE,
    correlation = FALSE,
    na.perc.threshold = 50
  )


  saveRDS(x, paste0(filename, ".rds"))
  DAP_write_excel(x, info.n = info.file, filename)
  # fpdf <- paste0(ordnr, "_ref_distance.pdf")
  # file.rename(fpdf, file.path("output/DAP", fpdf))
}
