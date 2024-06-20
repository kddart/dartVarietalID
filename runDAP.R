library(dartVarietalID)

ff <-
  list.files("DAP_input",
             pattern = "Counts.csv$",
             recursive = TRUE,
             full = TRUE)

for (i in 1:length(ff)) {

  counts.file <- ff[i]
  ordnr <- gsub("_Counts.csv", "", basename(counts.file))
  filename <- file.path(gsub("input", "output/DAP", counts.file))
  filename <- gsub("_Counts.csv$", "_DAP", filename)
  print(filename)

  info.file <- file.path("DAP_input",
                         gsub(
                           "Counts.csv$",
                           "variety-info-refined.csv",
                           basename(counts.file)
                         ))

  x <- runSampleAnalysis(
    counts.file = counts.file,
    info.file = info.file ,
    ncores = parallel::detectCores(),
    maf = 0.01,
    pop.size = 10,
    dis.mat = F,
    plot.ref = F,
    gen_dif = F,
    purity = F,
    overlap = F,
    correlation = F,
    na.perc.threshold = 80
  )
  saveRDS(x, "DMz24-9199_variety.rds")
  DAP_write_excel(x, info.n = info.file, filename)
}

