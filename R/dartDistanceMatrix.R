
dartDistanceMatrix <- function(genotypic_data,
                               sampleWiseAnalysis = FALSE,
                               ncores = parallel::detectCores() - 1) {

  if (class(genotypic_data)[1] != "matrix") {
    stop("genotypic_data must be a matrix")
  }

  tm = system.time({
    distances = dartDistanceMatrixCpp(genotypic_data,
                                      sampleWiseAnalysis = sampleWiseAnalysis,
                                      nThreads = ncores)
  })
  cat("dartDistanceMatrix took", tm["elapsed"], "seconds\n")
  return(distances)
}
