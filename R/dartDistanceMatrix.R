#' dartDistanceMatrix
#'
#' Calls the DArT implementation of distance matrix
#' @export
dartDistanceMatrix <- function(genotypic_data,
                               sampleWiseAnalysis = F,
                               ncores = parallel::detectCores()) {
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
