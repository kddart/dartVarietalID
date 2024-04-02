#' @name dart.differentiation
#' @title Population differentiation
#' @description
#' This function takes one filed sample and estimates its genetic
#' differentiation against each reference
#
#' @param ref A list of genlight objects each containing one reference
#'  [required].
#' @param unknown A genlight object containing a field sample [required].
#' @details
#' @return A \code{data.frame} consisting of differentiation measures for each
#'  population.
#' @export

dart.differentiation <- function(ref,
                                 sam) {

  pop.list <- lapply(ref,rbind,sam)

  dif_res <- lapply(pop.list,dartR:::utils.basic.stats)
  dif_res <- lapply(dif_res,"[[","overall")
  dif_res <- lapply(dif_res,"[",c("Fst","Fstp","Dest","Gst_H"))
  dif_res <- Reduce(rbind,dif_res)



  ret <- data.frame(
    # TargetID = NA,
                    # sample = NA,
                    variety = NA)

  # for each reference
  for (popx in 1:length(ref)) {
    # ret[popx, "sample"] <- as.character(ref[[popx]]$other$ind.metrics$sample[1])
    ret[popx, "variety"] <- as.character(ref[[popx]]$other$ind.metrics$variety[1])
    # ret[popx, "TargetID"] <- as.character(ref[[popx]]$other$ind.metrics$TargetID[1])
  }

  ret <- cbind(ret,dif_res)

  return(ret)

}

