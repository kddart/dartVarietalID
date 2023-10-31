#' @name dart.assignment
#' @title Population assignment probabilities
#' @description
#' This function takes one sample and estimates its probability of coming from
#' a reference
#
#' @param ref Genlight object containing the reference data [required].
#' @param unknown Genlight object containing the individual to be assigned to a
#' reference [required].
#' @details
#' This function is a re-implementation of the function multilocus_assignment
#'  from package gstudio.
#'  Description of the method used in this function can be found at:
#' https://dyerlab.github.io/applied_population_genetics/population-assignment.html
#' @return A \code{data.frame} consisting of assignment probabilities for each
#'  population.
#' @export

dart.assignment <- function(ref,
                            unknown) {

  require(dartR)

  gl2alleles <- function (gl) {
    x <- as.matrix(gl)
    homs1 <-
      paste(substr(gl@loc.all, 1, 1), "/", substr(gl@loc.all, 1, 1), sep = "")
    hets <- gl@loc.all
    homs2 <-
      paste(substr(gl@loc.all, 3, 3), "/", substr(gl@loc.all, 3, 3), sep = "")
    xx <- matrix(NA, ncol = ncol(x), nrow = nrow(x))
    for (i in 1:nrow(x)) {
      for (ii in 1:ncol(x)) {
        inp <- x[i, ii]
        if (!is.na(inp)) {
          if (inp == 0)
            xx[i, ii] <- homs1[ii]
          else if (inp == 1)
            xx[i, ii] <- hets[ii]
          else if (inp == 2)
            xx[i, ii] <- homs2[ii]
        } else{
          xx[i, ii] <- NA
        }
      }
    }
    xx <- gsub("/", ":", xx)
    return(xx)
  }

  unknown_pop <- data.frame(gl2alleles(unknown))

  pop_names <- names(ref)

  pop_list <- ref
  gl_alleles <- do.call(rbind, strsplit(ref[[1]]$loc.all, "/"))

  frequencies <- lapply(pop_list, function(y) {
    freq_allele <- gl.alf(y)
    freqs_gl <-
      data.frame(
        Allele1 = gl_alleles[, 1],
        Allele2 = gl_alleles[, 2],
        Frequency1 = freq_allele[, 1],
        Frequency2 = freq_allele[, 2]
      )
    return(freqs_gl)
  })

  ret <- data.frame(TargetID = pop_names,
                    Genotype = NA,
                    RefType = NA,
                    NumLoci = NA,
                    Probability = 0)

  for (popx in 1:length(ref)) {
    prob <- 1
    popfreq <- frequencies[[popx]]

    loc <-
      as.data.frame(do.call(rbind, strsplit(unname(
        unlist(unknown_pop)
      ), ":")))
    colnames(loc) <- c("a1", "a2")

    df_assign <- cbind(loc, popfreq)
    df_assign$hom1 <- df_assign$Frequency1 ^ 2
    df_assign$hom2 <- df_assign$Frequency2 ^ 2
    df_assign$het <- 2 * df_assign$Frequency1 * df_assign$Frequency2

    df_assign[which(df_assign$a1 == df_assign$a2 &
                      df_assign$a1 == df_assign$Allele1), c("hom2", "het")] <-
      0

    df_assign[which(df_assign$a1 == df_assign$a2 &
                      df_assign$a1 == df_assign$Allele2), c("hom1", "het")] <-
      0

    df_assign[which(df_assign$a1 != df_assign$a2), c("hom1", "hom2")] <-
      0

    df_assign$prob <- df_assign$hom1 + df_assign$hom2 + df_assign$het
    df_assign[which(is.na(df_assign$a1)), "prob"] <- NA
    df_assign[which(df_assign$prob == 0), "prob"] <- -1
    n_loc_tmp <- rowSums(cbind(is.na(loc$a1),is.na(popfreq$Frequency1)))
    n_loc <- sum(n_loc_tmp==0)
    df_assign$prob <- df_assign$prob / n_loc

    # assign probability
    ret[popx, "Genotype"] <-
      as.character(ref[[popx]]$other$ind.metrics$Genotype[1])
    ret[popx, "RefType"] <-
      as.character(ref[[popx]]$other$ind.metrics$RefType[1])
    ret[popx, "TargetID"] <-
      as.character(ref[[popx]]$other$ind.metrics$TargetID[1])
    ret[popx, "NumLoci"] <- n_loc
    ret[popx, "Probability"] <- sum(df_assign$prob, na.rm = TRUE)
  }

  ret <- ret[order(-ret$Probability), ]
  ret$Posterior <- ret$Probability / sum(ret$Probability)
  ret$Posterior <- round(ret$Posterior, 5)
  ret$Probability <- round(ret$Probability, 5) * 100

  return(ret)

}
