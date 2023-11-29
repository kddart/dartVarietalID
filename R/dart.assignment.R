#' @name dart.assignment
#' @title Population assignment probabilities
#' @description
#' This function takes one sample and estimates its probability of coming from
#' a reference
#
#' @param ref A list of genlight objects each containing one reference
#'  [required].
#' @param unknown A genlight object containing the sample to be assigned to a
#' reference [required].
#' @details
#' @return A \code{data.frame} consisting of assignment probabilities for each
#'  population.
#' @export

dart.assignment <- function(ref,
                            unknown
                            ) {

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
                    Probability = NA)

  # for each reference
  for (popx in 1:length(ref)) {
    popfreq <- frequencies[[popx]]

    loc <-
      as.data.frame(do.call(rbind, strsplit(unname(
        unlist(unknown_pop)
      ), ":")))
    colnames(loc) <- c("a1", "a2")

df_assign <- cbind(loc, popfreq)
# probability of being homozygote for the reference allele (p^2)
df_assign$hom1 <- df_assign$Frequency1 ^ 2
# probability of being homozygote for the alternative allele (q^2)
df_assign$hom2 <- df_assign$Frequency2 ^ 2
# probability of being heterozygote (2pq)
df_assign$het <- 2 * df_assign$Frequency1 * df_assign$Frequency2
# if sample is homozygote for the reference allele, set probability of
# homozygote for the alternative allele and heterozygote to 0
df_assign[which(df_assign$a1 == df_assign$a2 &
                  df_assign$a1 == df_assign$Allele1), c("hom2", "het")] <- 0
# if sample is homozygote for the alternative allele, set probability of
# homozygote for the reference allele and heterozygote to 0
df_assign[which(df_assign$a1 == df_assign$a2 &
                  df_assign$a1 == df_assign$Allele2), c("hom1", "het")] <- 0
# if sample is heterozygote, set probability of homozygote for the reference
# allele and homozygote for the alternative allele to 0
df_assign[which(df_assign$a1 != df_assign$a2), c("hom1", "hom2")] <- 0
# get the probability for each locus
df_assign$prob <- df_assign$hom1 + df_assign$hom2 + df_assign$het
# set NA if the loci of the sample is missing
df_assign[which(is.na(df_assign$a1)), "prob"] <- NA
# set -1 to loci that have probability of 0, i.e. if the reference is fixed
# for one allele and the sample is homozygote for the the other allele
df_assign[which(df_assign$prob == 0), "prob"] <- -1
# get the number of loci that do not have missing data in both, the sample and
# the reference
n_loc_tmp <- rowSums(cbind(is.na(loc$a1),is.na(popfreq$Frequency1)))
n_loc <- sum(n_loc_tmp==0)
# assign probability
ret[popx, "Genotype"] <- as.character(ref[[popx]]$other$ind.metrics$Genotype[1])
ret[popx, "RefType"] <- as.character(ref[[popx]]$other$ind.metrics$RefType[1])
ret[popx, "TargetID"] <- as.character(ref[[popx]]$other$ind.metrics$TargetID[1])
ret[popx, "NumLoci"] <- n_loc
# get the mean probability across all the loci
ret[popx, "Probability"] <- sum(df_assign$prob, na.rm = TRUE)/ n_loc
  }
# order references by probability
ret <- ret[order(-ret$Probability), ]
ret$Probability <- round(ret$Probability, 4) * 100

  return(ret)

}
