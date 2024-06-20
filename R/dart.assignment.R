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
                            sam,
                            unknown) {

  # sam_pop <- sam[[which(names(sam)==unknown$other$ind.metrics$variety[1])]]
  sam_pop <- sam[[which(names(sam)==unknown$other$ind.metrics$TargetID[1])]]

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

  pop.het_fun <- function(df) {
    df <- as.matrix(df)
    # n_loc <- apply(df, 2, function(y) {
    #   sum(!is.na(y))
    # })
    q_freq <- colMeans(df, na.rm = TRUE) / 2
    p_freq <- 1 - q_freq
    He.loc <- 2 * p_freq * q_freq
    # n_ind <- apply(df, 2, function(y) {
    #   sum(!is.na(y))
    # })
    # ### CP ### Unbiased He (i.e. corrected for sample size) hard
    # # coded for diploid
    # uHe.loc <-
    #   (2 * as.numeric(n_ind) / (2 * as.numeric(n_ind) - 1)) * He.loc
    # uHe <- mean(uHe.loc,na.rm=TRUE)
    He <- mean(He.loc,na.rm =TRUE)
    return(He)
  }

  # pop.het_fun <- function(df) {
  #   out <- colMeans(as.matrix(df) == 1, na.rm = T)
  #   He <- mean(out,na.rm =TRUE)
  #   return(He)
  # }

  Ho_sam <- pop.het_fun(sam_pop)

  # Ho <- lapply(1:length(ref),function(y){
  #   colMeans(as.matrix(ref[[y]]) == 1, na.rm = TRUE)
  # })
  #
  # HoA1 <- lapply(1:length(ref),function(y){
  #   colMeans(as.matrix(ref[[y]]) == 0, na.rm = TRUE)
  # })
  #
  # HoA2 <- lapply(1:length(ref),function(y){
  #   colMeans(as.matrix(ref[[y]]) == 2, na.rm = TRUE)
  # })

  ret <- data.frame(TargetID = pop_names,
                    sample = NA,
                    variety = NA,
                    NumLoci = NA,
                    Probability = NA,
                    Probability_corr = NA,
                    Het = NA)

  # for each reference
  for (popx in 1:length(ref)) {

    Ho_pop <- pop.het_fun(df=ref[[popx]])

    popfreq <- frequencies[[popx]]

    loc <-
      as.data.frame(do.call(rbind, strsplit(unname(
        unlist(unknown_pop)
      ), ":")))
    colnames(loc) <- c("a1", "a2")

df_assign <- cbind(loc, popfreq)
# Ho_pop <- mean(Ho[[popx]],na.rm = TRUE)
# df_assign$HoA1 <- HoA1[[popx]]
# df_assign$HoA2 <- HoA2[[popx]]
# df_assign$Ho2 <- (1-df_assign$Ho)/df_assign$Ho
# probability of being homozygote for the reference allele (p^2)
# df_assign$hom1 <- (1-(df_assign$Frequency1 ^ 2)) / df_assign$Ho
df_assign$hom1 <- df_assign$Frequency1 ^ 2

# / df_assign$Ho
# probability of being homozygote for the alternative allele (q^2)
# df_assign$hom2 <- (1-(df_assign$Frequency2 ^ 2)) / df_assign$Ho
df_assign$hom2 <- df_assign$Frequency2 ^ 2

# / df_assign$Ho
# probability of being heterozygote (2pq)
df_assign$het <- (2 * df_assign$Frequency1 * df_assign$Frequency2) / 0.5
# df_assign$het <- (2 * df_assign$Frequency1 * df_assign$Frequency2)

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
 # df_assign[which(df_assign$prob == 0), "prob"] <- -1
 # df_assign[which(df_assign$Ho>0.05 & df_assign$Ho<0.3), "prob"] <- NA

# get the number of loci that do not have missing data in both, the sample and
# the reference
n_loc_tmp <- rowSums(cbind(is.na(loc$a1),is.na(popfreq$Frequency1)))
n_loc <- sum(n_loc_tmp==0)
# assign probability
ret[popx, "sample"] <- as.character(ref[[popx]]$other$ind.metrics$sample[1])
ret[popx, "variety"] <- as.character(ref[[popx]]$other$ind.metrics$variety[1])
ret[popx, "TargetID"] <- as.character(ref[[popx]]$other$ind.metrics$TargetID[1])
ret[popx, "NumLoci"] <- n_loc
# get the mean probability across all the loci
het_correction <- 1 - ((Ho_pop+Ho_sam)/2)
# het_correction <- 1 - (Ho_pop)

# ret[popx, "Probability"] <- (sum(df_assign$prob, na.rm = TRUE)/ n_loc) / het_correction
ret[popx, "Probability"] <- (sum(df_assign$prob, na.rm = TRUE)/ n_loc)
ret[popx, "Probability_corr"] <- (sum(df_assign$prob, na.rm = TRUE)/ n_loc) / het_correction
ret[popx, "Het"] <- ((Ho_pop+Ho_sam)/2)

# ret[popx, "Probability2"] <- ret[popx, "Probability"]/(1-Ho_pop)
  }
# order references by probability
ret <- ret[order(-ret$Probability_corr), ]
ret$Probability <- round(ret$Probability * 100, 2)
ret$Probability_corr <- round(ret$Probability_corr * 100, 2)

# head(ret)
  return(ret)

}

