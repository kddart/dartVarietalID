#' @name counts2geno
#' @title Convert count data to genotypes
#' @param count.data An object created via the function read.dart.counts
#' [required].
#' @param pop.size Number of individuals to be simulated [default = 10].
#' @return A genlight object.
#' @export

counts2geno <- function(count.data,
                        pop.size = 10) {
  cdata <- count.data$counts
  cdata_pops <- count.data$ind_metrics$TargetID
  nsnp <- nrow(cdata) / 2
  npop <- ncol(cdata)
  esl <- seq(2, nrow(cdata), 2)
  loc_metrics <- count.data$loc_metrics[esl,]
  loc_names <- loc_metrics$uid
  ind_metrics <- count.data$ind_metrics

  # calculating allele frequencies

  allele_1 <- cdata[esl,]
  allele_2 <- cdata[-esl,]
  p_freq <- allele_1 / (allele_1 + allele_2)
  # assigning 0 to missing data (when counts are 0/0)
  freq_0_0 <- is.na(p_freq)
  p_freq[freq_0_0] <- 0

  # make chromosomes
  #to hack package checking...
  make_chr <- function() {
  }

  Rcpp::cppFunction(
    plugins = "cpp11",

    'StringVector make_chr(int j, NumericVector q) {
    StringVector out(j);
    int size = 1;
    IntegerVector x = IntegerVector::create(1,0);
    bool rep = false;
for (int i = 0; i < j; i++) {
std::ostringstream temp;
for (int z = 0; z < q.length(); z++) {
NumericVector p = NumericVector::create(q[z],1-q[z]);
   temp << sample(x, size, rep, p);
  }
      out[i] = temp.str();
    }
    return out;
  }'
  )

  pop_list <- as.list(1:npop)
  for (pop_n in 1:length(pop_list)) {
    chr_temp <- make_chr(j = pop.size * 2, q = p_freq[, pop_n])
    chr_pops <- split(chr_temp, 1:2)

    pop <- as.data.frame(matrix(ncol = 4, nrow = pop.size))
    # first column stores reference
    pop[, 1] <- ind_metrics$reference[pop_n]
    # second column stores TargetID
    pop[, 2] <- ind_metrics$TargetID[pop_n]
    # third and fourth columns stores chromosomes
    pop[, 3] <- chr_pops[1]
    pop[, 4] <- chr_pops[2]
    # fifth column stores sample
    pop[, 5] <- ind_metrics$sample[pop_n]
    # six column stores variety
    pop[, 6] <- ind_metrics$variety[pop_n]

    pop_list[[pop_n]] <- pop
  }

  df_genotypes <- data.table::rbindlist(pop_list)
  df_genotypes$id <-
    paste0(unlist(unname(df_genotypes[, 2])), "_", unlist(lapply(pop.size, function(x) {
      1:x
    })))

  # make genotypes
  #to hack package checking...
  make_geno <- function() {
  }

  Rcpp::cppFunction(
    plugins = "cpp11",
    'List make_geno(StringMatrix mat) {
    int ind = mat.nrow();
    int loc = strlen(mat(0,0));
    List out(ind);
for (int i = 0; i < ind; i++) {
 std::string chr1 (mat(i,0));
 std::string chr2 (mat(i,1));
 StringVector temp(loc);
for (int j = 0; j < loc; j++) {
 StringVector geno = StringVector::create(chr1[j],chr2[j]);
    temp[j] = collapse(geno);
  }
      out[i] = temp;
    }
    return out;
  }'
  )

  plink_temp <- as.matrix(df_genotypes[, 3:4])
  plink_ped <- make_geno(plink_temp)
  plink_ped_2 <- lapply(plink_ped, function(x) {
    x[x == "11"] <- 2
    x[x == "00"] <- 0
    x[x == "01"] <- 1
    x[x == "10"] <- 1

    return(x)

  })

  loc.names <- loc_names
  n.loc <- length(loc.names)

  misc.info <- lapply(1:5, function(i) {
    NULL
  })
  names(misc.info) <-
    c("reference", "TargetID", "sample", "variety", "ID")
  res <- list()
  temp <- as.data.frame(df_genotypes[, c("V1", "V2", "V5", "V6", "id")])

  for (i in 1:length(misc.info)) {
    misc.info[[i]] <- temp[, i]
  }

  txt <-
    lapply(plink_ped_2, function(e)
      suppressWarnings(as.integer(e)))

  res <-
    c(res, lapply(txt, function(e)
      new(
        "SNPbin", snp = e, ploidy = 2L
      )))

  res <- new("genlight", res, ploidy = 2L)

  indNames(res) <- misc.info$ID
  pop(res) <- misc.info$TargetID

  sep_res <- seppop(res)
  sep_res <- sep_res[order(match(names(sep_res), cdata_pops))]
  freq_0_0 <- t(freq_0_0)

  # setting loci with counts 0/0 as NA
  for (row_n in 1:nrow(freq_0_0)) {
    pop_tmp <- sep_res[[row_n]]
    matrix_tmp <- as.matrix(pop_tmp)
    NA_mat <- unname(which(freq_0_0[row_n, ] == TRUE))
    matrix_tmp[, NA_mat] <- NA
    sep_res[[row_n]]@gen <- lapply(1:nrow(matrix_tmp), function(i) {
      new("SNPbin", as.integer(matrix_tmp[i,]))
    })
  }

  res <- Reduce(rbind, sep_res)
  ploidy(res) <- 2
  locNames(res) <- loc.names
  res$other$ind.metrics <- as.data.frame(misc.info)
  res$other$loc.metrics <- loc_metrics
  # using "G/C" as dummy alleles
  res$loc.all <- rep("G/C", nLoc(res))
  res <- reset.flags(res)
  class(res) <- "dartR"

  return(res)

}
