library(dartR)
rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
                         pattern = ".rds",
                         full.names = T)

ncores <- 12
tot_sum <- as.data.frame(matrix(ncol = 8))
colnames(tot_sum) <- c(
  "TargetID.sample"  ,
  "sample.sample"   ,
  "TargetID.reference"  ,
  "sample.reference" ,
  "variety.reference" ,
  "NA.percentage"   ,
  "Probability.reference",
  "overlap"
)
overlap_res <- as.list(1:10)
message_parallel <- function(...) {
  system(sprintf('echo "\n%s\n"', paste0(..., collapse = "")))
}
# library(tictoc)
library(dartVarietalID)
n.varieties <- 5
for (i in 3:length(rds_files)) {
  # tic()
  print(i)
  t1 <- readRDS(rds_files[i])
  TargetID.sample <- t1$res_summary$TargetID.sample
  res_tmp <- t1$res_full
  test_pop_ref <- t1$gl.references
  test_pop_sam <- t1$gl.samples
  res_sum <- t1$res_summary

  # Separating populations
  sam_pops_sep <- seppop(test_pop_sam)
  ref_pops_sep <- seppop(test_pop_ref)

  #identify samples with all missing data
  NAs <-
    lapply(sam_pops_sep, function(x) {
      sum(sapply(x@gen, function(e) {
        length(e@NA.posi)
      }))
    })
  #total number of genotypes
  total_geno <- nInd(sam_pops_sep[[1]]) * nLoc(sam_pops_sep[[1]])
  all_NAs <-  which(NAs == total_geno)

  # remove sample if all data is missing
  if (length(all_NAs) > 0) {
    sam_pops_sep <- sam_pops_sep[-(all_NAs)]
  }
  print(paste(length(sam_pops_sep), "samples"))

  # fil_mono <- function(x) {
  #   mat <- as.matrix(x)
  #   hom_ref <- which(apply(mat, 2, \(y)all(y == 0, na.rm = T)) == T)
  #   hom_alt <- which(apply(mat, 2, \(y)all(y == 2, na.rm = T)) == T)
  #   all_na <- which(apply(mat, 2, \(y)all(is.na(y))) == T)
  #   loc.list <- unique(c(hom_ref, hom_alt, all_na))
  #   x <- x[, -loc.list]
  #   return(x)
  # }

  test_gl <- lapply(1:length(sam_pops_sep), function(z) {
    sam_name <- names(sam_pops_sep[z])
    # get the sample to test
    sam_gl <- sam_pops_sep[[z]]
    # get the closest references
    ref_pops <- res_tmp[[sam_name]][1:n.varieties, "variety"]
    ref_gl <-
      dartR::gl.keep.pop(test_pop_ref, pop.list = ref_pops, verbose = 0)
    tog <- rbind(ref_gl, sam_gl)
    # tog$other$loc.metrics <- ref_gl$other$loc.metrics
    # tog@other$loc.metrics.flags$monomorphs <- FALSE
    # tog <- dartR::gl.filter.callrate(tog,
    #                              threshold = 1,
    #                              verbose = 0)
    tog <- gl.filter.monomorphs(tog,verbose = 0)

    return(list(
      sam = sam_name,
      ref = ref_pops[1],
      tog_gl = tog
    ))
  })

  # # res <- parallel::mclapply(X = 1:length(test_gl),
  res <- parallel::mclapply(X = 1:10,

                            FUN = tryCatch({
                              function(x) {
                                overlap_proportion(
                                  tog_gl = test_gl[[x]][[3]],
                                  test.sample = test_gl[[x]][[1]],
                                  test.ref = test_gl[[x]][[2]],

                                  plot = FALSE
                                )
                              }
                            },
                            # error handling
                            error = function(e) {
                              return(NULL)
                            },
                            # warning handling
                            warning = function(w) {
                              return(NULL)
                            }),
                            mc.cores = ncores)
# tic()
  # res <- lapply(1:length(test_gl), function(x) {
  #   tmp <- tryCatch({
  #     overlap_proportion(
  #       tog_gl = test_gl[[x]][[3]],
  #       test.sample = test_gl[[x]][[1]],
  #       test.ref = test_gl[[x]][[2]],
  #       # test.sample = TargetID.sample[x],
  #       # full.report = res_tmp,
  #       # ref = test_pop_ref,
  #       # sam = test_pop_sam,
  #       # n.varieties=10,
  #       plot = F
  #     )
  #   },
  #   # error handling
  #   error = function(e) {
  #     return(NULL)
  #   },
  #   # warning handling
  #   warning = function(w) {
  #     return(NULL)
  #   })
  #
  #   return(tmp)
  #
  # })
  # toc()

  res <- lapply(res, function(x) {
    if (is.null(x)) {
      x <- matrix(nrow = 1, ncol = 2)
      return(x)
    } else{
      return(x)
    }
  })

  res2 <- lapply(res, "[", 2)
  res2 <- unlist(res2)
  overlap_res[[i]] <- as.numeric(res2)
  # toc()

}

overlap_res2 <- saveRDS(overlap_res,"overlap_res.rds")
lapply(overlap_res2,function(x){sum(is.na(x))})
tot_sum <- tot_sum[-1, ]

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
  list.files("/Users/mijangos/DAP_input",
             pattern = "Counts.csv$",
             recursive = TRUE,
             full = TRUE)
# ff <- ff[-6]

rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
                         pattern = ".rds",
                         full.names = T)

corr_res <- as.list(1:10)

for (i in 1:length(ff)) {
  counts.file <- ff[i]
  ordnr <- gsub("_Counts.csv", "", basename(counts.file))
  filename <- file.path(gsub("input", "output/DAP", counts.file))
  filename <- gsub("_Counts.csv$", "_DAP", filename)
  print(filename)

  #	if (file.exists(paste0(filename, ".rds"))) next
  # dir.create(dirname(filename), FALSE, TRUE)

  # ojo: using the fixed references
  info.file <- file.path("/Users/mijangos/DAP_input",
                         gsub(
                           "Counts.csv$",
                           "variety-info-refined.csv",
                           basename(counts.file)
                         ))

  # read in counts file and info file
  ref_sam <- read.dart.counts(counts.file = counts.file,
                              info.file = info.file)
  t1 <- readRDS(rds_files[i])
  res_summary <-  t1$res_summary

  counts <- ref_sam$counts
  cor_df <- res_summary[, c("TargetID.sample", "TargetID.reference")]
  r1 <- apply(cor_df, 1, function(x) {
    summary(lm(counts[, x[1]] ~
                 counts[, x[2]]))$r.squared
  })

  corr_res[[i]] <- r1

}
saveRDS(corr_res,"corr_res.rds")
tot_sum$corr <- r1_fin

saveRDS(tot_sum, "tot_sum.rds")
tot_sum <- readRDS("tot_sum.rds")

more100 <- tot_sum[which(tot_sum$Probability.reference > 100), ]
more100 <-
  more100[order(more100$Probability.reference, decreasing = T), ]
more100_pur <- more100[which(more100$purity)]

library(dartVarietalID)
rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
                         pattern = ".rds",
                         full.names = T)
crop_names <- lapply(basename(rds_files), function(x) {
  gsub("_DAP.rds$",
       "", x)
})
crop_names <- unlist(crop_names)
p1_fin <- NULL
ncores <- 12
purity_list <- as.list(1:10)
names(purity_list) <- crop_names
for (i in 8:length(ff)) {
  counts.file <- ff[i]
  ordnr <- gsub("_Counts.csv", "", basename(counts.file))
  filename <- file.path(gsub("input", "output/DAP", counts.file))
  filename <- gsub("_Counts.csv$", "_DAP", filename)
  print(filename)

  #	if (file.exists(paste0(filename, ".rds"))) next
  # dir.create(dirname(filename), FALSE, TRUE)

  # ojo: using the fixed references
  info.file <- file.path("/Users/mijangos/DAP_input",
                         gsub(
                           "Counts.csv$",
                           "variety-info-refined.csv",
                           basename(counts.file)
                         ))

  # read in counts file and info file
  ref_sam <- read.dart.counts(counts.file = counts.file,
                              info.file = info.file)
  t1 <- readRDS(rds_files[i])
  res_summary <-  t1$res_summary

  genotypic_counts <- ds14.genotypic(ds14.read(counts.file))
  infoFile <- readTargetInfoFile(file = info.file)
  assigned_test_reference <- res_summary$variety.reference
  names(assigned_test_reference) <- res_summary$TargetID.sample
  res_purity <- calculatePurity(genotypic_counts,
                                infoFile,
                                assigned_test_reference,
                                ncores)
  purity_list[[i]] <- unlist(unname(res_purity))

}

saveRDS(purity_list, "purity_list.rds")

tot_sum$purity <- p1_fin

library(ggplot2)
library(ggpointdensity)
library(viridis)
res_summary$overlap <- res_sum$overlap
res_summary$corr <- res_sum$corr

# tot_sum$Probability.reference <- scales::rescale(tot_sum$Probability.reference, to = c(0,100) )
ggplot(res_summary,
       aes(x = res_summary$Probability.reference,
           y = res_summary$corr)) +
  geom_pointdensity() +
  scale_color_viridis() +
  geom_smooth()

ggplot(tot_sum,
       aes(x = tot_sum$Probability.reference,
           y = tot_sum$corr)) +
  geom_pointdensity() +
  scale_color_viridis() +
  geom_smooth()

ggplot(tot_sum, aes(x = tot_sum$corr,
                    y = tot_sum$overlap)) +
  geom_point() +
  geom_smooth()


ggplot(more100,
       aes(x = more100$Probability.reference,
           y = more100$corr)) +
  geom_pointdensity() +
  # geom_point() + +
  scale_color_viridis() +
  geom_smooth()

ggplot(more100,
       aes(x = more100$Probability.reference,
           y = more100$purity)) +
  geom_pointdensity() +
  # geom_point() + +
  scale_color_viridis()
# +
#   geom_smooth()

ggplot(res_summary,
       aes(x = res_summary$corr,
           y = 1 - res_summary$purityPercent)) +
  geom_point() +
  geom_smooth()




r1 <- res_summary[which(res_summary$overlap < 5 &
                          res_summary$Probability.reference > 90), ]
r1 <- r1[order(r1$Probability.reference, decreasing = T), ]




dart.differentiation <- function(ref,
                                 sam) {
  pop.list <- lapply(ref, rbind, sam)

  dif_res <- lapply(pop.list, dartR:::utils.basic.stats)
  dif_res <- lapply(dif_res, "[[", "overall")
  dif_res <- lapply(dif_res, "[", c("Fst", "Fstp", "Dest", "Gst_H"))
  dif_res <- Reduce(rbind, dif_res)



  ret <- data.frame(# TargetID = NA,
    # sample = NA,
    variety = NA)

  # for each reference
  for (popx in 1:length(ref)) {
    # ret[popx, "sample"] <- as.character(ref[[popx]]$other$ind.metrics$sample[1])
    ret[popx, "variety"] <-
      as.character(ref[[popx]]$other$ind.metrics$variety[1])
    # ret[popx, "TargetID"] <- as.character(ref[[popx]]$other$ind.metrics$TargetID[1])
  }

  ret <- cbind(ret, dif_res)

  return(ret)

}
