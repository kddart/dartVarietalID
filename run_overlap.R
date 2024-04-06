rds_files <-  list.files("/Users/mijangos/DAP_output/DAP",
           pattern = ".rds",
           full.names = T)

tot_sum <- as.data.frame(matrix(ncol=8))
colnames(tot_sum) <- c("TargetID.sample"  ,     "sample.sample"   ,
                        "TargetID.reference"  ,  "sample.reference" ,
                        "variety.reference" ,   "NA.percentage"   ,
                        "Probability.reference","overlap")
library(dartVarietalID)
for(i in 1:length(rds_files)){
  t1 <- readRDS(rds_files[i])
  TargetID.sample <- t1$res_summary$TargetID.sample
  res_tmp <- t1$res_full
  test_pop_ref <- t1$gl.references
  test_pop_sam <- t1$gl.samples
  res_sum <- t1$res_summary

  # x <- which(TargetID.sample=="3043235")

  res <- lapply(1:length(TargetID.sample),function(x){
    # res <- lapply(900:914,function(x){
    tmp <- tryCatch({
      overlap_proportion(test.sample= TargetID.sample[x],
                         full.report = res_tmp,
                         ref = test_pop_ref,
                         sam = test_pop_sam,
                         n.varieties=5,
                         plot = T)
    },
    # error handling
    error = function(e) {
      return(NULL)
    },
    # warning handling
    warning = function(w) {
      return(NULL)
    })

    return(tmp)

  })

  res <- lapply(res,function(x){
    if(is.null(x)){
      x <- matrix(nrow = 1,ncol=2)
      return(x)
    }else{
      return(x)
    }
  })

  res2 <- lapply(res,"[",2)
  res2 <- unlist(res2)
  res_sum$overlap <- as.numeric(res2)

  tot_sum <- rbind(tot_sum,res_sum)


}
