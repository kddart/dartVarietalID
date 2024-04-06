ncores = parallel::detectCores() -1
pop.size = 10
dis.mat = F
plot.ref = F
gen_dif = F
purity = T
overlap = T
correlation =T
na.perc.threshold = 50


which(names(sam_pops_sep)=="3215600")
which(names(top_ind)=="3215600")
which(names(ref_pops_sep)=="Selian 13")
t_sam <- sam_pops_sep[[57]]
t_ref <- ref_pops_sep[[32]]
t_tog <- rbind(t_sam,t_ref)
t_ref$other$loc.metrics
t_tog$other$loc.metrics <- t_ref$other$loc.metrics
t_tog@other$loc.metrics.flags$monomorphs <- TRUE
res_het <- gl.report.heterozygosity(t_tog)
res_het$He


counts.file = counts.file
info.file = info.file
ncores = parallel::detectCores()
pop.size = 10
dis.mat = FALSE
plot.ref = FALSE
gen_dif = FALSE
purity = FALSE
overlap = FALSE
correlation = FALSE
na.perc.threshold = 50

plot(res_summary$Probability.reference,res_summary$overlap)
plot(res_summary$Probability.reference,res_summary$corr)
plot(res_summary$Probability.reference,res_summary$purityPercent)

library(dartVarietalID)
t1 <- readRDS("/Users/mijangos/DAP_output/DAP/DCob23-7823_DAP.rds")
TargetID.sample <- t1$res_summary$TargetID.sample
res_tmp <- t1$res_full
test_pop_ref <- t1$gl.references
test_pop_sam <- t1$gl.samples
res_sum <- t1$res_summary

x <- which(TargetID.sample=="3043235")

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

 res_sum[which(res_sum$overlap < 5 & res_sum$Probability.reference >95 ),]

plot(res_sum$Probability.reference,res_sum$overlap)

DCas23_7954_DAP <- res_sum
DCob22_7521_DAP <- res_sum
DCob23_7823_DAP <- res_sum

DCob23_7823_DAP$corr <- res_summary$corr
DCob23_7823_DAP[which(DCob23_7823_DAP$overlap < 2 & DCob23_7823_DAP$Probability.reference >98 ),]


library(ggplot2)
ggplot(DCas23_7954_DAP,aes(x= DCas23_7954_DAP$Probability.reference,
                           y= DCas23_7954_DAP$overlap))+
geom_point() +
  geom_smooth()

ggplot(DCob22_7521_DAP,aes(x= DCob22_7521_DAP$Probability.reference,
                           y= DCob22_7521_DAP$overlap))+
  geom_point() +
  geom_smooth()

ggplot(DCob23_7823_DAP,aes(x= DCob23_7823_DAP$Probability.reference,
                           y= DCob23_7823_DAP$corr))+
  geom_point() +
  geom_smooth()

res_tog <- rbind(DCas23_7954_DAP,DCob22_7521_DAP,DCob23_7823_DAP)

ggplot(res_tog,aes(y= res_tog$Probability.reference,
                           x= res_tog$overlap))+
  geom_point() +
  geom_smooth() +
  guides(fill=guide_legend(title="New Legend Title"))

library(dartVarietalID)
correlation <- TRUE
counts.file <- "/Users/mijangos/DAP_input/DCob22-7521_moreOrders_Counts.csv"
info.file <- "/Users/mijangos/DAP_input/DCob22-7521_moreOrders_variety-info-refined.csv"
# read in counts file and info file
ref_sam <- read.dart.counts(counts.file = counts.file,
                            info.file = info.file)


t1 <- readRDS("/Users/mijangos/DAP_output/DAP/DCob23-7823_DAP.rds")
TargetID.sample <- t1$res_summary$TargetID.sample
res_tmp <- t1$res_full
test_pop_ref <- t1$gl.references
test_pop_sam <- t1$gl.samples
res_summary <- t1$res_summary



if(correlation){

  counts <- ref_sam$counts
  cor_df <- res_summary[,c("TargetID.sample","TargetID.reference")]
  r1 <-apply(cor_df,1,function(x){
    summary(lm(counts[,x[1]]~
                 counts[,x[2]]))$r.squared
  })

  res_summary$corr <- r1


}

library(ggplot2)
ggplot(res_summary,aes(x= res_summary$Probability.reference,
                       y= res_summary$overlap))+
  geom_point() +
  geom_smooth()

ggplot(res_summary,aes(x= res_summary$corr,
                       y= res_summary$overlap))+
  geom_point() +
  geom_smooth()


ggplot(res_summary,aes(x= res_summary$Probability.reference,
                       y= res_summary$corr))+
  geom_point() +
  geom_smooth()

ggplot(res_summary,aes(x= res_summary$Probability.reference,
                       y= 1-res_summary$purityPercent))+
  geom_point() +
  geom_smooth()

ggplot(res_summary,aes(x= res_summary$corr,
                       y= 1-res_summary$purityPercent))+
  geom_point() +
  geom_smooth()




r1 <- res_summary[which(res_summary$overlap < 5 &
                          res_summary$Probability.reference >95),]
r1 <- r1[order(r1$Probability.reference,decreasing = T),]
r1






