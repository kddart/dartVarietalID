
plot.stats <- function(rds.file,
                       corr.data = NULL,
                       diff.data = NULL,
                       overlap.data = NULL,
                       purity.data = NULL,
                       n.ref = 10,
                       title_plot = NULL){
if(!is.null(title_plot)){
  crop_name <- title_plot
}else{
  crop_name <- gsub(".rds$","", basename(rds.file))
}


  t1 <- readRDS(rds.file)
  res_sum <-  t1$res_full
  names_sam <- names(res_sum)
  res_sum <- lapply(1:length(res_sum), function(x){
    tmp <- res_sum[[x]]
    tmp <- tmp[1:n.ref,]
    tmp$TargetID.sample <- names_sam[x]
    tmp$TargetID.reference <- tmp$variety
    return(tmp)
  })
  res_sum <-data.table::rbindlist(res_sum)

  p1 <- p2 <- p3 <- p4 <- NULL
 if(!is.null(overlap.data)){
   res_sum$overlap <- overlap.data
   p1 <-   ggplot(res_sum,aes(x= Probability_corr_scaled,
                              y= overlap/100))+
     geom_pointdensity()+
     scale_color_viridis() +
     geom_smooth()  +
     labs(title=crop_name) +
     xlim(0.60,1)+
     ylim(0,1) +
     xlab("Probability score") +
     ylab("Overlap") +
     theme_bw() +
     theme(plot.title = element_text(size=14)) +
     theme(legend.position="none")
 }

  if(!is.null(corr.data)){
    res_sum$corr <- corr.data
  p2 <-   ggplot(res_sum,aes(x= Probability_corr_scaled,
                             y= corr))+
    geom_pointdensity()+
    scale_color_viridis() +
    geom_smooth()  +
    xlim(0.25,1)+
    ylim(0,1) +
    xlab("Probability score") +
    ylab("Correlation") +
    theme_bw() +
    theme(legend.position="none") +
    labs(title=crop_name)
  }

  if(!is.null(diff.data)){
  res_sum$diff <- diff.data
  p3 <- ggplot(res_sum,aes(x= Probability_corr_scaled,
                           y= diff))+
    geom_pointdensity()+
    scale_color_viridis() +
    geom_smooth()  +
    xlim(0.25,1)+
    ylim(0,0.3) +
    xlab("Probability score") +
    ylab("Nei's distance") +
    theme_bw() +
    theme(legend.position="none")
  }

  if(!is.null(purity.data)){
  res_sum$purity <- purity.data
  p4 <- ggplot(res_sum,aes(x= Probability_corr_scaled,
                           y= purity))+
    geom_pointdensity()+
    scale_color_viridis() +
    geom_smooth()  +
    xlim(0.50,1)+
    ylim(0.7,1.45) +
    xlab("Probability score") +
    ylab("Purity score") +
    theme_bw()
  }

  # p5 <- p1 + p2 + p3 + p4
  p5 <- p2 + p3 + p4


  return(p5)

}

