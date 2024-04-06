
overlap_proportion <- function(test.sample,
                               full.report,
                               ref,
                               sam,
                               n.varieties=10,
                               plot = TRUE){


  # get the sample to test from the full report
  full_sample <- full.report[[test.sample]]
  # get the n closest references from the sample to test
  ref_pops <- full_sample[1:n.varieties,"variety"]
  # ref_pops <- full_sample[1,"variety"]

  # get the closest references
  ref_test <- dartR::gl.keep.pop(ref,pop.list = ref_pops,verbose=0)
  # get the sample to test
  sam_test <- dartR::gl.keep.pop(sam,pop.list = test.sample,verbose=0)
  # merging sample to test and closest references
  tog <- rbind(ref_test,sam_test)
  tog$other$loc.metrics <- ref_test$other$loc.metrics
  tog@other$loc.metrics.flags$monomorphs <- FALSE

  tog <- dartR::gl.filter.callrate(tog,
                            threshold = 1,
                            verbose = 0)
    # if unix
    if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
      pca <- adegenet::glPca(tog,
                              nf = 3,
                              parallel = FALSE,
                              loadings = FALSE)
    }

    ## if windows
    if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
      pca <- adegenet::glPca(tog,
                              nf = 3,
                              parallel = TRUE,
                              loadings = FALSE)
    }

pcoa_scores <- as.data.frame(pca$scores)
pcoa_scores$pop <- as.character(pop(tog))
sam_p <- which(pcoa_scores$pop==test.sample)
ref_p <- which(pcoa_scores$pop==full_sample[1,"variety"])
# get the center of ellipsoids
center_s <- unlist(unname(colMeans(pcoa_scores[sam_p,1:3])))
center_r <- unlist(unname(colMeans(pcoa_scores[ref_p,1:3])))
#generating the ellipsoids
ellipseSam <- rgl::ellipse3d(cov(pcoa_scores[sam_p,1:3]),
                        centre = center_s,
                        level = 0.999999)
ellipseRef <- rgl::ellipse3d(cov(pcoa_scores[ref_p,1:3]),
                        centre = center_r,
                        level = 0.999999)
# get the semi-axes of an ellipsoid
axes_ref <- estimateSemiAxes(mesh=ellipseRef,
                             centroid = center_r)
axes_sam <- estimateSemiAxes(mesh=ellipseSam,
                             centroid = center_s)
# get a box coordinates to put random points
box_r <- getBoundingBox(center_r,axes_ref)
box_s <- getBoundingBox(center_s,axes_sam)
# Number of points to generate
n <- 10000
# random points within a box containing the ellipsoid of the sample
lower_s <- box_s$min
upper_s <- box_s$max
points_s <- matrix(runif(n * 3, lower_s, upper_s), ncol = 3, byrow = TRUE)
# random points within a box containing the ellipsoid of the reference
lower_r <- box_r$min
upper_r <- box_r$max
points_r <- matrix(runif(n * 3, lower_r, upper_r), ncol = 3, byrow = TRUE)

means_r <- colMeans(t(ellipseRef$vb[1:3,]))
covariance_r <- cov(t(ellipseRef$vb[1:3,]))

means_s <- colMeans(t(ellipseSam$vb[1:3,]))
covariance_s <- cov(t(ellipseSam$vb[1:3,]))

# Getting the points of the box that are inside the sample ellipse
Z_s <- SIBER::pointsToEllipsoid(X=points_s,
                         Sigma= covariance_s,
                         mu=means_s )
test_s <- SIBER::ellipseInOut(Z_s, p = 0.7)
points_s <- points_s[test_s,]

# Getting the points of the box that are inside the reference ellipse
Z_r <- SIBER::pointsToEllipsoid(X=points_r,
                                Sigma= covariance_r,
                                mu=means_r)
test_r <- SIBER::ellipseInOut(Z_r, p = 0.7)
points_r <- points_r[test_r,]

#get the points of the sample that are inside the reference ellipse
Z_sINr <- SIBER::pointsToEllipsoid(X=points_s,
                         Sigma= covariance_r,
                         mu=means_r )
test_sINr <- SIBER::ellipseInOut(Z_sINr, p = 0.7)
points_sINr <- points_s[test_sINr,]

#get the points of the reference that are inside the sample ellipse
Z_rINs <- SIBER::pointsToEllipsoid(X=points_r,
                                   Sigma= covariance_s,
                                   mu=means_s )
test_rINs <- SIBER::ellipseInOut(Z_rINs, p = 0.7)
points_rINs <- points_r[test_rINs,]

if(nrow(points_sINr) > 0 | nrow(points_rINs) > 0 ){
  prop_sINr <- round((nrow(points_sINr)/nrow(points_s))*100,2)
  prop_rINs <- round((nrow(points_rINs)/nrow(points_r))*100,2)
  prop <- max(c(prop_sINr,prop_rINs))
}else{
  prop <- 0
}

if(plot){
p <- plot_ly() %>%
  add_trace(data = pcoa_scores,
            x = pcoa_scores[,1],
            y = pcoa_scores[,2],
            z = pcoa_scores[,3],
            type = "scatter3d",
            mode = 'markers'
            ,
            marker = list(size = 6),
            color = pcoa_scores[,"pop"],
            colors = c("black" ,polychrome(nPop(tog)))
  )%>%
  add_trace(x = ellipseSam$vb[1,],
            y = ellipseSam$vb[2,],
            z = ellipseSam$vb[3,],
            type = 'mesh3d',
            alphahull = 0,
            opacity = 0.4) %>%
  add_trace(x = ellipseRef$vb[1,],
            y = ellipseRef$vb[2,],
            z = ellipseRef$vb[3,],
            type = 'mesh3d',
            alphahull = 0,
            opacity = 0.4) %>%
  add_trace(data = points_sINr,
            x = points_sINr[,1],
            y = points_sINr[,2],
            z = points_sINr[,3],
            type = "scatter3d",
            mode = 'markers'
  ) %>%
  add_trace(data = prop_rINs,
            x = points_rINs[,1],
            y = points_rINs[,2],
            z = points_rINs[,3],
            type = "scatter3d",
            mode = 'markers'
  )
# %>%
#   add_trace(data = points_s,
#             x = points_s[,1],
#             y = points_s[,2],
#             z = points_s[,3],
#             type = "scatter3d",
#             mode = 'markers'
#   )

print(p)

}

return(cbind(test.sample,prop))

}
