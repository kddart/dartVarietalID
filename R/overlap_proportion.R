
overlap_proportion <- function(tog_gl,
                               test.sample,
                               test.ref,
                               # ref,
                               # sam,
                               # n.varieties=10,
                               plot = TRUE){

# print(test.sample)
  # tog_gl = test_gl[[5]][[3]]
          pcoa <- adegenet::glPca(tog_gl,
                                  nf = 3,
                                  loadings = FALSE)

pcoa_scores <- as.data.frame(pcoa$scores)
# pcoa_scores <- as.data.frame(apply(pcoa_scores,2,scales::rescale))
pcoa_scores$pop <- as.character(pop(tog_gl))
sam_p <- which(pcoa_scores$pop==test.sample)
ref_p <- which(pcoa_scores$pop==test.ref)
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
points_left <- sum(nrow(points_sINr), nrow(points_rINs),na.rm = T)
if(points_left > 0 ){
  prop_sINr <- round((nrow(points_sINr)/nrow(points_s))*100,2)
  prop_rINs <- round((nrow(points_rINs)/nrow(points_r))*100,2)
  prop <- max(c(prop_sINr,prop_rINs),na.rm = T)
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
            colors = c("black" ,polychrome(nPop(tog_gl)))
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
  add_trace(data = points_rINs,
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
