run.heatmap <- function(rds.file){

  t1 <- readRDS(rds.file)
  test_pop_ref <- t1$gl.references
  test_pop_ref_2 <- test_pop_ref
  pop(test_pop_ref_2) <-
    paste0(
      test_pop_ref_2$other$ind.metrics$TargetID,
      "_",
      test_pop_ref_2$other$ind.metrics$variety
    )

  test_pop_ref_2$other$ind.metrics$variety <- as.factor(test_pop_ref_2$other$ind.metrics$variety)
  t1 <- dartR::gl.dist.pop(test_pop_ref_2, method = "nei",
                           plot.out = FALSE,
                           verbose = 0)
  t1 <- as.matrix(t1)

    colors_pops <-
      polychrome(length(levels(
        test_pop_ref_2$other$ind.metrics$variety
      )))
    names(colors_pops) <-
      as.character(levels(test_pop_ref_2$other$ind.metrics$variety))

    df_colors_temp_1 <-
      as.data.frame(cbind(
        as.character(pop(test_pop_ref_2)),
        as.character(test_pop_ref_2$other$ind.metrics$variety)
      ))
    df_colors_temp_1 <- unique(df_colors_temp_1)
    df_colors_temp_1$order <- 1:nPop(test_pop_ref_2)
    colnames(df_colors_temp_1) <- c("ind", "pop", "order")

    df_colors_temp_2 <- as.data.frame(cbind(names(colors_pops), colors_pops))
    colnames(df_colors_temp_2) <- c("pop", "color")
    df_colors <- merge(df_colors_temp_1, df_colors_temp_2, by = "pop")
    df_colors$order <- as.numeric(df_colors$order)
    df_colors <- df_colors[order(df_colors$order),]

    df_colors_2 <- merge(data.frame(ind=colnames(t1)),
                         df_colors,by="ind" )

    palette.divergent <- dartR.base::gl.colors("div")

    pdf(
      paste0(strsplit(basename(rds.file), "_")[[1]][1], "_ref_distance.pdf"),
      width = nPop(test_pop_ref_2) / 5,
      height = nPop(test_pop_ref_2) / 5
    )
    heatmap.3(
      t1,
      margins = c(10, 10),
      ColSideColors = df_colors_2$color,
      RowSideColors = df_colors_2$color,
      sepcolor = "black",
      dendrogram = "column",
      trace = "none",
      col = viridis::turbo(255),
      colRow = df_colors_2$color,
      colCol = df_colors_2$color,
      density.info = "none",
      reorderfun = function(d, w){reorder(d, w, agglo.FUN = mean, na.rm = TRUE)},
      main = "Genetic distance (Nei's distance) between references",
      na.rm = TRUE,
      na.color = "grey"
    )
    # Close device
    dev.off()

    return(t1)

}
