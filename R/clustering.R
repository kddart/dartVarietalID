getShuffledPalette <- function(n, colorFunc = polychrome) {
  palette = colorFunc(n)
  return(sample(palette[1:n]))
}

getMergedBinsClusters <-
  function(distances,
           method,
           mapping = NULL,
           coloring = NULL,
           showLine = TRUE,
           useAutoPar = TRUE,
           text = NULL,
           colorFunc = polychrome) {
    set.seed(123)
    hc = hclust(as.dist(distances), method = method)

    library(dendextend)
    distancesByAxis = as.matrix(distances)

    if (!is.null(mapping)) {
      if (!all(rownames(distances) %in% names(mapping))) {
        stop("Missing mappings for: ",
             paste(
               rownames(distances)[rownames(distances) %in% names(mapping)],
               collapse = ", ",
               sep = ", "
             ))
      }
      rownames(distancesByAxis) = mapping[rownames(distancesByAxis)]
      colnames(distancesByAxis) = mapping[colnames(distancesByAxis)]
      if (!is.null(coloring)) {
        if (!all(rownames(distances) %in% names(coloring))) {
          stop("Missing mappings for: ",
               paste(
                 rownames(distances)[rownames(distances) %in% names(mapping)],
                 collapse = ", ",
                 sep = ", "
               ))
        }
        names(coloring) = mapping[names(coloring)]
      }
    }

    hc1 = hclust(as.dist(distancesByAxis), method = method)
    dend <- as.dendrogram(hc1)

    if (!is.null(coloring)) {
      map  = hc1$labels[hc1$order]
      unique_colors = length(unique(coloring))
      colorMap = getShuffledPalette(n = unique_colors, colorFunc = colorFunc)
      names(colorMap) = unique(coloring)
      labels_colors(dend) = colorMap[coloring[map]]
    } else{
      map  = as.integer(factor(hc1$labels[hc1$order]))
      map2Col = getShuffledPalette(n = length(unique(map)))
      names(map2Col) = unique(map)
      labels_colors(dend) = map2Col[as.character(map)]
    }

    return(list(hc = hc1, dend = dend))
  }

polychrome = function (n)
{
  pal <- c(
    "#3283FE",
    "#FEAF16",
    "#B00068",
    "#1CFFCE",
    "#90AD1C",
    "#2ED9FF",
    "#DEA0FD",
    "#AA0DFE",
    "#F8A19F",
    "#325A9B",
    "#C4451C",
    "#1C8356",
    "#85660D",
    "#B10DA1",
    "#FBE426",
    "#1CBE4F",
    "#FA0087",
    "#FC1CBF",
    "#F7E1A0",
    "#C075A6",
    "#782AB6",
    "#AAF400",
    "#BDCDFF",
    "#822E1C",
    "#B5EFB5",
    "#7ED7D1",
    "#1C7F93",
    "#D85FF7",
    "#683B79",
    "#66B0FF",
    "#5A5156",
    "#E4E1E3",
    "#F6222E",
    "#FE00FA",
    "#16FF32",
    "#3B00FB"
  )
  names(pal) <-
    c(
      "Dark_Purplish_Gray",
      "Purplish_White",
      "Vivid_Red",
      "Vivid_Purple",
      "Vivid_Yellowish_Green",
      "Strong_Purplish_Blue",
      "Vivid_Orange_Yellow",
      "Vivid_Purplish_Red",
      "Brilliant_Green",
      "Vivid_Yellow_Green",
      "Vivid_Blue",
      "Brilliant_Purple",
      "Vivid_Violet",
      "Strong_Pink",
      "Strong_Blue",
      "Strong_Reddish_Orange",
      "Vivid_Green",
      "Light_Olive_Brown",
      "Vivid_Reddish_Purple",
      "Vivid_Greenish_Yellow",
      "Vivid_Yellowish_Green",
      "Vivid_Red",
      "Vivid_Purplish_Red",
      "Pale_Yellow",
      "Strong_Reddish_Purple",
      "Vivid_Violet",
      "Vivid_Yellow_Green",
      "Very_Light_Blue",
      "Strong_Reddish_Brown",
      "Very_Light_Yellowish_Green",
      "Very_Light_Bluish_Green",
      "Deep_Greenish_Blue",
      "Vivid_Purple",
      "Deep_Purple",
      "Brilliant_Blue",
      "Vivid_Violet"
    )
  pal2 = colorspace::darken(pal, amount = 0.2)
  # pie(rep(1,length(pal2)),col = pal2)
  return(rep(pal2, ceiling(n / length(pal2)))[1:n])
}
