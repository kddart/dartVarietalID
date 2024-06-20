
# Function to estimate the semi-axes of an ellipsoid from a mesh3d object
estimateSemiAxes <- function(mesh,centroid) {
  # Calculate the distances of vertices from the centroid along each axis
  distances <- t(mesh$vb[1:3,]-centroid)
  # The semi-axes are approximated by the maximum absolute distance along each axis
  semiAxes <- apply(abs(distances), 2, max)
  return(semiAxes)
}

#get a box coordinates to put random points
getBoundingBox <- function(center, axes) {
  # Unpack the center coordinates and axes lengths
  x0 <- center[1]
  y0 <- center[2]
  z0 <- center[3]
  a <- axes[1]
  b <- axes[2]
  c <- axes[3]

  # Calculate the bounding box coordinates
  minX <- x0 - a
  maxX <- x0 + a
  minY <- y0 - b
  maxY <- y0 + b
  minZ <- z0 - c
  maxZ <- z0 + c

  return(list(min = c(minX, minY, minZ), max = c(maxX, maxY, maxZ)))
}


na_check <- function(x){
  temp <- sapply(x@gen, function(e){
    length(e@NA.posi)
  })
  na_tmp <- round((sum(temp) / (
    nInd(x) * nLoc(x)
  )) * 100, 2)
  return(na_tmp)
}

readTargetInfoFile <- function(file
                               # ,
                               # filterFields = c("Genotype",
                               #                  "ExPlateBarcode",
                               #                  "ExPlateWell")
                               ) {
  obj <- {
  }

  if (is.null(dim(file))) {
    obj$table <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  } else {
    obj$table <- file
  }

  # fix for no visible binding for global variable from CRAN checks
  TargetID <- SampleType <- RefType <- NULL

  refField <- "variety"
  # barcodeField <- "ExPlateBarcode"
  # wellField <- "ExPlateWell"
  # expectedFields <- c("TargetID", "SampleType", refField)
  expectedFields <- c("TargetID", "reference", refField)


  if (length(intersect(expectedFields, names(obj$table))) !=
      length(expectedFields)) {
    stop(
      paste(
        "Missing fields in Target Information File:",
        file,
        "\nExpected:",
        expectedFields,
        "\nFound:",
        names(obj$table)
      )
    )

  }
  rownames(obj$table) <- obj$table$TargetID
  # obj$table <- subset(obj$table, select = -c(TargetID))
  obj$table <-
    # obj$table[, names(obj$table) %in% union(expectedFields,
    #                                         filterFields)]
  obj$table[, expectedFields]

  obj$getReferences <- function() {
    refTable <-
      obj$table[!is.na(obj$table$variety) & obj$table$variety !=
                  "", ]
    # table <- subset(refTable, select = -c(SampleType))
    table <- subset(refTable, select = -c(reference))
    # table$ExPlateBarcode <-
    #   sapply(table$ExPlateBarcode, as.character)
    # return(sortByRefTypeBarcodeWell_2(table))
    return(table)

  }

  obj$getSamples <- function() {
    refTable <-
      obj$table[is.na(obj$table$variety) | obj$table$variety ==
                  "", ]
    table <- subset(refTable, select = -c(variety))
    # table$ExPlateBarcode <-
    #   sapply(table$ExPlateBarcode, as.character)
    # return(sortByRefTypeBarcodeWell_2(table, refField = "SampleType"))
    return(table)

  }
  return(obj)
}

sortByRefTypeBarcodeWell_2 <- function(table,
                                       refField = "variety",
                                       barcodeField = "ExPlateBarcode",
                                       wellField = "ExPlateWell",
                                       priorOrderFields = c()) {
  if (!is.null(barcodeField) &&
      !any(names(table) == barcodeField)) {
    stop(paste("Barcode field '", refField, "' is missing in table",
               sep = ""))
  }

  if (!is.null(wellField) && !any(names(table) == wellField)) {
    stop(paste("Well field '", refField, "' is missing in table",
               sep = ""))
  }

  table[, refField][table[, refField] == ""] <- NA

  if (!is.null(refField) &&
      !is.null(barcodeField) && !is.null(wellField)) {
    wellRows <- sapply(as.character(table[, wellField]), function(x) {
      return(strsplit(gsub("([0-9]+)", ",\\1", x), ",")[[1]][1])
    })
    wellColumns <-
      as.integer(sapply(as.character(table[, wellField]),
                        function(x) {
                          return(strsplit(gsub("([0-9]+)", ",\\1", x), ",")[[1]][2])
                        }))

    tableWithWell <- cbind(table, wellRows, wellColumns)
    if (length(priorOrderFields) == 0) {
      tr <- tableWithWell[order(table[, refField], table[,
                                                         barcodeField], wellRows, wellColumns, na.last = T), ]
    } else {
      for (f in priorOrderFields) {
        if (!any(f == colnames(table))) {
          stop("Field ", f, " does not exist in table")
        }
      }
      tr <- tableWithWell[order(table[, priorOrderFields],
                                table[, refField], table[, barcodeField], wellRows,
                                wellColumns, na.last = T), ]
    }
    return(tr[, !colnames(tr) %in% c("wellRows", "wellColumns")])
  }
}

gl2alleles <- function (gl) {
  x <- as.matrix(gl)
  homs1 <-
    paste(substr(gl@loc.all, 1, 1), "/", substr(gl@loc.all, 1, 1), sep = "")
  hets <- gl@loc.all
  homs2 <-
    paste(substr(gl@loc.all, 3, 3), "/", substr(gl@loc.all, 3, 3), sep = "")
  xx <- matrix(NA, ncol = ncol(x), nrow = nrow(x))
  for (i in 1:nrow(x)) {
    for (ii in 1:ncol(x)) {
      inp <- x[i, ii]
      if (!is.na(inp)) {
        if (inp == 0)
          xx[i, ii] <- homs1[ii]
        else if (inp == 1)
          xx[i, ii] <- hets[ii]
        else if (inp == 2)
          xx[i, ii] <- homs2[ii]
      } else{
        xx[i, ii] <- NA
      }
    }
  }
  xx <- gsub("/", ":", xx)
  return(xx)
}

# selecting the representative individual from the sample using PCA
rep_ind <- function(x){
  # removing missing data for PCA
  pop_test <- dartR::gl.filter.callrate(x,
                                        threshold = 1,
                                        verbose = 0)
  # test whether all individuals are the same (when all loci are monomorphic)
  pop_test_mat <- as.matrix(pop_test)
  test_var <- sum(apply(pop_test_mat, 2, function(x) {
    var(x) != 0
  }), na.rm = TRUE)
  # if individuals are different
  if (test_var > 0) {
    # if unix
    if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
      pcoa <- adegenet::glPca(pop_test,
                              nf = 3,
                              parallel = FALSE,
                              loadings = FALSE)
    }

    ## if windows
    if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
      pcoa <- adegenet::glPca(pop_test,
                              nf = 3,
                              parallel = TRUE,
                              loadings = FALSE)
    }

    pcoa_scores <- pcoa$scores
    means <- colMeans(pcoa_scores)
    covariance <- stats::cov(pcoa_scores)
    D <- stats::mahalanobis(x = pcoa_scores,
                            center = means,
                            cov = covariance,
                            toll = 1e-20)
    return(x[which.min(D),])
    # if individuals are the same, get the first individual
  } else{
    return(x[1,])
  }
}

heatmap.3 <-
  function (x,
            Rowv = TRUE,
            Colv = if (symm)
              "Rowv"
            else
              TRUE,
            distfun = dist,
            hclustfun = hclust,
            dendrogram = c("both",
                           "row", "column", "none"),
            reorderfun = function(d, w)
              reorder(d,
                      w),
            symm = FALSE,
            scale = c("none", "row", "column"),
            na.rm = TRUE,
            revC = identical(Colv, "Rowv"),
            add.expr,
            breaks,
            symbreaks = any(x < 0, na.rm = TRUE) ||
              scale != "none",
            col = "heat.colors",
            colsep,
            rowsep,
            sepcolor = "white",
            sepwidth = c(0.05, 0.05),
            cellnote,
            notecex = 1,
            notecol = "cyan",
            na.color = par("bg"),
            trace = c("column", "row", "both",
                      "none"),
            tracecol = "cyan",
            hline = median(breaks),
            vline = median(breaks),
            linecol = tracecol,
            margins = c(5, 5),
            ColSideColors,
            RowSideColors,
            cexRow = 0.2 + 1 / log10(nr),
            cexCol = 0.2 + 1 / log10(nc),
            labRow = NULL,
            labCol = NULL,
            srtRow = NULL,
            srtCol = NULL,
            adjRow = c(0,
                       NA),
            adjCol = c(NA, 0),
            offsetRow = 0.5,
            offsetCol = 0.5,
            colRow = NULL,
            colCol = NULL,
            key = TRUE,
            keysize = 1.5,
            density.info = c("histogram", "density", "none"),
            denscol = tracecol,
            symkey = any(x < 0, na.rm = TRUE) ||
              symbreaks,
            densadj = 0.25,
            key.title = NULL,
            key.xlab = NULL,
            key.ylab = NULL,
            key.xtickfun = NULL,
            key.ytickfun = NULL,
            key.par = list(),
            main = NULL,
            xlab = NULL,
            ylab = NULL,
            lmat = NULL,
            lhei = NULL,
            lwid = NULL,
            extrafun = NULL,
            ...)
  {
    scale01 <- function(x,
                        low = min(x),
                        high = max(x)) {
      x <- (x - low) / (high - low)
      x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
      "none"
    else
      match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
      col <- get(col, mode = "function")
    if (!missing(breaks) && any(duplicated(breaks)))
      stop("breaks may not contain duplicate values")
    if (!missing(breaks) && (scale != "none"))
      warning(
        "Using scale=\"row\" or scale=\"column\" when breaks are",
        "specified can produce unpredictable results.",
        "Please consider using only one or the other."
      )
    if (is.null(Rowv) || any(is.na(Rowv)))
      Rowv <- FALSE
    if (is.null(Colv) || any(is.na(Colv)))
      Colv <- FALSE
    else if (all(Colv == "Rowv"))
      Colv <- Rowv
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) &&
          (dendrogram %in% c("both", "row"))) {
        warning(
          "Discrepancy: Rowv is FALSE, while dendrogram is `",
          dendrogram,
          "'. Omitting row dendogram."
        )
        if (dendrogram == "both")
          dendrogram <- "column"
        else
          dendrogram <- "none"
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) &&
          (dendrogram %in% c("both", "column"))) {
        warning(
          "Discrepancy: Colv is FALSE, while dendrogram is `",
          dendrogram,
          "'. Omitting column dendogram."
        )
        if (dendrogram == "both")
          dendrogram <- "row"
        else
          dendrogram <- "none"
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
      if (length(rowInd) > nr || any(rowInd < 1 | rowInd >
                                     nr))
        stop("Rowv dendrogram doesn't match size of x")
      if (length(rowInd) < nr)
        nr <- length(rowInd)
    }
    else if (is.integer(Rowv)) {
      distr <- distfun(x)
      hcr <- hclustfun(distr)
      ddr <- as.dendrogram(hcr)
      ddr <- reorderfun(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      distr <- distfun(x)
      hcr <- hclustfun(distr)
      ddr <- as.dendrogram(hcr)
      ddr <- reorderfun(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (!isTRUE(Rowv)) {
      rowInd <- nr:1
      ddr <- as.dendrogram(hclust(dist(diag(nr))))
    }
    else {
      rowInd <- nr:1
      ddr <- as.dendrogram(Rowv)
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
      if (length(colInd) > nc || any(colInd < 1 | colInd >
                                     nc))
        stop("Colv dendrogram doesn't match size of x")
      if (length(colInd) < nc)
        nc <- length(colInd)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- order.dendrogram(ddc)
      }
      else
        colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      distc <- distfun(if (symm)
        x
        else
          t(x))
      hcc <- hclustfun(distc)
      ddc <- as.dendrogram(hcc)
      ddc <- reorderfun(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      distc <- distfun(if (symm)
        x
        else
          t(x))
      hcc <- hclustfun(distc)
      ddc <- as.dendrogram(hcc)
      ddc <- reorderfun(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (!isTRUE(Colv)) {
      colInd <- 1:nc
      ddc <- as.dendrogram(hclust(dist(diag(nc))))
    }
    else {
      colInd <- 1:nc
      ddc <- as.dendrogram(Colv)
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
      labRow <- if (is.null(rownames(x)))
        (1:nr)[rowInd]
    else
      rownames(x)
    else
      labRow <- labRow[rowInd]
    if (is.null(labCol))
      labCol <- if (is.null(colnames(x)))
        (1:nc)[colInd]
    else
      colnames(x)
    else
      labCol <- labCol[colInd]
    if (!is.null(colRow))
      colRow <- colRow[rowInd]
    if (!is.null(colCol))
      colCol <- colCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) <
        1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else
        breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (is(col, "function"))
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
      lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
      lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)
      if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) !=
            nc)
          stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1,] + 1, c(NA, 1), lmat[2,] +
                        1)
        lhei <- c(lhei[1], 0.2, lhei[2])
      }
      if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) !=
            nr)
          stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
                                             1), 1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat))
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat,
           widths = lwid,
           heights = lhei,
           respect = FALSE)
    plot.index <- 1
    if (!missing(RowSideColors)) {
      par(mar = c(margins[1], 0, 0, 0))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      plot.index <- plot.index + 1
    }
    if (!missing(ColSideColors)) {
      par(mar = c(0, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      plot.index <- plot.index + 1
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr"))
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else
      iy <- 1:nr
    image(
      1:nc,
      1:nr,
      x,
      xlim = 0.5 + c(0, nc),
      ylim = 0.5 +
        c(0, nr),
      axes = FALSE,
      xlab = "",
      ylab = "",
      col = col,
      breaks = breaks,
      ...
    )
    retval$carpet <- x
    if (exists("ddr"))
      retval$rowDendrogram <- ddr
    if (exists("ddc"))
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!gtools::invalid(na.color) & any(is.na(x))) {
      mmat <- ifelse(is.na(x), 1, NA)
      image(
        1:nc,
        1:nr,
        mmat,
        axes = FALSE,
        xlab = "",
        ylab = "",
        col = na.color,
        add = TRUE
      )
    }
    if (is.null(srtCol) && is.null(colCol))
      axis(
        1,
        1:nc,
        labels = labCol,
        las = 2,
        line = -0.5 +
          offsetCol,
        tick = 0,
        cex.axis = cexCol,
        hadj = adjCol[1],
        padj = adjCol[2]
      )
    else {
      if (is.null(srtCol) || is.numeric(srtCol)) {
        if (missing(adjCol) || is.null(adjCol))
          adjCol = c(1, NA)
        if (is.null(srtCol))
          srtCol <- 90
        xpd.orig <- par("xpd")
        par(xpd = NA)
        xpos <- axis(1,
                     1:nc,
                     labels = rep("", nc),
                     las = 2,
                     tick = 0)
        text(
          x = xpos,
          y = par("usr")[3] - (1 + offsetCol) *
            strheight("M"),
          labels = labCol,
          adj = adjCol,
          cex = cexCol,
          srt = srtCol,
          col = colCol
        )
        par(xpd = xpd.orig)
      }
      else
        warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow) && is.null(colRow)) {
      axis(
        4,
        iy,
        labels = labRow,
        las = 2,
        line = -0.5 + offsetRow,
        tick = 0,
        cex.axis = cexRow,
        hadj = adjRow[1],
        padj = adjRow[2]
      )
    }
    else {
      if (is.null(srtRow) || is.numeric(srtRow)) {
        xpd.orig <- par("xpd")
        par(xpd = NA)
        ypos <- axis(
          4,
          iy,
          labels = rep("", nr),
          las = 2,
          line = -0.5,
          tick = 0
        )
        text(
          x = par("usr")[2] + (1 + offsetRow) * strwidth("M"),
          y = ypos,
          labels = labRow,
          adj = adjRow,
          cex = cexRow,
          srt = srtRow,
          col = colRow
        )
        par(xpd = xpd.orig)
      }
      else
        warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab))
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab))
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
      eval(substitute(add.expr))
    if (!missing(colsep))
      for (csep in colsep)
        rect(
          xleft = csep + 0.5,
          ybottom = 0,
          xright = csep + 0.5 + sepwidth[1],
          ytop = ncol(x) +
            1,
          lty = 1,
          lwd = 1,
          col = sepcolor,
          border = sepcolor
        )
    if (!missing(rowsep))
      for (rsep in rowsep)
        rect(
          xleft = 0,
          ybottom = (ncol(x) +
                       1 - rsep) - 0.5,
          xright = nrow(x) + 1,
          ytop = (ncol(x) +
                    1 - rsep) - 0.5 - sepwidth[2],
          lty = 1,
          lwd = 1,
          col = sepcolor,
          border = sepcolor
        )
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in 1:length(colInd)) {
        if (!is.null(vline)) {
          abline(v = i - 0.5 + vline.vals,
                 col = linecol,
                 lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        lines(
          x = xv,
          y = yv,
          lwd = 1,
          col = tracecol,
          type = "s"
        )
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in 1:length(rowInd)) {
        if (!is.null(hline)) {
          abline(h = i - 0.5 + hline.vals,
                 col = linecol,
                 lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i,] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        lines(
          x = xv,
          y = yv,
          lwd = 1,
          col = tracecol,
          type = "s"
        )
      }
    }
    if (!missing(cellnote))
      text(
        x = c(row(cellnote)),
        y = c(col(cellnote)),
        labels = c(cellnote),
        col = notecol,
        cex = notecex
      )
    plot.index <- plot.index + 1
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
      flag <- try(plot(ddr,
                       horiz = TRUE,
                       axes = TRUE,
                       yaxs = "i"))
      if ("try-error" %in% class(flag)) {
        cond <- attr(flag, "condition")
        if (!is.null(cond) &&
            conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
          stop(
            "Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...)."
          )
      }
    }
    else
      plot.new()
    par(mar = c(10, 0, if (!is.null(main))
      5
      else
        0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      flag <- try(# plot(ddc, xaxs = "i")
        ddc %>%
          # Custom branches
          # set("branches_col", "grey") %>% set("branches_lwd", 3) %>%
          # Custom labels
          set("labels_col", colCol) %>%
          plot(xaxs = "i"))
      if ("try-error" %in% class(flag)) {
        cond <- attr(flag, "condition")
        if (!is.null(cond) &&
            conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
          stop(
            "Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...)."
          )
      }
    }
    else
      plot.new()
    if (!is.null(main))
      title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
      mar <- c(5, 4, 2, 1)
      if (!is.null(key.xlab) && is.na(key.xlab))
        mar[1] <- 2
      if (!is.null(key.ylab) && is.na(key.ylab))
        mar[2] <- 2
      if (!is.null(key.title) && is.na(key.title))
        mar[3] <- 1
      par(mar = mar,
          cex = 0.75,
          mgp = c(2, 1, 0))
      if (length(key.par) > 0)
        do.call(par, key.par)
      tmpbreaks <- breaks
      if (symkey) {
        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        min.raw <- -max.raw
        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
        min.raw <- min.breaks
        max.raw <- max.breaks
      }
      z <- seq(min.raw, max.raw, by = min(diff(breaks) / 100))
      image(
        z = matrix(z, ncol = 1),
        col = col,
        breaks = tmpbreaks,
        xaxt = "n",
        yaxt = "n"
      )
      par(usr = c(0, 1, 0, 1))
      if (is.null(key.xtickfun)) {
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        xargs <- list(at = xv, labels = lv)
      }
      else {
        xargs <- key.xtickfun()
      }
      xargs$side <- 1
      do.call(axis, xargs)
      if (is.null(key.xlab)) {
        if (scale == "row")
          key.xlab <- "Row Z-Score"
        else if (scale == "column")
          key.xlab <- "Column Z-Score"
        else
          key.xlab <- "Value"
      }
      if (!is.na(key.xlab)) {
        mtext(
          side = 1,
          key.xlab,
          line = par("mgp")[1],
          padj = 0.5,
          cex = par("cex") * par("cex.lab")
        )
      }
      if (density.info == "density") {
        dens <- density(
          x,
          adjust = densadj,
          na.rm = TRUE,
          from = min.scale,
          to = max.scale
        )
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[!omit]
        dens$y <- dens$y[!omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        lines(dens$x,
              dens$y / max(dens$y) * 0.95,
              col = denscol,
              lwd = 1)
        if (is.null(key.ytickfun)) {
          yargs <- list(at = pretty(dens$y) / max(dens$y) *
                          0.95,
                        labels = pretty(dens$y))
        }
        else {
          yargs <- key.ytickfun()
        }
        yargs$side <- 2
        do.call(axis, yargs)
        if (is.null(key.title))
          key.title <- "Color Key\nand Density Plot"
        if (!is.na(key.title))
          title(key.title)
        par(cex = 0.5)
        if (is.null(key.ylab))
          key.ylab <- "Density"
        if (!is.na(key.ylab))
          mtext(
            side = 2,
            key.ylab,
            line = par("mgp")[1],
            padj = 0.5,
            cex = par("cex") * par("cex.lab")
          )
      }
      else if (density.info == "histogram") {
        h <- hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        lines(
          hx,
          hy / max(hy) * 0.95,
          lwd = 1,
          type = "s",
          col = denscol
        )
        if (is.null(key.ytickfun)) {
          yargs <- list(at = pretty(hy) / max(hy) * 0.95,
                        labels = pretty(hy))
        }
        else {
          yargs <- key.ytickfun()
        }
        yargs$side <- 2
        do.call(axis, yargs)
        if (is.null(key.title))
          key.title <- "Color Key\nand Histogram"
        if (!is.na(key.title))
          title(key.title)
        par(cex = 0.5)
        if (is.null(key.ylab))
          key.ylab <- "Count"
        if (!is.na(key.ylab))
          mtext(
            side = 2,
            key.ylab,
            line = par("mgp")[1],
            padj = 0.5,
            cex = par("cex") * par("cex.lab")
          )
      }
      else {
        if (is.null(key.title))
          key.title <- "Color Key"
        if (!is.na(key.title))
          title(key.title)
      }
      if (trace %in% c("both", "column")) {
        vline.vals <- scale01(vline, min.raw, max.raw)
        if (!is.null(vline)) {
          abline(v = vline.vals,
                 col = linecol,
                 lty = 2)
        }
      }
      if (trace %in% c("both", "row")) {
        hline.vals <- scale01(hline, min.raw, max.raw)
        if (!is.null(hline)) {
          abline(v = hline.vals,
                 col = linecol,
                 lty = 2)
        }
      }
    }
    else {
      par(mar = c(0, 0, 0, 0))
      plot.new()
    }
    retval$colorTable <-
      data.frame(
        low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1],
        color = retval$col
      )
    retval$layout <- list(lmat = lmat,
                          lhei = lhei,
                          lwid = lwid)
    if (!is.null(extrafun))
      extrafun()
    invisible(retval)
  }


DAP_write_excel <- function(x, info.n, filename) {
  info <- read.csv(info.n)
  nms <- names(x$res_full)
  resfull <- do.call(rbind, lapply(
    1:length(x$res_full),
    \(i) data.frame(field_TID = nms[i], x$res_full[[i]])
  ))
  colnames(resfull)[1:4] <-
    c("field_Tid", "ref_Tid", "ref_id", "variety")
  info <- info[, c("TargetID", "sample")]
  colnames(info)[2] <- "field_id"
  full <- merge(resfull, info, by = 1, all.x = TRUE)
  full$var_rank <-
    with(full, stats::ave(Probability_corr_scaled, field_id, FUN = \(x) rank(1000 - x, ties.method =
                                                                               "min")))
  full$Probability <- full$Probability / 100
  x[[2]] <- full
  colnames(x$res_summary)[1:7] <-
    c("field_Tid",
      "field_id",
      "ref_Tid",
      "ref_id",
      "variety",
      "NA.perc",
      "Probability")
  x$res_summary$Probability <- x$res_summary$Probability / 100
  # x[[3]] <- as.data.frame(x$ref_distance, check.names = FALSE)
  writexl::write_xlsx(x[1:2], paste0(filename, ".xlsx"))
}
