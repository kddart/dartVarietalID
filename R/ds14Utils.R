ds14.verify <- function(report) {
  if (report[1, 1] != "*") {
    stop("Report does not contain expected '*' in first row/column")
  }
  startGenotypicColumn = which(report[1, ] != "*")[1]
  startGenotypicRow = which(report[, 1] != "*")[1] + 1
  if (startGenotypicRow > 1) {
    endMetaRow = startGenotypicRow - 1
  } else{
    stop(
      "Could not identify any sample meta rows. Use '*' at beginning of row to indicate a sample meta row"
    )
  }

  if (startGenotypicColumn > 1) {
    endMetaColumn = startGenotypicColumn - 1
  } else{
    stop(
      "Could not identify any marker meta columns. Use '*' in the first row to indicate a marker meta column"
    )
  }

  return(
    list(
      startGenotypicRow = startGenotypicRow,
      startGenotypicColumn = startGenotypicColumn
    )
  )
}

ds14.read <- function(file) {
  if (class(file) == "character") {
    if (file.exists(file)) {
      report = as.matrix(fread(file, header = FALSE, stringsAsFactors = FALSE))
    } else{
      stop("report is of type character, but does not appear to point to a valid file")
    }
  } else if (class(file) != "matrix" && class(file) != "data.frame") {
    stop("report must be either a matrix, data.frame or a file path to a report")
  } else{
    report = file
  }

  start = ds14.verify(report)

  startGenotypicColumn = start$startGenotypicColumn
  startGenotypicRow = start$startGenotypicRow


  rownames(report) = c(sapply(1:(startGenotypicRow - 1), function(x)
    paste("SampleMeta", x, collapse = "", sep = "")),
    report[-c(1:(startGenotypicRow - 1)), 1])
  colnames(report) = as.vector(report[startGenotypicRow - 1, ])

  message(
    "Report contains ",
    nrow(report) - startGenotypicRow + 1,
    " markers, ",
    ncol(report) - startGenotypicColumn + 1,
    " samples. ",
    startGenotypicRow - 1,
    " sample meta fields.\nMarker fields:",
    paste(c("", report[startGenotypicRow - 1, 1:(startGenotypicColumn -
                                                   1)]), collapse = "\n\t- ")
  )
  return(report)
}

ds14.subMeta <- function(report, field, data, sampleField = F) {
  if (sampleField) {
    meta = ds14.sampleMeta(report)
    if (any(field %in% colnames(meta))) {
      rowIdx = which(field == colnames(meta))[1]
      report[rowIdx, -c(1:(ncol(report) - nrow(meta)))] = data
    } else{
      stop("Meta field ",
           field,
           " does not exist in report: ",
           paste(colnames(meta), collapse = ", "))
    }
  } else{
    meta = ds14.markerMeta(report)
    if (any(field %in% colnames(meta))) {
      colIdx = which(field == colnames(meta))[1]
      report[-c(1:(nrow(report) - nrow(meta))), colIdx] = data
    } else{
      stop("Meta field ",
           field,
           " does not exist in report: ",
           paste(colnames(meta), collapse = ", "))
    }
  }
  return(report)
}


parseDf <- function(df) {
  for (i in 1:ncol(df)) {
    df[, i] = as.character(df[, i])
    if (suppressWarnings(all(!is.na(as.numeric(df[, i]))))) {
      df[, i] = as.numeric(df[, i])
    }
  }
  return(df)
}

ds14.genotypic <- function(report, sampleHeader = F) {
  start = ds14.verify(report)
  if (sampleHeader) {
    start$startGenotypicColumn = start$startGenotypicColumn + 1
  }
  return(
    extractGenotypicData(
      report,
      sampleIDRow = start$startGenotypicRow - 1,
      markerIDColumn = 1,
      startGenotypicRow = start$startGenotypicRow,
      startGenotypicColumn = start$startGenotypicColumn
    )
  )
}

ds14.markerMeta <- function(report, sampleHeader = F) {
  start = ds14.verify(report)
  if (sampleHeader) {
    start$startGenotypicColumn = start$startGenotypicColumn + 1
  }

  df = as.data.frame(report[-c(1:(start$startGenotypicRow - 1)), 1:(start$startGenotypicColumn -
                                                                      1)], stringsAsFactors = F)
  rownames(df) = 1:nrow(df)
  return(parseDf(df))
}

ds14.sampleMeta <- function(report, sampleHeader = F) {
  start = ds14.verify(report)
  if (sampleHeader) {
    start$startGenotypicColumn = start$startGenotypicColumn + 1
  }
  df = parseDf(data.frame(t(report[c(1:(start$startGenotypicRow - 1)), -c(1:(start$startGenotypicColumn -
                                                                               1))]), stringsAsFactors = F))
  if (sampleHeader) {
    colnames(df) = c(report[c(1:(start$startGenotypicRow - 2)), start$startGenotypicColumn -
                              1], "ID")
  }
  return(df)
}

ds14.addMeta <-
  function(report,
           field,
           data,
           afterField,
           sampleField = F) {
    report = as.matrix(report)
    if (sampleField) {
      meta = ds14.sampleMeta(report)
      if (any(afterField %in% colnames(meta))) {
        rowIdx = which(afterField == colnames(meta))[1]
        nf = matrix(c(rep("*", ncol(report) - nrow(meta)), data), nrow = 1)
        report = rbind(report[1:rowIdx, ], nf, report[-c(1:rowIdx), ])
      } else{
        stop(
          "Meta field ",
          afterField,
          " does not exist in report: ",
          paste(colnames(meta), collapse = ", ")
        )
      }
    } else{
      meta = ds14.markerMeta(report)
      if (any(afterField %in% colnames(meta))) {
        colIdx = which(afterField == colnames(meta))[1]
        nf = matrix(c(rep("*", nrow(report) - nrow(meta) - 1), field, data), ncol = 1)
        colnames(nf) = field
        report = cbind(report[, 1:colIdx], nf, report[, -c(1:colIdx)])
      } else{
        stop(
          "Meta field ",
          afterField,
          " does not exist in report: ",
          paste(colnames(meta), collapse = ", ")
        )
      }
    }
    return(report)
  }

ds14.removeMeta <- function(report, field, sampleField = F) {
  if (sampleField) {
    meta = ds14.sampleMeta(report)
    if (any(field %in% colnames(meta))) {
      rowIdx = which(field == colnames(meta))[1]
      report = rbind(report[1:(rowIdx - 1), ], report[-c(1:rowIdx), ])
    } else{
      stop("Meta field ",
           field,
           " does not exist in report: ",
           paste(colnames(meta), collapse = ", "))
    }
  } else{
    meta = ds14.markerMeta(report)
    if (any(field %in% colnames(meta))) {
      colIdx = which(field == colnames(meta))[1]
      report = cbind(report[, 1:(colIdx - 1)], report[, -c(1:colIdx)])
    } else{
      stop("Meta field ",
           field,
           " does not exist in report: ",
           paste(colnames(meta), collapse = ", "))
    }
  }
  return(report)
}

asStringMatrix <- function(m, na = "") {
  m = apply(m, 2, as.character)
  m[is.na(m)] = na
  return(m)
}

ds14.write <- function(report, file) {
  verify = ds14.verify(report)
  markerMeta = asStringMatrix(report[, 1:(verify$startGenotypicColumn -
                                            1)], na = "")
  genotypes = asStringMatrix(report[, -c(1:(verify$startGenotypicColumn -
                                              1))], na = "")
  reportFormatted = cbind(markerMeta, genotypes)
  write.table(
    reportFormatted,
    file,
    quote = F,
    sep = ",",
    row.names = F,
    col.names = F
  )
}

ds14.sort <- function(report, names, sampleField = F) {
  start = ds14.verify(report)
  if (sampleField) {
    if (all(names %in% (colnames(report)[-c(1:(start$startGenotypicColumn -
                                               1))]))) {
      nonHeaderIdx = -c(1:(start$startGenotypicColumn - 1))
      report = cbind(report[, 1:(start$startGenotypicColumn - 1)], report[, nonHeaderIdx][, which(colnames(report)[nonHeaderIdx] %in%
                                                                                                    names)])
    } else{
      stop(paste(sum(!names %in% (
        colnames(report)[-c(1:(start$startGenotypicColumn - 1))]
      )), " samples not found in report"))
    }
  } else{
    if (all(names %in% (rownames(report)[-c(1:(start$startGenotypicRow - 1))]))) {
      nonHeaderIdx = -c(1:(start$startGenotypicRow - 1))
      report = rbind(report[1:(start$startGenotypicRow - 1), ], report[nonHeaderIdx, ][which(rownames(report)[nonHeaderIdx] %in%
                                                                                               names), ])
    } else{
      stop(paste(sum(!names %in% (
        rownames(report)[-c(1:(start$startGenotypicRow - 1))]
      )), " markers not found in report"))
    }
  }
  return(report)
}

ds14.write.genotypic <- function(genotypic_data, file) {
  write.csv(
    genotypic_data,
    file,
    row.names = T,
    quote = F,
    na = "-"
  )
}

ds14.read.genotypic <- function(file) {
  return(read.csv(
    file,
    row.names = 1,
    stringsAsFactors = 1,
    na = "-",
    as.is = T,
    check.names = F
  ))
}

ds14.subGeno <- function(report, genotypes) {
  genotypesReport = ds14.genotypic(report)
  if (!base::all(colnames(genotypesReport) == colnames(genotypes), na.rm = T)) {
    stop("Mismatch of samples")
  }
  if (!base::all(rownames(genotypesReport) == rownames(genotypes), na.rm = T)) {
    stop("Mismatch of markers")
  }
  genotypes = apply(genotypes, 2, as.character)
  genotypes[is.na(genotypes)] = "-"


  priorGenotypicColumn = ncol(report) - ncol(genotypes)
  priorGenotypicRow = nrow(report) - nrow(genotypes)
  report[-c(1:priorGenotypicRow), -c(1:priorGenotypicColumn)] = genotypes

  return(report)
}

extractGenotypicData <-
  function(report,
           sampleIDRow,
           markerIDColumn,
           startGenotypicRow,
           startGenotypicColumn,
           addTag = F) {
    report_dim = dim(report)

    genotypic_data = suppressWarnings(apply(report[-c(1:(startGenotypicRow -
                                                           1)), -c(1:(startGenotypicColumn - 1))], 2, as.numeric))

    cnames = as.vector(sapply(report[sampleIDRow, -c(1:(startGenotypicColumn -
                                                          1))], as.character))
    if (addTag) {
      cnames = addTag(cnames)

    }
    colnames(genotypic_data) = cnames


    rnames = sapply(report[-c(1:(startGenotypicRow - 1)), markerIDColumn], as.character)
    if (addTag) {
      rnames = addTag(rnames)

    }
    rownames(genotypic_data) = as.vector(rnames)


    return(genotypic_data)
  }
