#' verify
#' @param report file name
#' @noRd
ds14.verify <- function(report) {

  if (report[1, 1] != "*") {
    stop("Report does not contain expected '*' in first row/column")
  }

  startGenotypicColumn <- which(report[1, ] != "*")[1]
  startGenotypicRow <- which(report[, 1] != "*")[1] + 1
  if (startGenotypicRow > 1) {
    endMetaRow <- startGenotypicRow - 1
  } else {
    stop(
      "Could not identify any sample meta rows. Use '*' at beginning of row to indicate a sample meta row"
    )
  }

  if (startGenotypicColumn > 1) {
    endMetaColumn <- startGenotypicColumn - 1
  } else {
    stop(
      "Could not identify any marker meta columns. Use '*' in the first row to indicate a marker meta column"
    )
  }

  return(
    list(
      startGenotypicRow = startGenotypicRow, startGenotypicColumn = startGenotypicColumn
    )
  )
}

#' read
#' @param file file name
#' @noRd
ds14.read <- function(file) {
  if (is.character(file)) {
    if (file.exists(file)) {
      report <- as.matrix(fread(file, header = FALSE, stringsAsFactors = FALSE))
    } else {
      stop("report is of type character, but does not appear to point to a valid file")
    }
  } else if (!is.matrix(file) &&
             !is.data.frame(file)) {
    stop("report must be either a matrix, data.frame or a file path to a report")
  } else {
    report <- file
  }

  start <- ds14.verify(report)

  startGenotypicColumn <- start$startGenotypicColumn
  startGenotypicRow <- start$startGenotypicRow


  rownames(report) <- c(
    sapply(
      1:(startGenotypicRow - 1), function(x) paste("SampleMeta", x, collapse = "", sep = "")
    ),
    report[-c(1:(startGenotypicRow - 1)),
           1]
  )
  colnames(report) <- as.vector(report[startGenotypicRow - 1, ])

  message(
    "Report contains ", nrow(report) -
      startGenotypicRow + 1, " markers, ", ncol(report) -
      startGenotypicColumn + 1, " samples. ", startGenotypicRow - 1, " sample meta fields.\nMarker fields:",
    paste(
      c("", report[startGenotypicRow - 1, 1:(startGenotypicColumn - 1)]),
      collapse = "\n\t- "
    )
  )
  return(report)
}

#' subMeta
#' @param report description
#' @param field description
#' @param data description
#' @param sampleField description
#' @noRd
ds14.subMeta <- function(report,
                         field,
                         data,
                         sampleField = FALSE) {
  if (sampleField) {
    meta <- ds14.sampleMeta(report)
    if (any(field %in% colnames(meta))) {
      rowIdx <- which(field == colnames(meta))[1]
      report[rowIdx, -c(
        1:(ncol(report) -
             nrow(meta))
      )] <- data
    } else {
      stop(
        "Meta field ", field, " does not exist in report: ", paste(
          colnames(meta),
          collapse = ", "
        )
      )
    }
  } else {
    meta <- ds14.markerMeta(report)
    if (any(field %in% colnames(meta))) {
      colIdx <- which(field == colnames(meta))[1]
      report[-c(
        1:(nrow(report) -
             nrow(meta))
      ),
      colIdx] <- data
    } else {
      stop(
        "Meta field ", field, " does not exist in report: ", paste(
          colnames(meta),
          collapse = ", "
        )
      )
    }
  }
  return(report)
}

#' parseDf
#' @param df description
#' @noRd
parseDf <- function(df) {
  for (i in 1:ncol(df)) {
    df[, i] <- as.character(df[, i])
    if (suppressWarnings(all(!is.na(as.numeric(df[, i]))))) {
      df[, i] <- as.numeric(df[, i])
    }
  }
  return(df)
}

#' genotypic
#' @param report description
#' @param sampleHeader description
#' @noRd
ds14.genotypic <- function(report,
                           sampleHeader = FALSE) {
  start <- ds14.verify(report)
  if (sampleHeader) {
    start$startGenotypicColumn <- start$startGenotypicColumn +
      1
  }
  return(
    extractGenotypicData(
      report, sampleIDRow = start$startGenotypicRow -
        1, markerIDColumn = 1,
      startGenotypicRow = start$startGenotypicRow,
      startGenotypicColumn = start$startGenotypicColumn
    )
  )
}

#' marker Meta
#' @param report description
#' @param sampleHeader description
#' @noRd
ds14.markerMeta <- function(report,
                            sampleHeader = FALSE) {
  start <- ds14.verify(report)
  if (sampleHeader) {
    start$startGenotypicColumn <- start$startGenotypicColumn + 1
  }

  df <- as.data.frame(
    report[-c(1:(start$startGenotypicRow - 1)),
           1:(start$startGenotypicColumn - 1)], stringsAsFactors = FALSE
  )
  rownames(df) <- 1:nrow(df)
  return(parseDf(df))
}

#' sample Meta
#' @param report description
#' @param sampleHeader description
#' @noRd
ds14.sampleMeta <- function(report,
                            sampleHeader = FALSE) {
  start <- ds14.verify(report)
  if (sampleHeader) {
    start$startGenotypicColumn <- start$startGenotypicColumn + 1
  }
  df <- parseDf(
    data.frame(
      t(
        report[c(1:(start$startGenotypicRow - 1)),
               -c(1:(start$startGenotypicColumn - 1))]
      ),
      stringsAsFactors = F
    )
  )
  if (sampleHeader) {
    colnames(df) <- c(
      report[c(1:(start$startGenotypicRow - 2)),
             start$startGenotypicColumn - 1], "ID"
    )
  }
  return(df)
}

#' addMeta
#' @param report description
#' @param field description
#' @param data description
#' @param afterField description
#' @param sampleField description
#' @noRd
ds14.addMeta <- function(report,
                         field,
                         data,
                         afterField,
                         sampleField = FALSE) {
  report <- as.matrix(report)
  if (sampleField) {
    meta <- ds14.sampleMeta(report)
    if (any(afterField %in% colnames(meta))) {
      rowIdx <- which(afterField == colnames(meta))[1]
      nf <- matrix(
        c(
          rep(
            "*", ncol(report) -
              nrow(meta)
          ),
          data
        ),
        nrow = 1
      )
      report <- rbind(
        report[1:rowIdx, ], nf, report[-c(1:rowIdx),
        ]
      )
    } else {
      stop(
        "Meta field ", afterField, " does not exist in report: ", paste(
          colnames(meta),
          collapse = ", "
        )
      )
    }
  } else {
    meta <- ds14.markerMeta(report)
    if (any(afterField %in% colnames(meta))) {
      colIdx <- which(afterField == colnames(meta))[1]
      nf <- matrix(
        c(
          rep(
            "*", nrow(report) -
              nrow(meta) -
              1
          ),
          field, data
        ),
        ncol = 1
      )
      colnames(nf) <- field
      report <- cbind(report[, 1:colIdx], nf, report[, -c(1:colIdx)])
    } else {
      stop(
        "Meta field ", afterField, " does not exist in report: ", paste(
          colnames(meta),
          collapse = ", "
        )
      )
    }
  }
  return(report)
}

#' removeMeta
#' @param report description
#' @param field description
#' @param sampleField description
#' @noRd
ds14.removeMeta <- function(report,
                            field,
                            sampleField = FALSE) {
  if (sampleField) {
    meta <- ds14.sampleMeta(report)
    if (any(field %in% colnames(meta))) {
      rowIdx <- which(field == colnames(meta))[1]
      report <- rbind(
        report[1:(rowIdx - 1), ], report[-c(1:rowIdx),
        ]
      )
    } else {
      stop(
        "Meta field ", field, " does not exist in report: ", paste(
          colnames(meta),
          collapse = ", "
        )
      )
    }
  } else {
    meta <- ds14.markerMeta(report)
    if (any(field %in% colnames(meta))) {
      colIdx <- which(field == colnames(meta))[1]
      report <- cbind(report[, 1:(colIdx - 1)], report[, -c(1:colIdx)])
    } else {
      stop(
        "Meta field ", field, " does not exist in report: ", paste(
          colnames(meta),
          collapse = ", "
        )
      )
    }
  }
  return(report)
}

#' as String Matrix
#' @param m description
#' @param na description
#' @noRd
asStringMatrix <- function(m,
                           na = "") {
  m <- apply(m, 2, as.character)
  m[is.na(m)] <- na
  return(m)
}

#' write
#' @param report description
#' @param file description
#' @noRd
ds14.write <- function(report,
                       file) {
  verify <- ds14.verify(report)
  markerMeta <- asStringMatrix(
    report[, 1:(verify$startGenotypicColumn - 1)],
    na = ""
  )
  genotypes <- asStringMatrix(
    report[, -c(1:(verify$startGenotypicColumn - 1))],
    na = ""
  )
  reportFormatted <- cbind(markerMeta, genotypes)
  write.table(
    reportFormatted, file, quote = FALSE, sep = ",", row.names = FALSE,
    col.names = FALSE
  )
}

#' sort
#' @param report description
#' @param names description
#' @param sampleField description
#' @noRd
ds14.sort <- function(report,
                      names,
                      sampleField = FALSE) {
  start <- ds14.verify(report)
  if (sampleField) {
    if (all(
      names %in% (colnames(report)[-c(
        1:(start$startGenotypicColumn -
           1)
      )])
    )) {
      nonHeaderIdx <- -c(
        1:(start$startGenotypicColumn -
             1)
      )
      report <- cbind(
        report[, 1:(start$startGenotypicColumn -
                      1)], report[, nonHeaderIdx][,
                                                  which(
                                                    colnames(report)[nonHeaderIdx] %in%
                                                      names
                                                  )]
      )
    } else {
      stop(
        paste(
          sum(
            !names %in% (colnames(report)[-c(
              1:(start$startGenotypicColumn -
                   1)
            )])
          ),
          " samples not found in report"
        )
      )
    }
  } else {
    if (all(
      names %in% (rownames(report)[-c(
        1:(start$startGenotypicRow -
           1)
      )])
    )) {
      nonHeaderIdx <- -c(
        1:(start$startGenotypicRow -
             1)
      )
      report <- rbind(
        report[1:(start$startGenotypicRow -
                    1), ], report[nonHeaderIdx,
                    ][which(
                      rownames(report)[nonHeaderIdx] %in%
                        names
                    ),
                    ]
      )
    } else {
      stop(
        paste(
          sum(
            !names %in% (rownames(report)[-c(
              1:(start$startGenotypicRow -
                   1)
            )])
          ),
          " markers not found in report"
        )
      )
    }
  }
  return(report)
}

#' extract Genotypic Data
#' @param report description
#' @param sampleIDRow description
#' @param markerIDColumn description
#' @param startGenotypicRow description
#' @param startGenotypicColumn description
#' @param addTag description
#' @noRd
extractGenotypicData <- function(report,
                                 sampleIDRow,
                                 markerIDColumn,
                                 startGenotypicRow,
                                 startGenotypicColumn,
                                 addTag = FALSE
) {
  report_dim <- dim(report)

  genotypic_data <- suppressWarnings(
    apply(
      report[-c(1:(startGenotypicRow - 1)),
             -c(1:(startGenotypicColumn - 1))],
      2, as.numeric
    )
  )

  cnames <- as.vector(
    sapply(
      report[sampleIDRow, -c(1:(startGenotypicColumn - 1))],
      as.character
    )
  )
  if (addTag) {
    cnames <- addTag(cnames)

  }
  colnames(genotypic_data) <- cnames


  rnames <- sapply(
    report[-c(1:(startGenotypicRow - 1)),
           markerIDColumn], as.character
  )
  if (addTag) {
    rnames <- addTag(rnames)

  }
  rownames(genotypic_data) <- as.vector(rnames)


  return(genotypic_data)
}
