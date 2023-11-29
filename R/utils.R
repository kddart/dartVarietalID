
readTargetInfoFile <- function(file,
                               filterFields = c("Genotype",
                                                "ExPlateBarcode",
                                                "ExPlateWell")) {
  obj <- {}

  if (is.null(dim(file))) {
    obj$table <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  } else {
    obj$table <- file
  }

  # fix for no visible binding for global variable from CRAN checks
  TargetID <- SampleType <- RefType <- NULL

  refField <- "RefType"
  barcodeField <- "ExPlateBarcode"
  wellField <- "ExPlateWell"
  expectedFields <- c("TargetID", "SampleType", refField)

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
  obj$table <- subset(obj$table, select = -c(TargetID))
  obj$table <-
    obj$table[, names(obj$table) %in% union(expectedFields,
                                            filterFields)]

  obj$getReferences <- function() {
    refTable <-
      obj$table[!is.na(obj$table$RefType) & obj$table$RefType !=
                  "",]
    table <- subset(refTable, select = -c(SampleType))
    table$ExPlateBarcode <-
      sapply(table$ExPlateBarcode, as.character)
    return(sortByRefTypeBarcodeWell_2(table))
  }

  obj$getSamples <- function() {
    refTable <-
      obj$table[is.na(obj$table$RefType) | obj$table$RefType ==
                  "",]
    table <- subset(refTable, select = -c(RefType))
    table$ExPlateBarcode <-
      sapply(table$ExPlateBarcode, as.character)
    return(sortByRefTypeBarcodeWell_2(table, refField = "SampleType"))
  }
  return(obj)
}

sortByRefTypeBarcodeWell_2 <- function(table,
                                       refField = "RefType",
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
                                                         barcodeField], wellRows, wellColumns, na.last = T),]
    } else {
      for (f in priorOrderFields) {
        if (!any(f == colnames(table))) {
          stop("Field ", f, " does not exist in table")
        }
      }
      tr <- tableWithWell[order(table[, priorOrderFields],
                                table[, refField], table[, barcodeField], wellRows,
                                wellColumns, na.last = T),]
    }
    return(tr[,!colnames(tr) %in% c("wellRows", "wellColumns")])
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
