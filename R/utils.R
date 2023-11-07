getFileViaRegex <- function(path, pattern) {
  available_files = list.files(path = path,
                               pattern = pattern,
                               full.names = TRUE)
  if (length(available_files) == 0) {
    stop("Could not find file in directory ", path, " Regex is ", pattern)
  } else if (length(available_files) > 1) {
    stop("Too many files matched in directory ",
         path,
         " Regex is ",
         pattern)
  }
  return(available_files)
}

readTargetInfoFile <-
  function(file,
           filterFields = c("Genotype", "ExPlateBarcode", "ExPlateWell")) {
    obj = {

    }
    if (is.null(dim(file))) {
      obj$table = read.csv(file, header = T, stringsAsFactors = F)
    } else{
      obj$table = file
    }

    refField = "RefType"
    barcodeField = "ExPlateBarcode"
    wellField = "ExPlateWell"
    expectedFields = c("TargetID", "SampleType", refField)

    if (length(intersect(expectedFields, names(obj$table))) != length(expectedFields)) {
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
    rownames(obj$table) = obj$table$TargetID
    obj$table = subset(obj$table, select = -c(TargetID))
    obj$table = obj$table[, names(obj$table) %in% union(expectedFields, filterFields)]

    obj$getReferences <- function() {
      refTable =  obj$table[!is.na(obj$table$RefType) &
                              obj$table$RefType != "",]
      table = subset(refTable, select = -c(SampleType))
      table$ExPlateBarcode = sapply(table$ExPlateBarcode, as.character)
      return(sortByRefTypeBarcodeWell_2(table))
    }

    obj$getSamples = function() {
      refTable =  obj$table[is.na(obj$table$RefType) |
                              obj$table$RefType == "",]
      table = subset(refTable, select = -c(RefType))
      table$ExPlateBarcode = sapply(table$ExPlateBarcode, as.character)
      return(sortByRefTypeBarcodeWell_2(table, refField = 'SampleType'))
    }
    return(obj)
  }

sortByRefTypeBarcodeWell_2 <-
  function(table,
           refField = "RefType",
           barcodeField = "ExPlateBarcode",
           wellField = "ExPlateWell",
           priorOrderFields = c()) {
    if (is.null(refField)) {
      stop("refField must be provided")
    }

    if (!is.null(barcodeField) &&
        !any(names(table) == barcodeField)) {
      stop(paste("Barcode field '", refField, "' is missing in table", sep = ""))

    }
    if (!is.null(wellField) && !any(names(table) == wellField)) {
      stop(paste("Well field '", refField, "' is missing in table", sep = ""))

    }

    table[, refField][table[, refField] == ""] = NA

    if (!is.null(refField) &&
        !is.null(barcodeField) && !is.null(wellField)) {
      wellRows = sapply(as.character(table[, wellField]), function(x) {
        return(strsplit(gsub("([0-9]+)", ",\\1", x), ",")[[1]][1])
      })
      wellColumns = as.integer(sapply(as.character(table[, wellField]), function(x) {
        return(strsplit(gsub("([0-9]+)", ",\\1", x), ",")[[1]][2])
      }))

      tableWithWell  = cbind(table, wellRows, wellColumns)
      if (length(priorOrderFields) == 0) {
        tr = tableWithWell[order(table[, refField], table[, barcodeField], wellRows, wellColumns, na.last = T),]
      } else{
        for (f in priorOrderFields) {
          if (!any(f == colnames(table))) {
            stop("Field ", f, " does not exist in table")
          }
        }
        tr = tableWithWell[order(table[, priorOrderFields], table[, refField], table[, barcodeField], wellRows, wellColumns, na.last = T),]
      }
      return(tr[,!colnames(tr) %in% c("wellRows", "wellColumns")])
    } else if (!is.null(refField) && !is.null(barcodeField)) {
      return(table[order(refFieldEmpty, table[, barcodeField]),])
    } else{
      return(table[order(refFieldEmpty),])
    }
  }
